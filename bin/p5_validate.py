#!/usr/bin/env python3
"""P5 validation: does the merged matrix recover lineage, without being told it?

Real isolates have no variant-level truth, so the pilot's endpoint is truth-free
(section 4.2). This is the strongest of those checks: reference selection was by
SNP-profile distance and never by lineage label, so the tb-profiler assignments
are information the pipeline never saw. If samples of the same lineage are each
other's nearest neighbours in the matrix, the merge preserved sample identity;
if they are not, the matrix has scrambled samples in a way no amount of staring
at it would reveal.

Distance is over ALT/REF states only. ABSENT and NOCALL cells are excluded
pairwise rather than imputed, because imputing them to REF is exactly the
encoding hazard T6 measured, and it would make divergent samples look artificially
close at precisely the positions where they differ most.
"""
import numpy as np
import argparse, collections, csv, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="",
                    help="the dense matrix; legacy input, for cohorts built "
                         "before the states array existed")
    ap.add_argument("--states-array", default="",
                    help="directory holding states.u8.npy (p5_matrix.py). "
                         "Needs --keys and --sites; replaces --matrix")
    ap.add_argument("--keys", default="")
    ap.add_argument("--sites", default="")
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--region", default="core",
                    help="restrict to one region, or 'all'")
    ap.add_argument("--out", default="refbias/p5/validation.tsv")
    a = ap.parse_args()

    meta = {r["sample"]: r for r in
            csv.DictReader(open(a.refmap, newline=""), delimiter="\t")}
    if a.states_array:
        # Rows of the states array, selected by sites.tsv. The array's codes
        # already put ALT at 0 and REF at 1; everything else becomes the
        # not-comparable code. At 10,000 isolates the selected core rows are
        # a few GB of uint8, against the dense text matrix of about 23 GB.
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from p5_states_io import open_array
        keys = list(csv.DictReader(open(a.keys, newline=""), delimiter="\t"))
        A, samples = open_array(a.states_array, keys)
        sel = [int(r["key_index"]) for r in
               csv.DictReader(open(a.sites, newline=""), delimiter="\t")
               if a.region == "all" or r["region"] == a.region]
        CODE = {"ALT": 0, "REF": 1}
        NA = 255
        n_used = len(sel)
        if n_used:
            sub = np.asarray(A[np.asarray(sel, dtype=np.int64)])
            sub = np.where(sub <= 1, sub, NA).astype(np.uint8)
            M = sub.T                                   # samples x sites
        else:
            M = np.zeros((len(samples), 0), dtype=np.uint8)
    else:
        with open(a.matrix, newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            hdr = next(rd)
            fixed = hdr.index("n_nocall") + 1
            # skip the artefact-flag columns: they sit between the counts and the
            # samples, and treating them as samples adds two pseudo-isolates
            while fixed < len(hdr) and hdr[fixed] in ("panel_af", "h37rv_minor"):
                fixed += 1
            samples = hdr[fixed:]
            ri = hdr.index("region")
            # ONE uint8 ARRAY, AND THE PAIRWISE PASS VECTORISED. This check is
            # O(samples^2 x sites) by nature, and as a Python triple loop over
            # lists of strings that is 997 x 996 x 173,685 comparisons at 997
            # isolates -- 1.7e11 of them. gwas1000 spent six and a half hours here
            # and had not finished one of the four steps after the matrix. The
            # same computation over a coded array, one sample against all the
            # others at a time, is minutes.
            CODE = {"ALT": 0, "REF": 1}          # everything else is not comparable
            NA = 255
            rows_used = []
            for row in rd:
                if a.region != "all" and row[ri] != a.region:
                    continue
                rows_used.append([CODE.get(v, NA) for v in row[fixed:]])
            n_used = len(rows_used)
            M = (np.asarray(rows_used, dtype=np.uint8).T if rows_used
                 else np.zeros((len(samples), 0), dtype=np.uint8))

    if n_used == 0:
        print("no sites selected", file=sys.stderr); return 1

    # THE WHOLE DISTANCE MATRIX AS THREE MATRIX PRODUCTS, not a pass per
    # sample. Vectorising the inner loop took this from hours to seconds at 200
    # isolates, but it was still O(samples^2 x sites) of numpy element work --
    # 5.4e13 at 10,000 isolates, which is a day. The same element count through
    # BLAS is two orders of magnitude faster, and it is exactly a matrix
    # product: with A the ALT indicator and R the REF indicator,
    #     comparable = (A+R)(A+R)^T      differing = A R^T + R A^T
    # Sites are processed in blocks so the indicators never have to exist for
    # the whole key set at once; only the two samples x samples accumulators do,
    # which at 10,000 isolates is 400 MB each.
    idx = {s: i for i, s in enumerate(samples)}
    ns = len(samples)
    both_n = np.zeros((ns, ns), dtype=np.float32)
    diff_n = np.zeros((ns, ns), dtype=np.float32)
    BLK = max(1, int(2e8 // max(ns, 1)))
    for s0 in range(0, M.shape[1], BLK):
        blk = M[:, s0:s0 + BLK]
        A = (blk == CODE["ALT"]).astype(np.float32)
        R = (blk == CODE["REF"]).astype(np.float32)
        C = A + R
        both_n += C @ C.T
        diff_n += A @ R.T + R @ A.T
    DIST = np.where(both_n > 0, diff_n / np.maximum(both_n, 1), 1.0)

    def dist_row(x):
        """(distance, comparable sites) from x to every sample."""
        i = idx[x]
        return DIST[i], both_n[i]

    print(f"  {n_used} sites in region '{a.region}', {len(samples)} samples\n")
    rows = []
    correct = 0
    for s in samples:
        dv, _ = dist_row(s)
        dv[idx[s]] = np.inf
        j = int(np.argmin(dv))
        nn_d, nn = float(dv[j]), samples[j]
        lin = (meta.get(s, {}).get("sample_lineage") or "").strip()
        nnlin = (meta.get(nn, {}).get("sample_lineage") or "").strip()
        same = (lin and nnlin and lin == nnlin)
        # a lineage with only one representative cannot have a same-lineage
        # neighbour, so it is not scored either way
        n_same_lin = sum(1 for t in samples if t != s
                         and (meta.get(t, {}).get("sample_lineage") or "").strip() == lin)
        scored = bool(lin) and n_same_lin > 0
        if scored and same:
            correct += 1
        rows.append(dict(sample=s, lineage=lin or "-", nearest=nn,
                         nearest_lineage=nnlin or "-",
                         distance=round(nn_d, 5), scorable=("yes" if scored else "no"),
                         nearest_same_lineage=("yes" if same else "no")))
        print(f"  {s:<17s}{(lin or '-'):<11s}-> {nn:<17s}"
              f"{(nnlin or '-'):<11s}d={nn_d:.4f}  "
              f"{'MATCH' if same else ('unscored' if not scored else 'MISMATCH')}")

    scorable = [r for r in rows if r["scorable"] == "yes"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"\n  nearest neighbour shares the lineage label: "
          f"{correct}/{len(scorable)} scorable samples"
          + (f" ({100*correct/len(scorable):.1f}%)" if scorable else ""))
    print(f"    {len(rows)-len(scorable)} unscorable (sole representative of "
          f"their lineage, or unlabelled)")
    print(f"\n  the pipeline never saw these labels: selection was by SNP")
    print(f"  distance, and stage 1 measured label-gating as 0.6 points worse")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
