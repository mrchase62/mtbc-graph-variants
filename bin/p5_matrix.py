#!/usr/bin/env python3
"""P5 step 3: assemble the cohort matrix and check it without truth.

Emits one row per key and one column per sample, with the four states kept
distinct rather than collapsed into a missing value.

The validation is section 4.2's: a phylogeny built from the matrix should recover
the tb-profiler lineage assignments. That information was never shown to the
pipeline -- selection was by SNP distance, not by label -- so agreement is a real
check, and it is the one most likely to catch a merge that silently scrambles
samples, which is the failure mode a matrix cannot reveal by inspection.
"""
import numpy as np
import argparse, bisect, collections, csv, gzip, os, sys


def rd(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--keys", required=True)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-carriers", type=int, default=1)
    ap.add_argument("--graph-vcf",
                    default="graphs/CX333.s10k.k23.K15/all_variants.nolab.vcf.gz",
                    help="source of panel allele frequency, for the "
                         "reference-artefact flag")
    ap.add_argument("--no-dense", action="store_true",
                    help="skip writing matrix.tsv. The dense text form is "
                         "keys x samples: 783 MB at 997 isolates and 26.5 GB "
                         "at 10,000, of which about 98% is REF or NOCALL and "
                         "recoverable from the sparse per-sample files. The "
                         "counts and the summary are still produced.")
    a = ap.parse_args()

    meta = {r["sample"]: r for r in rd(a.refmap)}
    keys = rd(a.keys)

    # ONE uint8 ARRAY, NOT 997 DICTS OF STRINGS. The previous version held
    # every sample's whole state table in memory at once as
    # {key: (state, allele)}. At 997 isolates and 173,685 keys that is 173
    # million Python dict entries, each a string key and a tuple of two
    # strings, and gwas1000's p5finish died OUT_OF_MEMORY at the 8 GB ceiling
    # 98 seconds in. The same information is 173,685 x 997 single bytes --
    # 173 MB -- once the four states are encoded, and the allele field is not
    # used by this script at all.
    # One reader for every consumer; the sparse format lives in p5_states_io.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from p5_states_io import CODE, NAME, load_states, samples_in
    samples = samples_in(a.dir, meta)
    if not samples:
        print("no state files found", file=sys.stderr); return 1
    M = load_states(a.dir, keys, samples)

    # Reference-artefact flag, as a general screen rather than a list of the 437
    # sites section 21.1 happened to find. A site where most of the 333-genome
    # panel carries the non-H37Rv allele is a position where H37Rv is
    # unrepresentative of the species, not a variant 23 isolates converged on --
    # and 95% of the ubiquitous sites were of exactly that kind. Flagging by
    # panel frequency catches the same class at any carrier count, including
    # sites carried by two isolates that would never have drawn attention.
    panel_af = {}
    if os.path.exists(a.graph_vcf):
        op = gzip.open if a.graph_vcf.endswith(".gz") else open
        with op(a.graph_vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 8)
                try:
                    pos = int(f[1])
                except (ValueError, IndexError):
                    continue
                for kv in f[7].split(";"):
                    if kv.startswith("AF="):
                        try:
                            v = max(float(x) for x in kv[3:].split(","))
                        except ValueError:
                            break
                        if v > panel_af.get(pos, 0):
                            panel_af[pos] = v
                        break
        print(f"  panel allele frequencies loaded for {len(panel_af)} positions")
    else:
        print(f"  WARNING: no graph VCF at {a.graph_vcf}; the "
              f"reference-artefact flag will be empty", file=sys.stderr)

    # Counts per row straight off the array, and no row accumulation: the
    # rows are written as they are formed, so peak memory is the array and
    # one row rather than the whole matrix twice over.
    n_alt_all = (M == CODE["ALT"]).sum(axis=1)
    keep_mask = n_alt_all >= a.min_carriers
    dropped_noalt = int((~keep_mask).sum())
    n_ref_all = (M == CODE["REF"]).sum(axis=1)
    n_abs_all = (M == CODE["ABSENT"]).sum(axis=1)
    n_noc_all = (M == CODE["NOCALL"]).sum(axis=1)

    import contextlib
    with (open(os.devnull, "w") if a.no_dense
          else open(a.out, "w", newline="")) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["key", "frame", "region", "kind", "h37rv_pos", "node",
                    "acc_locus", "canonical_ref", "n_alt", "n_ref",
                    "n_absent", "n_nocall", "panel_af", "h37rv_minor"] + samples)
        n_flag = 0
        n_kept = 0
        for i, k in enumerate(keys):
            if not keep_mask[i]:
                continue
            n_kept += 1
            col = [NAME[v] for v in M[i]]
            n_alt = int(n_alt_all[i])
            af = ""
            minor = ""
            if k["frame"] == "h37rv" and k["h37rv_pos"]:
                v = panel_af.get(int(k["h37rv_pos"]))
                if v is not None:
                    af = round(v, 4)
                    minor = 1 if v >= 0.5 else 0
                    n_flag += minor
            w.writerow([k["key"], k["frame"], k["region"], k["kind"],
                        k["h37rv_pos"], k["node"], k["acc_locus"],
                        k["canonical_ref"], n_alt, int(n_ref_all[i]),
                        int(n_abs_all[i]), int(n_noc_all[i]), af, minor] + col)

    print(f"  {n_kept} sites x {len(samples)} samples")
    print(f"    flagged h37rv_minor (panel AF >= 0.5, H37Rv unrepresentative): "
          f"{n_flag}")
    print(f"    keys dropped with no ALT carrier: {dropped_noalt}")
    tot = collections.Counter()
    for st, code in CODE.items():
        tot[st] = int((M[keep_mask] == code).sum())
    n = sum(tot.values())
    print(f"\n  cell states across the matrix ({n} cells)")
    for st in ("ALT", "REF", "ABSENT", "NOCALL"):
        print(f"    {st:<8s}{tot[st]:>9d}  {100*tot[st]/n:>5.1f}%")
    print(f"\n  ABSENT is the state a plain VCF cannot express: "
          f"{tot['ABSENT']} cells")
    print(f"    T6 measured 508,039 such cells being written as REF or missing")

    byreg = collections.Counter(k["region"] for i, k in enumerate(keys)
                                if keep_mask[i])
    print(f"\n  sites by region")
    for r, c in byreg.most_common():
        print(f"    {r:<22s}{c:>7d}")
    sing = int(((n_alt_all == 1) & keep_mask).sum())
    ubiq = int(((n_alt_all == len(samples)) & keep_mask).sum())
    print(f"\n  {sing} singleton sites, {ubiq} carried by all {len(samples)} "
          f"samples")
    print(f"    a site carried by every sample is a reference artefact candidate,")
    print(f"    not a shared variant, since H37Rv is the canonical frame")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
