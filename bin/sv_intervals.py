#!/usr/bin/env python3
"""Define the structural-variant interval set ONCE, for genotyping.

WHY THE ARM IS BEING REBUILT AROUND THIS. The SV block used to take the union
of per-sample caller breakpoints, cluster them, and treat the samples that
contributed a call as the carriers. Four tests in assoc/SV_SCATTER_RESULTS.md
showed what that produces: tightening the clustering tolerance to 10 bp splits
each row into sub-clusters that are 80.4% singletons, so two isolates carrying
the same deletion rarely report the same breakpoint even within 10 bp. The row
was a bag of individually placed calls rather than one event genotyped across
samples, and on a tree it looked like one independent origin per carrier.

The same tests showed the way out. Carriers earned by measuring depth over a
fixed interval were markedly more clade-consistent than caller-called ones at
matched carrier counts, 0.75-0.80 against 0.91-1.00, because a fixed interval
measured identically in every sample is reproducible where a per-sample
breakpoint is not. So the interval set is defined first, here, and every
sample is then genotyped at every interval.

TWO SOURCES, AND THE FIRST IS THE ONE THE CALLERS CANNOT SEE.

  graph   The CX333 deconstruct VCF already carries 4,805 deletions of 50 bp
          or more against H37Rv, genotyped across all 332 panel genomes. Those
          are exactly the deletions absorbed into the panel -- an isolate
          aligned to a reference that carries one sees contiguous sequence and
          no caller ever emits a call, however sensitive. Taking them from the
          graph also gives, per interval, the exact set of references that
          carry it, so an isolate's inherited state is a lookup rather than an
          inference from a projection failing.

  caller  Deletions some isolate reported that the panel does not contain.
          These are kept, but as INTERVALS to genotype, not as carrier sets.

An interval from the caller set is dropped when a graph interval already
covers it under this project's tolerance -- max(200 bp, 20% of length) in
position, 1.5x in length -- so the two sources do not produce duplicate rows.

INSERTIONS ARE NOT HERE, and that is a scope decision rather than an
oversight. Depth across an insertion point is the same whether or not the
sample carries the inserted sequence, so the instrument this catalogue exists
to feed cannot genotype them. They stay in the caller-derived block, which is
presence-only and labelled as such.
"""
import argparse, bisect, collections, csv, gzip, os, sys

TOL_CAP = 20000


def tol(length):
    return min(max(200, int(0.2 * max(length, 1))), TOL_CAP)


def same_event(p1, l1, p2, l2):
    if abs(p1 - p2) > tol(max(l1, l2)):
        return False
    lo, hi = sorted((max(l1, 1), max(l2, 1)))
    return hi / lo <= 1.5


def overlap_frac(a, b):
    """Reciprocal overlap of two deletions {start, end}: shared bases over the
    longer one's length."""
    shared = min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
    longer = max(a["end"] - a["start"] + 1, b["end"] - b["start"] + 1)
    return max(0, shared) / max(1, longer)


MIN_OVERLAP = 0.5


def same_deletion(s1, l1, s2, l2):
    """THE ONE RULE for "these two deletions are the same event", used to
    cluster the graph's deletions, to decide whether a caller deletion is
    already in the catalogue, and (merge_cohort_vcf.py) whether the catalogue
    supersedes a caller row. Starts are FIRST DELETED BASES, 1-based.

    Near in position and length (same_event), AND sharing at least half of the
    longer deletion's bases. The position tolerance alone -- max(200 bp, 20%
    of length) -- let two 58 bp deletions 150 bp apart match without sharing a
    base; the clustering had the overlap test (audit P4P5-5) but the two
    matching steps did not, so 65 scale200 caller deletions, 25 of them
    sharing no base with the interval, were dropped as covered (review 2,
    R2-GENO-2)."""
    if not same_event(s1, l1, s2, l2):
        return False
    return overlap_frac({"start": s1, "end": s1 + max(l1, 1) - 1},
                        {"start": s2, "end": s2 + max(l2, 1) - 1}) >= MIN_OVERLAP


def main():
    ap = argparse.ArgumentParser()
    # no default: it was CX333's; p5_svgt.sh passes the build's
    ap.add_argument("--graph-vcf", required=True,
                    help="<build>/assets/graph_collapsed.vcf.gz")
    ap.add_argument("--sv-matrix", default="",
                    help="a cohort's sv_matrix.tsv, for caller-derived "
                         "intervals the panel does not contain")
    ap.add_argument("--is6110-gff",
                    default="data/annotation/H37Rv_IS6110.pansn.gff",
                    help="H37Rv's own element copies, as breakpoint landmarks")
    ap.add_argument("--is6110-keys", default="",
                    help="a cohort's IS6110 key table, so sites the cohort "
                         "carries but H37Rv does not are landmarks too")
    ap.add_argument("--is6110-window", type=int, default=500,
                    help="how close a breakpoint must be to a landmark to be "
                         "flagged IS6110-proximal")
    ap.add_argument("--min-caller-carriers", type=int, default=1,
                    help="take a caller deletion into the catalogue only if at "
                         "least this many isolates called it (the matrix's "
                         "n_alt). Below it, the deletion stays in the VCF "
                         "presence-only, from the caller matrix, flagged "
                         "UNCATALOGUED; it is not lost")
    ap.add_argument("--min-len", type=int, default=50)
    ap.add_argument("--max-len", type=int, default=100000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    op = gzip.open if a.graph_vcf.endswith(".gz") else open
    panel, rows = [], []
    with op(a.graph_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                panel = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            if "," in alt or alt.startswith("<"):
                continue
            d = len(ref) - len(alt)
            if d < a.min_len or d > a.max_len:
                continue
            pos = int(f[1])
            # VCF deletions are left-anchored: the deleted bases are the REF
            # allele minus its first base.
            start, end = pos + len(alt), pos + len(ref) - 1
            carriers = [panel[i] for i, g in enumerate(f[9:])
                        if g.split(":")[0] == "1"]
            rows.append(dict(source="graph", svtype="DEL", start=start,
                             end=end, svlen=d, n_ref_carriers=len(carriers),
                             ref_carriers=",".join(sorted(carriers))))
    print(f"  graph: {len(rows):,} deletion records {a.min_len}-{a.max_len} "
          f"bp over {len(panel)} panel genomes")
    # COLLAPSE FIRST. vg deconstruct emits one record per bubble and nested
    # bubbles describe the same deletion several times, so the raw records
    # carry thousands of duplicates -- 3,243 of them shared a start and a
    # length outright. Cluster them on the project's own tolerance and union
    # the carrier sets, or the catalogue counts one event many times and the
    # genotyper probes the same interval repeatedly.
    rows.sort(key=lambda r: (r["start"], r["svlen"]))
    merged = []
    for r in rows:
        hit = None
        for c in reversed(merged[-50:]):
            if same_deletion(c["start"], c["svlen"], r["start"], r["svlen"]):
                hit = c
                break
        if hit is None:
            r["_car"] = set(r.pop("ref_carriers").split(",")) - {""}
            merged.append(r)
        else:
            # THE CLUSTER KEEPS ITS FIRST MEMBER'S COORDINATES (audit P4P5-5).
            # It used to widen to the union of its members and compare the next
            # deletion with that widened span, so a chain of overlapping losses
            # in a tandem repeat grew an interval longer than most of its
            # members -- 125 of 278 multi-member intervals beyond the 1.5x
            # length tolerance -- and a reference carrying a 116 bp loss was
            # listed as a carrier of a 174 bp interval it only half covers.
            # Every member is matched against the fixed representative and must
            # overlap it by half of each one's length, so a carrier's own
            # deletion covers at least half of the interval it is listed on.
            # The position tolerance alone lets two 58 bp deletions 150 bp
            # apart merge without sharing a base. start, end and svlen stay
            # consistent (review 5.7).
            hit["_car"] |= set(r["ref_carriers"].split(",")) - {""}
    for r in merged:
        r["ref_carriers"] = ",".join(sorted(r["_car"]))
        r["n_ref_carriers"] = len(r.pop("_car"))
    print(f"    collapsed to {len(merged):,} distinct intervals")
    rows = merged
    if rows:
        nc = [r["n_ref_carriers"] for r in rows]
        print(f"    references carrying one: median {sorted(nc)[len(nc)//2]}, "
              f"{sum(1 for x in nc if x == 0):,} carried by none, "
              f"{sum(1 for x in nc if x == 1):,} by exactly one")

    n_add = n_dup = n_few = 0
    if a.sv_matrix:
        idx = sorted((r["start"], r["svlen"], i) for i, r in enumerate(rows))
        starts = [x[0] for x in idx]
        for r in csv.DictReader(open(a.sv_matrix, newline=""), delimiter="\t"):
            if r.get("svtype") != "DEL":
                continue
            try:
                pos, L = int(r["h37rv_pos"]), int(r["svlen"])
            except (ValueError, KeyError):
                continue
            if L < a.min_len or L > a.max_len:
                continue
            if int(r.get("n_alt") or 0) < a.min_caller_carriers:
                n_few += 1
                continue
            lo = bisect.bisect_left(starts, pos - TOL_CAP)
            hi = bisect.bisect_right(starts, pos + TOL_CAP)
            # the matrix's h37rv_pos is the anchor; the first deleted base
            # is the next one, which is what the catalogue's start records
            if any(same_deletion(pos + 1, L, s, l) for s, l, _ in idx[lo:hi]):
                n_dup += 1
                continue
            rows.append(dict(source="caller", svtype="DEL", start=pos + 1,
                             end=pos + L, svlen=L, n_ref_carriers=0,
                             ref_carriers=""))
            n_add += 1
        print(f"  caller: {n_add:,} intervals added, {n_dup:,} already covered "
              f"by a graph interval"
              + (f", {n_few:,} with fewer than {a.min_caller_carriers} carriers "
                 f"left presence-only" if n_few else ""))

    # ---- IS6110 proximity.
    # WHY IT IS ON EVERY ROW. At ppe38 the project measured 204 SV records in
    # 14.8 kb of which 103 are IS6110-derived, 13.7 per kb against a background
    # median of 5. One element insertion enters the callset as an INS where the
    # sample carries it, a DEL where the reference does, and a COMPLEX block
    # substitution where there is flanking rearrangement -- three
    # representations of one event, beside the dedicated IS6110 arm's fourth.
    # Flagging the intervals whose ends sit on a landmark is what lets the two
    # arms be cross-checked instead of double-counting.
    marks = []
    if a.is6110_gff and os.path.exists(a.is6110_gff):
        for line in open(a.is6110_gff):
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) > 4 and f[3].isdigit():
                marks += [int(f[3]), int(f[4])]
    n_ref_marks = len(marks) // 2
    if a.is6110_keys and os.path.exists(a.is6110_keys):
        for r in csv.DictReader(open(a.is6110_keys), delimiter="\t"):
            if r.get("frame") == "h37rv" and (r.get("h37rv_pos") or "").isdigit():
                p0 = int(r["h37rv_pos"])
                marks += [p0, p0 + 1355]
    marks = sorted(set(marks))

    def near(p):
        i = bisect.bisect_left(marks, p)
        return any(0 <= j < len(marks) and abs(marks[j] - p) <= a.is6110_window
                   for j in (i - 1, i))

    # ---- support tier. Complete assemblies are a different instrument from
    # short-read SV calling, and measured against the cohort tree they behave
    # like one: median independent origins per carrier falls monotonically with
    # the number of panel assemblies that carry the event -- 1.00 where none
    # do, 0.77 at one, 0.67 at two to five, 0.40 at six or more -- and the
    # difference survives matching on carrier count. So the tier is written on
    # every row rather than left to be recomputed.
    for r in rows:
        n = int(r["n_ref_carriers"])
        r["support_tier"] = ("A_multi_assembly" if n >= 2 else
                             "B_single_assembly" if n == 1 else
                             "C_caller_only")
    tc = collections.Counter(r["support_tier"] for r in rows)
    for k in sorted(tc):
        print(f"    tier {k:<20}{tc[k]:>6,}")

    n_prox = 0
    for r in rows:
        r["is6110_prox"] = int(bool(marks) and (near(r["start"]) or
                                                near(r["end"])))
        n_prox += r["is6110_prox"]
    if marks:
        print(f"  IS6110 landmarks: {n_ref_marks} H37Rv copies -> "
              f"{len(marks):,} coordinates; {n_prox:,} of {len(rows):,} "
              f"intervals have an end within {a.is6110_window} bp "
              f"({n_prox / max(1, len(rows)):.1%})")

    rows.sort(key=lambda r: (r["start"], r["svlen"]))
    seen = collections.Counter()
    for r in rows:
        base = f"svi:{r['svtype']}:{r['start']}:{r['svlen']}"
        seen[base] += 1
        # A suffix only where two intervals genuinely share start and length
        # after collapsing; the id has to be unique because everything
        # downstream keys on it.
        r["interval"] = base if seen[base] == 1 else f"{base}.{seen[base]}"
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "interval", "source", "svtype", "start", "end", "svlen",
            "support_tier", "is6110_prox", "n_ref_carriers", "ref_carriers"])
        w.writeheader()
        w.writerows(rows)
    dup = len(rows) - len({r["interval"] for r in rows})
    if dup:
        sys.exit(f"FATAL: {dup} duplicate interval ids")
    # The UNION of the intervals, not the sum of their lengths: they overlap
    # heavily, and summing makes the catalogue look as though it covers most
    # of the genome.
    spans = sorted((r["start"], r["end"]) for r in rows)
    union, cs, ce = 0, None, None
    for s0, e0 in spans:
        if cs is None:
            cs, ce = s0, e0
        elif s0 <= ce + 1:
            ce = max(ce, e0)
        else:
            union += ce - cs + 1; cs, ce = s0, e0
    if cs is not None:
        union += ce - cs + 1
    print(f"  {len(rows):,} intervals total, covering {union:,} bp of H37Rv "
          f"({union / 4411532:.1%}) once overlaps are merged")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
