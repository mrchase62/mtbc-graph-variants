#!/usr/bin/env python3
"""Give each catalogue locus the SET of anchor positions seen across the panel.

is6110/bin/is6110_copy_matrix.py groups copies by single-linkage chaining at a 200 bp
window and assigns each locus the MEDIAN of its cluster. Single linkage is
unbounded -- 40 of 432 clusters exceed the 200 bp window and the widest spans
868 bp -- so an individual genome's junction can sit far from the median. 122 of
432 clusters are wider than 30 bp, which is the clip test's tolerance, and the
dispersion is strongly lineage-dependent: 33.8% of lineage-1 anchors are more
than 30 bp from their locus against 5.7% for lineage 2. That ordering tracks
recall (0.670 vs 0.924) and is the cause of the lineage-1 gap.

WHY NOT JUST USE EACH GENOME'S OWN ANCHOR
Because it is circular. The anchor comes from that genome's assembly, and a new
short-read sample has no assembly -- the detector would be handed the answer.

WHAT IS LEGITIMATE
The catalogue is built from the panel, so the positions where panel genomes place
their junctions at a locus are catalogue knowledge, available for any new sample.
Each locus therefore carries a small set of candidate positions -- the distinct
anchors observed, merged at fine resolution -- and the detector tests each. That
is bounded in a way the +/-400 bp peak search was not: a handful of specific
candidates rather than every offset, which is why it should not carry the same
precision cost.
"""
import argparse, collections, csv, glob, os, statistics, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--copies", default="loci/is6110_cx333/copies/*.tsv")
    ap.add_argument("--loci", default="refbias/t11/loci_clean.tsv")
    ap.add_argument("--member-window", type=int, default=200,
                    help="anchors within this of the locus belong to it "
                         "(matches is6110/bin/is6110_copy_matrix.py)")
    ap.add_argument("--merge", type=int, default=10,
                    help="anchors within this of each other are one candidate")
    ap.add_argument("--min-support", type=int, default=2,
                    help="genomes needed before a candidate position is kept")
    ap.add_argument("--out", default="refbias/t11/anchor_sets.tsv")
    a = ap.parse_args()

    anchors = collections.defaultdict(list)
    for p in glob.glob(a.copies):
        for r in csv.DictReader(open(p), delimiter="\t"):
            v = r.get("left_ref")
            if v in (None, "", "NA"):
                continue
            anchors[r["sample"]].append(int(v))
    allpos = sorted(p for s in anchors for p in anchors[s])
    loci = [int(r["pos"]) for r in csv.DictReader(open(a.loci), delimiter="\t")]

    rows = []
    for L in loci:
        near = [p for p in allpos if abs(p - L) <= a.member_window]
        # merge nearby anchors into candidate positions
        cands, cur = [], []
        for p in sorted(near):
            if cur and p - cur[-1] > a.merge:
                cands.append(cur); cur = []
            cur.append(p)
        if cur:
            cands.append(cur)
        keep = [(int(statistics.median(c)), len(c)) for c in cands
                if len(c) >= a.min_support]
        keep.sort(key=lambda t: -t[1])
        if not keep:
            keep = [(L, 0)]
        rows.append(dict(locus=f"IS_{L}", pos=L, n_candidates=len(keep),
                         candidates=";".join(f"{p}:{n}" for p, n in keep),
                         offsets=";".join(str(p - L) for p, _ in keep),
                         max_offset=max(abs(p - L) for p, _ in keep)))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    nc = [r["n_candidates"] for r in rows]
    mo = sorted(r["max_offset"] for r in rows)
    print(f"  {len(rows)} loci")
    print(f"    candidates per locus: median {sorted(nc)[len(nc)//2]}, "
          f"max {max(nc)}, mean {sum(nc)/len(nc):.1f}")
    print(f"    loci with >1 candidate: {sum(1 for x in nc if x > 1)}")
    print(f"    furthest candidate from the catalogued position: "
          f"median {mo[len(mo)//2]} bp, 90th {mo[int(.9*len(mo))]} bp, max {mo[-1]} bp")
    print(f"    loci where a candidate sits >30 bp away (invisible before): "
          f"{sum(1 for x in mo if x > 30)}")
    print(f"\n  written: {a.out}")


if __name__ == "__main__":
    sys.exit(main())
