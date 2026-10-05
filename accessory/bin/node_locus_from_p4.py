#!/usr/bin/env python3
"""node -> accessory locus, taken from p4's own anchor, for LEVEL 2.

WHY NOT SEQUENCE CONTAINMENT. The census matches a node to a locus by
containing the node's sequence in the locus's representative allele, because a
node-frame record carries no H37Rv coordinate in the merged VCF. That works for
the census's purpose but it cannot support level 2: of 20,299 nodes bearing a
node-frame record, only 626 are >= 20 bp -- a SNP node is one base -- and of
those only 12.7% match exactly one of 802 loci, the rest matching several loci
that share repeat content. 1,264 of 22,193 records placed, most ambiguously.

WHAT THIS USES INSTEAD. p4_place.py already did this work per sample and kept
the answer. When it emits an off-path record it has the record's H37Rv anchor in
hand, and for a record beyond --near-tol of the path it calls acc_name(anchor),
which is the nearest catalogued accessory locus within --acc-tol. That is a
coordinate link, it is unambiguous, and it is already written to column
acc_locus of every <sample>.placed.tsv. This script does nothing but collect
those assignments across the cohort and check they agree.

A record labelled off_path_near gets no locus by design: p4 measured 96.3% of
off-path calls within 50 bp of the path, which are small-indel interiors sitting
inside the same reading frame as their anchor. Those are ordinary small variants
and level 2 does not apply to them. Only off_path_accessory records are inside
an insert, and only they can be conditioned on carrying it.
"""
import argparse, collections, csv, glob, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p4dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.p4dir, "*.placed.tsv")))
    if not files:
        sys.exit(f"FATAL: no *.placed.tsv under {a.p4dir}")

    # node -> Counter of loci it was assigned, across samples
    assign = collections.defaultdict(collections.Counter)
    anchor = {}
    n_acc = n_near = 0
    for i, f in enumerate(files, 1):
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("frame") != "node":
                    continue
                reg = r.get("region", "")
                if reg == "off_path_near":
                    n_near += 1
                    continue
                if reg != "off_path_accessory":
                    continue
                n_acc += 1
                nd = r.get("node") or ""
                if not nd:
                    continue
                assign[nd][r.get("acc_locus") or ""] += 1
                anchor.setdefault(nd, r.get("h37rv_pos") or "")
        if i % 200 == 0:
            print(f"  {i}/{len(files)} samples, {len(assign):,} nodes")

    rows = []
    disagree = 0
    for nd, c in sorted(assign.items(), key=lambda t: int(t[0])):
        named = {k: v for k, v in c.items() if k}
        if len(named) > 1:
            disagree += 1
        top = max(named.items(), key=lambda t: t[1])[0] if named else ""
        rows.append(dict(node=nd, locus=top, n_samples=sum(c.values()),
                         n_named=sum(named.values()),
                         n_distinct_loci=len(named),
                         h37rv_anchor=anchor.get(nd, "")))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    placed = sum(1 for r in rows if r["locus"])
    print(f"\n  off_path_accessory cells {n_acc:,}   off_path_near (not level 2) "
          f"{n_near:,}")
    print(f"  {len(rows):,} accessory nodes, {placed:,} carry a locus assignment"
          f" ({100*placed/max(1,len(rows)):.1f}%)")
    print(f"  nodes assigned different loci in different samples: {disagree:,}")
    print(f"  distinct loci referenced: "
          f"{len({r['locus'] for r in rows if r['locus']}):,}")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
