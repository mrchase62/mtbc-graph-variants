#!/usr/bin/env python3
"""Pick the nearest panel genome for a real isolate, from its own H37Rv calls.

Stage 1 established that selecting by SNP-profile distance beats selecting by
tb-profiler label, and that a label-gated search is actively worse because in
16.1% of genomes the nearest reference lies outside the label. So this ignores
the label and compares profiles directly.

The isolate's profile is its called SNPs in H37Rv coordinates -- the thing a
standard pipeline already produces. The panel's profiles come from the graph's
SNP matrix in the same coordinates. Distance is the count of sites where the two
disagree, over sites the panel has genotyped.

No self-exclusion and no near-clone exclusion here, unlike stage 1. Those existed
because a panel genome is its own nearest neighbour and the panel contains
near-clones, which would have made the leave-one-out measurement meaningless. A
real isolate is not in the panel, so its nearest genuine relative is exactly what
should be chosen.
"""
import argparse, collections, csv, gzip, os, sys
import numpy as np


def open_maybe_gz(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="isolate calls in H37Rv coordinates")
    ap.add_argument("--panel-snps", default="graphs/CX333.s10k.k23.K15/snps.vcf.gz")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # panel matrix: position -> (ref, alt, per-genome genotype)
    hdr = None
    for line in open_maybe_gz(a.panel_snps):
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t"); break
    samples = hdr[9:]
    pos_idx, rows = {}, []
    for line in open_maybe_gz(a.panel_snps):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        alts = f[4].upper().split(",")
        if len(alts) != 1 or len(f[3]) != 1 or len(alts[0]) != 1:
            continue
        pos_idx[int(f[1])] = len(rows)
        rows.append([int(g) if g.isdigit() else -1 for g in f[9:]])
    G = np.array(rows, dtype=np.int8)

    # isolate profile: 1 where it calls ALT at a panel site, 0 where it does not
    called = set()
    for line in open_maybe_gz(a.vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[6] not in (".", "PASS"):
            continue
        if len(f[3]) == 1 and len(f[4]) == 1 and f[4] not in (".", "*"):
            called.add(int(f[1]))
    q = np.zeros(G.shape[0], dtype=np.int8)
    for p, i in pos_idx.items():
        if p in called:
            q[i] = 1
    ok = G >= 0
    d = ((G != q[:, None]) & ok).sum(axis=0)

    order = np.argsort(d)
    with open(a.out, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rank", "reference", "snp_distance"])
        for r, i in enumerate(order[: a.top], 1):
            w.writerow([r, samples[i], int(d[i])])
    best = samples[order[0]]
    print(f"  {os.path.basename(a.vcf)}: {len(called)} SNP calls, "
          f"{G.shape[0]} panel sites; nearest {best} at {int(d[order[0]])} "
          f"(next {samples[order[1]]} at {int(d[order[1]])})")
    print(best)
    return 0


if __name__ == "__main__":
    sys.exit(main())
