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

A site is an allele, (pos, ref, alt) trimmed to its minimal form, not a
position (audit GRAPHVCF-4 / P0P2-5). The panel VCF is split biallelic, so a
position with C>A and C>G is two rows; keyed by position the second overwrote the
first, the first was scored as if no isolate ever carried it, and an isolate's
C>G matched the C>A carriers. A padded SNP such as `CG>TG` is trimmed to `C>T`
and kept; it used to be dropped. The isolate's profile is the alleles its GT
calls.

Two things are deliberately NOT done (audit P0P2-6 and P0P2-7, left as
decisions):
  - the distance is a raw count, not normalised by the number of sites each
    genome has genotyped, and an isolate site without a call counts as REF;
  - H37Rv is never a candidate. It is the deconstruct reference path, so it has
    no column in the panel VCF; the matched arm is a reference other than H37Rv.

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


def trim(pos, ref, alt):
    """Minimal allele: shared suffix, then shared prefix down to one base."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="isolate calls in H37Rv coordinates")
    ap.add_argument("--panel-snps", default="graphs/CX333.s10k.k23.K15/snps.vcf.gz")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # panel matrix: (pos, ref, alt) -> per-genome genotype
    hdr = None
    for line in open_maybe_gz(a.panel_snps):
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t"); break
    samples = hdr[9:]
    key_idx, rows = {}, []
    for line in open_maybe_gz(a.panel_snps):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        alts = f[4].upper().split(",")
        if len(alts) != 1:
            continue
        k = trim(int(f[1]), f[3].upper(), alts[0])
        if len(k[1]) != 1 or len(k[2]) != 1:
            continue
        g = [int(x) if x.isdigit() else -1 for x in f[9:]]
        if k in key_idx:
            # a repeated key (an uncollapsed panel VCF) is one allele: union
            # its genotypes, ALT over REF over missing, rather than overwrite
            rows[key_idx[k]] = [max(x, y) for x, y in zip(rows[key_idx[k]], g)]
            continue
        key_idx[k] = len(rows)
        rows.append(g)
    G = np.array(rows, dtype=np.int8)

    # isolate profile: 1 where it calls that ALLELE at a panel site, else 0
    called = set()
    for line in open_maybe_gz(a.vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[6] not in (".", "PASS") or len(f) < 10:
            continue
        alts = f[4].upper().split(",")
        gt = f[9].split(":")[0].replace("|", "/").split("/")
        for i in {int(x) for x in gt if x.isdigit() and int(x) > 0}:
            if i > len(alts) or alts[i - 1] in (".", "*"):
                continue
            k = trim(int(f[1]), f[3].upper(), alts[i - 1])
            if len(k[1]) == 1 and len(k[2]) == 1:
                called.add(k)
    q = np.zeros(G.shape[0], dtype=np.int8)
    for k, i in key_idx.items():
        if k in called:
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
