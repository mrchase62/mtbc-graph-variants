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

ONLY SITES BOTH SIDES CALLED ARE COMPARED, and the ranking is per site (the
user's decisions D20 and D24, 2026-10-07). Before, the distance was a raw
count over the sites each GENOME had genotyped, and an isolate site with no
call counted as REF. So a genome with more missing panel cells had fewer
chances to mismatch (72,932-75,249 comparable sites per genome, against
typical distances of 100-450), and a site the isolate never covered matched
every REF genome. An isolate site is called where it has at least --min-depth
reads (the isolate's own H37Rv BAM, or a precomputed `samtools depth` table)
and no filtered, indel or complex record overlaps it. `snp_distance` is the
mismatch count over those sites, `n_compared` their number, and the ranking is
by `distance_per_site`. A genome compared over fewer than --min-compared-frac
of the best-covered genome's sites is not a candidate, so a sparse genome
cannot win on a small denominator.

H37Rv is never a candidate here (audit P0P2-7). It is the deconstruct
reference path, so it has no column in the panel VCF.

No self-exclusion and no near-clone exclusion here, unlike stage 1. Those existed
because a panel genome is its own nearest neighbour and the panel contains
near-clones, which would have made the leave-one-out measurement meaningless. A
real isolate is not in the panel, so its nearest genuine relative is exactly what
should be chosen.
"""
import argparse, collections, csv, gzip, os, subprocess, sys, tempfile
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


def read_depth(a, positions):
    """{pos: depth} at `positions` (1-based, H37Rv), from --depth or by
    running samtools depth on --bam over a BED of exactly those sites."""
    want = set(positions)
    if a.depth:
        lines = open_maybe_gz(a.depth)
    else:
        hdr = subprocess.run([a.samtools, "view", "-H", a.bam], check=True,
                             capture_output=True, text=True).stdout
        sq = [l.split("SN:", 1)[1].split("\t")[0] for l in hdr.splitlines()
              if l.startswith("@SQ")]
        if len(sq) != 1:
            sys.exit(f"FATAL: {a.bam} has {len(sq)} contigs; expected H37Rv's one")
        tmp = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
        for p in positions:
            tmp.write(f"{sq[0]}\t{p - 1}\t{p}\n")
        tmp.close()
        try:
            lines = subprocess.run(
                [a.samtools, "depth", "-a", "-Q", str(a.min_mapq), "-b",
                 tmp.name, a.bam], check=True, capture_output=True,
                text=True).stdout.splitlines()
        finally:
            os.unlink(tmp.name)
    out = {}
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[1].isdigit() and int(f[1]) in want:
            out[int(f[1])] = int(f[2])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="isolate calls in H37Rv coordinates")
    # no default: it was CX333's; p1_select_reference.sh passes the build's
    ap.add_argument("--panel-snps", required=True,
                    help="<build>/assets/panel_snps.vcf.gz")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--out", required=True)
    ap.add_argument("--bam", default="",
                    help="the isolate's H37Rv BAM: its depth at panel sites "
                         "decides which sites it called (D24)")
    ap.add_argument("--samtools", default=os.environ.get("MTB_SAMTOOLS", "samtools"))
    ap.add_argument("--depth", default="",
                    help="instead of --bam: `samtools depth -a` output "
                         "(chrom, pos, depth) at the panel sites")
    ap.add_argument("--min-depth", type=int, default=5,
                    help="reads needed for an isolate site to count as called "
                         "(P5's --min-dp)")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-compared-frac", type=float, default=0.5)
    a = ap.parse_args()
    if bool(a.bam) == bool(a.depth):
        sys.exit("FATAL: pass exactly one of --bam or --depth. Without the "
                 "isolate's coverage an uncovered site would count as REF "
                 "(D24)")

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

    # isolate profile: 1 where it calls that ALLELE at a panel site, else 0.
    # Positions under a filtered record, or under an indel or complex
    # record's REF, are not REF evidence: not called.
    called, masked = set(), set()
    for line in open_maybe_gz(a.vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        pos, ref = int(f[1]), f[3].upper()
        alts = f[4].upper().split(",")
        gt = f[9].split(":")[0].replace("|", "/").split("/")
        carried = [alts[int(x) - 1] for x in gt if x.isdigit()
                   and 0 < int(x) <= len(alts) and alts[int(x) - 1] not in (".", "*")]
        if f[6] not in (".", "PASS"):
            masked.update(range(pos, pos + len(ref)))
            continue
        for alt in carried:
            k = trim(pos, ref, alt)
            if len(k[1]) == 1 and len(k[2]) == 1:
                called.add(k)
            else:
                masked.update(range(pos, pos + len(ref)))
    q = np.zeros(G.shape[0], dtype=np.int8)
    for k, i in key_idx.items():
        if k in called:
            q[i] = 1

    # which panel sites the isolate covers
    want = sorted({k[0] for k in key_idx})
    depth = read_depth(a, want)
    cov = np.array([depth.get(k[0], 0) >= a.min_depth and k[0] not in masked
                    for k in sorted(key_idx, key=key_idx.get)], dtype=bool)
    if not cov.any():
        sys.exit(f"FATAL: the isolate covers none of {G.shape[0]:,} panel "
                 f"sites at {a.min_depth}x; wrong BAM or contig?")
    ok = (G >= 0) & cov[:, None]
    d = ((G != q[:, None]) & ok).sum(axis=0)
    n = ok.sum(axis=0)
    rate = d / np.maximum(n, 1)
    eligible = n >= a.min_compared_frac * n.max()
    rate = np.where(eligible, rate, np.inf)

    # by rate; a tie by count, then by name, so the order is reproducible
    order = sorted(range(len(samples)), key=lambda i: (rate[i], d[i], samples[i]))
    with open(a.out, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rank", "reference", "snp_distance", "n_compared",
                    "distance_per_site"])
        for r, i in enumerate(order[: a.top], 1):
            if not eligible[i]:
                break
            w.writerow([r, samples[i], int(d[i]), int(n[i]), f"{rate[i]:.6g}"])
    best, nxt = order[0], order[1] if len(order) > 1 else order[0]
    print(f"  {os.path.basename(a.vcf)}: {len(called)} SNP calls; "
          f"{int(cov.sum()):,} of {G.shape[0]:,} panel sites covered at "
          f">= {a.min_depth}x; {int((~eligible).sum())} genomes below "
          f"{a.min_compared_frac:.0%} of the best comparison; nearest "
          f"{samples[best]} at {int(d[best])}/{int(n[best])} "
          f"(next {samples[nxt]} at {int(d[nxt])}/{int(n[nxt])})")
    print(samples[best])
    return 0


if __name__ == "__main__":
    sys.exit(main())
