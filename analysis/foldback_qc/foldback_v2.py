#!/usr/bin/env python3
"""Fold-back chimera scan, version 2 (approved 2026-10-11; plan in
analysis/qc_harmonize/PLAN.md, part A). One sample per call.

  foldback_v2.py --bam <bam|cram> --ref <fasta> --sample <s> --out <tsv>
                 [--sites <tsv.gz>] [--frac 1.0]

A fold-back molecule reads forward and then back along the reverse
complement of nearby sequence. Aligned, the read either continues as a second
piece on the opposite strand close by (split read), or, when that piece is
too short to align, ends in a soft clip whose sequence is the reverse
complement of nearby reference.

Reads counted (the denominator N): primary, mapped, MAPQ >= 20, not marked
duplicate. With --frac < 1, reads are kept by a hash of their name (both
mates together), not by position, so the subsample is spread over the genome.

Per read, at most one class:
  split  any SA entry with MAPQ >= 20 on the same contig, the opposite strand,
         starting within 1 kb of the read's start
  clip   otherwise, a soft clip >= 15 bp after quality trimming (outward from
         the alignment end, stopping at the first 5-base window with mean
         quality below 20, as breakpoint_caller.py) whose 15 bases next to
         the alignment end, reverse-complemented, occur exactly in the
         reference within 500 bp of that end
  clip_forward / clip_other (reported, not fold-back): the clip's 15 bases
         occur forward within 500 bp (tandem-like), or nowhere nearby
The junction position is the alignment end next to the clip (for split
reads, the end on the side of the larger clip).

The old count (controls/split_reads.awk: first SA entry, same contig,
opposite strand, within 1 kb, any MAPQ, all primary alignments) is computed
in the same pass, for comparison.

--sites writes junction positions with their read counts (split, clip), for
the cross-sample recurrence step (sites that fold back in many samples are a
property of the reference, not of one library).
"""
import argparse
import collections
import gzip
import statistics
import zlib

import pysam

MINQ = 20
MINCLIP = 15
WIN = 500
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def trimmed_len(quals):
    """Clip length kept by the sliding-window trim, quals read outward."""
    if quals is None:
        return 10 ** 9
    for i in range(len(quals) - 4):
        if sum(quals[i:i + 5]) < 5 * MINQ:
            return i
    return len(quals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sites")
    ap.add_argument("--frac", type=float, default=1.0)
    a = ap.parse_args()
    bam = pysam.AlignmentFile(a.bam, reference_filename=a.ref)
    fa = pysam.FastaFile(a.ref)
    seqs = {}
    keep_below = int(a.frac * 1000)
    n = collections.Counter()
    sites = collections.Counter()
    rlen, isize = [], []
    for r in bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary:
            continue
        if a.frac < 1 and zlib.crc32(r.query_name.encode()) % 1000 >= keep_below:
            continue
        if r.is_unmapped:
            n["old_denominator"] += 1
            continue
        n["old_denominator"] += 1
        ctg, st = r.reference_name, "-" if r.is_reverse else "+"
        sa = []
        if r.has_tag("SA"):
            for x in r.get_tag("SA").rstrip(";").split(";"):
                c, p, s, cg, mq, _ = x.split(",")
                sa.append((c, int(p), s, int(mq)))
        if sa and sa[0][0] == ctg and sa[0][2] != st and abs(sa[0][1] - (r.reference_start + 1)) < 1000:
            n["old_foldback"] += 1
        if r.is_duplicate:
            n["dup"] += 1
            continue
        if r.mapping_quality < MINQ:
            continue
        n["reads"] += 1
        if len(rlen) < 20000:
            rlen.append(r.query_length)
        if r.is_proper_pair and r.template_length > 0 and len(isize) < 20000:
            isize.append(r.template_length)
        ct = r.cigartuples
        lc = ct[0][1] if ct[0][0] == 4 else 0
        rcl = ct[-1][1] if ct[-1][0] == 4 else 0
        if any(c == ctg and s != st and q >= MINQ and abs(p - (r.reference_start + 1)) < 1000 for c, p, s, q in sa):
            n["split"] += 1
            j = r.reference_end if rcl >= lc else r.reference_start
            sites[(ctg, j, "split")] += 1
            continue
        if max(lc, rcl) < MINCLIP:
            continue
        seq, q = r.query_sequence, r.query_qualities
        right = rcl >= lc
        if right:
            k = trimmed_len(q[-rcl:] if q is not None else None)
            near = seq[-rcl:][:min(k, rcl)][:MINCLIP]  # outward from the end
            j = r.reference_end
        else:
            k = trimmed_len(q[:lc][::-1] if q is not None else None)
            near = seq[:lc][::-1][:min(k, lc)][:MINCLIP][::-1]  # the bases next to the start, read order
            j = r.reference_start
        if len(near) < MINCLIP:
            continue
        if ctg not in seqs:
            seqs[ctg] = fa.fetch(ctg).upper()
        win = seqs[ctg][max(0, j - WIN):j + WIN]
        if rc(near) in win:
            n["clip"] += 1
            sites[(ctg, j, "clip")] += 1
        elif near in win:
            n["clip_forward"] += 1
        else:
            n["clip_other"] += 1
    N = n["reads"] or 1
    row = dict(sample=a.sample, reads=n["reads"],
               foldback_per_1k=round(1000 * (n["split"] + n["clip"]) / N, 4),
               split_per_1k=round(1000 * n["split"] / N, 4),
               clip_per_1k=round(1000 * n["clip"] / N, 4),
               clip_forward_per_1k=round(1000 * n["clip_forward"] / N, 4),
               clip_other_per_1k=round(1000 * n["clip_other"] / N, 4),
               old_per_1k=round(1000 * n["old_foldback"] / max(1, n["old_denominator"]), 4),
               dup_frac=round(n["dup"] / max(1, n["old_denominator"]), 4),
               read_len_median=statistics.median(rlen) if rlen else "",
               insert_median=statistics.median(isize) if isize else "",
               frac=a.frac)
    with open(a.out, "w") as fo:
        fo.write("\t".join(row) + "\n" + "\t".join(str(v) for v in row.values()) + "\n")
    if a.sites:
        with gzip.open(a.sites, "wt") as fo:
            fo.write("contig\tpos\tkind\treads\n")
            for (c, p, kind), v in sorted(sites.items()):
                fo.write(f"{c}\t{p}\t{kind}\t{v}\n")


if __name__ == "__main__":
    main()
