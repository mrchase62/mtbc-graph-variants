#!/usr/bin/env python3
"""Is an inversion called from short reads between two repeat copies seen in
an isolate's Nanopore reads? (report only; run from runroot with $MTB_PY
after 02_align.sbatch)

  lr_inversion.py <isolate> <A> <B>

A and B are the called inversion ends (1-based, matched reference R). At
each end E, reads (primary or supplementary, MAPQ >= 20 for the counted
alignment) are classified:
  reference  one alignment (no deletion >= 500 bp) runs continuously from
             2 kb before E to 2 kb after it, i.e. the junction at E is in R's
             order
  inverted   the read has one alignment within 3 kb of A and another within
             3 kb of B, on opposite strands (the junction an inversion makes)
  other      reads touching E that are neither
Also listed per end: insertions and deletions of >= 300 bp inside alignments
within 6 kb of E (rounded to 100 bp), with the number of reads carrying them.
An "inverted" read can also be a read whose inserted repeat copy was aligned
to a copy of the same repeat near the other end; the indel list tells these
apart (2026-10-11, the 383 kb INV at 1,980,173-2,363,063 on GCF_001275565).
A read crossing E continuously with 2 kb on each side has unique sequence
on at least one side (an IS6110 copy is 1.36 kb), so it cannot be a read
from another copy.
"""
import argparse
import collections

import pysam

BAMS = "/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_ont/bam"
NEAR = 3000


def alns(bam, ctg, lo, hi):
    out = collections.defaultdict(list)
    for r in bam.fetch(ctg, lo, hi):
        if r.is_secondary or r.is_unmapped:
            continue
        out[r.query_name].append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("isolate")
    ap.add_argument("A", type=int)
    ap.add_argument("B", type=int)
    a = ap.parse_args()
    bam = pysam.AlignmentFile(f"{BAMS}/{a.isolate}.bam")
    ctg = bam.references[0]
    A, B = a.A - 1, a.B - 1
    at_a = alns(bam, ctg, A - NEAR, A + NEAR)
    at_b = alns(bam, ctg, B - NEAR, B + NEAR)
    inv = set()
    for q in set(at_a) & set(at_b):
        sa = {r.is_reverse for r in at_a[q] if r.mapping_quality >= 20}
        sb = {r.is_reverse for r in at_b[q] if r.mapping_quality >= 20}
        if any(x != y for x in sa for y in sb):
            inv.add(q)
    res = {}
    for end, (pos, reads) in {"A": (A, at_a), "B": (B, at_b)}.items():
        lo, hi = pos - 2000, pos + 2000
        n = collections.Counter()
        for q, rs in reads.items():
            if not any(r.reference_start <= pos <= r.reference_end for r in rs):
                continue
            if q in inv:
                n["inverted"] += 1
            elif any(r.mapping_quality >= 20 and r.reference_start <= lo and r.reference_end >= hi
                     and not any(op == 2 and k >= 500 for op, k in r.cigartuples) for r in rs):
                n["reference"] += 1
            else:
                n["other"] += 1
        big = collections.Counter()
        for rs in reads.values():
            for r in rs:
                if r.mapping_quality < 20 or not r.reference_start <= pos <= r.reference_end:
                    continue
                p = r.reference_start
                for op, k in r.cigartuples:
                    if op in (1, 2) and k >= 300 and abs(p - pos) <= 6000:
                        big[("INS" if op == 1 else "DEL", round(p + 1, -2), round(k, -2))] += 1
                    if op in (0, 2, 3, 7, 8):
                        p += k
        n["indels"] = "; ".join(f"{t} {q}+{k}bp x{c}" for (t, q, k), c in big.most_common(3) if c >= 3)
        res[end] = n
    print("\t".join([a.isolate] + [f"{e}:ref={res[e]['reference']},inv={res[e]['inverted']},other={res[e]['other']} [{res[e]['indels']}]"
                                    for e in ("A", "B")]))


if __name__ == "__main__":
    main()
