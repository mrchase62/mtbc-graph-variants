#!/usr/bin/env python3
"""Is an insertion seen in short reads present in an isolate's Nanopore
reads? (report only; run from runroot with $MTB_PY after 02_align.sbatch)

  lr_site.py <isolate> <pos> [--min-ins 1000] [--motif SEQ]

Reads (primary or supplementary, MAPQ >= 20) covering pos with >= 200 bp of
reference on both sides are counted:
  inserted   an insertion >= --min-ins within 50 bp of pos
  through    no insertion or deletion >= 50 bp within 50 bp of pos
  other
Reads clipped >= --min-ins within 50 bp of pos (the rest aligned elsewhere)
are counted separately. With --motif (e.g. an IS6110 terminal inverted
repeat), each inserted or clipped sequence is tested for it on either strand
(exact, 26 bp).
"""
import argparse

import pysam

BAMS = "/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_ont/bam"
COMP = str.maketrans("ACGT", "TGCA")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("isolate")
    ap.add_argument("pos", type=int, help="1-based position in the matched reference")
    ap.add_argument("--min-ins", type=int, default=1000)
    ap.add_argument("--motif", default="")
    a = ap.parse_args()
    m = a.motif.upper()
    mrc = m.translate(COMP)[::-1]
    has = (lambda s: bool(m) and (m in s.upper() or mrc in s.upper()))
    b = pysam.AlignmentFile(f"{BAMS}/{a.isolate}.bam")
    ctg = b.references[0]
    p0 = a.pos - 1
    n = dict(inserted=0, inserted_motif=0, through=0, other=0, clipped=0, clipped_motif=0)
    sizes = []
    for r in b.fetch(ctg, p0 - 1, p0 + 1):
        if r.is_secondary or r.mapping_quality < 20:
            continue
        ct = r.cigartuples
        lc = ct[0][1] if ct[0][0] == 4 else 0
        rc = ct[-1][1] if ct[-1][0] == 4 else 0
        q = r.query_sequence or ""
        if rc >= a.min_ins and abs(r.reference_end - p0) <= 50:
            n["clipped"] += 1
            n["clipped_motif"] += has(q[-rc:])
            continue
        if lc >= a.min_ins and abs(r.reference_start - p0) <= 50:
            n["clipped"] += 1
            n["clipped_motif"] += has(q[:lc])
            continue
        if r.reference_start > p0 - 200 or r.reference_end < p0 + 200:
            continue
        p, qp, big, sv = r.reference_start, 0, None, False
        for op, k in ct:
            if op == 1 and abs(p - p0) <= 50:
                if k >= a.min_ins:
                    big = q[qp:qp + k]
                if k >= 50:
                    sv = True
            if op == 2 and k >= 50 and p - 50 <= p0 <= p + k + 50:
                sv = True
            if op in (0, 2, 3, 7, 8):
                p += k
            if op in (0, 1, 4, 7, 8):
                qp += k
        if big is not None:
            n["inserted"] += 1
            n["inserted_motif"] += has(big)
            sizes.append(len(big))
        elif not sv:
            n["through"] += 1
        else:
            n["other"] += 1
    span = n["inserted"] + n["through"]
    print(f"{a.isolate} {ctg}:{a.pos}  {n}  inserted sizes {sorted(sizes)}")
    if span:
        print(f"fraction with insertion among spanning reads: {n['inserted'] / span:.3f} ({n['inserted']}/{span})")


if __name__ == "__main__":
    main()
