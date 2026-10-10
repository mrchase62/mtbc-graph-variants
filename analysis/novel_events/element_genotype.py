#!/usr/bin/env python3
"""Genotype the tatC-helY element (68 bp hairpin, H37Rv 2,352,065) per
isolate from its short-read alignment to the matched reference R, so calls
at the site can be judged by genotype instead of masked (user, 2026-10-10).

Why: in isolates that carry the element, short reads show a low-fraction
precise 68 bp deletion (callers report DEL at AF 0.1-0.4). Nanopore reads of
12 setE carriers have the element in 839/843 spanning reads, so that signal
is a short-read library/sequencing artifact (longread/lr_element.py). A real
loss leaves no reads inside the element.

Steps, per isolate:
  1. Locate the site in R by the two 30-mer flanks of artifact_mask.py.
     R has the element when the flanks are 68 +- 10 bp apart, lacks it when
     they are adjacent (gap < 20 bp); otherwise unresolved.
  2. If R has it: mean depth inside the element (10 bp in from each end;
     MAPQ >= 20 primary reads) over mean depth of the 200 bp outside each
     flank.
  3. genotype  absent   ratio < ABSENT_MAX (0.05: real losses measure 0.00,
                        carriers 0.08-0.66 on 46 Phase B isolates)
               present  otherwise
               ref_lacks R lacks the element (an INS call there may be real)
               unresolved
Rule applied by score.py --element-gt: a DEL call overlapping the site
(+-50 bp) is dropped when the genotype is present. INS calls and DEL calls in
absent isolates are kept. The ratio is the record of the artifact's strength.

Writes one row per isolate: sample, reference, site_start, site_end (0-based
half-open, the sequence between the flanks), ref_has_element, inside_depth,
flank_depth, ratio, genotype.
"""
import argparse
import csv
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from artifact_mask import SITES, find, read_fasta  # noqa: E402

ABSENT_MAX = 0.05
FLANK = 200
ELEMENT = 68


def depth(b, ctg, s, e):
    if e <= s:
        return 0.0
    cov = b.count_coverage(ctg, s, e, quality_threshold=0,
                           read_callback=lambda r: not (r.is_secondary or r.is_supplementary
                                                        or r.mapping_quality < 20))
    return sum(sum(x[i] for x in cov) for i in range(e - s)) / (e - s)


def locate(seq):
    _, left, right = SITES[0]
    L, R = find(seq, left), find(seq, right)
    if len(L) != 1 or len(R) != 1 or abs(L[0] - R[0]) > 1000:
        return None
    lo, hi = min(L[0], R[0]) + 30, max(L[0], R[0])
    return lo, hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refmap", required=True)
    ap.add_argument("--bam-dir", required=True, help="<sample>.bam aligned to the matched reference")
    ap.add_argument("--refs", default="refbias/build/7713a8d71d8e-fix1/refs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    seqs = {}
    rows = []
    for r in csv.DictReader(open(a.refmap), delimiter="\t"):
        s, g = r["sample"], r["reference"]
        if g not in seqs:
            seqs[g] = read_fasta(f"{a.refs}/{g}.fasta")[1]
        row = dict(sample=s, reference=g, site_start="", site_end="", ref_has_element="", inside_depth="",
                   flank_depth="", ratio="", genotype="unresolved")
        iv = locate(seqs[g])
        bam_path = f"{a.bam_dir}/{s}.bam"
        if iv:
            lo, hi = iv
            row.update(site_start=lo, site_end=hi)
            if hi - lo < 20:
                row.update(ref_has_element=0, genotype="ref_lacks")
            elif abs(hi - lo - ELEMENT) <= 10 and os.path.exists(bam_path):
                bam = pysam.AlignmentFile(bam_path)
                ctg = bam.references[0]
                ins = depth(bam, ctg, lo + 10, hi - 10)
                fl = (depth(bam, ctg, lo - FLANK, lo) + depth(bam, ctg, hi, hi + FLANK)) / 2
                ratio = ins / fl if fl else float("nan")
                row.update(ref_has_element=1, inside_depth=round(ins, 1), flank_depth=round(fl, 1),
                           ratio=round(ratio, 3), genotype="absent" if ratio < ABSENT_MAX else "present")
        rows.append(row)
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n = {}
    for x in rows:
        n[x["genotype"]] = n.get(x["genotype"], 0) + 1
    print(n)


if __name__ == "__main__":
    main()
