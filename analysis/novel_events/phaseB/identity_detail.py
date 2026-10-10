#!/usr/bin/env python3
"""Why do an isolate's reads and assembly disagree? (report only; run from
runroot with $MTB_PY after ours_false.py)

For each isolate given (default: the three identity-check failures plus
TB3334 as a clean comparison), the SNP sets of ours_false.py (assembly vs
reads, callable sequence) are split into shared / assembly_only /
reads_only. At every site the clean-read BAM is piled up (MAPQ >= 20,
base quality >= 20) and the alt fraction = reads with the SNP base / all
reads is recorded.
  mixed sample          intermediate alt fractions (0.2-0.8) at many sites
  assembly error        assembly_only sites with reads ~0 alt, often clustered
  reads drifted         reads_only sites with alt ~1 (extra fixed mutations)
Also reports, per class, how many sites fall in PE/PPE-like clusters
(another discordant site within 1 kb). Writes out/qc51/identity_detail.tsv.
"""
import collections
import csv
import sys

import pysam

sys.path.insert(0, "../analysis/novel_events/phaseB")
sys.path.insert(0, "../analysis/novel_events")
import ours_false as of  # noqa: E402
import score  # noqa: E402

DEFAULT = ["mar_QC_7", "mar_TB1612", "mar_TB3054", "mar_TB3334"]


def alt_frac(bam, ctg, pos, alt):
    n = collections.Counter()
    for col in bam.pileup(ctg, pos - 1, pos, truncate=True, min_mapping_quality=20, min_base_quality=20,
                          stepper="nofilter", ignore_orphans=False):
        for pr in col.pileups:
            if not pr.is_del and not pr.is_refskip:
                n[pr.alignment.query_sequence[pr.query_position]] += 1
    d = sum(n.values())
    return d, (n[alt] / d if d else float("nan"))


def main():
    samples = sys.argv[1:] or DEFAULT
    ref = {r["sample"]: r["reference"] for r in csv.DictReader(open(of.RM), delimiter="\t")}
    masks = score.load_masks(of.RM, "../analysis/novel_events/out/masks")
    rows = []
    for s in samples:
        a = {x for x in of.asm_snps(s) if not score.masked(masks[s], x[0], x[0])}
        r = {x for x in of.read_snps(s) if not score.masked(masks[s], x[0], x[0])}
        bam = pysam.AlignmentFile(f"{of.BAM}/{s}.bam")
        ctg = bam.references[0]
        disc = sorted(p for p, _ in a ^ r)
        for cls, sites in (("shared", a & r), ("assembly_only", a - r), ("reads_only", r - a)):
            for p, alt in sorted(sites):
                d, f = alt_frac(bam, ctg, p, alt)
                near = sum(1 for q in disc if q != p and abs(q - p) <= 1000)
                rows.append(dict(sample=s, reference=ref[s], cls=cls, pos=p, alt=alt, depth=d,
                                 alt_frac=round(f, 3), discordant_within_1kb=near))
    with open("../analysis/novel_events/phaseB/out/qc51/identity_detail.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("sample  class  n | alt fraction: <0.2  0.2-0.8  >0.8 | median depth | clustered(>=1 other within 1kb)")
    for s in samples:
        for cls in ("shared", "assembly_only", "reads_only"):
            v = [x for x in rows if x["sample"] == s and x["cls"] == cls]
            if not v:
                continue
            f = [x["alt_frac"] for x in v]
            dp = sorted(x["depth"] for x in v)[len(v) // 2]
            print(f"{s:15s} {cls:14s} {len(v):4d} | {sum(x < 0.2 for x in f):4d} {sum(0.2 <= x <= 0.8 for x in f):6d} "
                  f"{sum(x > 0.8 for x in f):5d} | {dp:4d} | {sum(x['discordant_within_1kb'] > 0 for x in v)}")


if __name__ == "__main__":
    main()
