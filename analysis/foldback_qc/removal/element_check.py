#!/usr/bin/env python3
"""Would removing fold-back reads remove evidence for the tatC-helY element
(68 bp hairpin, H37Rv 2,352,065)? Report only (user, 2026-10-12). Run from
runroot with $MTB_PY.

The element is a palindrome, so a read through it could in principle look
like a fold-back (forward, then back along the reverse complement nearby).
Two checks, on existing outputs only:

1. Junctions at the element. The full-read fold-back scans (foldback_v2.py
   --sites) list every fold-back junction with its read count. Per sample,
   the element is located in the scanned reference by its two 30-mer flanks
   (element_genotype.locate), and the fold-back reads with a junction within
   the element +-50 bp are counted, against chance: the sample's fold-back
   reads genome-wide times the window's share of the genome.
   Samples: selfcal54 (each on its own genome) and the 9 chimeric controls
   (fbremove, all reads, own genome).
2. Genotype before and after removal. element_genotype's measure (depth
   inside the element over depth of the 200 bp flanks, MAPQ >= 20 primary
   reads; absent below 0.05) on the 9 chimeric controls, all reads against
   every fold-back pair removed:
     fbremove   own genome (_all vs _fb100)
     fbrecall   nearest other panel genome (fbrecall vs fbrecall100)
   Plus the reads lost inside the element and in the flanks.
Writes out/element_check.tsv and prints a summary.
"""
import csv
import gzip
import os
import sys

import pysam

sys.path.insert(0, "../analysis/novel_events")
from element_genotype import ABSENT_MAX, ELEMENT, FLANK, depth, locate  # noqa: E402
from artifact_mask import read_fasta  # noqa: E402

D = "../analysis/foldback_qc"
REFS = "refbias/build/7713a8d71d8e-fix1/refs"
PAD = 50


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


SEQ = {}


def site(g):
    if g not in SEQ:
        seq = read_fasta(f"{REFS}/{g}.fasta")[1]
        SEQ[g] = (locate(seq), len(seq))
    return SEQ[g]


def junctions(path, lo, hi):
    tot = near = 0
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n = int(r["reads"])
            tot += n
            if lo - PAD <= int(r["pos"]) <= hi + PAD:
                near += n
    return tot, near


def geno(bam_path, lo, hi):
    bam = pysam.AlignmentFile(bam_path)
    ctg = bam.references[0]
    ins = depth(bam, ctg, lo + 10, hi - 10)
    fl = (depth(bam, ctg, lo - FLANK, lo) + depth(bam, ctg, hi, hi + FLANK)) / 2
    ratio = ins / fl if fl else float("nan")
    return ins, fl, ratio, "absent" if ratio < ABSENT_MAX else "present"


def main():
    rows = []
    # 1. junctions at the element
    pts = []
    for r in tsv("refbias/selfcal54/run/p1/refmap.tsv"):
        pts.append((r["sample"], r["reference"], f"{D}/selfcal/out/sites/{r['sample']}.tsv.gz", "selfcal54"))
    for r in tsv("refbias/fbremove/run/p1/refmap.tsv"):
        if r["sample"].endswith("_all"):
            pts.append((r["sample"], r["reference"], f"{D}/removal/out/sites/{r['sample']}.tsv.gz", "fbremove"))
    print("1. Fold-back reads with a junction at the element (+-50 bp)")
    print("sample\tcohort\treference\telement_in_ref\tfoldback_reads_genome\tat_element\texpected")
    for s, g, p, coh in pts:
        if not os.path.exists(p):
            continue
        iv, L = site(g)
        has = "unresolved" if not iv else "yes" if abs(iv[1] - iv[0] - ELEMENT) <= 10 else \
            "no" if iv[1] - iv[0] < 20 else "unresolved"
        if not iv:
            continue
        tot, near = junctions(p, *iv)
        exp = tot * (iv[1] - iv[0] + 2 * PAD) / L
        rows.append(dict(check="junctions", sample=s, cohort=coh, reference=g, element_in_ref=has,
                         foldback_reads=tot, at_element=near, expected=round(exp, 2)))
        print(f"{s}\t{coh}\t{g}\t{has}\t{tot}\t{near}\t{exp:.2f}")
    j = [r for r in rows if r["check"] == "junctions"]
    jy = [r for r in j if r["element_in_ref"] == "yes"]
    print(f"  samples {len(j)} ({len(jy)} with the element in the reference): "
          f"{sum(r['at_element'] for r in jy)} fold-back reads at the element, "
          f"{sum(r['expected'] for r in jy):.1f} expected by chance; "
          f"samples with any: {sum(r['at_element'] > 0 for r in jy)}")

    # 2. genotype before and after removal
    print("\n2. Element genotype, all reads vs fold-back pairs removed")
    print("sample\treference\tinside_all\tinside_fb100\tflank_all\tflank_fb100\tratio_all\tratio_fb100\tgt_all\tgt_fb100")
    arms = [("fbremove", r["sample"][:-4], r["reference"], f"refbias/work/fbremove_p2/{r['sample']}.bam",
             f"refbias/work/fbremove_p2/{r['sample'][:-4]}_fb100.bam")
            for r in tsv("refbias/fbremove/run/p1/refmap.tsv") if r["sample"].endswith("_all")]
    arms += [("fbrecall", r["sample"], r["reference"], f"refbias/work/fbrecall_p2/{r['sample']}.bam",
              f"refbias/work/fbrecall100_p2/{r['sample']}_fb100.bam")
             for r in tsv("refbias/fbrecall/run/p1/refmap.tsv") if r["sample"].startswith("mar_")]
    for arm, s, g, b0, b1 in arms:
        iv, _ = site(g)
        if not iv or abs(iv[1] - iv[0] - ELEMENT) > 10 or not (os.path.exists(b0) and os.path.exists(b1)):
            print(f"{arm}:{s}\t{g}\tno element in reference or BAM missing")
            continue
        a, b = geno(b0, *iv), geno(b1, *iv)
        rows.append(dict(check=f"genotype_{arm}", sample=s, reference=g, element_in_ref="yes",
                         inside_all=round(a[0], 1), inside_fb100=round(b[0], 1),
                         flank_all=round(a[1], 1), flank_fb100=round(b[1], 1),
                         ratio_all=round(a[2], 3), ratio_fb100=round(b[2], 3), gt_all=a[3], gt_fb100=b[3]))
        print(f"{arm}:{s}\t{g}\t{a[0]:.1f}\t{b[0]:.1f}\t{a[1]:.1f}\t{b[1]:.1f}\t{a[2]:.3f}\t{b[2]:.3f}\t{a[3]}\t{b[3]}")
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(f"{D}/removal/out/element_check.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
