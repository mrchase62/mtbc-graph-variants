#!/usr/bin/env python3
"""Ancestral allele from the outgroup's own genotype, for every variant class.

WHY THIS EXISTS. The project's `AA` annotation is Fitch parsimony over the
333-genome tree, computed from a SNP ALIGNMENT -- so it covers SNPs and nothing
else. Every indel in the merged VCF therefore comes out `unpolarised`, and the
event writer falls back to assuming ALT is derived. Where the reference is the
odd genome out that assumption is exactly backwards, and the gains and losses
counted on the variant are inverted with it.

The case that exposed it: a 68 bp insertion at 2,352,065, in the promoter
region of helY. M. canettii carries it and so do 229 of the 332 complete panel
genomes; H37Rv does not. The reconstruction read 11 gains and 31 losses of an
insertion, when the truth is 31 independent deletions and 11 reversions.

The graph's deconstruct VCF genotypes the outgroup at every variant it
contains, indels included, so the ancestral allele can simply be read off.

PRECEDENCE, AND WHY THIS DOES NOT OVERRIDE `AA`. Fitch over 333 genomes and a
tree is a stronger inference than one genome's allele: it uses the whole panel
and the topology, and it can resolve a site where the outgroup itself carries
a derived state by homoplasy. So this table is consumed only where `AA` did not
resolve. `panel_af` is carried alongside, because an outgroup allele that is
also common in the panel is a much safer ancestral call than one that is rare
in it, and the reader should be able to see which they have.

    $MTB_PY bin/panel_polarity.py --out refbias/assets/panel_polarity.tsv
"""
import argparse, os, subprocess, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-vcf",
                    default="graphs/CX333.s10k.k23.K15/all_variants.decomposed.vcf.gz")
    ap.add_argument("--outgroup", default="GCF_035581225")
    ap.add_argument("--bcftools",
                    default=os.environ.get("MTB_BCFTOOLS", "bcftools"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    have = subprocess.run([a.bcftools, "query", "-l", a.panel_vcf],
                          capture_output=True, text=True).stdout.split()
    if a.outgroup not in have:
        sys.exit(f"FATAL: {a.outgroup} is not a sample of {a.panel_vcf}")
    q = subprocess.run(
        [a.bcftools, "query", "-s", a.outgroup,
         "-f", "%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN[\t%GT]\n", a.panel_vcf],
        capture_output=True, text=True, check=True).stdout

    n = kept = anc_alt = 0
    with open(a.out, "w") as fh:
        fh.write("pos\tref\talt\toutgroup_gt\tancestral\tpanel_af\n")
        for line in q.splitlines():
            f = line.split("\t")
            if len(f) < 6:
                continue
            n += 1
            pos, ref, alt, ac, an, gt = f[0], f[1].upper(), f[2].upper(), \
                f[3], f[4], f[5]
            if "," in alt or gt not in ("0", "1"):
                continue          # multiallelic or the outgroup has no call
            try:
                af = int(ac) / int(an)
            except (ValueError, ZeroDivisionError):
                af = float("nan")
            anc = "ALT" if gt == "1" else "REF"
            anc_alt += gt == "1"
            kept += 1
            fh.write(f"{pos}\t{ref}\t{alt}\t{gt}\t{anc}\t{af:.4f}\n")
    print(f"  {n:,} panel records read, {kept:,} usable "
          f"(biallelic with an outgroup call)")
    print(f"  the outgroup carries the ALT at {anc_alt:,} of them "
          f"({anc_alt / max(1, kept):.1%}), so ALT is the ancestral allele there")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
