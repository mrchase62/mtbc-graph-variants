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

INPUT: THE COLLAPSED VCF (audit PGB-9). all_variants.collapsed.vcf.gz has one
record per trimmed (chrom, pos, ref, alt), genotypes unioned and multiallelic
records split, so every allele gets a row -- on CX333, 628 SNP alleles that sit
in multiallelic decomposed records had none.

ONE ROW PER (chrom, pos, ref, alt), NOT PER RECORD. The decomposed VCF repeats
a key once per allele path through a bubble, and a genome carrying the ALT has
GT 1 in only one of the copies. So if it is given, the copies are pooled: the outgroup is ALT
if it is 1 in any copy, REF if it is 0 in some copy and 1 in none; and
`panel_af` counts each genome once, as a carrier if it is 1 in any copy, over
the genomes with a call in any copy. Taking one record's genotype and its own
AC/AN, as this did, let the last copy decide the outgroup's allele (audit
GRAPHVCF-3) and understated the panel frequency of a split allele.

    $MTB_PY bin/panel_polarity.py --panel-vcf <build>/assets/graph_collapsed.vcf.gz \
        --outgroup <outgroup> --out <build>/assets/panel_polarity.tsv

P0 step panel_polarity runs it; the outgroup is config MTB_OUTGROUP.
"""
import argparse, collections, itertools, os, subprocess, sys


MISSING = {".", "./.", ".|."}


def main():
    ap = argparse.ArgumentParser()
    # no defaults: they were CX333's collapsed VCF and its canettii; P0 step
    # panel_polarity passes the build's VCF and its outgroup
    ap.add_argument("--panel-vcf", required=True,
                    help="<build>/assets/graph_collapsed.vcf.gz")
    ap.add_argument("--outgroup", required=True)
    ap.add_argument("--bcftools",
                    default=os.environ.get("MTB_BCFTOOLS", "bcftools"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    have = subprocess.run([a.bcftools, "query", "-l", a.panel_vcf],
                          capture_output=True, text=True).stdout.split()
    if a.outgroup not in have:
        sys.exit(f"FATAL: {a.outgroup} is not a sample of {a.panel_vcf}")
    og = have.index(a.outgroup)
    proc = subprocess.Popen(
        [a.bcftools, "query",
         "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", a.panel_vcf],
        stdout=subprocess.PIPE, text=True)

    def records():
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            yield f[0], f[1], f[2].upper(), f[3].upper(), f[4:]

    n = kept = anc_alt = n_dup = 0
    with open(a.out, "w") as fh:
        fh.write("chrom\tpos\tref\talt\toutgroup_gt\tancestral\tpanel_af"
                 "\tn_records\n")
        # duplicates share a position, and the VCF is sorted, so pooling one
        # position at a time sees every copy of a key
        for _, grp in itertools.groupby(records(), key=lambda r: r[:2]):
            by_key = collections.OrderedDict()
            for chrom, pos, ref, alt, gts in grp:
                n += 1
                if "," in alt:
                    continue          # multiallelic
                by_key.setdefault((chrom, pos, ref, alt), []).append(gts)
            for (chrom, pos, ref, alt), copies in by_key.items():
                n_dup += len(copies) > 1
                og_gts = {c[og] for c in copies}
                gt = "1" if "1" in og_gts else "0" if "0" in og_gts else None
                if gt is None:
                    continue          # the outgroup has no call in any copy
                # per genome: a carrier if 1 in any copy, called if any call
                col = ["1" if "1" in c else "0" if c - MISSING else "."
                       for c in map(set, zip(*copies))] \
                    if len(copies) > 1 else copies[0]
                carrier = sum(g == "1" for g in col)
                called = sum(g not in MISSING for g in col)
                af = carrier / called if called else float("nan")
                anc = "ALT" if gt == "1" else "REF"
                anc_alt += gt == "1"
                kept += 1
                fh.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gt}\t{anc}\t"
                         f"{af:.4f}\t{len(copies)}\n")
    if proc.wait() != 0:
        sys.exit(f"FATAL: bcftools query failed on {a.panel_vcf}")
    print(f"  {n:,} panel records read, {kept:,} usable keys "
          f"(biallelic with an outgroup call); {n_dup:,} keys pooled over "
          f"duplicate records")
    print(f"  the outgroup carries the ALT at {anc_alt:,} of them "
          f"({anc_alt / max(1, kept):.1%}), so ALT is the ancestral allele there")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
