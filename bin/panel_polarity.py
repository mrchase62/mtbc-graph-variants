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

KEYS ARE NORMALISED THE WAY THE COHORT VCF'S ARE (review 2, R2-INT-1). The
collapsed VCF trims alleles but does not left-align them (decision D22), while
every H37Rv-frame key in a merged cohort VCF goes through
mtb_norm.normalise, which does. The event writer looks this table up by exact
(chrom, pos, ref, alt), so an indel in a repeat was missed whenever the two
placed it differently: on scale200, 2,951 of 8,711 H37Rv-frame indels matched
only after left-alignment, 103 of them with the outgroup carrying the ALT.
Each key on the H37Rv path is therefore passed through the same function with
the same H37Rv sequence, and keys that normalise to the same allele are pooled
by the rule above.

    $MTB_PY bin/panel_polarity.py --panel-vcf <build>/assets/graph_collapsed.vcf.gz \
        --h37rv <build>/refs/GCF_000195955.fasta \
        --outgroup <outgroup> --out <build>/assets/panel_polarity.tsv

P0 step panel_polarity runs it; the outgroup is config MTB_OUTGROUP.
"""
import argparse, collections, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mtb_norm


MISSING = {".", "./.", ".|."}


def main():
    ap = argparse.ArgumentParser()
    # no defaults: they were CX333's collapsed VCF and its canettii; P0 step
    # panel_polarity passes the build's VCF and its outgroup
    ap.add_argument("--panel-vcf", required=True,
                    help="<build>/assets/graph_collapsed.vcf.gz")
    ap.add_argument("--outgroup", required=True)
    ap.add_argument("--h37rv", required=True,
                    help="H37Rv FASTA (<build>/refs/<acc>.fasta); keys are "
                         "left-aligned against it as the cohort VCF's are")
    ap.add_argument("--h37rv-contig", default="NC_000962.3",
                    help="H37Rv contig name; panel contigs are matched on the "
                         "part after the last '#' (PanSN)")
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

    seq = "".join(l.strip() for l in open(a.h37rv) if not l.startswith(">"))
    if not seq:
        sys.exit(f"FATAL: no sequence in {a.h37rv}")

    # pool every copy of a normalised key. Normalising moves an indel left, so
    # copies of one allele need not share a position any more; the table is
    # small (about 90,000 keys), so a dict over the whole VCF is fine
    pooled = collections.OrderedDict()
    n = n_multi = n_moved = n_offpath = 0
    for chrom, pos, ref, alt, gts in records():
        n += 1
        if "," in alt:
            n_multi += 1
            continue
        pos = int(pos)
        if chrom.split("#")[-1] == a.h37rv_contig:
            npos, nref, nalt = mtb_norm.normalise(seq, pos, ref, alt)
            n_moved += (npos, nref, nalt) != (pos, ref, alt)
            pos, ref, alt = npos, nref, nalt
        else:
            n_offpath += 1
        pooled.setdefault((chrom, pos, ref, alt), []).append(gts)
    if proc.wait() != 0:
        sys.exit(f"FATAL: bcftools query failed on {a.panel_vcf}")

    kept = anc_alt = n_dup = 0
    with open(a.out, "w") as fh:
        fh.write("chrom\tpos\tref\talt\toutgroup_gt\tancestral\tpanel_af"
                 "\tn_records\n")
        for (chrom, pos, ref, alt), copies in sorted(pooled.items(),
                key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3])):
            n_dup += len(copies) > 1
            og_gts = {c[og] for c in copies}
            gt = "1" if "1" in og_gts else "0" if "0" in og_gts else None
            if gt is None:
                continue              # the outgroup has no call in any copy
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
    print(f"  {n:,} panel records read ({n_multi:,} multiallelic skipped); "
          f"{n_moved:,} keys changed by left-alignment"
          + (f"; {n_offpath:,} off the H37Rv path, not normalised" if n_offpath else ""))
    print(f"  {kept:,} usable keys (biallelic with an outgroup call); "
          f"{n_dup:,} keys pooled over several records")
    print(f"  the outgroup carries the ALT at {anc_alt:,} of them "
          f"({anc_alt / max(1, kept):.1%}), so ALT is the ancestral allele there")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
