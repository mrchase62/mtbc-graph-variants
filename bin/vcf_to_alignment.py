#!/usr/bin/env python3
"""Turn a biallelic SNP VCF into a FASTA alignment for tree building.

Each sample gets the REF base where GT=0, the ALT base where GT=1, and N where
the call is missing. Sites that end up constant after filtering are dropped:
an alignment of only variable sites needs an ascertainment-bias correction
(`+ASC`) and IQ-TREE refuses to run it if constant sites are present.
"""
import argparse
import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-missing", type=float, default=0.10,
                    help="drop a site if more than this fraction of samples are missing")
    ap.add_argument("--sites-out", default="",
                    help="TSV of the site behind each alignment column: "
                         "chrom, pos, ref, alt, in column order")
    ap.add_argument("--ref-sample", default="",
                    help="emit the VCF's reference genome as an all-REF "
                         "sequence under this name; it has no sample column "
                         "and would otherwise be absent from the tree")
    ap.add_argument("--min-minor", type=int, default=1,
                    help="drop a site whose minor allele count falls below this")
    args = ap.parse_args()

    v = pysam.VariantFile(args.vcf)
    samples = list(v.header.samples)
    cols = {s: [] for s in samples}
    sites = []
    if args.ref_sample:
        cols[args.ref_sample] = []
    kept = dropped_missing = dropped_const = 0

    ACGT = set("ACGT")
    dropped_ambig = 0
    for rec in v:
        if len(rec.alts) != 1 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
        # Reject ambiguous alleles. This VCF contains records whose ALT is
        # literally "N" (an ambiguous base in the assembly). They pass an
        # allele-count filter -- two distinct allele indices are present -- but
        # produce an alignment column with one real state plus N, which IQ-TREE
        # then counts as invariant and refuses to run +ASC on.
        if rec.ref.upper() not in ACGT or rec.alts[0].upper() not in ACGT:
            dropped_ambig += 1
            continue
        gts = {}
        miss = 0
        for s in samples:
            a = rec.samples[s].allele_indices
            g = a[0] if a and a[0] is not None else None
            if g is None:
                miss += 1
            gts[s] = g
        if miss / len(samples) > args.max_missing:
            dropped_missing += 1
            continue
        called = [g for g in gts.values() if g is not None]
        n1 = sum(called)
        n0 = len(called) - n1
        if min(n0, n1) < args.min_minor:
            dropped_const += 1
            continue
        ref, alt = rec.ref.upper(), rec.alts[0].upper()
        for s in samples:
            g = gts[s]
            cols[s].append("N" if g is None else (ref if g == 0 else alt))
        # The VCF's reference genome has no sample column, so a tree built
        # straight from this alignment omits it -- which for CX333 means H37Rv,
        # the panel's best-annotated genome and the frame every coordinate in
        # this project is reported in, missing from its own phylogeny. It is
        # the all-REF sequence by definition, so it is emitted as one when
        # --ref-sample is given.
        if args.ref_sample:
            cols[args.ref_sample].append(ref)
        # Record which site each alignment column came from. Without this the
        # alignment is unusable for anything that has to join back to a
        # coordinate -- notably the ancestral-allele reconstruction, whose
        # whole output is a per-position statement.
        sites.append((rec.chrom, rec.pos, ref, alt))
        kept += 1

    with open(args.out, "w") as fh:
        # iterate cols, not samples: the reference sequence added by
        # --ref-sample is a key in cols and not in samples, so writing from
        # `samples` silently dropped it and produced an alignment that looked
        # correct at 332 sequences.
        for s in cols:
            fh.write(f">{s}\n")
            seq = "".join(cols[s])
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

    print(f"[vcf_to_alignment] {len(samples)} samples")
    print(f"  sites kept                  : {kept}")
    print(f"  dropped, >{args.max_missing:.0%} missing        : {dropped_missing}")
    print(f"  dropped, constant/too rare  : {dropped_const}")
    print(f"  dropped, ambiguous allele   : {dropped_ambig}")
    if args.sites_out:
        with open(args.sites_out, "w") as fh:
            fh.write("column\tchrom\tpos\tref\talt\n")
            for i, (c, pos, r, a) in enumerate(sites):
                fh.write(f"{i}\t{c}\t{pos}\t{r}\t{a}\n")
        print(f"  -> {args.sites_out}  ({len(sites)} sites)")
    print(f"  -> {args.out}")


if __name__ == "__main__":
    main()
