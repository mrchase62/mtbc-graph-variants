#!/usr/bin/env python3
"""Append an outgroup sequence to a cohort SNP alignment, from the panel VCF.

WHY. A cohort of isolates carries no outgroup, so a tree built from it alone
has no information about where its root sits. Every branch below a node whose
state parsimony cannot resolve is written undetermined by
assoc/bin/write_event_matrix.py, so an unrooted-in-practice cohort tree throws
away exactly the deep events an association test would most like to count. The
panel VCF already holds `M. canettii` (GCF_035581225) genotyped against the
same H37Rv coordinates, so the outgroup can be read off rather than sequenced.

The outgroup's only job is to place the root. Prune it before the writer runs
(`--prune-tips` there), or its tip will be asked for a genotype it does not
have in the cohort VCF.

THE ALLELE RULE, and where it says N rather than guessing:

  exact (pos, ref, alt) record in the panel VCF   -> its genotype, 0 -> REF,
                                                     1 -> ALT, missing -> N
  a record at that position, but a different ALT  -> REF if the outgroup is
                                                     reference there, else N,
                                                     since a third allele is
                                                     not representable in a
                                                     biallelic column
  no record at that position at all               -> REF. The panel VCF holds
                                                     every site that varies
                                                     across the panel, and the
                                                     outgroup is in the panel,
                                                     so silence means it
                                                     matches the reference
  site not in the panel frame (node_* contigs)    -> N

    $MTB_PY assoc/bin/add_outgroup.py \
        --alignment data/trees/<cohort>.snps.fasta \
        --sites data/trees/<cohort>.sites.tsv \
        --out data/trees/<cohort>.outgroup.fasta
"""
import argparse, collections, csv, os, subprocess, sys, tempfile

PANEL_VCF = "graphs/CX333.s10k.k23.K15/all_variants.decomposed.vcf.gz"
PANEL_CONTIG = "GCF_000195955#1#NC_000962.3"


def read_fasta(path):
    seqs, name, buf = collections.OrderedDict(), None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--panel-vcf", default=PANEL_VCF)
    ap.add_argument("--panel-contig", default=PANEL_CONTIG)
    ap.add_argument("--outgroup", default="GCF_035581225",
                    help="the panel sample to emit (default M. canettii)")
    ap.add_argument("--cohort-contig", default="NC_000962.3",
                    help="the cohort VCF's name for the reference contig; "
                         "sites on any other contig get N")
    a = ap.parse_args()

    seqs = read_fasta(a.alignment)
    sites = list(csv.DictReader(open(a.sites), delimiter="\t"))
    width = len(next(iter(seqs.values())))
    if width != len(sites):
        sys.exit(f"FATAL: alignment is {width} columns but the sites table has "
                 f"{len(sites)} rows; they must be the same order and length")
    if a.outgroup in seqs:
        sys.exit(f"FATAL: {a.outgroup} is already in the alignment")

    bcftools = os.environ.get("MTB_BCFTOOLS", "bcftools")
    have = subprocess.run([bcftools, "query", "-l", a.panel_vcf],
                          capture_output=True, text=True).stdout.split()
    if a.outgroup not in have:
        sys.exit(f"FATAL: {a.outgroup} is not a sample of {a.panel_vcf}")

    # ONE INDEXED SWEEP, NOT 46,099 FETCHES. Asking pysam for each site in turn
    # decodes all 332 genotype columns of every record it touches; restricting
    # bcftools to one sample and a regions file does the same work once. The
    # first version of this script was still running after three minutes and
    # this one finishes in seconds.
    want = collections.defaultdict(list)
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as rf:
        for s in sites:
            if s["chrom"] == a.cohort_contig:
                rf.write(f"{a.panel_contig}\t{s['pos']}\n")
        regions = rf.name
    try:
        q = subprocess.run(
            [bcftools, "query", "-s", a.outgroup, "-R", regions,
             "-f", "%POS\t%REF\t%ALT\t[%GT]\n", a.panel_vcf],
            capture_output=True, text=True, check=True)
    finally:
        os.unlink(regions)
    for line in q.stdout.splitlines():
        pos, ref, alt, gt = line.split("\t")
        want[int(pos)].append((ref.upper(), alt.upper(), gt))

    out, tally = [], collections.Counter()
    for s in sites:
        if s["chrom"] != a.cohort_contig:
            out.append("N"); tally["off_panel_contig"] += 1
            continue
        pos, ref, alt = int(s["pos"]), s["ref"].upper(), s["alt"].upper()
        exact, other_nonref, any_rec = None, False, False
        for rref, ralt, gt in want.get(pos, ()):
            any_rec = True
            g = None if gt in (".", "./.", ".|.") else int(gt.split("/")[0]
                                                           .split("|")[0])
            if rref == ref and "," not in ralt and ralt == alt:
                exact = g
            elif g not in (0, None):
                other_nonref = True
        if exact is not None:
            out.append(alt if exact == 1 else ref)
            tally["alt" if exact == 1 else "ref_exact"] += 1
        elif exact is None and any_rec:
            if other_nonref:
                out.append("N"); tally["third_allele"] += 1
            else:
                out.append(ref); tally["ref_other_record"] += 1
        else:
            out.append(ref); tally["ref_no_record"] += 1

    with open(a.out, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n{seq}\n")
        fh.write(f">{a.outgroup}\n{''.join(out)}\n")

    n = len(out)
    print(f"  alignment    {len(seqs)} sequences x {width} sites")
    print(f"  outgroup     {a.outgroup}")
    for k in ("ref_exact", "alt", "ref_other_record", "ref_no_record",
              "third_allele", "off_panel_contig"):
        if tally[k]:
            print(f"    {k:<20s}{tally[k]:>8,}  {tally[k] / n:6.2%}")
    amb = tally["third_allele"] + tally["off_panel_contig"]
    print(f"  N in the outgroup: {amb:,} of {n:,} ({amb / n:.2%})")
    print(f"  -> {a.out}  ({len(seqs) + 1} sequences)")


if __name__ == "__main__":
    main()
