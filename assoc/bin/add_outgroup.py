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

THE ALLELE RULE, and where it says N rather than guessing. Every panel record
whose REF span covers the site is consulted, not only those starting at it: an
outgroup MNP, deletion or complex record that starts upstream says what the
outgroup has at the site too.

  an exact (pos, ref, alt) record with GT 1       -> ALT. The collapsed VCF
                                                     has one record per key;
                                                     the decomposed one repeats
                                                     a key once per allele path
                                                     with the outgroup's ALT in
                                                     any one copy, so every
                                                     copy is consulted
  exact records present, every one missing        -> N
  a covering record the outgroup carries (GT>=1)  -> the base that allele puts
                                                     at the site: an SNP or MNP
                                                     gives its own base, an
                                                     indel's unchanged anchor
                                                     base gives REF, and any
                                                     other position inside an
                                                     indel or complex allele
                                                     (a deletion included) is
                                                     N. A base that is neither
                                                     REF nor ALT, or two
                                                     records that disagree, is
                                                     N, since a third allele is
                                                     not representable in a
                                                     biallelic column
  an exact record with GT 0 and nothing above     -> REF
  covering records, all missing or GT 0, at least
  one missing and no exact GT 0                   -> N: the outgroup has no
                                                     call over that stretch
  no covering record at all                       -> REF. The panel VCF holds
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
import argparse, bisect, collections, csv, os, subprocess, sys

# THE COLLAPSED FILE, NOT THE DECOMPOSED ONE (audit PGB-9): one record per
# trimmed (pos, ref, alt) with each genome's genotype unioned over vcfwave's
# duplicates, and multiallelic records split, so an allele is found under its
# minimal key. The rule above still tolerates duplicate keys, so the decomposed
# file gives the same answer at every site it holds (CX333: identical outgroup
# rows for scale200 and gwas1000).
#
# NO DEFAULT FILE: it was CX333's, so a new build's chain that omitted the
# argument read the old graph's outgroup genotypes. The chain passes
# <build>/assets/graph_collapsed.vcf.gz, P0's copy of the graph's
# all_variants.collapsed.vcf.gz (step assets).
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


def base_at(rpos, rref, allele, site):
    """The base an allele of a record starting at rpos puts at `site`, or None
    where it cannot be said: inside a deletion, or past the anchor of an indel
    or complex allele whose bases do not line up with REF's."""
    off = site - rpos
    if len(allele) == len(rref):             # SNP or MNP: base for base
        return allele[off]
    if off == 0 and allele[:1] == rref[:1]:  # an indel's unchanged anchor
        return rref[0]
    return None


def outgroup_allele(pos, ref, alt, recs):
    """REF, ALT or N for the outgroup at one biallelic site; see THE ALLELE
    RULE above. recs: (POS, REF, [ALTs], GT or None) for every panel record
    whose REF span covers pos. Returns (call, tally key)."""
    def is_exact(rpos, rref, ralts, g):
        return (rpos == pos and rref == ref and alt in ralts
                and (g is None or g == 0 or ralts[g - 1] == alt))
    exact = [r[3] for r in recs if is_exact(*r)]
    if any(g is not None and g > 0 for g in exact):
        return "ALT", "alt"
    if exact and all(g is None for g in exact):
        return "N", "missing_exact"
    bases, unknown, missing = set(), False, False
    if exact:                                 # at least one GT 0 here
        bases.add(ref)
    for rpos, rref, ralts, g in recs:
        if is_exact(rpos, rref, ralts, g):
            continue                          # the exact records, above
        if g is None:
            missing = True
            continue
        if g == 0:
            continue
        b = base_at(rpos, rref, ralts[g - 1], pos)
        if b is None:
            unknown = True
        else:
            bases.add(b)
    if unknown:
        return "N", "in_outgroup_indel"
    if len(bases) > 1:
        return "N", "third_allele"
    if bases:
        b = next(iter(bases))
        if b == alt:
            return "ALT", "alt_covering_record"
        if b == ref:
            return "REF", "ref_exact" if exact else "ref_covering_record"
        return "N", "third_allele"
    if missing:
        return "N", "missing_covering"
    return "REF", "ref_covering_record" if recs else "ref_no_record"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--panel-vcf", required=True,
                    help="the build's collapsed graph VCF")
    ap.add_argument("--panel-contig", default=PANEL_CONTIG)
    # no default: the outgroup is the build's (build_info.tsv `outgroup`,
    # from config MTB_OUTGROUP), which cohort_assoc_tail.sh passes
    ap.add_argument("--outgroup", required=True,
                    help="the panel sample to emit: the build's outgroup")
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

    # ONE SWEEP, NOT 46,099 FETCHES. Asking pysam for each site in turn
    # decodes all 332 genotype columns of every record it touches; restricting
    # bcftools to one sample and streaming the contig once does the work in
    # seconds. (A regions file of every site, which this used before, makes
    # bcftools seek once per site and took over six minutes on scale200.)
    #
    # Each record is filed under EVERY site its REF span covers, not only under
    # its own POS: keyed on POS alone, an outgroup MNP or deletion starting
    # before the site was never consulted and the site fell through to "no
    # record -> REF" (audit TP-1).
    site_pos = sorted({int(s["pos"]) for s in sites
                       if s["chrom"] == a.cohort_contig})
    want = collections.defaultdict(list)
    q = subprocess.run(
        [bcftools, "query", "-s", a.outgroup, "-r", a.panel_contig,
         "-f", "%POS\t%REF\t%ALT\t[%GT]\n", a.panel_vcf],
        capture_output=True, text=True, check=True)
    for line in q.stdout.splitlines():
        pos, ref, alt, gt = line.split("\t")
        pos, ref = int(pos), ref.upper()
        g = None if gt in (".", "./.", ".|.") else int(gt.split("/")[0]
                                                       .split("|")[0])
        rec = (pos, ref, alt.upper().split(","), g)
        i = bisect.bisect_left(site_pos, pos)
        while i < len(site_pos) and site_pos[i] <= pos + len(ref) - 1:
            want[site_pos[i]].append(rec)
            i += 1

    out, tally = [], collections.Counter()
    for s in sites:
        if s["chrom"] != a.cohort_contig:
            out.append("N"); tally["off_panel_contig"] += 1
            continue
        pos, ref, alt = int(s["pos"]), s["ref"].upper(), s["alt"].upper()
        call, why = outgroup_allele(pos, ref, alt, want.get(pos, ()))
        out.append({"REF": ref, "ALT": alt}.get(call, "N"))
        tally[why] += 1

    with open(a.out, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n{seq}\n")
        fh.write(f">{a.outgroup}\n{''.join(out)}\n")

    n = len(out)
    print(f"  alignment    {len(seqs)} sequences x {width} sites")
    print(f"  outgroup     {a.outgroup}")
    for k in ("ref_exact", "alt", "ref_covering_record", "alt_covering_record",
              "ref_no_record", "missing_exact", "missing_covering",
              "in_outgroup_indel", "third_allele", "off_panel_contig"):
        if tally[k]:
            print(f"    {k:<20s}{tally[k]:>8,}  {tally[k] / n:6.2%}")
    amb = sum(tally[k] for k in ("missing_exact", "missing_covering",
                                 "in_outgroup_indel", "third_allele",
                                 "off_panel_contig"))
    print(f"  N in the outgroup: {amb:,} of {n:,} ({amb / n:.2%})")
    print(f"  -> {a.out}  ({len(seqs) + 1} sequences)")


if __name__ == "__main__":
    main()
