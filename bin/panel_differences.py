#!/usr/bin/env python3
"""P0: each panel genome's own differences from H37Rv, as one table.

WHY THIS EXISTS (HANDOFF fix list, ABS-1 and ABS-2). P5 states a site ABSENT
when the sample's matched reference R lacks it, and the SV block reports R's
own large differences from H37Rv as the sample's ("inherited"). Neither was
checked against the sample's reads, and measured on scale200_fix that was
wrong often enough to matter: 2,279 + 161 core ABSENT cells contradicted by
the reads and by local assembly (analysis/inherited_absent/README.md), and
every inherited insertion reported present in every sample matched to its R
(analysis/inherited_matched_check/RESULTS.md).

The checks that replace the assumption need, per sample, WHICH of R's
differences explain a call: the ends of R's deletion over an ABSENT site (to
look for reads that cross them), R's allele for an inherited SV (to place its
junction on R's own sequence), and the sites where R has no genotype at all
(its path does not cross the site in place -- an inversion, a rearrangement,
an assembly break -- which no single junction can test). All of that is a
property of the panel, not of a cohort, so it is read once per build from the
collapsed graph VCF and written here, rather than re-parsed from the VCF by
every per-sample task.

WHAT A ROW IS. One per (accession, graph record) where the accession's allele
differs from H37Rv in one of three ways:

  DEL    the accession's allele is shorter than REF. The shared prefix is
         trimmed and what remains of REF is deleted: `start`..`end` are the
         deleted H37Rv bases, 1-based and inclusive. Every size, not only
         >= 50 bp, because a small-variant ABSENT can sit in a 5 bp deletion
         of R (analysis/inherited_absent/junction_check.py, r_deletions(),
         the reference version of this rule).
  INS    the accession's allele is longer than REF by at least --min-ins
         (50 bp, the SV threshold P4b uses). `start` = `end` = the VCF POS,
         the anchor base the insertion follows.
  NOGT   the accession has no genotype (`.`) at the record: `start`..`end`
         span REF.

`vcf_pos` and `allele` (the accession's ALT index) identify the record, so a
consumer that needs the alleles themselves (R's sequence at an inherited SV)
reads that one record from the same graph VCF by tabix rather than from a
copy here. The table and the VCF are one build's assets, and the header line
`#graph_vcf_sha256` records which VCF the table was read from, so a consumer
refuses a table made from another graph VCF (ABS-2 item 6: one graph VCF per
run; truth7's P4b rows came from an older VCF whose records sat a few bp from
the current one's).

H37Rv itself is not a column of the graph VCF -- it is the VCF's reference
path -- and has no row: a sample matched to H37Rv has no inherited
difference. Every accession that IS a column is listed in the `#accessions`
header line, so a consumer can tell "this genome has no differences" from
"this genome is not in the table", and refuses the latter.

    panel_differences.py --graph-vcf <build>/assets/graph_collapsed.vcf.gz \\
        --out <build>/assets/panel_differences.tsv.gz
"""
import argparse, collections, gzip, hashlib, os, sys

COLUMNS = ["accession", "kind", "start", "end", "svlen", "vcf_pos", "ref_len",
           "allele"]


def record_differences(pos, ref, alts, gts, panel, min_ins):
    """(accession, kind, start, end, svlen, allele) for one graph record.

    `gts` are the genotype strings in panel order. A deletion's deleted bases
    are REF after the prefix it shares with the allele; the trailing part of a
    complex allele is counted as deleted, as junction_check.py did."""
    out = []
    ref = ref.upper()
    for acc, g in zip(panel, gts):
        g = g.split(":", 1)[0]
        if g == ".":
            out.append((acc, "NOGT", pos, pos + len(ref) - 1, len(ref), 0))
            continue
        if not g.isdigit() or g == "0":
            continue
        i = int(g)
        if i > len(alts):
            sys.exit(f"FATAL: genotype {g} at {pos} for {acc} names allele {i} "
                     f"of {len(alts)}; the graph VCF is malformed")
        alt = alts[i - 1].upper()
        if len(alt) < len(ref):
            k = 0
            while k < len(alt) and ref[k] == alt[k]:
                k += 1
            s, e = pos + k, pos + len(ref) - 1 - (len(alt) - k)
            if e >= s:
                out.append((acc, "DEL", s, e, e - s + 1, i))
        elif len(alt) - len(ref) >= min_ins:
            out.append((acc, "INS", pos, pos, len(alt) - len(ref), i))
    return out


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph-vcf", required=True,
                    help="<build>/assets/graph_collapsed.vcf.gz")
    ap.add_argument("--min-ins", type=int, default=50,
                    help="insertions shorter than this are not listed; the "
                         "inherited SV half (P4b) starts at 50 bp")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if not os.path.exists(a.graph_vcf):
        sys.exit(f"FATAL: no graph VCF at {a.graph_vcf}")

    op = gzip.open if a.graph_vcf.endswith(".gz") else open
    panel, rows, n_rec = None, [], 0
    with op(a.graph_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                panel = line.rstrip("\n").split("\t")[9:]
                continue
            if panel is None:
                sys.exit(f"FATAL: {a.graph_vcf} has a record before its #CHROM line")
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 + len(panel):
                sys.exit(f"FATAL: record at {f[1] if len(f) > 1 else '?'} has "
                         f"{len(f) - 9} genotypes for {len(panel)} panel genomes")
            n_rec += 1
            pos, ref, alts = int(f[1]), f[3], f[4].split(",")
            for acc, kind, s, e, ln, i in record_differences(
                    pos, ref, alts, f[9:], panel, a.min_ins):
                rows.append((acc, kind, s, e, ln, pos, len(ref), i))
    if panel is None:
        sys.exit(f"FATAL: {a.graph_vcf} has no #CHROM line")
    if not n_rec:
        sys.exit(f"FATAL: {a.graph_vcf} has no records")
    rows.sort(key=lambda r: (r[0], r[2], r[3], r[1]))
    tmp = a.out + ".tmp"
    with (gzip.open(tmp, "wt") if a.out.endswith(".gz") else open(tmp, "w")) as fo:
        fo.write("#accessions\t" + ",".join(panel) + "\n")
        fo.write("#graph_vcf_sha256\t" + file_sha256(a.graph_vcf) + "\n")
        fo.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")
    os.replace(tmp, a.out)
    n = collections.Counter(r[1] for r in rows)
    print(f"  {n_rec:,} graph records over {len(panel)} panel genomes: "
          + ", ".join(f"{k} {n[k]:,}" for k in ("DEL", "INS", "NOGT")))
    print(f"  -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
