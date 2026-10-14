#!/usr/bin/env python3
"""Add the outgroup kit's genomes to the panel SNP alignment, for the panel tree.

WHY. The panel tree (docs/PANEL_TREE.md section 1) is built from the graph
VCF's SNP alignment (bin/vcf_to_alignment.py), which has a row only for the
genomes in the graph. The P200 panel keeps lineage 8 (GCF_012923765) and
canettii ET1291 (GCF_035581225) OUT of the graph, as the "outgroup kit", but
the ancestral reconstruction (bin/panel_ancestral.py) needs them on the tree:
canettii roots it, and both take part in Fitch at SNP keys. Their rows come
from their direct alignment to H37Rv instead.

THE ROW. For each alignment column (sites table: chrom, pos, ref, alt), the
kit genome's base at that H37Rv position from its states file
(analysis/panel_checks/bin/genome_states.py: npz array 'state', one uint8 per
H37Rv base, 0 = N, 1-4 = ACGT): the base if it is the column's REF or ALT,
N if uncalled or a third base. Graph genomes carrying a third base are N in
vcf_to_alignment.py's rows too (audit TP-7), so every column stays the
biallelic site it was made as, and panel_ancestral.py reads the kit's states
by the same rule.

WHAT IS LOST. The columns are the panel's variable sites only, so a kit
genome's own private SNPs are not in the alignment and its terminal branch is
short. The topology, and the rooting on canettii, are what the
reconstruction uses; branch lengths are not.

    $MTB_PY bin/append_kit_rows.py --alignment data/trees/<panel>.snps.fasta \\
        --sites data/trees/<panel>.sites.tsv \\
        --kit GCF_012923765=<states>/GCF_012923765.npz \\
        --kit GCF_035581225=<states>/GCF_035581225.npz \\
        --out data/trees/<panel>.kit.snps.fasta
"""
import argparse
import csv
import sys

import numpy as np

BASES = "NACGT"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alignment", required=True,
                    help="vcf_to_alignment.py's FASTA")
    ap.add_argument("--sites", required=True,
                    help="vcf_to_alignment.py's --sites-out table")
    ap.add_argument("--kit", action="append", required=True,
                    metavar="ACC=STATES.npz",
                    help="a kit genome and its genome_states.py npz; repeatable")
    ap.add_argument("--h37rv-contig", default="NC_000962.3",
                    help="sites on other contigs get N")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    names, seqs = [], {}
    for line in open(a.alignment):
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            seqs[names[-1]] = []
        elif names:
            seqs[names[-1]].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    if not seqs:
        sys.exit(f"FATAL: no sequences in {a.alignment}")
    width = {len(s) for s in seqs.values()}
    if len(width) != 1:
        sys.exit(f"FATAL: {a.alignment} rows differ in length: {sorted(width)}")
    width = width.pop()
    sites = list(csv.DictReader(open(a.sites), delimiter="\t"))
    if len(sites) != width:
        sys.exit(f"FATAL: {a.sites} has {len(sites)} sites for {width} columns")

    kit = []
    for arg in a.kit:
        if "=" not in arg:
            sys.exit(f"FATAL: --kit takes ACC=STATES.npz, not {arg!r}")
        acc, path = arg.split("=", 1)
        if acc in seqs or acc in dict(kit):
            sys.exit(f"FATAL: {acc} already has a row; a kit genome must be "
                     f"outside the graph")
        d = np.load(path)
        if "state" not in d.files:
            sys.exit(f"FATAL: {path} has no 'state' array")
        kit.append((acc, d["state"]))

    with open(a.out, "w") as fh:
        for n in names:
            fh.write(f">{n}\n")
            for i in range(0, width, 60):
                fh.write(seqs[n][i:i + 60] + "\n")
        for acc, st in kit:
            row, called = [], 0
            for s in sites:
                pos, ref, alt = int(s["pos"]), s["ref"].upper(), s["alt"].upper()
                b = "N"
                if s["chrom"].split("#")[-1] == a.h37rv_contig:
                    if not 1 <= pos <= len(st):
                        sys.exit(f"FATAL: site {pos} is outside {acc}'s states "
                                 f"({len(st):,} bases)")
                    b = BASES[st[pos - 1]] if st[pos - 1] < 5 else "N"
                    if b not in (ref, alt):
                        b = "N"
                called += b != "N"
                row.append(b)
            seq = "".join(row)
            fh.write(f">{acc}\n")
            for i in range(0, width, 60):
                fh.write(seq[i:i + 60] + "\n")
            print(f"  {acc}: {called:,} of {width:,} columns called")
    print(f"  -> {a.out}  ({len(names) + len(kit)} rows)")


if __name__ == "__main__":
    main()
