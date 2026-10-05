#!/usr/bin/env python3
"""Per-assembly length, ambiguity and contig stats -- the divergence-independent QC.

The 2025 collinearity QC keyed on MUMmer block count against H37Rv, which in MTBC
is dominated by real biology (IS6110 positions, PE_PGRS, RDs) rather than assembly
error: the median is 93 blocks against a threshold of 8, so 75% of genomes were
flagged. The statistics here do not depend on divergence from the reference, so
they separate assembly quality from lineage.
"""
import argparse, glob, gzip, os, sys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for p in sorted(glob.glob(os.path.join(a.dir, "*.fna.gz"))):
        acc = os.path.basename(p).replace(".fna.gz", "")
        seqs, cur = [], []
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur: seqs.append("".join(cur)); cur = []
                else:
                    cur.append(line.strip().upper())
        if cur: seqs.append("".join(cur))
        total = sum(len(s) for s in seqs)
        n = sum(s.count("N") for s in seqs)
        amb = sum(1 for s in seqs for c in s if c not in "ACGTN")
        runs = 0; longest = 0
        for s in seqs:
            i = 0
            while i < len(s):
                if s[i] == "N":
                    j = i
                    while j < len(s) and s[j] == "N": j += 1
                    runs += 1; longest = max(longest, j - i); i = j
                else:
                    i += 1
        rows.append((acc, len(seqs), total, n, round(100.0 * n / total, 5) if total else 0,
                     runs, longest, amb))
    with open(a.out, "w") as fh:
        fh.write("accession\tcontigs\tlength\tN_count\tN_pct\tN_runs\tlongest_N_run\tother_ambiguous\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"[assembly_qc_stats] {len(rows)} assemblies -> {a.out}")


if __name__ == "__main__":
    main()
