#!/usr/bin/env python3
"""Per-assembly length, ambiguity and contig stats -- the divergence-independent QC.

The 2025 collinearity QC keyed on MUMmer block count against H37Rv, which in MTBC
is dominated by real biology (IS6110 positions, PE_PGRS, RDs) rather than assembly
error: the median is 93 blocks against a threshold of 8, so 75% of genomes were
flagged. The statistics here do not depend on divergence from the reference, so
they separate assembly quality from lineage.

SINGLE-BASE RUNS (audit PGB-7, HANDOFF 0f). A pure run of one base of
--max-base-run (100) bp or more is flagged (`flag` = single_base_run). No MTBC
genome carries one: in CX333 the other 330 genomes have no run of even 50 bp.
The three that do are assembly artifacts -- GCF_050259585 (1,015 bp poly-G,
the two-colour Illumina no-signal artifact), GCF_045348265 (484 bp poly-A
inside a 1,243 bp A-rich stretch, ONT) and GCF_039770655 (377 bp poly-T in
(CCATT)n, the lineage 9 reference). Flagged genomes need a recorded decision
(mask, replace or exclude) before build_panel.py will write a panel.
"""
import argparse, glob, gzip, itertools, os, sys


def longest_base_run(seqs):
    """(length, base, contig index, 1-based start) of the longest run of one
    A/C/G/T base, and the number of such runs >= each threshold is left to the
    caller via base_runs()"""
    best = (0, "", -1, 0)
    for ci, s in enumerate(seqs):
        pos = 0
        for b, g in itertools.groupby(s):
            n = sum(1 for _ in g)
            if b in "ACGT" and n > best[0]:
                best = (n, b, ci, pos + 1)
            pos += n
    return best


def base_runs(seqs, min_len):
    """[(contig index, 1-based start, length, base)] for runs >= min_len"""
    out = []
    for ci, s in enumerate(seqs):
        pos = 0
        for b, g in itertools.groupby(s):
            n = sum(1 for _ in g)
            if b in "ACGT" and n >= min_len:
                out.append((ci, pos + 1, n, b))
            pos += n
    return out


def read_seqs(p):
    op = gzip.open if p.endswith(".gz") else open
    seqs, cur = [], []
    with op(p, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur: seqs.append("".join(cur)); cur = []
            else:
                cur.append(line.strip().upper())
    if cur: seqs.append("".join(cur))
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--suffix", default=".fna.gz",
                    help="assembly file suffix; .dnaA_rotated.fasta for data/rotated")
    ap.add_argument("--max-base-run", type=int, default=100,
                    help="flag a pure single-base run of at least this length")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    paths = sorted(glob.glob(os.path.join(a.dir, "*" + a.suffix)))
    if not paths:
        sys.exit(f"no *{a.suffix} files in {a.dir}")
    for p in paths:
        acc = os.path.basename(p)[:-len(a.suffix)]
        seqs = read_seqs(p)
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
        bl, bb, _, bp = longest_base_run(seqs)
        nbr = len(base_runs(seqs, a.max_base_run))
        rows.append((acc, len(seqs), total, n, round(100.0 * n / total, 5) if total else 0,
                     runs, longest, amb, bl, bb, bp, nbr,
                     "single_base_run" if nbr else ""))
    with open(a.out, "w") as fh:
        fh.write("accession\tcontigs\tlength\tN_count\tN_pct\tN_runs\tlongest_N_run\t"
                 "other_ambiguous\tlongest_base_run\tlongest_base_run_base\t"
                 f"longest_base_run_pos\tbase_runs_ge{a.max_base_run}\tflag\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    fl = [r for r in rows if r[-1]]
    print(f"[assembly_qc_stats] {len(rows)} assemblies -> {a.out}")
    print(f"  single-base run >= {a.max_base_run} bp: {len(fl)}")
    for r in fl:
        print(f"    {r[0]}  {r[8]} bp of {r[9]} at {r[10]:,}")


if __name__ == "__main__":
    main()
