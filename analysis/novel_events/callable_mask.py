#!/usr/bin/env python3
"""Uncallable sequence in a matched reference R's own coordinates, for
novel-event calling (the user's decision, 2026-10-09: call novel events only
in mappable sequence and flag the rest as uncallable).

Novel events are called in R's frame, so the H37Rv repeat mask cannot be used
as it is. Two sources, unioned:

1. Measured on R, with the same tests as bin/build_repeat_mask.py (imported,
   not copied):
   - paralogy: a position inside a 50-mer that occurs more than once in R,
     counting both strands (reads from one copy pile onto another);
   - tandem: a 300 bp window whose 9-mers recurring 3 or more times cover
     10% or more of it (the alignment slides, indel positions are ambiguous).
   Intervals of 100 bp or more, gaps of 50 bp or less closed.
   Dropped from build_repeat_mask.py, because they need an H37Rv gene
   annotation that R does not have: the gene extension (a gene 25% masked is
   masked in full) and the name-based list. Step 2 brings both in where R
   shares the sequence with H37Rv.
2. The H37Rv repeat mask (the build's assets/repeat_mask.bed), lifted to R:
   H37Rv aligned to R with minimap2 -cx asm5, mask intervals carried over
   with paftools liftover. An interval whose ends do not both lift (it spans
   a difference between the genomes) is lost; source 1 still covers R's
   own repeats there.

Output BED (R contig, 0-based start, end, reason: paralog, tandem, h37rv_mask
or a combination), and one summary line on stderr.
"""
import argparse
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin"))
from build_repeat_mask import nonunique_positions, tandem_positions, to_intervals  # noqa: E402


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    break
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    return name, "".join(seq).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="matched reference FASTA (one contig)")
    ap.add_argument("--h37rv", required=True, help="H37Rv FASTA")
    ap.add_argument("--h37rv-mask", required=True, help="the build's repeat_mask.bed")
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--k8", required=True)
    ap.add_argument("--paftools", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ctg, seq = read_fasta(a.ref)
    N = len(seq)
    par, _ = nonunique_positions(seq, 50)
    tan, _ = tandem_positions(seq, 300, 50, 0.10)
    why = [0] * N  # bit 1 paralog, 2 tandem, 4 lifted
    for s, e in to_intervals(par, 100, 50):
        for j in range(s, e):
            why[j] |= 1
    for s, e in to_intervals(tan, 100, 50):
        for j in range(s, e):
            why[j] |= 2
    h_ctg, _ = read_fasta(a.h37rv)
    with tempfile.TemporaryDirectory() as td:
        bed = os.path.join(td, "m.bed")
        with open(a.h37rv_mask) as fi, open(bed, "w") as fo:
            for line in fi:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and f[1].isdigit():
                    fo.write(f"{h_ctg}\t{f[1]}\t{f[2]}\n")
        paf = os.path.join(td, "h.paf")
        with open(paf, "w") as fo:
            subprocess.run([a.minimap2, "-cx", "asm5", "--secondary=no", a.ref, a.h37rv],
                           stdout=fo, stderr=subprocess.DEVNULL, check=True)
        lifted = subprocess.run([a.k8, a.paftools, "liftover", paf, bed],
                                capture_output=True, text=True, check=True).stdout
    n_lift = 0
    for line in lifted.splitlines():
        f = line.split("\t")
        if len(f) >= 3 and f[0] == ctg:
            n_lift += 1
            for j in range(int(f[1]), min(N, int(f[2]))):
                why[j] |= 4
    names = {1: "paralog", 2: "tandem", 4: "h37rv_mask"}
    rows, s = [], None
    for j in range(N + 1):
        v = why[j] if j < N else 0
        if s is not None and v != why[s]:
            rows.append((s, j, why[s]))
            s = None
        if v and s is None:
            s = j
    with open(a.out, "w") as fo:
        for s, e, v in rows:
            fo.write(f"{ctg}\t{s}\t{e}\t{','.join(n for b, n in names.items() if v & b)}\n")
    tot = sum(1 for v in why if v)
    print(f"{os.path.basename(a.ref)}: {N:,} bp, uncallable {tot:,} bp ({tot / N:.1%}); "
          f"measured {sum(1 for v in why if v & 3):,}, lifted {sum(1 for v in why if v & 4):,} "
          f"({n_lift} intervals lifted)", file=sys.stderr)


if __name__ == "__main__":
    main()
