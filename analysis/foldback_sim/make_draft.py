#!/usr/bin/env python3
"""A synthetic long-read draft: the truth genome with known errors of the
kind short-read polishing exists to fix. 300 one-base indels in homopolymers
of 3 bp or more (half insertions, half deletions; the typical ONT/PacBio CLR
error) and 30 substitutions, at random positions at least 200 bp apart.
Writes the draft FASTA and a table of the errors in draft coordinates."""
import argparse
import random


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--errors", required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    name, seq = None, []
    for line in open(a.truth):
        if line.startswith(">"):
            name = name or line[1:].split()[0]
        else:
            seq.append(line.strip())
    g = "".join(seq).upper()
    homo = [i for i in range(1000, len(g) - 1000) if g[i] == g[i + 1] == g[i + 2] and g[i - 1] != g[i]]
    picks = sorted(rng.sample(homo, 3000))
    sites, last = [], -1000
    for i in picks:
        if i - last >= 200:
            sites.append(i)
            last = i
    sites = rng.sample(sites, 300)
    other = [i for i in rng.sample(range(1000, len(g) - 1000), 200)
             if all(abs(i - s) >= 200 for s in sites)][:30]
    edits = [(i, "ins" if k % 2 == 0 else "del") for k, i in enumerate(sorted(sites))] + \
            [(i, "snp") for i in other]
    edits.sort(reverse=True)
    s = list(g)
    for i, t in edits:
        if t == "ins":
            s.insert(i, g[i])
        elif t == "del":
            del s[i]
        else:
            s[i] = rng.choice([b for b in "ACGT" if b != g[i]])
    with open(a.out, "w") as fo:
        fo.write(f">{name}_draft\n")
        d = "".join(s)
        for x in range(0, len(d), 80):
            fo.write(d[x:x + 80] + "\n")
    with open(a.errors, "w") as fo:
        fo.write("truth_pos\ttype\n")
        for i, t in sorted(edits):
            fo.write(f"{i + 1}\t{t}\n")


if __name__ == "__main__":
    main()
