#!/usr/bin/env python3
"""Recurrent read-artifact sites, masked as uncallable in each matched
reference (user's decision, 2026-10-10). Written beside callable_mask.py's
BEDs as <reference>.artifact.bed; score.py unions the two only with
--artifact-mask (off by default: under evaluation, user 2026-10-10).

One site so far, found in phaseB/dysgu_false.py: a 67 bp palindromic element
(GGGCTGGCGAGCAGACGCAAAATCCCCCGCACGCCCGGCGTGTCGGGGGATTTTGCGTCTGCTCGCC, a
strong hairpin between short direct repeats) near 2.35 Mb. dysgu and our
caller call its deletion at AF 0.2-0.4 in 12 setE isolates and in 3 clean
controls whose own genome is the reference, so it is not a fixed event. The
element is present in 41 of the 51 matched references and absent from 10,
so the site also differs between strains; masking gives up any true event
there.

A site is located by its two flanking 30-mers, exact match on either strand,
and the sequence between them (the element, or its absence) is masked.
Where only one flank is found (rearranged), 100 bp either side of it is
masked. Neither found: nothing masked, reported.

Run from runroot. Reads the references in the build used for the masks
(07_masks.sbatch), for every <reference>.bed in --mask-dir.
"""
import argparse
import glob
import os
import re

SITES = [("palindrome_67bp", "CGGCCAGCTCAGTCACGTCGCCGCCGCCTC", "TTGACCGCGCCCGCTCGCGGCTAGCGGGCC")]
COMP = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                break
            name = line[1:].split()[0]
        else:
            seq.append(line.strip())
    return name, "".join(seq).upper()


def find(seq, q):
    rc = q.translate(COMP)[::-1]
    return [m.start() for m in re.finditer(q, seq)] + [m.start() for m in re.finditer(rc, seq)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refs", default="refbias/build/7713a8d71d8e-fix1/refs")
    ap.add_argument("--mask-dir", default="../analysis/novel_events/out/masks")
    a = ap.parse_args()
    refs = sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{a.mask_dir}/*.bed")
                  if not p.endswith(".artifact.bed"))
    for g in refs:
        name, seq = read_fasta(f"{a.refs}/{g}.fasta")
        rows = []
        for site, left, right in SITES:
            L, R = find(seq, left), find(seq, right)
            if len(L) > 1 or len(R) > 1:
                print(f"{g}\t{site}\tflank found more than once, not masked")
                continue
            if L and R and abs(L[0] - R[0]) < 1000:
                lo, hi = min(L[0], R[0]), max(L[0], R[0]) + 30
                how = "between flanks"
            elif L or R:
                p = (L or R)[0]
                lo, hi = max(0, p - 100), p + 130
                how = "one flank, +-100 bp"
            else:
                print(f"{g}\t{site}\tnot found")
                continue
            rows.append((name, lo, hi, site))
            print(f"{g}\t{site}\t{lo}-{hi}\t{hi - lo} bp\t{how}")
        with open(f"{a.mask_dir}/{g}.artifact.bed", "w") as fo:
            for r in rows:
                fo.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
