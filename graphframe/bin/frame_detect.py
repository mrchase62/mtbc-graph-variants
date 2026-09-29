#!/usr/bin/env python3
"""Decide, from the data, which frame a VCF's POS column is in.

WHY THIS IS NEEDED AND WHY ASSUMING IS NOT ENOUGH. The frames are not uniform
across this project. P1 and P2 align to refbias/build/<id>/refs/<R>.fasta and
their VCFs are in the REFS frame; the stage 2 simulation aligned to the
dnaA-rotated panel sequence and its VCFs are in the PANEL frame. Both name the
same contig and both have the same length, so nothing in a VCF header
distinguishes them. Reading one as the other moves every coordinate by the
rotation, which is what GRAPH_FRAME_RESOLUTION.md is about -- and assuming the
refs frame for stage 2 is the same mistake pointed the other way.

The test needs no metadata: take single-base REF alleles and ask which frame
reproduces them from the reference FASTA. The right frame gives ~100% and the
wrong one gives the base composition floor, around 25%, so there is no
judgement call.

    frame_detect.py --vcf <f.vcf.gz> --accession GCF_x [--refs <dir>]

Exit status is 0 when a frame is decided, 2 when neither frame fits, which
means something other than a rotation is wrong.
"""
import argparse, gzip, importlib.util, os, sys

_s = importlib.util.spec_from_file_location(
    "graph_frame", os.path.join(os.path.dirname(os.path.abspath(__file__)), "graph_frame.py"))
gf = importlib.util.module_from_spec(_s); _s.loader.exec_module(gf)

DEFAULT_REFS = "refbias/build/7713a8d71d8e/refs"


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def detect(vcf, accession, frames, refs_dir=DEFAULT_REFS, limit=2000,
           decide_at=0.90):
    """(frame, n, frac_refs, frac_panel). frame is 'refs', 'panel', 'either'
    when the accession is stored identically so the question is empty, or None
    when neither fits."""
    if frames.identity(accession):
        return "either", 0, 1.0, 1.0
    seq = load_seq(os.path.join(refs_dir, accession + ".fasta"))
    L = len(seq)
    flip = frames.flipped(accession)
    n = k_refs = k_panel = 0
    for line in _open(vcf):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        ref = f[3].upper()
        if len(ref) != 1 or ref not in "ACGT":
            continue
        try:
            p = int(f[1]) - 1
        except ValueError:
            continue
        n += 1
        if n > limit:
            break
        if seq[p % L] == ref:
            k_refs += 1
        # Read as panel: convert the index back, and COMPLEMENT the base when
        # the panel stores this accession reverse complemented. Without the
        # complement the panel hypothesis scores exactly 0.00 for every flipped
        # accession -- a complemented base is never equal to itself -- which
        # reads as "neither frame fits" and silently drops the sample.
        b = seq[frames.to_refs(accession, p) % L]
        if flip:
            b = b.translate(gf.COMP)
        if b == ref:
            k_panel += 1
    if not n:
        return None, 0, 0.0, 0.0
    fr, fp = k_refs / n, k_panel / n
    if fr >= decide_at and fr > fp:
        return "refs", n, fr, fp
    if fp >= decide_at and fp > fr:
        return "panel", n, fr, fp
    return None, n, fr, fp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--refs", default=DEFAULT_REFS)
    ap.add_argument("--table", default=None)
    ap.add_argument("--limit", type=int, default=2000)
    a = ap.parse_args()
    frame, n, fr, fp = detect(a.vcf, a.accession, gf.Frames(a.table), a.refs,
                              a.limit)
    print(f"{a.accession}\t{frame}\tn={n}\trefs={fr:.3f}\tpanel={fp:.3f}")
    return 0 if frame else 2


if __name__ == "__main__":
    sys.exit(main())
