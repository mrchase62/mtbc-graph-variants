#!/usr/bin/env python3
"""Build the ORF panel and locus window used to profile the PPE38 / plcC-glyS locus.

Step 1 of the PPE38 pipeline. Produces, under data/annotation/ppe38/:

  orf_panel.fasta          13 ORFs used to profile every genome's locus
  locus_window_h37rv.fasta 28 kb H37Rv window used to LOCATE the locus in an assembly

The panel is assembled from three sources, deliberately:

  * esxN.2 and esxJ.3 come from the supplied anchors_user.fasta. They are NOT in
    H37Rv -- neither maps to it even at asm20 -- because H37Rv is a derived
    haplotype at this locus. They have to come from outside the reference.
  * plcC, plcB, plcA, PPE38 and IS6110 come from the 2025 archive ORF set.
  * PPE39, esxP, Rv2348c, Rv2354, PPE40 and glyS are cut from H37Rv by
    coordinate. PPE39 in particular is absent from the 2025 set and is central
    to the locus, sitting between PPE38 and the IS6110-8 insertion.
"""
import argparse, os, pysam

# H37Rv gene intervals, 1-based inclusive, with strand (from the RefSeq annotation)
H37RV_GENES = [
    ("esxP",    2626222, 2626519, "-"),
    ("Rv2348c", 2626653, 2626980, "-"),
    ("PPE39",   2634527, 2635592, "-"),
    ("Rv2354",  2635627, 2635954, "+"),
    ("PPE40",   2637687, 2639535, "-"),
    ("glyS",    2639672, 2641064, "-"),
]
FROM_ARCHIVE = ["plcC", "plcB", "plcA", "PPE38", "IS6110"]
FROM_ANCHORS = ["esxN.2", "esxJ.3"]
WINDOW = (2620000, 2648000)   # generous flank around plcC-glyS
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="PanSN H37Rv FASTA (indexed)")
    ap.add_argument("--archive-orfs", required=True,
                    help="2025 PPE38_esxN.2_esxJ.3.fasta")
    ap.add_argument("--anchors", required=True,
                    help="anchors_user.fasta with esxN.2 and esxJ.3")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ref = pysam.FastaFile(args.ref); chrom = ref.references[0]
    arc = pysam.FastaFile(args.archive_orfs)
    anc = pysam.FastaFile(args.anchors)

    panel = []
    for n in FROM_ARCHIVE:
        panel.append((n, arc.fetch(n).upper()))
    for n in FROM_ANCHORS:
        s = anc.fetch(n).upper()
        # store on the same strand as the archive copies so hits are comparable
        if n in arc.references and arc.fetch(n).upper() == rc(s):
            s = rc(s)
        panel.append((n, s))
    for n, s0, e0, strand in H37RV_GENES:
        seq = ref.fetch(chrom, s0 - 1, e0).upper()
        panel.append((n, rc(seq) if strand == "-" else seq))

    p = os.path.join(args.outdir, "orf_panel.fasta")
    with open(p, "w") as fh:
        for n, s in panel:
            fh.write(f">{n}\n{s}\n")
    pysam.faidx(p)

    w = ref.fetch(chrom, WINDOW[0] - 1, WINDOW[1]).upper()
    wp = os.path.join(args.outdir, "locus_window_h37rv.fasta")
    with open(wp, "w") as fh:
        fh.write(f">H37Rv_plcC_glyS_window {WINDOW[0]}-{WINDOW[1]}\n{w}\n")
    pysam.faidx(wp)

    print(f"  orf_panel.fasta          {len(panel)} ORFs")
    for n, s in panel:
        print(f"    {n:<9s} {len(s):>5d} bp")
    print(f"  locus_window_h37rv.fasta {len(w)} bp  (H37Rv {WINDOW[0]}-{WINDOW[1]})")


if __name__ == "__main__":
    main()
