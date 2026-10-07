#!/usr/bin/env python3
"""The one place that knows which coordinate frame a position is in.

THE PROBLEM THIS EXISTS TO STOP. The pangenome graph was built from
`data/fastas/mtb.complex333.fasta.gz`, whose sequences are dnaA-rotated. Read
alignment uses `refbias/build/<id>/refs/<accession>.fasta`, the sequence as
deposited. For 110 of the 333 accessions those two are the same sequence written
from a different origin, and 22 of those are also written on the opposite
strand. `odgi position` reads and writes PANEL coordinates. Every BAM, VCF, gVCF,
GFF and depth file in this project is in REFS coordinates. Handing one to the
other without converting silently queries a base hundreds of kilobases away, and
`graphframe/docs/GRAPH_FRAME_RESOLUTION.md` records four call sites that did.

The conversion is exact and lossless -- the two strings have equal length and
differ only by origin and strand -- so nothing is approximated here.

USE IT AT THE BOUNDARY, NOWHERE ELSE. Convert when writing an `odgi position`
input file and when reading its output, and keep every other coordinate in the
project in the refs frame. A coordinate that travels any distance in the panel
frame is a coordinate whose frame will eventually be forgotten.

    import importlib.util, os
    _s = importlib.util.spec_from_file_location(
        "graph_frame", os.path.join("graphframe", "bin", "graph_frame.py"))
    graph_frame = importlib.util.module_from_spec(_s); _s.loader.exec_module(graph_frame)
    fr = graph_frame.Frames()

    fr.to_panel(acc, refs_pos0)     -> panel_pos0
    fr.to_refs(acc, panel_pos0)     -> refs_pos0
    fr.to_refs_interval(acc, lo, hi)-> (lo, hi) in refs, re-sorted if flipped
    fr.flipped(acc)                 -> True when the panel path is the reverse
                                       complement of the refs FASTA
    fr.complement_needed(acc, odgi_strand)
                                    -> whether a base read at the projected
                                       position must be complemented before it
                                       is compared with a refs-frame allele

WHY `complement_needed` IS NOT JUST odgi's OWN COLUMN. `odgi position` reports
`strand.vs.ref` for the relation between the source and target STEPS inside the
graph. It knows nothing about the refs FASTA, so for a globally flipped
accession it reports `+` while the base still needs complementing. Measured: 13
of 13 projected rows on GCF_965124535 match only after complementing, with odgi
reporting `+` for all 13.
"""
import csv, os

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# The build's frame table (P0 step frames). It was graphframe/results/, the
# CX333 graph's table, read whenever MTB_GRAPH_FRAMES was unset -- so a script
# run on a new build outside refbias_run.sh would have converted every
# projection with the old graph's frames. Now MTB_GRAPH_FRAMES, else
# <MTB_BUILD_DIR>/assets/graph_frame_offsets.tsv, else an error.
DEFAULT_TABLE = (os.path.join(os.environ["MTB_BUILD_DIR"], "assets",
                              "graph_frame_offsets.tsv")
                 if os.environ.get("MTB_BUILD_DIR") else "")


def rc(s):
    return s.translate(COMP)[::-1]


class Frames:
    def __init__(self, table=None):
        self.table = table or os.environ.get("MTB_GRAPH_FRAMES") or DEFAULT_TABLE
        if not self.table:
            raise SystemExit("FATAL: no frame table: set MTB_GRAPH_FRAMES or "
                             "MTB_BUILD_DIR (<build>/assets/graph_frame_offsets.tsv)")
        self.f = {}
        with open(self.table, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if not r["strand"]:
                    continue          # unresolved; treated as unknown below
                self.f[r["accession"]] = (r["strand"], int(r["offset"]),
                                          int(r["panel_len"]))

    # -- queries ----------------------------------------------------------
    def known(self, acc):
        return acc in self.f

    def flipped(self, acc):
        return self._get(acc)[0] == "-"

    def length(self, acc):
        return self._get(acc)[2]

    def identity(self, acc):
        """True when panel and refs are the same string, so no conversion is
        needed. 223 of 333 accessions."""
        st, off, _ = self._get(acc)
        return st == "+" and off == 0

    # -- conversion -------------------------------------------------------
    def to_panel(self, acc, refs_pos0):
        st, off, L = self._get(acc)
        return ((refs_pos0 + off) % L) if st == "+" else ((off - refs_pos0) % L)

    def to_refs(self, acc, panel_pos0):
        st, off, L = self._get(acc)
        return ((panel_pos0 - off) % L) if st == "+" else ((off - panel_pos0) % L)

    def to_refs_interval(self, acc, panel_lo, panel_hi):
        a = self.to_refs(acc, panel_lo)
        b = self.to_refs(acc, panel_hi)
        return (a, b) if a <= b else (b, a)

    def complement_needed(self, acc, odgi_strand):
        return self.flipped(acc) != (odgi_strand == "-")

    # -- odgi glue --------------------------------------------------------
    def query_line(self, path_name, acc, refs_pos0):
        """One line of an `odgi position -F` input file, in the panel frame.
        The strand stays `+`: a flipped accession is handled by the coordinate,
        because odgi's own strand flag does not mean the same thing."""
        return f"{path_name},{self.to_panel(acc, refs_pos0)},+\n"

    def _get(self, acc):
        try:
            return self.f[acc]
        except KeyError:
            raise KeyError(
                f"{acc} has no measured frame in {self.table}; "
                f"run graphframe/bin/graph_frame_offsets.py") from None


def main():
    """Print the conversion for one accession and position, for spot checks."""
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("accession")
    ap.add_argument("pos", type=int, help="0-based")
    ap.add_argument("--from-panel", action="store_true")
    ap.add_argument("--table", default=None)
    a = ap.parse_args()
    fr = Frames(a.table)
    st, off, L = fr.f[a.accession]
    q = fr.to_refs(a.accession, a.pos) if a.from_panel else fr.to_panel(a.accession, a.pos)
    print(f"{a.accession}\tstrand {st}\toffset {off}\tlen {L}\t"
          f"{'panel->refs' if a.from_panel else 'refs->panel'}\t{a.pos} -> {q}")


if __name__ == "__main__":
    main()
