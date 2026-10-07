#!/usr/bin/env python3
"""Convert `odgi position` input and output between the refs and panel frames.

A filter, so a call site changes by two pipes and nothing else has to learn
about frames. See graph_frame.py for why the two frames differ and
graphframe/docs/GRAPH_FRAME_RESOLUTION.md for what went wrong without this.

    to-panel    stdin: `path,pos0,strand` lines destined for `odgi position -F`
                stdout: the same lines with pos0 converted refs -> panel

    from-panel  stdin: `odgi position` output, either the `-r` form
                (source, target, dist.to.ref, strand.vs.ref) or the `-v` form
                (source, `node,offset,strand`)
                stdout: the same lines with every PATH coordinate converted
                panel -> refs, and strand.vs.ref restated as the relation
                between the two REFS sequences

Node ids and node offsets are not converted: a node is not in either frame.

STRAND. `to-panel` always writes `+`. A globally flipped accession is handled by
the coordinate, because odgi's strand flag describes steps inside the graph and
not the relation to the refs FASTA.

`from-panel` REPLACES column 4 with the relation between the two REFS
sequences, which is `-` when exactly one of the two accessions is stored
reverse complemented in the panel and `+` otherwise. odgi's own value is kept,
appended as a fifth column, so nothing is lost.

Column 4 does NOT incorporate odgi's flag, so a reader must apply column 5
itself. odgi does NOT already account for a locally inverted step. Where
odgi's own flag is `-` (source and target walk the node in opposite
directions) the reported target is one base past the homolog and the source
reads complemented relative to column 4. Measured by reading the source's
31-mer against the target's (bin/p4_place.py, parse_pos_file):

    col 4  col 5  homolog              GCF_000193185    scale200 P4
    +      -      t-1, complemented    8,684 / 8,694    167 / 228
    -      -      t-1, same strand          -              7 / 7
    +      +      t, same strand      34,298 / 34,465

bin/p4_place.py, bin/p4b_place_sv.py and is6110/bin/is6110_project_sites.py
shift and toggle there; bin/p5_states.py does the same for the H37Rv -> R
direction. The earlier reading of this note, that odgi accounts for the
inverted step and column 5 must be ignored, rested on 4 pilot SNPs (75% "as
is") and was wrong.

The `-v` form's node offset is likewise counted along the source path's
walking direction: one base of a node walked `+` by one path and `-` by
another has offsets off and L-1-off. Node keys are written with the node's
forward offset (bin/mtb_norm.py, forward_offset).

A PATH PROJECTED ONTO ITSELF IS THE IDENTITY (D21). Where the source and the
target are the same path, `from-panel` writes the source position as the
target, distance 0, odgi flag `+`. odgi does not: where the path passes a
node more than once it answers with one of the copies. With H37Rv as the
matched reference, scale200's pinned arm had 729 of 51,139 composed P4
records (all in PE/PPE tandem repeats, 552 at 3.93-3.95 Mb) moved by up to
1 kb, and P5 projects the other way through the same nodes.
"""
import argparse, importlib.util, os, sys

_s = importlib.util.spec_from_file_location(
    "graph_frame", os.path.join(os.path.dirname(os.path.abspath(__file__)), "graph_frame.py"))
gf = importlib.util.module_from_spec(_s); _s.loader.exec_module(gf)


def split_pathpos(tok):
    """`GCF_x#1#NZ_y,12345,+` -> (path, pos, strand); None if not that shape."""
    p = tok.rsplit(",", 2)
    if len(p) != 3 or "#" not in p[0]:
        return None
    try:
        return p[0], int(p[1]), p[2]
    except ValueError:
        return None


def to_panel(fr, fh_in, fh_out, strict):
    n = c = 0
    for line in fh_in:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            fh_out.write(line + "\n"); continue
        t = split_pathpos(line)
        if t is None:
            if strict:
                sys.exit(f"frame_convert: cannot parse {line!r}")
            fh_out.write(line + "\n"); continue
        path, pos, _ = t
        acc = path.split("#")[0]
        n += 1
        if not fr.known(acc):
            if strict:
                sys.exit(f"frame_convert: {acc} has no measured frame")
            fh_out.write(line + "\n"); continue
        if not fr.identity(acc):
            c += 1
        fh_out.write(f"{path},{fr.to_panel(acc, pos)},+\n")
    print(f"[frame_convert] to-panel: {n} positions, {c} converted", file=sys.stderr)


def from_panel(fr, fh_in, fh_out, strict):
    n = c = s = 0
    for line in fh_in:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            fh_out.write(line + "\n"); continue
        f = line.split("\t")
        flip = False
        same = (len(f) >= 3 and "#" in f[1].split(",", 1)[0]
                and f[0].rsplit(",", 2)[0] == f[1].rsplit(",", 2)[0])
        if same:
            # the identity, before conversion: one path, one frame
            f[1] = f[0].rsplit(",", 1)[0] + ",+"
            f[2] = "0"
            if len(f) >= 4 and f[3] in ("+", "-"):
                f[3] = "+"
            s += 1
        for i in (0, 1):
            if i >= len(f):
                break
            t = split_pathpos(f[i])
            if t is None:
                continue                      # `node,offset,strand`: not a frame
            path, pos, st = t
            acc = path.split("#")[0]
            if not fr.known(acc):
                if strict:
                    sys.exit(f"frame_convert: {acc} has no measured frame")
                continue
            if fr.flipped(acc):
                flip = not flip
            if not fr.identity(acc):
                c += 1
            f[i] = f"{path},{fr.to_refs(acc, pos)},{st}"
        if len(f) >= 4 and f[3] in ("+", "-"):
            f = f[:4] + [f[3]] + f[4:]        # keep odgi's own flag alongside
            f[3] = "-" if flip else "+"
        n += 1
        fh_out.write("\t".join(f) + "\n")
    print(f"[frame_convert] from-panel: {n} records, {c} coordinates converted"
          + (f", {s} projected onto their own path (identity)" if s else ""),
          file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("direction", choices=["to-panel", "from-panel"])
    ap.add_argument("--table", default=None)
    ap.add_argument("--lenient", action="store_true",
                    help="pass through what cannot be parsed or has no measured "
                         "frame, instead of failing. Off by default: a silently "
                         "unconverted coordinate is the whole problem.")
    a = ap.parse_args()
    fr = gf.Frames(a.table)
    (to_panel if a.direction == "to-panel" else from_panel)(
        fr, sys.stdin, sys.stdout, not a.lenient)


if __name__ == "__main__":
    main()
