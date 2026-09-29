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

That column 4 does NOT incorporate odgi's flag is a measured decision, not an
oversight. Scoring P4's projected SNPs against H37Rv's own base over the
23-isolate pilot:

    cause of `-`            n     allele matches as-is    complemented
    accession stored rc    31              0.0%              100.0%
    odgi's own local `-`    4             75.0%                0.0%
    neither (`+`)         547            100.0%                0.0%

So `odgi position` already accounts for a locally inverted step when it reports
the target position, and the only thing that still needs complementing is the
storage flip this file knows about. Folding odgi's flag in as well double-counts
it and corrupts the allele at exactly the inverted sites.
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
    n = c = 0
    for line in fh_in:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            fh_out.write(line + "\n"); continue
        f = line.split("\t")
        flip = False
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
    print(f"[frame_convert] from-panel: {n} records, {c} coordinates converted",
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
