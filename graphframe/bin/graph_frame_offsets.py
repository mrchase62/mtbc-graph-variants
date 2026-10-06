#!/usr/bin/env python3
"""Measure the coordinate frame of every graph path against the refs/ FASTA of
the same accession, by sequence.

WHY THIS EXISTS. The pangenome graph was built from data/fastas/mtb.complex333
.fasta.gz, whose sequences are dnaA-rotated. Read alignment in this project uses
refbias/build/<id>/refs/<accession>.fasta, which is the sequence as deposited.
For most accessions the two are the same string, so a coordinate means the same
base in both; for a minority they are the same sequence written from a different
origin, and for at least one they are on opposite strands. A coordinate handed to
`odgi position` is interpreted in the PANEL frame, so for those accessions a
refs-frame coordinate silently queries the wrong base and the projection lands
tens to hundreds of kb away from the truth.

The measurement is a string search, not an alignment: take probes from the panel
sequence, find them in the refs sequence, and require every probe to give the
same offset. That leaves nothing to infer.

Output columns:
  accession, panel_len, refs_len, strand, offset, agree
where `strand` is + when the panel sequence reads the same way as refs and - when
it is the reverse complement, and `offset` converts between the two:

  strand +   panel_index = (refs_index + offset) mod L
  strand -   panel_index = (offset - refs_index) mod L , base complemented

`agree` is ok when every probe gave the same answer.
"""
import argparse, glob, gzip, os, re, subprocess, sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta_seq(path):
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                out.append(line.strip())
    return "".join(out).upper()


def panel_index(fai):
    names = {}
    for line in open(fai):
        f = line.split("\t")
        names[f[0].split("#")[0]] = (f[0], int(f[1]))
    return names


def probe(samtools, panel, name, start, n):
    r = subprocess.run([samtools, "faidx", panel, f"{name}:{start + 1}-{start + n}"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return ""
    return "".join(r.stdout.split("\n")[1:]).upper()


def main():
    ap = argparse.ArgumentParser()
    # no defaults: they were CX333's panel FASTA and build 7713a8d71d8e's
    # refs; P0 step frames passes the build's
    ap.add_argument("--panel", required=True,
                    help="the FASTA the graph was built from")
    ap.add_argument("--refs", required=True, help="<build>/refs")
    ap.add_argument("--samtools", default=os.environ.get("MTB_SAMTOOLS", "samtools"))
    ap.add_argument("--probe-len", type=int, default=200)
    ap.add_argument("--probes", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    idx = panel_index(a.panel + ".fai")
    rows = []
    for acc in sorted(idx):
        name, plen = idx[acc]
        rpath = os.path.join(a.refs, acc + ".fasta")
        if not os.path.exists(rpath):
            rows.append((acc, plen, "", "", "", "no refs fasta"))
            continue
        s = read_fasta_seq(rpath)
        seen = set()
        note = "ok"
        for k in range(a.probes):
            start = (plen // (a.probes + 1)) * (k + 1)
            p = probe(a.samtools, a.panel, name, start, a.probe_len)
            if len(p) < a.probe_len:
                note = "short probe"
                continue
            fwd = [m.start() for m in re.finditer("(?=" + p + ")", s)]
            rev = [m.start() for m in re.finditer("(?=" + rc(p) + ")", s)]
            if len(fwd) == 1 and not rev:
                seen.add(("+", (start - fwd[0]) % len(s)))
            elif len(rev) == 1 and not fwd:
                # panel[start .. start+n-1] == rc(refs[r .. r+n-1]) so
                # panel[start] pairs with refs[r+n-1]
                seen.add(("-", (start + rev[0] + a.probe_len - 1) % len(s)))
            else:
                note = f"{len(fwd)} fwd / {len(rev)} rc hits"
        if len(seen) == 1:
            st, off = seen.pop()
            rows.append((acc, plen, len(s), st, off, note))
        else:
            rows.append((acc, plen, len(s), "", "",
                         note if note != "ok" else f"{len(seen)} distinct answers"))

    with open(a.out, "w") as fh:
        fh.write("accession\tpanel_len\trefs_len\tstrand\toffset\tagree\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    same = sum(1 for r in rows if r[3] == "+" and r[4] == 0)
    rot = sum(1 for r in rows if r[3] == "+" and r[4] not in ("", 0))
    flip = sum(1 for r in rows if r[3] == "-")
    bad = sum(1 for r in rows if r[3] == "")
    print(f"{len(rows)} accessions: {same} same frame, {rot} rotated, "
          f"{flip} reverse complemented, {bad} unresolved")
    return 0


if __name__ == "__main__":
    sys.exit(main())
