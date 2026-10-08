#!/usr/bin/env python3
"""Phase A truth set: the structural differences between an assembly and
the matched reference R that P1 chose for its simulated reads.

The reads were simulated from the assembly with no mutations, so every
difference between the assembly and R is a real event the caller should
find in the P2 alignment to R.

Method:
1. minimap2 -cx asm5 --cs, R as target and the assembly as query.
2. Indels of 50 bp or more inside alignments, from paftools.js call
   (-L 10000 -l 1000).
3. Events between alignments: consecutive primary alignments in query
   order, on the same strand and in order on R:
   - R skips more than the assembly, by >= 50: DEL;
   - the assembly skips more than R, by >= 50: INS;
   - both skip by >= 50: REPL, a replacement;
   - a strand change: INV;
   - out of order on R: REARR.
   One jump wrapping from R's end to its start is the circular origin and
   is dropped.

Output, in R coordinates (1-based): sample, reference, type, r_start,
r_end, ref_len, alt_len, source (cs or break).
"""
import argparse
import csv
import gzip
import os
import subprocess
import sys
import tempfile

MIN = 50


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--refmap", required=True)
    ap.add_argument("--refs", required=True)
    ap.add_argument("--asm-dir", required=True)
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--k8", required=True)
    ap.add_argument("--paftools", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ref_of = {r["sample"]: r["reference"] for r in
              csv.DictReader(open(a.refmap), delimiter="\t")}
    os.makedirs(a.workdir, exist_ok=True)
    rows = []
    for g in csv.DictReader(open(a.genomes), delimiter="\t"):
        s, R = g["sample"], ref_of.get(g["sample"])
        if not R:
            print(f"{s}: no reference in refmap", file=sys.stderr)
            continue
        rfa = os.path.join(a.refs, f"{R}.fasta")
        paf = os.path.join(a.workdir, f"{s}.paf")
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as t:
            t.write(gzip.open(os.path.join(a.asm_dir, g["file"]), "rt").read())
        with open(paf, "w") as fo:
            subprocess.run([a.minimap2, "-cx", "asm5", "--cs", "-t", "2", rfa, t.name],
                           stdout=fo, stderr=subprocess.DEVNULL, check=True)
        os.remove(t.name)
        # indels inside alignments
        srt = subprocess.run(f"sort -k6,6 -k8,8n {paf}", shell=True,
                             capture_output=True, text=True, check=True).stdout
        call = subprocess.run([a.k8, a.paftools, "call", "-L", "10000", "-l", "1000", "-"],
                              input=srt, capture_output=True, text=True, check=True).stdout
        n_cs = 0
        for line in call.splitlines():
            c = line.split("\t")
            if c[0] != "V":
                continue
            ref = c[6] if c[6] != "-" else ""
            alt = c[7] if c[7] != "-" else ""
            if abs(len(ref) - len(alt)) < MIN and min(len(ref), len(alt)) < MIN:
                continue
            kind = "DEL" if len(ref) > len(alt) else "INS"
            st, en = int(c[2]) + 1, max(int(c[3]), int(c[2]) + 1)
            rows.append((s, R, kind, st, en, len(ref), len(alt), "cs"))
            n_cs += 1
        # events between alignments
        rlen = None
        alns = []
        for line in open(paf):
            c = line.split("\t")
            if "tp:A:P" not in line:
                continue
            rlen = int(c[6])
            alns.append(dict(qs=int(c[2]), qe=int(c[3]), strand=c[4],
                             ts=int(c[7]), te=int(c[8])))
        alns.sort(key=lambda x: x["qs"])
        n_br = 0
        for x, y in zip(alns, alns[1:]):
            qgap = y["qs"] - x["qe"]
            if x["strand"] != y["strand"]:
                kind, st, en = "INV", min(x["te"], y["ts"]) + 1, max(x["te"], y["ts"])
                rg = 0
            else:
                if x["strand"] == "+":
                    rg = y["ts"] - x["te"]
                    st, en = x["te"] + 1, y["ts"]
                else:
                    rg = x["ts"] - y["te"]
                    st, en = y["te"] + 1, x["ts"]
                if rlen and abs(rg) > rlen / 2:  # wrap at the circular origin
                    continue
                if rg < -MIN:
                    kind = "REARR"
                elif rg - qgap >= MIN and qgap < MIN:
                    kind = "DEL"
                elif qgap - rg >= MIN and rg < MIN:
                    kind = "INS"
                elif qgap >= MIN and rg >= MIN:
                    kind = "REPL"
                else:
                    continue
            rows.append((s, R, kind, st, max(st, en), max(rg, 0), max(qgap, 0), "break"))
            n_br += 1
        print(f"{s} vs {R}: {len(alns)} alignments, {n_cs} indels >= {MIN} bp, "
              f"{n_br} events between alignments", file=sys.stderr)
    with open(a.out, "w") as fo:
        fo.write("sample\treference\ttype\tr_start\tr_end\tref_len\talt_len\tsource\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
