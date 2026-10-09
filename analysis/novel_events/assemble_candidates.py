#!/usr/bin/env python3
"""Step 1 of the Phase A follow-up: local assembly of each prototype
candidate, to type and size it.

For one sample:
1. The prototype's events (breakpoint_caller.py) are grouped into windows:
   candidates within 500 bp share one window.
2. For each window, reads are collected from the P2 BAM (reads aligned to the
   matched reference R) within 1.5 kb, with their mates and unmapped mates,
   and assembled with SPAdes --isolate (local_assembly.collect_reads /
   write_fastq, as in the ABS-1 audit).
3. Contigs of 300 bp or more are aligned to R with minimap2 -cx asm5.
4. Events are read off the contigs with the truth-set rules (truth.py):
   - indels of 50 bp or more inside an alignment (from the CIGAR);
   - gaps between consecutive alignments of one contig: DEL, INS, REPL,
     INV or REARR.
   Only events within the window are kept.
5. A window whose contigs give events contributes those, typed and sized.
   A window whose contigs give none keeps the prototype's original calls,
   noted "not assembled".
Window selection (Phase B, real reads): a window is assembled only if it has
a typed candidate (DEL, INS, INV, REARR) or a one-sided cluster (BND) with
support of at least --bnd-min-frac of the sample's median depth; at most
--max-windows are assembled, those with the most support first. The rest keep
the prototype's calls, noted "not assembled (weak BND)" or "(cap)".
Output: an events table in breakpoint_caller.py's format, for score.py.
"""
import argparse
import collections
import csv
import os
import re
import shutil
import statistics
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "inherited_absent"))
from local_assembly import collect_reads, write_fastq  # noqa: E402

PAD = 1500
MIN = 50


def contig_events(paf, lo, hi):
    by = {}
    for line in open(paf):
        c = line.rstrip("\n").split("\t")
        cg = next((x[5:] for x in c[12:] if x.startswith("cg:Z:")), "")
        by.setdefault(c[0], []).append(dict(qs=int(c[2]), qe=int(c[3]), strand=c[4],
                                            ts=int(c[7]), te=int(c[8]), cg=cg))
    ev = []
    for name, al in by.items():
        for x in al:  # indels inside an alignment
            p = x["ts"]
            for n, op in re.findall(r"(\d+)([MIDNSHP=X])", x["cg"]):
                n = int(n)
                if op == "D" and n >= MIN and lo <= p <= hi:
                    ev.append(("DEL", p + 1, p + n, n))
                if op == "I" and n >= MIN and lo <= p <= hi:
                    ev.append(("INS", p, p, n))
                if op in "M=XDN":
                    p += n
        al.sort(key=lambda a: a["qs"])
        for x, y in zip(al, al[1:]):  # gaps between alignments
            qgap = y["qs"] - x["qe"]
            if x["strand"] != y["strand"]:
                st, en = min(x["te"], y["ts"]) + 1, max(x["te"], y["ts"])
                if lo <= st <= hi:
                    ev.append(("INV", st, en, abs(en - st)))
                continue
            if x["strand"] == "+":
                rg, st, en = y["ts"] - x["te"], x["te"] + 1, y["ts"]
            else:
                rg, st, en = x["ts"] - y["te"], y["te"] + 1, x["ts"]
            if not lo <= st <= hi:
                continue
            if rg < -MIN:
                ev.append(("INS", st, st, -rg))  # tandem duplication
            elif rg - qgap >= MIN and qgap < MIN:
                ev.append(("DEL", st, max(st, en), rg))
            elif qgap - rg >= MIN and rg < MIN:
                ev.append(("INS", st, st, qgap))
            elif qgap >= MIN and rg >= MIN:
                ev.append(("REPL", st, max(st, en), max(rg, qgap)))
    return sorted(set(ev), key=lambda e: e[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True, help="prototype events table")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True, help="matched reference FASTA")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--spades", required=True)
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--bnd-min-frac", type=float, default=0.0)
    ap.add_argument("--max-windows", type=int, default=0, help="0 = no cap")
    a = ap.parse_args()
    ev = list(csv.DictReader(open(a.events), delimiter="\t"))
    bam = pysam.AlignmentFile(a.bam)
    ctg = bam.references[0]
    os.makedirs(a.workdir, exist_ok=True)
    mmi = os.path.join(a.workdir, "ref.mmi")
    subprocess.run([a.minimap2, "-x", "asm5", "-d", mmi, a.ref],
                   stderr=subprocess.DEVNULL, check=True)
    # windows
    wins = []
    for e in sorted(ev, key=lambda x: int(x["start"])):
        s, t = int(e["start"]), int(e["end"])
        if wins and s - wins[-1]["hi"] <= 500:
            wins[-1]["hi"] = max(wins[-1]["hi"], t)
            wins[-1]["ev"].append(e)
        else:
            wins.append(dict(lo=s, hi=t, ev=[e]))
    L = bam.lengths[0]
    dep = [bam.count(ctg, x, x + 1) for x in range(1000, L - 1000, max(1, L // 2000))]
    med = statistics.median(dep) if dep else 0
    for w in wins:
        w["sup"] = max(int(e["support"] or 0) for e in w["ev"])
        w["skip"] = None
        if all(e["type"] == "BND" for e in w["ev"]) and w["sup"] < a.bnd_min_frac * med:
            w["skip"] = "weak BND"
    if a.max_windows:
        cand = sorted((w for w in wins if not w["skip"]), key=lambda w: -w["sup"])
        for w in cand[a.max_windows:]:
            w["skip"] = "cap"
    rows, n_asm = [], 0
    for i, w in enumerate(wins):
        lo, hi = w["lo"], w["hi"]
        if w["skip"]:
            rows += [dict(r, source=r["source"] + f";not assembled ({w['skip']})") for r in w["ev"]]
            continue
        if hi - lo > 20000:
            rows += [dict(r, source=r["source"] + ";not assembled (window > 20 kb)") for r in w["ev"]]
            continue
        d = os.path.join(a.workdir, f"w{i:04d}")
        os.makedirs(d, exist_ok=True)
        reads = collect_reads(bam, ctg, lo - PAD, hi + PAD)
        p1, p2, ps, npair, nsing = write_fastq(reads, d)
        sp = os.path.join(d, "spades")
        cmd = a.spades.split() + ["--isolate", "-k", "21,33,55,77", "-t", "1", "-m", "4", "-o", sp]
        if npair:
            cmd += ["-1", p1, "-2", p2]
        if nsing:
            cmd += ["-s", ps]
        r = subprocess.run(cmd, capture_output=True, text=True)
        found = []
        contigs = os.path.join(sp, "contigs.fasta")
        if r.returncode == 0 and os.path.exists(contigs):
            keep = os.path.join(d, "contigs.fasta")
            with open(keep, "w") as fo:
                for c in pysam.FastxFile(contigs):
                    if len(c.sequence) >= 300:
                        fo.write(f">{c.name}\n{c.sequence}\n")
            paf = os.path.join(d, "contigs.paf")
            with open(paf, "w") as fo:
                subprocess.run([a.minimap2, "-cx", "asm5", "--secondary=no", mmi, keep],
                               stdout=fo, stderr=subprocess.DEVNULL, check=True)
            found = contig_events(paf, lo - 100, hi + 100)
        if os.path.realpath(sp).startswith(os.path.realpath(a.workdir) + os.sep) and os.path.isdir(sp):
            shutil.rmtree(sp)
        for f in (p1, p2, ps):
            if os.path.exists(f):
                os.remove(f)
        if found:
            n_asm += 1
            for t, s, e, n in found:
                rows.append(dict(sample=a.sample, type=t, start=s, end=e, size=n, support="",
                                 left="", right="", depth_ratio="", source="assembled"))
        else:
            rows += [dict(r, source=r["source"] + ";not assembled") for r in w["ev"]]
    cols = ["sample", "type", "start", "end", "size", "support", "left", "right",
            "depth_ratio", "source"]
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    sk = collections.Counter(w["skip"] for w in wins)
    print(f"{a.sample}: {len(wins)} windows, {sk[None]} assembled, {sk['weak BND']} weak BND, "
          f"{sk['cap']} over the cap; {n_asm} gave events by assembly", file=sys.stderr)


if __name__ == "__main__":
    main()
