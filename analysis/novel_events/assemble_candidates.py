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
Long-event guard (2026-10-10, as breakpoint_caller.py v5): a contig event
(DEL, INV, REARR or REPL) spanning more than LONG_EVENT (5 kb) of R is
kept only if it is a DEL with depth inside below 30% of the median. A
contig holding an inserted repeat element (IS6110) aligns to another copy
anywhere in R and otherwise gave DEL or INV calls of 0.2-3 Mb, replacing a
correct INS from the caller (N1177, 1,312). If nothing is left, the window
keeps the caller's calls, noted "long contig event dropped".
Contig-length guard (2026-10-11, user): a tandem duplication read off one
contig (the next alignment starts before the previous one ended on R) is
kept only if the duplicated length is at most the contig's length, since
the contig must hold the segment twice. A contig through an inserted repeat
copy aligns its copy part to another copy far back in R, which was read as
duplications of 77 kb-1.7 Mb at true IS6110 insertions (QC_1, 3,189,869).
Window selection (Phase B, real reads): a window is assembled only if it has
a typed candidate (DEL, INS, INV, REARR) or a one-sided cluster (BND) with
support of at least --bnd-min-frac of the sample's median depth; at most
--max-windows are assembled, those with the most support first. The rest keep
the prototype's calls, noted "not assembled (weak BND)" or "(cap)".
AF floor (2026-10-12, user; foldback_qc/recall/affloor.py): with --min-af, a
caller event with read support is kept only if support / reads overlapping
its start in the P2 BAM is at least --min-af. Assembly-typed events (no read
count) are always kept; dropping them halved recall. 0.1 is the production
setting: no combined recall cost in part D, precision in moderate fold-back
libraries 0.35 -> 0.62. Default 0 (off), as in all runs before 2026-10-12.
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
LONG_EVENT = 5000


def contig_events(paf, lo, hi):
    by = {}
    for line in open(paf):
        c = line.rstrip("\n").split("\t")
        cg = next((x[5:] for x in c[12:] if x.startswith("cg:Z:")), "")
        by.setdefault(c[0], []).append(dict(ql=int(c[1]), qs=int(c[2]), qe=int(c[3]), strand=c[4],
                                            ts=int(c[7]), te=int(c[8]), cg=cg))
    ev, n_guard = [], 0
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
                if -rg <= x["ql"]:  # contig-length guard
                    ev.append(("INS", st, st, -rg))  # tandem duplication
                else:
                    n_guard += 1
            elif rg - qgap >= MIN and qgap < MIN:
                ev.append(("DEL", st, max(st, en), rg))
            elif qgap - rg >= MIN and rg < MIN:
                ev.append(("INS", st, st, qgap))
            elif qgap >= MIN and rg >= MIN:
                ev.append(("REPL", st, max(st, en), max(rg, qgap)))
    return sorted(set(ev), key=lambda e: e[1]), n_guard


def af_floor(rows, bam, ctg, min_af):
    """Drop caller events with support / depth at start < min_af; keep events
    without a read count (assembly-typed, support empty or 0)."""
    if not min_af:
        return rows
    out = []
    for r in rows:
        sup = int(r["support"]) if str(r.get("support", "")).isdigit() else 0
        if sup:
            st = int(r["start"])
            dp = bam.count(ctg, max(0, st - 1), st)
            if dp and min(1.0, sup / dp) < min_af:
                continue
        out.append(r)
    return out


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
    ap.add_argument("--min-af", type=float, default=0.0,
                    help="AF floor on caller events with read support; 0 = off (production: 0.1)")
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
    def depth_ratio(s, e):
        pts = range(s + 50, e - 50, max(1, (e - s - 100) // 20))
        return statistics.mean(bam.count(ctg, x - 1, x) for x in pts) / med if med else 1.0

    def long_unsupported(f):
        t, s, e, n = f
        span = max(abs(e - s), n if t in ("DEL", "REPL") else 0)
        if t not in ("DEL", "INV", "REARR", "REPL") or span <= LONG_EVENT:
            return False
        return not (t == "DEL" and depth_ratio(s, e) < 0.3)

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
        dropped_long = False
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
            found, n_guard = contig_events(paf, lo - 100, hi + 100)
            n0 = len(found)
            found = [f for f in found if not long_unsupported(f)]
            if (n0 or n_guard) and not found:
                dropped_long = True
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
            note = ";long contig event dropped" if dropped_long else ";not assembled"
            rows += [dict(r, source=r["source"] + note) for r in w["ev"]]
    n_pre = len(rows)
    rows = af_floor(rows, bam, ctg, a.min_af)
    cols = ["sample", "type", "start", "end", "size", "support", "left", "right",
            "depth_ratio", "source"]
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    sk = collections.Counter(w["skip"] for w in wins)
    print(f"{a.sample}: {len(wins)} windows, {sk[None]} assembled, {sk['weak BND']} weak BND, "
          f"{sk['cap']} over the cap; {n_asm} gave events by assembly; "
          f"AF floor {a.min_af} dropped {n_pre - len(rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()
