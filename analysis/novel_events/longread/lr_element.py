#!/usr/bin/env python3
"""Do the setE isolates' Nanopore reads contain molecules lacking the
tatC-helY element (and the other repeat-family members)? Run from runroot with
the mtb_pangenome_qc python (pysam) after 02_align.sbatch.

Members are located in each isolate's matched reference R as in
phaseB/repeat_family.py (H37Rv flanks; tatC-helY and fadE22 by sequence).

A read counts when its alignment (primary or supplementary, MAPQ >= 20)
covers the member with >= 200 bp of reference on both sides. Per read:
  present  >= 70% of the member's bases aligned, and no deletion >= 30 bp
           overlapping it
  absent   a deletion within +-20 bp of the member's length overlapping it,
           or < 20% of the member's bases aligned
  other    anything else (partial; reported, not interpreted)
Nanopore errors are mostly short indels, so a member-sized deletion inside a
read anchored >= 200 bp on both sides is a molecule without the member.

Compared, per isolate, with the short-read deletion-junction rate at the same
member (phaseB/out/qc51/repeat_family.isolates.tsv, junction reads per 100x
flank depth; QC_10 is not in that table).

Writes longread/out/lr_element.tsv (isolate x member) and prints tatC-helY
and a per-member summary.
"""
import collections
import csv
import os
import re
import sys

import pysam

N = "../analysis/novel_events"
sys.path.insert(0, f"{N}/phaseB")
from repeat_family import EXTRA, FL, LOCI, REFS, locate, locate_seq, rc, read_fasta  # noqa: E402

BAMS = "/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_ont/bam"
RUNS = f"{N}/longread/runs.tsv"
SR = f"{N}/phaseB/out/qc51/repeat_family.isolates.tsv"
OUT = f"{N}/longread/out"
ANCHOR = 200


def members_for(h37, refseqs):
    ms = []
    for r in csv.DictReader(open(LOCI), delimiter="\t"):
        a = int(r["pos"]) - 1
        ms.append(dict(name=f"{r['pos']}_{r['gene']}", left=h37[a - FL:a], right=h37[a + 68:a + 68 + FL]))
    core = "GATTTTGAGGCGATTCTGCG"  # the fadE22 copy, as in repeat_family.py
    fade = None
    for sq in refseqs.values():
        hits = [m.start() for q in (core, rc(core)) for m in re.finditer(q, sq)]
        if len(hits) == 1:
            p = hits[0]
            fade = sq[p - 30:p + 38] if core in sq else rc(sq[p - 18:p + 50])
            break
    ms.append(dict(name="tatC-helY", seq=EXTRA["tatC-helY"]))
    if fade:
        ms.append(dict(name="fadE22", seq=fade))
    return ms


def classify(bam, ctg, a, b):
    ln = b - a
    n = collections.Counter()
    for r in bam.fetch(ctg, a, b):
        if r.is_secondary or r.is_unmapped or r.mapping_quality < 20:
            continue
        if r.reference_start > a - ANCHOR or r.reference_end < b + ANCHOR:
            continue
        p, aligned, dels = r.reference_start, 0, []
        for op, k in r.cigartuples:
            if op in (0, 7, 8):
                aligned += max(0, min(p + k, b) - max(p, a))
            if op == 2 and p < b and p + k > a:
                dels.append(k)
            if op in (0, 2, 3, 7, 8):
                p += k
        if any(abs(d - ln) <= 20 for d in dels) or aligned < 0.2 * ln:
            n["absent"] += 1
        elif aligned >= 0.7 * ln and not any(d >= 30 for d in dels):
            n["present"] += 1
        else:
            n["other"] += 1
    return n


def main():
    h37 = read_fasta(os.environ["MTB_H37RV"])
    runs = list(csv.DictReader(open(RUNS), delimiter="\t"))
    refseqs = {r["reference"]: read_fasta(f"{REFS}/{r['reference']}.fasta") for r in runs}
    members = members_for(h37, refseqs)
    sr = {(r["sample"], r["member"]): r for r in csv.DictReader(open(SR), delimiter="\t") if r["found"] == "1"}
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for run in runs:
        s, g = run["isolate"], run["reference"]
        path = f"{BAMS}/{s}.bam"
        if not os.path.exists(path + ".bai"):
            print(f"{s}: no BAM yet", file=sys.stderr)
            continue
        bam = pysam.AlignmentFile(path)
        ctg = bam.references[0]
        for m in members:
            iv = locate_seq(refseqs[g], m["seq"]) if "seq" in m else locate(refseqs[g], m["left"], m["right"])
            if iv is None or iv[1] - iv[0] < 30:
                continue
            a, b = iv
            n = classify(bam, ctg, a, b)
            tot = sum(n.values())
            x = sr.get((s, m["name"]))
            jr = ""
            if x and float(x["flank_depth"] or 0):
                jr = round(100 * int(x["junction_reads"]) / float(x["flank_depth"]), 2)
            rows.append(dict(isolate=s, reference=g, member=m["name"], r_start=a + 1, r_end=b, length=b - a,
                             spanning_reads=tot, present=n["present"], absent=n["absent"], other=n["other"],
                             absent_frac=round(n["absent"] / tot, 4) if tot else "",
                             sr_junction_per_100x=jr, sr_dysgu_af=x["dysgu_af"] if x else ""))
    with open(f"{OUT}/lr_element.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("tatC-helY per isolate: spanning, present, absent, other, absent_frac | short-read junction/100x, dysgu AF")
    for r in rows:
        if r["member"] == "tatC-helY":
            print(f"  {r['isolate']:12s} {r['spanning_reads']:5d} {r['present']:5d} {r['absent']:4d} {r['other']:4d} "
                  f"{r['absent_frac']!s:7s} | {r['sr_junction_per_100x']!s:6s} {r['sr_dysgu_af']}")
    print("per member, all isolates: spanning, absent, isolates with any absent read")
    agg = collections.defaultdict(lambda: [0, 0, 0])
    for r in rows:
        v = agg[r["member"]]
        v[0] += r["spanning_reads"]
        v[1] += r["absent"]
        v[2] += r["absent"] > 0
    for k, v in sorted(agg.items(), key=lambda x: -x[1][1] / max(1, x[1][0])):
        print(f"  {k:20s} {v[0]:6d} {v[1]:5d} {v[2]:3d}  {v[1] / max(1, v[0]):.4f}")


if __name__ == "__main__":
    main()
