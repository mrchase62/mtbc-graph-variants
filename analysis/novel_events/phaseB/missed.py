#!/usr/bin/env python3
"""Why are true events missed? Report only, existing outputs (run from
runroot with $MTB_PY after rerun_compare.py).

Isolates and rules as rerun_compare.py "after": marinQC51 clean reads, P1
reference unchanged, phaseB/exclude.tsv removed; calls are our caller +
local assembly (proto_asm) and dysgu PASS AF >= 0.2, after the tatC-helY
genotype rule and the HaplotypeCaller hand-off. A miss is a callable true
event (novel or setE_polished) with no compatible matching call (score.py's
typed recall).

Cause, first that applies:
  sized_or_placed  a compatible-type call within 500 bp failed score.py's
                   position / size rules
  wrong_type       a call within 50 bp of the start or end, incompatible type
  bnd_only         our caller has an untyped breakpoint (BND) within 50 bp
  dysgu_filtered   a dysgu record matches (compatible type) but is not PASS
                   or has AF < 0.2
  dropped_by_rule  a matching call removed by the tatC-helY or HC rule
  cluster_rejected our caller's clusters table has a clip cluster within
                   50 bp that it did not keep
  no_call          none of these
For every miss, clipped reads (clip >= 10 bp, within 10 bp of either
breakpoint) are counted at MAPQ >= 20 and MAPQ < 20 in the clean P2 BAM:
evidence only at low MAPQ means the reads cannot be placed (repeat).
Writes out/qc51/missed.tsv and prints a summary.
"""
import collections
import csv
import os
import sys

import pysam

N = "../analysis/novel_events"
sys.path.insert(0, N)
import score  # noqa: E402

Q = f"{N}/phaseB/out/qc51"
P2 = "refbias/marinQC51/run/p2"
BAM = "refbias/work/marinQC51_p2"
RM = "refbias/marinQC51/run/p1/refmap.tsv"
RB = "refbias/marinB63/run/p1/refmap.tsv"


def dysgu_all(path):
    out = []
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        fmt = dict(zip(c[8].split(":"), c[9].split(":"))) if len(c) > 9 else {}
        try:
            af = float(fmt.get("AF", "nan"))
        except ValueError:
            af = float("nan")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in c[7].split(";"))
        t = {"DUP": "INS", "BND": "REARR", "TRA": "REARR"}.get(info.get("SVTYPE", ""), info.get("SVTYPE", ""))
        pos = int(c[1])
        end = int(info["END"]) if info.get("END", "").lstrip("-").isdigit() else pos
        sv = info.get("SVLEN", "").lstrip("-")
        size = int(sv) if sv.isdigit() and int(sv) > 0 else (max(pos, end) - pos + 1 if t in ("DEL", "INV") and end > pos else None)
        out.append(dict(type=t, start=pos, end=max(pos, end), size=size, filter=c[6], af=af))
    return out


def clips(bam, ctg, t):
    hi = lo = 0
    for q in {t["start"], t["end"]}:
        for a in bam.fetch(ctg, max(0, q - 10), q + 10):
            if a.is_secondary or a.is_supplementary or a.is_unmapped:
                continue
            ct = a.cigartuples
            if (ct[0][0] == 4 and ct[0][1] >= 10 and abs(a.reference_start + 1 - q) <= 10) or \
                    (ct[-1][0] == 4 and ct[-1][1] >= 10 and abs(a.reference_end - q) <= 10):
                if a.mapping_quality >= 20:
                    hi += 1
                else:
                    lo += 1
    return hi, lo


def main():
    ref = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM), delimiter="\t")}
    refb = {r["sample"]: r["reference"] for r in csv.DictReader(open(RB), delimiter="\t")}
    excl = {r["sample"] for r in csv.DictReader(open(f"{N}/phaseB/exclude.tsv"), delimiter="\t")}
    groups = {r["sample"]: r["group"] for r in csv.DictReader(open(f"{N}/phaseB/out/samples.tsv"), delimiter="\t")}
    keep = [s for s in ref if ref[s] == refb.get(s) and s not in excl and groups.get(s) in ("novel", "setE_polished")]
    masks = score.load_masks(RM, f"{N}/out/masks")
    egt = score.load_element_gt(f"{N}/phaseB/out/element_gt.marinQC51.tsv")
    rows = []
    for g in ("novel", "setE_polished"):
        for r in csv.DictReader(open(f"{Q}/truth.{g}.tsv"), delimiter="\t"):
            s = r["sample"]
            if s not in keep:
                continue
            t = dict(type=r["type"], start=int(r["r_start"]), end=int(r["r_end"]),
                     size=max(int(r["ref_len"]), int(r["alt_len"])), ref_len=int(r["ref_len"]), alt_len=int(r["alt_len"]))
            if score.masked(masks[s], t["start"], t["end"]):
                continue
            rows.append((g, s, t))
    out = []
    cache = {}
    for g, s, t in rows:
        if s not in cache:
            hc = score.load_hc_small(f"{P2}/{s}.vcf.gz")
            raw_ours = score.read_proto(f"{Q}/asm/{s}.events.tsv")
            raw_dys = dysgu_all(f"{P2}/{s}.dysgu.vcf")
            keep_c = lambda c: not score.masked(masks[s], c["start"], c["end"]) and \
                not score.element_artifact(egt, s, c) and not score.hc_small_indel(hc, c)
            ours = [c for c in raw_ours if keep_c(c)]
            dys = [c for c in raw_dys if c["filter"] in ("PASS", ".") and c["af"] >= 0.2 and keep_c(c)]
            cl = list(csv.DictReader(open(f"{Q}/proto/{s}.events.tsv.clusters.tsv"), delimiter="\t"))
            cache[s] = (raw_ours, raw_dys, ours, dys, cl, pysam.AlignmentFile(f"{BAM}/{s}.bam"))
        raw_ours, raw_dys, ours, dys, cl, bam = cache[s]
        calls = ours + dys
        comp = score.COMPAT[t["type"]]
        if any(c["type"] in comp and score.matches(c, t) for c in calls):
            continue
        near = lambda c: abs(c["start"] - t["start"]) <= 50 or abs(c["start"] - t["end"]) <= 50
        if any(c["type"] in comp and abs(c["start"] - t["start"]) <= 500 for c in calls if c["type"] != "BND"):
            cause = "sized_or_placed"
        elif any(near(c) and c["type"] not in comp for c in calls if c["type"] != "BND"):
            cause = "wrong_type"
        elif any(near(c) for c in ours if c["type"] == "BND"):
            cause = "bnd_only"
        elif any(c["type"] in comp and score.matches(c, t) for c in raw_dys):
            cause = "dysgu_filtered"
        elif any(c["type"] in comp and score.matches(c, t) for c in raw_ours):
            cause = "dropped_by_rule"
        elif any(abs(int(x["pos"]) - q) <= 50 and x["kept"] == "0" for x in cl for q in (t["start"], t["end"])):
            cause = "cluster_rejected"
        else:
            cause = "no_call"
        hi, lo = clips(bam, bam.references[0], t)
        d = next((c for c in raw_dys if c["type"] in comp and score.matches(c, t)), None)
        out.append(dict(group=g, sample=s, type=t["type"], start=t["start"], end=t["end"], ref_len=t["ref_len"],
                        alt_len=t["alt_len"], is6110_sized=int(1340 <= t["alt_len"] <= 1370 or 1340 <= t["ref_len"] <= 1370),
                        cause=cause, clipped_mapq20=hi, clipped_lowmapq=lo,
                        dysgu=f"{d['filter']} AF={d['af']}" if d else ""))
    with open(f"{Q}/missed.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    tot = collections.Counter(g for g, _, _ in rows)
    for g in ("novel", "setE_polished"):
        o = [r for r in out if r["group"] == g]
        print(f"\n{g}: callable true events {tot[g]}, missed {len(o)}")
        for c, n in collections.Counter(r["cause"] for r in o).most_common():
            x = [r for r in o if r["cause"] == c]
            print(f"  {c:17s} {n:3d}  types {dict(collections.Counter(r['type'] for r in x))}  "
                  f"IS6110-sized {sum(r['is6110_sized'] for r in x)}")


if __name__ == "__main__":
    main()
