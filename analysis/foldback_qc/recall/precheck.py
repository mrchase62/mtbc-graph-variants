#!/usr/bin/env python3
"""Why is precision low in part D's kept group? (report only). Run from
runroot with $MTB_PY.

False calls: callable typed calls (our caller + assembly; dysgu PASS AF >=
0.2) of the kept isolates (selfcal54, v2 <= 4.9) that match no truth event
under score.py's strict rule (matches() with a compatible type).

For each false call, the reference window around it (call +-500 bp, at most
6 kb) is aligned to the isolate's own genome (minimap2 -cx asm5 --cs). What
the two genomes show there:
  truth_gap   an indel >= 50 bp inside the window: a real difference that
              truth.py did not report (or reported elsewhere)
  indel_20_49 largest indel 20-49 bp: real, below truth's 50 bp cut-off
  no_sv       largest indel < 20 bp: no SV-scale difference, a caller false call
  complex     the window does not align as one piece covering >= 90% of it
              (split, inverted or missing in the own genome)
Also: SNP density in the window (divergence), distance to the nearest truth
event, and call type and size. Writes out/precheck.tsv and prints a summary.
"""
import collections
import csv
import os
import re
import subprocess
import sys
import tempfile

import pysam

sys.path.insert(0, "../analysis/novel_events")
from score import COMPAT, load_masks, masked, matches, read_proto, read_vcf  # noqa: E402

D = "../analysis/foldback_qc/recall"
RM = "refbias/fbrecall/run/p1/refmap.tsv"
REFS = "refbias/build/7713a8d71d8e-fix1/refs"
M2 = "/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin/minimap2"
CUT = 4.9
PAD, MAXW = 500, 6000


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def cs_stats(cs):
    snps = len(re.findall(r"\*[acgtn][acgtn]", cs))
    indels = [len(x) - 1 for x in re.findall(r"[+-][acgtn]+", cs)]
    return snps, max(indels, default=0)


def main():
    rm = {r["sample"]: r["reference"] for r in tsv(RM)}
    own = {r["sample"]: r["own_genome"] for r in tsv(f"{D}/genomes.tsv")}
    rate = {p["sample"]: float(p["rate"]) for p in tsv("../analysis/foldback_qc/selfcal/out/calib_points.tsv")}
    truth = collections.defaultdict(list)
    for r in tsv(f"{D}/out/truth.tsv"):
        truth[r["sample"]].append(dict(type=r["type"], start=int(r["r_start"]), end=int(r["r_end"]),
                                       size=max(int(r["ref_len"]), int(r["alt_len"]))))
    masks = load_masks(RM, "../analysis/novel_events/out/masks")
    kept = sorted(s for s in rm if s.startswith("ps_") and rate[s] <= CUT)
    rows = []
    for s in kept:
        calls = [("ours", c) for c in read_proto(f"{D}/out/asm/{s}.events.tsv")] + \
                [("dysgu_af20", c) for c in read_vcf(f"refbias/fbrecall/run/p2/{s}.dysgu.vcf", True, 0.2)]
        false = []
        for caller, c in calls:
            if c["type"] == "BND" or masked(masks[s], c["start"], c["end"]):
                continue
            if any(c["type"] in COMPAT[t["type"]] and matches(c, t) for t in truth[s]):
                continue
            false.append((caller, c))
        if not false:
            continue
        ref = pysam.FastaFile(f"{REFS}/{rm[s]}.fasta")
        ctg = ref.references[0]
        with tempfile.TemporaryDirectory() as td:
            fa = os.path.join(td, "w.fa")
            with open(fa, "w") as fo:
                for i, (caller, c) in enumerate(false):
                    a = max(0, min(c["start"], c["end"]) - PAD)
                    b = min(ref.get_reference_length(ctg), max(c["start"], c["end"]) + PAD)
                    if b - a > MAXW:
                        b = a + MAXW
                    fo.write(f">w{i}\n{ref.fetch(ctg, a, b)}\n")
            paf = subprocess.run([M2, "-cx", "asm5", "--cs", "-t", "2", f"{REFS}/{own[s]}.fasta", fa],
                                 capture_output=True, text=True, check=True).stdout
        alns = collections.defaultdict(list)
        for line in paf.splitlines():
            c = line.split("\t")
            cs = next((x[5:] for x in c[12:] if x.startswith("cs:Z:")), "")
            alns[c[0]].append(dict(qlen=int(c[1]), qs=int(c[2]), qe=int(c[3]),
                                   prim="tp:A:P" in line, cs=cs))
        for i, (caller, c) in enumerate(false):
            al = [x for x in alns.get(f"w{i}", []) if x["prim"]]
            dist = min((min(abs(c["start"] - t["start"]), abs(c["start"] - t["end"])) for t in truth[s]),
                       default=10 ** 9)
            if len(al) == 1 and (al[0]["qe"] - al[0]["qs"]) >= 0.9 * al[0]["qlen"]:
                snps, mx = cs_stats(al[0]["cs"])
                cls = "truth_gap" if mx >= 50 else "indel_20_49" if mx >= 20 else "no_sv"
                dens = round(1000 * snps / al[0]["qlen"], 2)
            else:
                cls, mx, dens = "complex", "", ""
            rows.append(dict(sample=s, reference=rm[s], caller=caller, type=c["type"], start=c["start"],
                             end=c["end"], size=c.get("size") or "", cls=cls, max_indel=mx,
                             snps_per_kb=dens, nearest_truth=dist))
    with open(f"{D}/out/precheck.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"kept isolates: {len(kept)}; false typed calls: {len(rows)}")
    for caller in ("ours", "dysgu_af20"):
        rr = [r for r in rows if r["caller"] == caller]
        print(f"\n{caller}: {len(rr)} false calls")
        cnt = collections.Counter(r["cls"] for r in rr)
        for k in ("truth_gap", "indel_20_49", "no_sv", "complex"):
            sub = [r for r in rr if r["cls"] == k]
            ty = collections.Counter(r["type"] for r in sub)
            near = sum(r["nearest_truth"] <= 500 for r in sub)
            print(f"  {k:12s} {cnt[k]:4d}  types {dict(ty)}  within 500 bp of a truth event {near}")
        nosv = [r for r in rr if r["cls"] == "no_sv"]
        if nosv:
            per = collections.Counter(r["sample"] for r in nosv)
            print(f"  no_sv per isolate: {sorted(per.values(), reverse=True)}")
    dd = [float(r["snps_per_kb"]) for r in rows if r["snps_per_kb"] != ""]
    if dd:
        dd.sort()
        print(f"\nSNP density in aligned windows (per kb): median {dd[len(dd)//2]}, 90th pct {dd[int(0.9*len(dd))]}")


if __name__ == "__main__":
    main()
