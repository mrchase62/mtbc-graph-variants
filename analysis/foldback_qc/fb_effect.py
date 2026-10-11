#!/usr/bin/env python3
"""What fold-back chimeras do to SV calls (FOLDBACK_SV_EFFECT_PLAN.md, parts
A-C; report only). Every isolate here is mapped to its own genome, so every
SV call is false. Run from runroot with $MTB_PY.

  fb_effect.py --refmap <refmap> --p2 <p2 dir> --bam-dir <P2 work dir>
               --asm-dir <caller+assembly dir> --sites-dir <fold-back sites>
               --out <prefix>

Per call (callable sequence only: score.py's mask, span +-50 bp):
  caller   ours (caller v6 + assembly; typed calls and BND), dysgu PASS,
           delly PASS
  type, size
  AF       dysgu FORMAT AF; delly (DV+RV)/(DR+DV+RR+RV); ours: supporting
           reads / reads overlapping the start position in the P2 BAM
  junction distance from either breakpoint to the nearest fold-back
           junction of the same sample (foldback_v2.py --sites, all reads)
Chance level for "within --near bp (default 20) of a junction": the share of
the sample's callable sequence within that distance of a junction (f); for a call with two
breakpoints 1-(1-f)^2. Attributable = observed share - chance share.
Writes <prefix>.calls.tsv (one row per call) and <prefix>.summary.tsv
(sample x caller x AF bin counts, observed and chance near-junction
counts), and prints totals.
"""
import argparse
import bisect
import collections
import csv
import gzip
import os
import sys

import pysam

sys.path.insert(0, "../analysis/novel_events")
from score import load_masks, masked, read_proto  # noqa: E402

NEAR = 20  # --near; at 150 bp a third of a heavily chimeric genome is "near" by chance
BINS = ((0, 0.1, "<0.1"), (0.1, 0.2, "0.1-0.2"), (0.2, 0.5, "0.2-0.5"), (0.5, 9, ">=0.5"))


def af_bin(x):
    if x is None:
        return "na"
    return next(lab for lo, hi, lab in BINS if lo <= x < hi)


def vcf_calls(path, kind):
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[6] not in ("PASS", "."):
            continue
        info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in c[7].split(";"))
        fmt = dict(zip(c[8].split(":"), c[9].split(":"))) if len(c) > 9 else {}
        pos = int(c[1])
        end = int(info["END"]) if info.get("END", "").lstrip("-").isdigit() else pos
        if info.get("CHR2") and info["CHR2"] != c[0]:
            end = pos
        af = None
        try:
            if kind == "dysgu":
                af = float(fmt["AF"])
            else:
                alt = int(fmt["DV"]) + int(fmt["RV"])
                tot = alt + int(fmt["DR"]) + int(fmt["RR"])
                af = alt / tot if tot else None
        except (KeyError, ValueError):
            pass
        svlen = info.get("SVLEN", "").lstrip("-")
        out.append(dict(type=info.get("SVTYPE", ""), start=pos, end=max(pos, end),
                        size=int(svlen) if svlen.isdigit() else abs(end - pos), af=af))
    return out


def load_sites(path):
    pos = []
    if os.path.exists(path):
        with gzip.open(path, "rt") as fh:
            next(fh)
            pos = sorted(int(x.split("\t")[1]) for x in fh)
    return pos


def dist(sites, p):
    if not sites:
        return 10 ** 9
    i = bisect.bisect_left(sites, p)
    return min(abs(sites[j] - p) for j in (i - 1, i) if 0 <= j < len(sites))


def chance(sites, m, glen):
    """Share of callable sequence within NEAR bp of a junction."""
    cov, last = 0, -1
    iv = []
    for p in sites:
        s, e = max(0, p - NEAR), min(glen, p + NEAR + 1)
        if iv and s <= iv[-1][1]:
            iv[-1][1] = max(iv[-1][1], e)
        else:
            iv.append([s, e])
    starts, mk = m
    callable_len = glen - sum(e - s for s, e in mk)
    for s, e in iv:
        n = e - s
        for ms, me in mk[max(0, bisect.bisect_right(starts, e) - 8):bisect.bisect_right(starts, e)]:
            n -= max(0, min(e, me) - max(s, ms))
        cov += max(0, n)
    return cov / callable_len if callable_len else 0


def main():
    global NEAR
    ap = argparse.ArgumentParser()
    for k in ("refmap", "p2", "bam_dir", "asm_dir", "sites_dir", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    ap.add_argument("--mask-dir", default="../analysis/novel_events/out/masks")
    ap.add_argument("--ref-dir", default="refbias/build/7713a8d71d8e-fix1/refs")
    ap.add_argument("--near", type=int, default=NEAR)
    a = ap.parse_args()
    NEAR = a.near
    rm = {r["sample"]: r["reference"] for r in csv.DictReader(open(a.refmap), delimiter="\t")}
    masks = load_masks(a.refmap, a.mask_dir)
    rows, summ = [], []
    for s in sorted(rm):
        m = masks[s]
        fai = [x.split("\t") for x in open(f"{a.ref_dir}/{rm[s]}.fasta.fai")]
        glen = int(fai[0][1])
        sites = load_sites(f"{a.sites_dir}/{s}.tsv.gz")
        f = chance(sites, m, glen)
        calls = {"dysgu": vcf_calls(f"{a.p2}/{s}.dysgu.vcf", "dysgu"),
                 "delly": vcf_calls(f"{a.p2}/{s}.delly.vcf", "delly")}
        ours = []
        bam = pysam.AlignmentFile(f"{a.bam_dir}/{s}.bam")
        ev = list(csv.DictReader(open(f"{a.asm_dir}/{s}.events.tsv"), delimiter="\t")) \
            if os.path.exists(f"{a.asm_dir}/{s}.events.tsv") else []
        for c in ev:
            st, en = int(c["start"]), int(c["end"])
            dp = bam.count(fai[0][0], max(0, st - 1), st)
            sup = int(c["support"]) if c.get("support", "").isdigit() else None
            ours.append(dict(type=c["type"], start=st, end=en, size=int(c["size"] or 0),
                             af=min(1.0, sup / dp) if sup is not None and dp else None))
        calls["ours"] = ours
        cnt = collections.Counter()
        for caller, cs in calls.items():
            for c in cs:
                if masked(m, c["start"], c["end"]):
                    continue
                d = min(dist(sites, c["start"]), dist(sites, c["end"]))
                two = c["end"] - c["start"] > NEAR
                b = af_bin(c["af"])
                rows.append(dict(sample=s, caller=caller, type=c["type"], start=c["start"], end=c["end"],
                                 size=c["size"], af="" if c["af"] is None else round(c["af"], 3),
                                 af_bin=b, junction_dist=d, near=int(d <= NEAR)))
                for key in (b, "all"):
                    cnt[(caller, key, "calls")] += 1
                    cnt[(caller, key, "near")] += d <= NEAR
                    cnt[(caller, key, "chance")] += 1 - (1 - f) ** (2 if two else 1)
        for caller in ("ours", "dysgu", "delly"):
            for b in ("all", "<0.1", "0.1-0.2", "0.2-0.5", ">=0.5", "na"):
                summ.append(dict(sample=s, junctions=len(sites), chance_frac=round(f, 5), caller=caller,
                                 af_bin=b, calls=cnt[(caller, b, "calls")], near=cnt[(caller, b, "near")],
                                 near_chance=round(cnt[(caller, b, "chance")], 2)))
    for name, data in (("calls", rows), ("summary", summ)):
        with open(f"{a.out}.{name}.tsv", "w") as fo:
            w = csv.DictWriter(fo, fieldnames=list(data[0]), delimiter="\t")
            w.writeheader()
            w.writerows(data)
    print(f"{len(rm)} samples, {len(rows)} callable calls")


if __name__ == "__main__":
    main()
