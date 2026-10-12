#!/usr/bin/env python3
"""Per-call SV rules shared by production (p4b_place_sv.py) and evaluation
(analysis/novel_events/score.py), so the two cannot drift apart. Moved out of
score.py unchanged on 2026-10-12 when our caller was wired into production.

A call is a dict with type (DEL, INS, REPL, INV, REARR, BND), start and end
(1-based, R coordinates).

tatC-helY genotype rule (user, 2026-10-10; element_genotype.py's table): a
DEL call overlapping the element site +-50 bp is dropped when the isolate's
genotype is present (reads inside the element). The low-fraction DEL there
is a short-read artifact (Nanopore: the element is in 839/843 spanning
reads, longread/lr_element.py). INS calls, and DEL calls in isolates whose
genotype is absent, are kept.

Small-indel hand-off (user, 2026-10-10; P2's HaplotypeCaller VCF): an SV
call is dropped when HaplotypeCaller has a PASS indel of 10-49 bp within 60
bp of its start that every read explains: depth >= 50% of the isolate's
median depth (median DP of its PASS SNPs) and alt reads >= 90% of AD. bwa
clips reads at a 20-30 bp indel, so the SV callers report it again as an
INS of unknown size; HaplotypeCaller calls it at its true size. A large
insertion near a small HC call is kept: reads cannot cross it, so HC's depth
there collapses (8129_04's 1,360 bp IS6110: DP 17 of ~100, AD 6,11).

Callable sequence (user, 2026-10-09; callable_mask.py's BED for R): a call
whose span +-50 bp touches the mask is uncallable.
"""
import bisect
import csv
import gzip
import os

TOL = 50


def load_element_gt(path):
    """element_genotype.py table -> sample: (1-based start, end, genotype)."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["site_start"] != "":
            out[r["sample"]] = (int(r["site_start"]) + 1, int(r["site_end"]), r["genotype"])
    return out


def element_artifact(gt, s, c):
    """A DEL call at the tatC-helY site in an isolate that carries the element."""
    x = gt.get(s) if gt else None
    return bool(x) and x[2] == "present" and c["type"] == "DEL" and \
        c["start"] <= x[1] + TOL and c["end"] >= x[0] - TOL


def load_hc_small(path):
    """HaplotypeCaller VCF -> [(pos, indel length, DP, alt fraction)] for PASS
    10-49 bp indels, and the median DP of PASS SNPs."""
    if not os.path.exists(path):
        return [], 0
    ind, snp_dp = [], []
    for line in gzip.open(path, "rt"):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[6] not in ("PASS", "."):
            continue
        fmt = dict(zip(c[8].split(":"), c[9].split(":")))
        try:
            dp = int(fmt.get("DP", "0"))
            ad = [int(x) for x in fmt.get("AD", "").split(",")]
        except ValueError:
            continue
        alts = c[4].split(",")
        if len(c[3]) == 1 and all(len(a) == 1 for a in alts):
            snp_dp.append(dp)
            continue
        for i, alt in enumerate(alts, 1):
            d = len(alt) - len(c[3])
            if 10 <= abs(d) < 50 and len(ad) > i and sum(ad):
                ind.append((int(c[1]), d, dp, ad[i] / sum(ad)))
    snp_dp.sort()
    return ind, (snp_dp[len(snp_dp) // 2] if snp_dp else 0)


def hc_small_indel(hc, c):
    """An SV call that is a HaplotypeCaller 10-49 bp indel explained by every read."""
    if not hc or c["type"] == "BND":
        return False
    ind, med = hc
    return any(abs(p - c["start"]) <= 60 and dp >= 0.5 * med and af >= 0.9 for p, _, dp, af in ind)


def load_bed(path):
    """A BED as masked()'s (sorted starts, sorted 0-based half-open intervals)."""
    iv = sorted((int(f[1]), int(f[2])) for f in (x.split("\t") for x in open(path)))
    return [x[0] for x in iv], iv


def masked(m, s, e):
    if m is None:
        return False
    starts, iv = m
    lo, hi = s - 1 - TOL, e + TOL  # 1-based event span +-50 bp, as 0-based half-open
    i = bisect.bisect_right(starts, hi)
    return any(iv[j][1] > lo for j in range(max(0, i - 8), i))
