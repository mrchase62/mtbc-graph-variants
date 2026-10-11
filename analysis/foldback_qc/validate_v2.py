#!/usr/bin/env python3
"""Validation of the fold-back scan v2 (report only; run from runroot with
$MTB_PY after 03_v2_validate.sbatch). Reads v2/rates.tsv and v2/sites/.

1. Synthetic (foldsim): true fold-back reads per 1,000 = junction-crossing
   reads (refbias/foldsim/junctions, r1_crosses + r2_crosses) over 2 x pairs
   (pairs = replaced pairs / fraction, conditions.tsv). Detected = the
   sample's rate minus the same genome's rate at 0. Sensitivity for v2
   (fold-back = split + clip) and for the old count.
2. Phase B (raw reads, marinB63): v2 vs old per group; Spearman between them.
3. False SV calls per isolate on the same raw reads, in callable sequence
   (score.py masks): dysgu PASS and our caller (v3b, after assembly), typed
   calls matching no true event (all calls in the in_panel controls, whose
   own genome is the reference... is not true for all; controls have no
   truth rows, so every callable call counts). Spearman of each fold-back
   metric with these counts.
4. Recurrent sites: within a reference shared by >= 3 Phase B isolates,
   junction positions (+-2 bp) carried by >= 3 of them; each isolate's v2
   rate is also given without them.
Writes v2/validation.tsv (per sample) and prints the summary.
"""
import collections
import csv
import gzip
import os
import sys

N = "../analysis/novel_events"
Q = "../analysis/foldback_qc"
sys.path.insert(0, N)
from score import COMPAT, load_masks, masked, matches, read_proto, read_vcf  # noqa: E402

RM = "refbias/marinB63/run/p1/refmap.tsv"


def spearman(x, y):
    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            for k in range(i, j + 1):
                r[o[k]] = (i + j) / 2
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def main():
    rates = {r["sample"]: r for r in csv.DictReader(open(f"{Q}/v2/rates.tsv"), delimiter="\t")}
    f = lambda r, k: float(r[k])  # noqa: E731
    # 1. synthetic
    cond = {r["sample"]: r for r in csv.DictReader(open(f"{Q}/../foldback_sim/conditions.tsv"), delimiter="\t")}
    print("1. Synthetic: true fold-back reads per 1,000 vs detected above the same genome at 0")
    print("   sample                true    v2    old   v2_sens old_sens")
    rows = []
    for s, c in sorted(cond.items()):
        if c["model"] == "none" or s not in rates:
            continue
        jt = f"refbias/foldsim/junctions/{s}.tsv"
        nj = cross = 0
        for r in csv.DictReader(open(jt), delimiter="\t"):
            nj += 1
            cross += int(r["r1_crosses"]) + int(r["r2_crosses"])
        pairs = nj / float(c["frac"])
        true = 1000 * cross / (2 * pairs)
        base = rates[s.rsplit("_", 2)[0] + "_0"] if s.rsplit("_", 2)[0] + "_0" in rates else None
        b2 = f(base, "foldback_per_1k") if base else 0
        bo = f(base, "old_per_1k") if base else 0
        v2 = f(rates[s], "foldback_per_1k") - b2
        old = f(rates[s], "old_per_1k") - bo
        rows.append((s, true, v2, old))
        print(f"   {s:20s} {true:6.2f} {v2:6.2f} {old:6.2f}   {v2 / true:5.2f}   {old / true:5.2f}")
    # 2-3. Phase B
    grp = {r["sample"]: r["group"] for r in csv.DictReader(open(f"{N}/phaseB/out/samples.tsv"), delimiter="\t")}
    masks = load_masks(RM, f"{N}/out/masks")
    truth = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{N}/phaseB/out/truth.tsv"), delimiter="\t"):
        truth[r["sample"]].append(dict(type=r["type"], start=int(r["r_start"]), end=int(r["r_end"]),
                                       size=max(int(r["ref_len"]), int(r["alt_len"]))))
    out = []
    for s, g in sorted(grp.items()):
        if s not in rates:
            continue
        m = masks.get(s)
        fc = {}
        for lab, cs in (("dysgu_pass", read_vcf(f"refbias/marinB63/run/p2/{s}.dysgu.vcf", True)),
                        ("ours", read_proto(f"{N}/phaseB/out/v3b/asm/{s}.events.tsv"))):
            cs = [c for c in cs if c["type"] != "BND" and not masked(m, c["start"], c["end"])]
            fc[lab] = sum(1 for c in cs if not any(c["type"] in COMPAT[t["type"]] and matches(c, t)
                                                     for t in truth.get(s, [])))
        out.append(dict(sample=s, group=g, v2=f(rates[s], "foldback_per_1k"), split=f(rates[s], "split_per_1k"),
                        clip=f(rates[s], "clip_per_1k"), clip_other=f(rates[s], "clip_other_per_1k"),
                        old=f(rates[s], "old_per_1k"), false_dysgu=fc["dysgu_pass"], false_ours=fc["ours"]))
    # 4. recurrent sites
    ref = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM), delimiter="\t")}
    by_ref = collections.defaultdict(list)
    for r in out:
        by_ref[ref[r["sample"]]].append(r["sample"])
    sites = {}
    for r in out:
        p = f"{Q}/v2/sites/{r['sample']}.tsv.gz"
        sites[r["sample"]] = collections.Counter()
        for x in csv.DictReader(gzip.open(p, "rt"), delimiter="\t"):
            sites[r["sample"]][int(x["pos"])] += int(x["reads"])
    for g, ss in by_ref.items():
        rec = set()
        if len(ss) >= 3:
            carriers = collections.Counter()
            for s in ss:
                for p in {q // 5 for q in sites[s]}:
                    carriers[p] += 1
            rec = {p for p, k in carriers.items() if k >= 3}
        for r in out:
            if r["sample"] in ss:
                n_all = sum(sites[r["sample"]].values())
                n_rec = sum(v for p, v in sites[r["sample"]].items() if p // 5 in rec)
                r["recurrent_frac"] = round(n_rec / n_all, 3) if n_all else 0
                r["v2_nonrecurrent"] = round(r["v2"] * (1 - (n_rec / n_all if n_all else 0)), 4)
                r["ref_group_size"] = len(ss)
    with open(f"{Q}/v2/validation.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(out[0]), delimiter="\t", restval="")
        w.writeheader()
        w.writerows(out)
    print("\n2. Phase B raw reads, per group: median v2 / old per 1,000 (min-max v2)")
    for g in ("novel", "setE_polished", "in_panel"):
        xs = [r for r in out if r["group"] == g]
        v = sorted(r["v2"] for r in xs)
        o = sorted(r["old"] for r in xs)
        print(f"   {g:14s} n={len(xs):2d}  v2 {v[len(v) // 2]:.3f}  old {o[len(o) // 2]:.3f}  (v2 {v[0]:.3f}-{v[-1]:.3f})")
    print(f"   Spearman v2 vs old, all 62: {spearman([r['v2'] for r in out], [r['old'] for r in out]):.3f}")
    print("\n3. Spearman with false SV calls in callable sequence (raw reads), all 62 isolates:")
    for k in ("v2", "old", "split", "clip", "clip_other"):
        print(f"   {k:10s} dysgu PASS {spearman([r[k] for r in out], [r['false_dysgu'] for r in out]):+.3f}"
              f"   ours {spearman([r[k] for r in out], [r['false_ours'] for r in out]):+.3f}")
    print("\n4. Recurrent junction sites (reference shared by >= 3 isolates):")
    for r in out:
        if r.get("ref_group_size", 0) >= 3:
            print(f"   {r['sample']:16s} v2 {r['v2']:.3f}  recurrent share {r['recurrent_frac']:.2f}  "
                  f"without {r['v2_nonrecurrent']:.3f}")


if __name__ == "__main__":
    main()
