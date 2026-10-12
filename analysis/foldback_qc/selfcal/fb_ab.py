#!/usr/bin/env python3
"""Parts A and B of ../FOLDBACK_SV_EFFECT_PLAN.md on selfcal54 (report only).
Run from runroot with $MTB_PY after ../fb_effect.py wrote out/effect.*.tsv.

Points: selfcal54 QC PASS isolates (out/effect.*) and the 9 chimeric controls,
all reads (../removal/out/effect.*, samples <isolate>_all). Every call is
false (own genome); callable sequence only.

Bands by v2 fold-back rate (calib_points.tsv): < 1, 1-4.9 (kept by the
cut-off), > 4.9 (excluded), and the chimeric controls.

A. Attribution: a call is "near" when a breakpoint lies within 20 bp of a
   fold-back junction. Chance: the per-call probability from fb_effect.py
   (share of callable sequence within 20 bp of a junction, 1-(1-f)^2 for two
   breakpoints). Attributable share = (near - expected) / (calls - expected),
   the fraction of calls explained beyond chance.
B. Dose-response:
   1. calls per isolate by AF bin and band;
   2. per band, the lowest AF floor at which the band's mean per-isolate
      false calls (ours typed + dysgu) falls to the clean level, the < 1 band's
      mean + 2 SD at the same floor;
   3. negative binomial regression, false calls (ours typed + dysgu AF >= 0.2)
      ~ log10 rate + library type + median insert + read length + depth,
      selfcal isolates only (library metadata); Poisson with robust SE if the
      NB fit fails.
"""
import collections
import csv
import math
import statistics as st

import numpy as np

S = "../analysis/foldback_qc/selfcal/out"
R = "../analysis/foldback_qc/removal/out"
CUT = 4.9
FLOORS = (0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def main():
    pts = {p["sample"]: p for p in tsv(f"{S}/calib_points.tsv")}
    calls = tsv(f"{S}/effect.calls.tsv") + [
        dict(r, sample=r["sample"][:-4]) for r in tsv(f"{R}/effect.calls.tsv") if r["sample"].endswith("_all")]
    summ = tsv(f"{S}/effect.summary.tsv") + [
        dict(r, sample=r["sample"][:-4]) for r in tsv(f"{R}/effect.summary.tsv") if r["sample"].endswith("_all")]
    chance = {r["sample"]: float(r["chance_frac"]) for r in summ}

    def band(s):
        p = pts[s]
        if p["group"] == "chimeric":
            return "chimeric"
        r = float(p["rate"])
        return "<1" if r < 1 else "1-4.9" if r <= CUT else ">4.9"

    samples = [s for s in pts if pts[s]["group"] in ("selfcal", "chimeric") and s in chance]
    bands = ("<1", "1-4.9", ">4.9", "chimeric")
    nb = collections.Counter(band(s) for s in samples)

    def kind(c):
        if c["caller"] == "ours":
            return "ours BND" if c["type"] == "BND" else "ours typed"
        return c["caller"]

    print("A. Attribution to fold-back junctions (within 20 bp; chance from callable share)")
    print("band\tisolates\tcaller\tcalls\tnear\texpected\tattributable\tchance share")
    for b in bands:
        for k in ("ours typed", "ours BND", "dysgu"):
            cs = [c for c in calls if c["sample"] in samples and band(c["sample"]) == b and kind(c) == k]
            if not cs:
                continue
            exp = 0.0
            for c in cs:
                f = chance[c["sample"]]
                exp += 1 - (1 - f) ** 2 if c["end"] != c["start"] else f
            near = sum(c["near"] == "1" for c in cs)
            att = (near - exp) / (len(cs) - exp) if len(cs) > exp else float("nan")
            print(f"{b}\t{nb[b]}\t{k}\t{len(cs)}\t{near}\t{exp:.1f}\t{att:.2f}\t{exp / len(cs):.2f}")
    print("\nA2. Attributable share by type and AF bin, all isolates with rate >= 1 (selfcal) pooled")
    print("caller\ttype\taf_bin\tcalls\tnear\texpected\tattributable")
    grp = collections.defaultdict(list)
    for c in calls:
        if c["sample"] in samples and band(c["sample"]) in ("1-4.9", ">4.9"):
            grp[(kind(c), c["type"], c["af_bin"])].append(c)
    for key in sorted(grp):
        cs = grp[key]
        if len(cs) < 5:
            continue
        exp = sum((1 - (1 - chance[c["sample"]]) ** 2) if c["end"] != c["start"] else chance[c["sample"]] for c in cs)
        near = sum(c["near"] == "1" for c in cs)
        att = (near - exp) / (len(cs) - exp) if len(cs) > exp else float("nan")
        print(f"{key[0]}\t{key[1]}\t{key[2]}\t{len(cs)}\t{near}\t{exp:.1f}\t{att:.2f}")

    print("\nB1. Mean false calls per isolate by AF bin (ours typed / dysgu PASS)")
    bins = ("<0.1", "0.1-0.2", "0.2-0.5", ">=0.5")
    print("band\tisolates\t" + "\t".join(f"ours {x}" for x in bins) + "\t" + "\t".join(f"dysgu {x}" for x in bins))
    per = {s: collections.Counter() for s in samples}
    for c in calls:
        if c["sample"] in per and kind(c) in ("ours typed", "dysgu"):
            per[c["sample"]][(kind(c), c["af_bin"])] += 1
    for b in bands:
        ss = [s for s in samples if band(s) == b]
        v = [st.mean(per[s][("ours typed", x)] for s in ss) for x in bins] + \
            [st.mean(per[s][("dysgu", x)] for s in ss) for x in bins]
        print(f"{b}\t{len(ss)}\t" + "\t".join(f"{x:.1f}" for x in v))

    print("\nB2. Mean false calls per isolate (ours typed + dysgu) with AF >= floor; clean = < 1 band mean + 2 SD")
    print("floor\tclean level\t" + "\t".join(bands[1:]))
    afs = {s: [float(c["af"]) for c in calls if c["sample"] == s and kind(c) in ("ours typed", "dysgu") and c["af"] != ""]
           for s in samples}
    reach = {}
    for f in FLOORS:
        cnt = {s: sum(a >= f for a in afs[s]) for s in samples}
        cl = [cnt[s] for s in samples if band(s) == "<1"]
        lvl = st.mean(cl) + 2 * st.stdev(cl)
        row = []
        for b in bands[1:]:
            m = st.mean(cnt[s] for s in samples if band(s) == b)
            row.append(f"{m:.1f}{'*' if m <= lvl else ''}")
            if m <= lvl and b not in reach:
                reach[b] = f
        print(f"{f}\t{lvl:.1f}\t" + "\t".join(row))
    print("(* at or below the clean level)")
    for b in bands[1:]:
        print(f"  {b}: lowest floor reaching the clean level: {reach.get(b, 'none up to 0.5')}")

    print("\nB3. Negative binomial: y (ours typed + dysgu AF >= 0.2) ~ log10 rate + library + insert + read length + depth (selfcal)")
    import statsmodels.api as sm
    sc = [pts[s] for s in samples if pts[s]["group"] == "selfcal"]
    lib = [p["adapter"] if p["adapter"] in ("truseq", "nextera_tn5") else "unclear" for p in sc]
    X = np.array([[math.log10(float(p["rate"]) + 0.01), l == "nextera_tn5", l == "unclear",
                   float(p["ins_med"]) / 100, float(p["read_len"]) / 100, float(p["depth"]) / 100]
                  for p, l in zip(sc, lib)], dtype=float)
    X = sm.add_constant(X)
    y = np.array([int(p["y"]) for p in sc])
    names = ["const", "log10 rate", "Tn5 vs TruSeq", "unclear vs TruSeq", "insert/100", "read len/100", "depth/100"]
    try:
        m = sm.NegativeBinomial(y, X).fit(disp=0, maxiter=500)
        if not m.mle_retvals.get("converged", False) or not np.all(np.isfinite(m.bse)):
            raise RuntimeError("NB did not converge")
        lab = "negative binomial"
        params, bse, pv = m.params[:-1], m.bse[:-1], m.pvalues[:-1]
        print(f"  dispersion alpha {m.params[-1]:.2f}")
    except Exception as e:  # noqa: BLE001
        print(f"  NB failed ({e}); Poisson GLM with robust (HC0) SE")
        m = sm.GLM(y, X, family=sm.families.Poisson()).fit(cov_type="HC0")
        lab = "Poisson (robust SE)"
        params, bse, pv = m.params, m.bse, m.pvalues
    print(f"  {lab}, n {len(y)}; rate ratios (exp coefficient)")
    for n, b, se, p in zip(names, params, bse, pv):
        print(f"  {n:18s} RR {math.exp(b):7.2f}  (95% CI {math.exp(b - 1.96 * se):.2f}-{math.exp(b + 1.96 * se):.2f})  p {p:.2g}")


if __name__ == "__main__":
    main()
