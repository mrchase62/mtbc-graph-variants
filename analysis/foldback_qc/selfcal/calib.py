#!/usr/bin/env python3
"""Fold-back cut-off calibration: rate alone vs rate plus insert size
(CALIBRATION_PLAN.md, "Added 2026-10-11 (user): rate plus insert size";
report only). Run from runroot with $MTB_PY.

Points (all QC'd reads mapped to their own genome, so every call is false):
  selfcal54 isolates with QC PASS; the 11 clean Phase B controls (both from
  out/selfcal_table.tsv); the 9 chimeric controls (../removal, <isolate>_all).
Outcome y: false calls in callable sequence, our typed calls + dysgu AF >= 0.2.
Predictors: fold-back v2 rate per 1,000; median insert and the share of pairs
  with insert < 2x read length (samtools stats IS histogram).
Threshold T = mean + 2 SD of y over the 11 clean controls. A rule's cut-offs
  are the ones keeping the most points while the kept points average y <= T.
  Rule 1: rate <= c.  Rule 2: rate <= c, or rate <= c2 when median insert >= i.
Leave-one-study-out: fit on the other studies, apply to the held-out one.
Also reported, not pre-specified: junctions with >= 5 fold-back reads
  (recurrence), Spearman correlations, a negative binomial fit (descriptive).
Writes out/calib_points.tsv and prints the report.
"""
import csv
import gzip
import math
import statistics as st
import sys

import numpy as np
from scipy.stats import spearmanr

S = "../analysis/foldback_qc/selfcal/out"
R = "../analysis/foldback_qc/removal"


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def insert_stats(path):
    """median insert, share of pairs with insert < 2 x mean read length"""
    hist, rl = [], None
    for line in open(path):
        if line.startswith("SN\taverage length:"):
            rl = float(line.split("\t")[2])
        elif line.startswith("IS\t"):
            c = line.split("\t")
            hist.append((int(c[1]), int(c[2])))
    tot = sum(n for _, n in hist)
    acc, med = 0, None
    for x, n in hist:
        acc += n
        if med is None and acc >= tot / 2:
            med = x
    short = sum(n for x, n in hist if x < 2 * rl) / tot
    return med, short, rl


def recurrent(path, k=5):
    with gzip.open(path, "rt") as fh:
        return sum(1 for r in csv.DictReader(fh, delimiter="\t") if int(r["reads"]) >= k)


def load():
    pts = []
    for r in tsv(f"{S}/selfcal_table.tsv"):
        if r["set"] == "selfcal54" and r["qc_status"] != "PASS":
            continue
        grp = "clean" if r["set"] == "phaseB_control" else "selfcal"
        med, short, rl = insert_stats(f"{S}/stats/{r['sample']}.full.txt")
        pts.append(dict(sample=r["sample"], group=grp,
                        study="marin2024" if grp == "clean" else r["study"],
                        instrument=r["instrument"], adapter=r["adapter_call"],
                        rate=float(r["foldback_per_1k"]), ins_med=med, ins_short=short,
                        read_len=rl, depth=float(r["final_depth"] or 0),
                        recur=recurrent(f"{S}/sites/{r['sample']}.tsv.gz"),
                        ours=int(r["noise_sv_typed"]), dysgu=int(r["noise_dysgu_af20"])))
    noise = {r["sample"]: r for r in tsv(f"{R}/out/self_noise.tsv")}
    for s in (x.split()[0] for x in open(f"{R}/isolates.txt") if x.strip()):
        fb = tsv(f"{R}/out/foldback/{s}_all.tsv")[-1]
        med, short, rl = insert_stats(f"{R}/out/stats/{s}_all.full.txt")
        n = noise[f"{s}_all"]
        pts.append(dict(sample=s, group="chimeric", study="marin2024", instrument="",
                        adapter="", rate=float(fb["foldback_per_1k"]), ins_med=med,
                        ins_short=short, read_len=rl, depth=100.0,
                        recur=recurrent(f"{R}/out/sites/{s}_all.tsv.gz"),
                        ours=int(n["sv_typed"]), dysgu=int(n["dysgu_af20"])))
    for p in pts:
        p["y"] = p["ours"] + p["dysgu"]
    return pts


def fit1(pts, T):
    best = None
    for c in sorted({p["rate"] for p in pts}):
        kept = [p for p in pts if p["rate"] <= c]
        if st.mean(p["y"] for p in kept) <= T:
            best = (c, len(kept))
    return best


def keep2(p, c, c2, i):
    return p["rate"] <= c or (p["rate"] <= c2 and p["ins_med"] >= i)


def fit2(pts, T):
    rates = sorted({p["rate"] for p in pts})
    ins = sorted({p["ins_med"] for p in pts})
    best = None
    for a, c in enumerate(rates):
        for c2 in rates[a:]:
            for i in ins:
                kept = [p for p in pts if keep2(p, c, c2, i)]
                if kept and st.mean(p["y"] for p in kept) <= T:
                    key = (len(kept), -c2, c, -i)   # most kept; then the narrowest extension
                    if best is None or key > best[0]:
                        best = (key, (c, c2, i))
    return best[1], best[0][0]


def explore(pts, T):
    """Not pre-specified: rule 3, recurrent junctions <= h (same criterion and
    leave-one-study-out), and false calls by rate band."""
    def fit3(ps, t):
        best = None
        for h in sorted({p["recur"] for p in ps}):
            kept = [p for p in ps if p["recur"] <= h]
            if st.mean(p["y"] for p in kept) <= t:
                best = (h, len(kept))
        return best
    h, n3 = fit3(pts, T)
    sc = [p for p in pts if p["group"] == "selfcal"]
    print(f"\nExploratory rule 3 (not pre-specified): junctions >= 5 reads <= {h}: keeps {n3} of "
          f"{len(pts)}, {sum(p['recur'] <= h for p in sc)} of {len(sc)} selfcal")
    tot1 = tot3 = 0
    for s in sorted({p["study"] for p in pts}):
        ho = [p for p in pts if p["study"] == s]
        tr = [p for p in pts if p["study"] != s]
        trc = [p["y"] for p in tr if p["group"] == "clean"]
        Ts = st.mean(trc) + 2 * st.stdev(trc) if len(trc) > 1 else T
        a1, a3 = fit1(tr, Ts)[0], fit3(tr, Ts)[0]
        k1 = [p for p in ho if p["rate"] <= a1]
        k3 = [p for p in ho if p["recur"] <= a3]
        tot1 += len(k1)
        tot3 += len(k3)
        if {p["sample"] for p in k1} != {p["sample"] for p in k3}:
            print(f"  LOSO {s}: rule 1 (c {a1:.2f}) keeps {len(k1)} y={[p['y'] for p in k1]}; "
                  f"rule 3 (h {a3}) keeps {len(k3)} y={[p['y'] for p in k3]}")
    print(f"  held-out kept: rule 1 {tot1}, rule 3 {tot3}")
    print("\nFalse calls (y) by fold-back rate band, all points:")
    print("band\tn\tmedian y\tmax y\tn with y > T\tn with y >= 10")
    for lo, hi in ((0, 1), (1, 3), (3, 6), (6, 13), (13, 1e9)):
        b = [p["y"] for p in pts if lo <= p["rate"] < hi]
        if b:
            print(f"{lo}-{hi if hi < 1e9 else ''}\t{len(b)}\t{st.median(b)}\t{max(b)}\t"
                  f"{sum(y > T for y in b)}\t{sum(y >= 10 for y in b)}")


def main():
    pts = load()
    with open(f"{S}/calib_points.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(pts[0]), delimiter="\t")
        w.writeheader()
        w.writerows(pts)
    clean = [p["y"] for p in pts if p["group"] == "clean"]
    T = st.mean(clean) + 2 * st.stdev(clean)
    print(f"points: {sum(p['group']=='selfcal' for p in pts)} selfcal, "
          f"{len(clean)} clean, {sum(p['group']=='chimeric' for p in pts)} chimeric")
    print(f"clean controls: y = ours typed + dysgu AF>=0.2, mean {st.mean(clean):.2f}, "
          f"SD {st.stdev(clean):.2f}; T = {T:.2f}\n")

    print("sample\tgroup\tstudy\tinstrument\tadapter\trate\tins_med\tins_short\tread_len\tdepth\trecur5\tours\tdysgu\ty")
    for p in sorted(pts, key=lambda p: p["rate"]):
        print(f"{p['sample']}\t{p['group']}\t{p['study']}\t{p['instrument']}\t{p['adapter']}\t"
              f"{p['rate']:.2f}\t{p['ins_med']}\t{p['ins_short']:.3f}\t{p['read_len']:.0f}\t"
              f"{p['depth']:.0f}\t{p['recur']}\t{p['ours']}\t{p['dysgu']}\t{p['y']}")

    c1, n1 = fit1(pts, T)
    (c, c2, i), n2 = fit2(pts, T)
    sc = [p for p in pts if p["group"] == "selfcal"]
    k1 = sum(p["rate"] <= c1 for p in sc)
    k2 = sum(keep2(p, c, c2, i) for p in sc)
    print(f"\nRule 1: rate <= {c1:.2f}: keeps {n1} of {len(pts)} points, {k1} of {len(sc)} selfcal")
    print(f"Rule 2: rate <= {c:.2f}, or rate <= {c2:.2f} with median insert >= {i}: "
          f"keeps {n2} of {len(pts)}, {k2} of {len(sc)} selfcal")
    extra = [p for p in pts if keep2(p, c, c2, i) and not p["rate"] <= c1]
    lost = [p for p in pts if p["rate"] <= c1 and not keep2(p, c, c2, i)]
    print("  rule 2 adds:", ", ".join(f"{p['sample']}({p['rate']:.1f}/{p['ins_med']}bp/y={p['y']})" for p in extra) or "none")
    print("  rule 2 drops:", ", ".join(f"{p['sample']}({p['rate']:.1f}/y={p['y']})" for p in lost) or "none")

    print("\nLeave-one-study-out (fit on the rest, apply to the held-out study):")
    print("study\tn\tc1\trule1 kept (y)\trule2 (c, c2, i)\trule2 kept (y)")
    tot1 = tot2 = 0
    ex_y = []
    for s in sorted({p["study"] for p in pts}):
        ho = [p for p in pts if p["study"] == s]
        tr = [p for p in pts if p["study"] != s]
        trc = [p["y"] for p in tr if p["group"] == "clean"]
        Ts = st.mean(trc) + 2 * st.stdev(trc) if len(trc) > 1 else T
        a1 = fit1(tr, Ts)[0]
        (a, a2, ai), _ = fit2(tr, Ts)
        h1 = [p for p in ho if p["rate"] <= a1]
        h2 = [p for p in ho if keep2(p, a, a2, ai)]
        ex_y += [p["y"] for p in h2 if p not in h1]
        tot1 += len(h1)
        tot2 += len(h2)
        print(f"{s}\t{len(ho)}\t{a1:.2f}\t{len(h1)} ({sum(p['y'] for p in h1)})\t"
              f"({a:.2f}, {a2:.2f}, {ai})\t{len(h2)} ({sum(p['y'] for p in h2)})")
    print(f"held-out kept: rule 1 {tot1}, rule 2 {tot2}; isolates only rule 2 kept: "
          f"{len(ex_y)}, y = {ex_y}")

    print("\nSpearman with y (not pre-specified except rate and insert):")
    for lab, sub in (("all points", pts), ("selfcal only", sc)):
        ys = [p["y"] for p in sub]
        for k in ("rate", "recur", "ins_med", "ins_short", "read_len"):
            r, pv = spearmanr([p[k] for p in sub], ys)
            print(f"  {lab:13s} {k:9s} rho {r:+.2f}  p {pv:.2g}")

    explore(pts, T)

    try:
        import statsmodels.api as sm
        X = np.array([[math.log10(p["rate"] + 0.01), p["ins_med"] / 100,
                       math.log10(p["rate"] + 0.01) * p["ins_med"] / 100] for p in pts])
        X = sm.add_constant(X)
        y = np.array([p["y"] for p in pts])
        m = sm.NegativeBinomial(y, X).fit(disp=0)
        print("\nNegative binomial (descriptive): y ~ log10 rate + insert/100 + interaction")
        for n, b, se in zip(["const", "log10 rate", "insert/100", "interaction", "alpha"], m.params, m.bse):
            print(f"  {n:12s} {b:+.3f}  SE {se:.3f}")
    except Exception as e:  # noqa: BLE001
        print("NB fit failed:", e, file=sys.stderr)


if __name__ == "__main__":
    main()
