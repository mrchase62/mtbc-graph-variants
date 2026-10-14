#!/usr/bin/env python3
"""Steps 2-5 of the fold-back vs mixed-call proposal. Analysis only.

Reads collect.py's outputs. "Mixed" = PASS calls with 0.1 <= AF < 0.9,
unmasked; trace (< 0.1) is reported, never modelled (user).

2. Distribution of mixed calls per isolate, by lineage.
3. Dose-response across fold-back bins within lineage; negative binomial
   model of mixed SNV count on log(1 + fold-back) adjusted for log depth,
   read length, error rate, duplicate rate, lineage and the TB-Profiler
   mixed-lineage flag.
4. Cutoff table: excess mixed calls per retained isolate over the clean
   baseline (fold-back < 0.5), and isolates kept, at 1.0 / 1.5 / 2.0 / 2.5.
5. Site level: recurrent mixed sites, their enrichment in high (> 2.5) vs
   clean (< 0.5) isolates, AF, strand balance, clustering (ECNT > 1),
   barcode status, inverted repeats nearby, DR-gene overlap.

    analyze.py --indir out --outdir out
"""
import argparse, collections, math, os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
H37RV = f"{REPO}/runroot/refbias/build/7713a8d71d8e-fix1/refs/GCF_000195955.fasta"
TBDB = "/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/tb-profiler/share/tbprofiler/tbdb.bed"
LIN = ["lineage1", "lineage2", "lineage3", "lineage4"]
BINS = [0, 0.25, 0.5, 1, 1.5, 2, 2.5, 5, 10, np.inf]
MIX = ("mixed_high", "mixed_mid", "mixed_low")


def ir_score(seq, i, w=150, k=12):
    """Number of k-mers within +-w of i whose reverse complement also occurs there."""
    s = seq[max(0, i - w): i + w].upper()
    rc = s[::-1].translate(str.maketrans("ACGT", "TGCA"))
    km = {s[j:j + k] for j in range(len(s) - k + 1)}
    return sum(1 for j in range(len(rc) - k + 1) if rc[j:j + k] in km)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    out = []
    P = lambda *x: out.append(" ".join(str(y) for y in x))

    d = pd.read_csv(f"{a.indir}/isolates.tsv", sep="\t", low_memory=False)
    for t in ("snv", "indel"):
        d[f"mix_{t}"] = sum(d[f"{t}_pass_{c}"] for c in MIX)
        d[f"trace_{t}"] = d[f"{t}_pass_trace"]
    d["fb"] = d["inv_chimera_per_1k"].astype(float)
    d["fbbin"] = pd.cut(d["fb"], BINS, right=False)
    d["mixlin"] = ((d["conflict_lineages"].fillna("").astype(str).str.len() > 0)
                   | (pd.to_numeric(d["n_mixed"], errors="coerce").fillna(0) > 0)).astype(int)
    d = d[d["lineage"].isin(LIN)].copy()
    P(f"ISOLATES: {len(d):,} lineage 1-4 standard-pass isolates with Mutect2 calls and a fold-back rate;"
      f" {int((d.fb > 2.5).sum()):,} above 2.5; TB-Profiler mixed-lineage flag on {int(d.mixlin.sum()):,}")

    # 2. distribution
    P("\n2. MIXED CALLS PER ISOLATE (PASS, unmasked, 0.1 <= AF < 0.9): median / p90 / p99 / max")
    for t in ("snv", "indel"):
        for l in LIN:
            x = d.loc[d.lineage == l, f"mix_{t}"]
            P(f"  {t:5s} {l}: {x.median():g} / {x.quantile(.9):g} / {x.quantile(.99):g} / {x.max():g}")
    P("  trace (AF < 0.1, reported only): SNV median "
      f"{d.trace_snv.median():g}, p90 {d.trace_snv.quantile(.9):g}; indel median {d.trace_indel.median():g}")

    # 3. dose-response
    P("\n3. DOSE-RESPONSE: mean (median) mixed SNV / mixed indel per isolate by fold-back bin")
    g = d.groupby(["lineage", "fbbin"], observed=True)
    tab = g.agg(n=("sample", "size"), snv_mean=("mix_snv", "mean"), snv_med=("mix_snv", "median"),
                indel_mean=("mix_indel", "mean"), indel_med=("mix_indel", "median")).reset_index()
    tab.to_csv(f"{a.outdir}/dose_response.tsv", sep="\t", index=False)
    for l in LIN:
        P(f"  {l}")
        for _, r in tab[tab.lineage == l].iterrows():
            P(f"    {str(r.fbbin):14s} n={r.n:6,d}  SNV {r.snv_mean:6.1f} ({r.snv_med:g})  indel {r.indel_mean:6.1f} ({r.indel_med:g})")
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    mids = [0.125, 0.375, 0.75, 1.25, 1.75, 2.25, 3.75, 7.5, 15]
    for k, t in enumerate(("snv", "indel")):
        for l in LIN:
            s = tab[tab.lineage == l]
            ax[k].plot([mids[BINS.index(b.left)] for b in s.fbbin], s[f"{t}_mean"], marker="o", label=l)
        ax[k].set_xscale("log"); ax[k].axvline(2.5, ls="--", c="grey")
        ax[k].set_xlabel("fold-back chimeras per 1,000 reads"); ax[k].set_ylabel(f"mean mixed {t.upper()} calls")
    ax[0].legend(); fig.tight_layout(); fig.savefig(f"{a.outdir}/dose_response.png", dpi=130)

    m = d.dropna(subset=["meandepth", "read_len", "error_rate", "duplicate_rate"]).copy()
    m["lfb"] = np.log1p(m.fb); m["ldepth"] = np.log(m.meandepth.astype(float))
    P("\n   Negative binomial, mixed count ~ log(1 + fold-back) + log depth + read length + error rate"
      " + duplicate rate + lineage + mixed-lineage flag")
    for t in ("snv", "indel"):
        f = (f"mix_{t} ~ lfb + ldepth + read_len + error_rate + duplicate_rate + C(lineage) + mixlin")
        fit = smf.glm(f, data=m, family=sm.families.NegativeBinomial(alpha=1.0)).fit()
        # alpha by method of moments on the Poisson residuals
        pfit = smf.glm(f, data=m, family=sm.families.Poisson()).fit()
        mu = pfit.mu; y = m[f"mix_{t}"]
        alpha = max(1e-6, (((y - mu) ** 2 - y) / mu ** 2).mean())
        fit = smf.glm(f, data=m, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
        ci = fit.conf_int().loc["lfb"]
        P(f"  {t}: n={len(m):,}, alpha={alpha:.3f}; per unit log(1+fold-back) the count changes x{math.exp(fit.params['lfb']):.3f}"
          f" (95% CI {math.exp(ci[0]):.3f}-{math.exp(ci[1]):.3f}); fold-back 0.25 -> 2.5 -> 10 predicts"
          f" x{math.exp(fit.params['lfb'] * (math.log1p(2.5) - math.log1p(0.25))):.2f}"
          f" and x{math.exp(fit.params['lfb'] * (math.log1p(10) - math.log1p(0.25))):.2f}")
        open(f"{a.outdir}/nb_{t}.txt", "w").write(fit.summary().as_text())

    # 4. cutoff table
    P("\n4. CUTOFFS: isolates kept, and mean mixed SNV / indel per kept isolate above the clean baseline (< 0.5)")
    for l in LIN:
        x = d[d.lineage == l]; base = x[x.fb < 0.5]
        P(f"  {l} (all {len(x):,}; baseline n={len(base):,}, SNV {base.mix_snv.mean():.1f}, indel {base.mix_indel.mean():.1f})")
        for c in (1.0, 1.5, 2.0, 2.5):
            k = x[x.fb <= c]
            P(f"    <= {c:3.1f}: keep {len(k):6,} ({len(k) / len(x):.0%}); excess SNV {k.mix_snv.mean() - base.mix_snv.mean():+.2f},"
              f" indel {k.mix_indel.mean() - base.mix_indel.mean():+.2f}")

    # 5. sites
    c = pd.read_csv(f"{a.indir}/mixed_calls.tsv.gz", sep="\t",
                    dtype={"ref": str, "alt": str})
    c = c[(c["pass"] == 1) & c["sample"].isin(d["sample"])]
    grp = d.set_index("sample").fb
    c["fb"] = c["sample"].map(grp)
    n_hi, n_lo = int((d.fb > 2.5).sum()), int((d.fb < 0.5).sum())
    c["sb"] = np.minimum(c.alt_fwd, c.alt_rev) / np.maximum(c.alt_fwd + c.alt_rev, 1)
    s = c.groupby(["pos", "ref", "alt", "type"]).agg(
        isolates=("sample", "size"), hi=("fb", lambda v: int((v > 2.5).sum())), lo=("fb", lambda v: int((v < 0.5).sum())),
        af_median=("af", "median"), strand_minor=("sb", "median"), clustered=("ecnt", lambda v: float((v > 1).mean())),
        barcode=("barcode", "max")).reset_index()
    s["rate_hi"] = s.hi / n_hi; s["rate_lo"] = s.lo / n_lo
    s["enrichment"] = (s.hi + 0.5) / n_hi / ((s.lo + 0.5) / n_lo)
    rec = s[s.isolates >= 20].copy()
    rec["fisher_p"] = [fisher_exact([[h, n_hi - h], [l_, n_lo - l_]], alternative="greater")[1] for h, l_ in zip(rec.hi, rec.lo)]
    rec["artifact_like"] = ((rec.enrichment >= 3) & (rec.fisher_p < 1e-6)).astype(int)
    seq = "".join(l.strip() for l in open(H37RV) if not l.startswith(">"))
    rec["ir_kmers"] = [ir_score(seq, p - 1) for p in rec.pos]
    dr = [(int(l.split("\t")[1]), int(l.split("\t")[2]), l.split("\t")[4]) for l in open(TBDB) if l.strip()]
    rec["dr_gene"] = [next((g for a_, b_, g in dr if a_ <= p <= b_), "") for p in rec.pos]
    rec.sort_values(["artifact_like", "enrichment"], ascending=False).to_csv(f"{a.outdir}/recurrent_sites.tsv", sep="\t", index=False)
    P(f"\n5. SITES: {len(s):,} distinct mixed PASS calls; {len(rec):,} recurrent (>= 20 isolates)")
    P(f"  high (> 2.5) {n_hi:,} isolates, clean (< 0.5) {n_lo:,}")
    for lab, sub in (("artifact-like (enrichment >= 3, p < 1e-6)", rec[rec.artifact_like == 1]),
                     ("other recurrent", rec[rec.artifact_like == 0])):
        if len(sub):
            P(f"  {lab}: {len(sub):,} sites; median AF {sub.af_median.median():.2f}; strand minor share {sub.strand_minor.median():.2f};"
              f" clustered {sub.clustered.median():.2f}; barcode {int(sub.barcode.sum())}; DR-gene {int((sub.dr_gene != '').sum())};"
              f" inverted-repeat 12-mers within 150 bp median {sub.ir_kmers.median():g}")
    # DR genes overall: mixed calls per isolate in DR genes by fold-back bin
    inside = np.zeros(4411533, bool)
    for a_, b_, _ in dr: inside[a_:b_ + 1] = True
    c["dr"] = inside[c.pos.values]
    per = c[c.dr].groupby("sample").size()
    d["mix_dr"] = d["sample"].map(per).fillna(0)
    P("\n   Mixed calls in DR genes per isolate by fold-back bin (expect flat if real):")
    for b, x in d.groupby("fbbin", observed=True):
        P(f"    {str(b):14s} n={len(x):6,d}  mean {x.mix_dr.mean():.3f}")
    open(f"{a.outdir}/summary.txt", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
