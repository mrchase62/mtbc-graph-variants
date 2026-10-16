#!/usr/bin/env python3
"""Sanity checks of the new1k_p200 merged VCF (P200 build 3893a3762ed1, all
fixes) before the 10K (user, 2026-10-16: "score the 1K VCF first").
Analysis only. There is no assembly truth for collection isolates, so three
independent comparisons, all on H37Rv-frame SNP records outside the build's
repeat mask:

A. P1 frame: each isolate's merged-VCF SNP genotypes against its own P1
   H37Rv-frame haploid HaplotypeCaller calls (refbias/work/new1k_p200_p1/
   <s>.h37rv.vcf.gz). The merged VCF comes from the matched-reference
   alignment projected through the graph, so agreement tests the whole chain.
   P1 has no record = REF (approximate: also no coverage).
B. Drug resistance (check-dr-positive-controls): every TB-Profiler dr_variant
   SNP at frequency >= 0.9 (collection pipeline JSON) should be ALT in the
   merged VCF; plus RRDR status per isolate against the cohort draw's label.
C. Old build: the 525 isolates shared with gwas1000 (CX333 build, merged VCF
   runroot/refbias/gwas1000/p5/merged.vcf.gz): genotype agreement at SNP keys
   present in both VCFs.

    score.py --outdir out
"""
import argparse, collections, gzip, json, os, subprocess, sys
import numpy as np

R = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
RR = f"{R}/runroot/refbias"
NEW = f"{RR}/new1k_p200/p5/merged.vcf.gz"
OLD = f"{RR}/gwas1000/p5/merged.vcf.gz"
P1W = f"{RR}/work/new1k_p200_p1"
MASK = f"{RR}/build/3893a3762ed1/assets/repeat_mask.bed"
COH = f"{R}/analysis/cohort_1k/cohort_1k.tsv"
TBP = ("/n/netscratch/sfortune_lab/Lab/pculviner/notebooks/260728_sv_exploration/"
       "pipeline_outputs/{0}/results/{0}.tbprofiler.json")
BC = os.environ["MTB_BCFTOOLS"]
L = 4411532


def snp_matrix(vcf, extra_filter=""):
    """(samples, keys [(pos, alt)], G int8: 1 ALT, 0 REF, -1 no call)."""
    s = subprocess.run([BC, "query", "-l", vcf], capture_output=True, text=True, check=True).stdout.split()
    flt = 'INFO/CLASS="small" && INFO/KIND="SNP" && INFO/FRAME="h37rv"' + extra_filter
    p = subprocess.Popen([BC, "query", "-i", flt, "-f", "%POS\t%REF\t%ALT[\t%GT]\n", vcf],
                         stdout=subprocess.PIPE, text=True)
    keys, rows = [], []
    for l in p.stdout:
        f = l.rstrip("\n").split("\t")
        alt = f[2].split(",")[0]
        if len(f[1]) != 1 or len(alt) != 1: continue
        keys.append((int(f[0]), alt))
        rows.append([1 if g == "1" else 0 if g == "0" else -1 for g in f[3:]])
    return s, keys, np.array(rows, np.int8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    mask = np.zeros(L + 200, bool)
    for l in open(MASK):
        c, s_, e = l.split("\t")[:3]; mask[int(s_):int(e)] = True
    out = []
    P = lambda *x: (out.append(" ".join(map(str, x))), print(*x, flush=True))

    samples, keys, G = snp_matrix(NEW)
    um = np.array([not mask[p - 1] for p, _ in keys])
    P(f"new1k_p200 merged VCF: {len(samples)} isolates, {len(keys):,} H37Rv-frame SNP keys "
      f"({int(um.sum()):,} outside the mask)")
    kidx = {k: i for i, k in enumerate(keys)}
    nocall = (G[um] < 0).mean(axis=0)
    P(f"  no-call share per isolate (unmasked SNP keys): median {np.median(nocall):.4f}, "
      f"p99 {np.quantile(nocall, .99):.4f}, max {nocall.max():.4f}")

    # A. P1 frame
    rows = []
    for j, s in enumerate(samples):
        p1 = {}
        with gzip.open(f"{P1W}/{s}.h37rv.vcf.gz", "rt") as fh:
            for l in fh:
                if l[0] == "#": continue
                f = l.split("\t", 10)
                gt = f[9].split(":")[0]
                if gt in ("0", "."): continue
                alt = f[4].split(",")[int(gt) - 1] if gt.isdigit() else ""
                if len(f[3]) == 1 and len(alt) == 1 and not mask[int(f[1]) - 1]:
                    p1[(int(f[1]), alt)] = 1
        col = G[:, j]
        m_alt = {keys[i] for i in np.where((col == 1) & um)[0]}
        both = len(m_alt & set(p1)); m_only = len(m_alt - set(p1)); p_only_keys = set(p1) - m_alt
        p_only_nocall = sum(1 for k in p_only_keys if k in kidx and col[kidx[k]] < 0)
        p_only_ref = sum(1 for k in p_only_keys if k in kidx and col[kidx[k]] == 0)
        p_only_absent = len(p_only_keys) - p_only_nocall - p_only_ref
        rows.append((s, both, m_only, len(p_only_keys), p_only_nocall, p_only_ref, p_only_absent))
    with open(f"{a.outdir}/p1_concordance.tsv", "w") as fo:
        fo.write("sample\tboth_alt\tmerged_only\tp1_only\tp1_only_nocall\tp1_only_ref\tp1_only_no_key\n")
        for r in rows: fo.write("\t".join(map(str, r)) + "\n")
    tot = np.array([r[1:] for r in rows]).sum(0)
    P("\nA. Merged VCF vs each isolate's own P1 H37Rv-frame calls (unmasked SNPs)")
    P(f"  ALT in both {tot[0]:,}; merged only {tot[1]:,} ({tot[1] / max(1, tot[0] + tot[1]):.3%} of merged ALT); "
      f"P1 only {tot[2]:,} ({tot[2] / max(1, tot[0] + tot[2]):.3%} of P1 ALT): of those NOCALL in merged {tot[3]:,}, "
      f"REF {tot[4]:,}, no merged key {tot[5]:,}")
    mo = np.array([r[2] for r in rows]); po = np.array([r[3] for r in rows])
    P(f"  per isolate: merged-only median {np.median(mo):g} (p99 {np.quantile(mo, .99):g}); "
      f"P1-only median {np.median(po):g} (p99 {np.quantile(po, .99):g})")
    worst = sorted(rows, key=lambda r: -(r[2] + r[3]))[:5]
    P("  worst isolates (sample, both, merged-only, P1-only):", [(r[0], r[1], r[2], r[3]) for r in worst])

    # B. drug resistance
    coh = {r.split("\t")[0]: r.rstrip("\n").split("\t") for r in open(COH)}
    hdr = coh.pop("sample")
    ri = hdr.index("rrdr")
    found, total, by_gene, rr = 0, 0, collections.Counter(), collections.Counter()
    miss_ex = []
    for j, s in enumerate(samples):
        try:
            jd = json.load(open(TBP.format(s)))
        except (OSError, ValueError):
            continue
        rr_vcf = 0
        for v in jd.get("dr_variants", []):
            ch = v.get("change", "") or ""
            ref, alt = v.get("ref", ""), v.get("alt", "")
            pos = int(v.get("pos", 0)); fr = float(v.get("freq", 1))
            if fr < 0.9 or len(ref) != 1 or len(alt) != 1: continue
            g = v.get("gene_name") or v.get("gene") or ""
            i = kidx.get((pos, alt))
            state = "noKEY" if i is None else ("ALT" if G[i, j] == 1 else "REF" if G[i, j] == 0 else "NOCALL")
            total += 1; found += state == "ALT"; by_gene[(g, state)] += 1
            if state != "ALT" and len(miss_ex) < 10: miss_ex.append((s, g, pos, alt, state, "masked" if mask[pos - 1] else ""))
        for (p_, alt_), i in kidx.items():
            if 761082 <= p_ <= 761162 and G[i, j] == 1: rr_vcf = 1
        if s in coh:
            rr[(int(coh[s][ri]), rr_vcf)] += 1
    P("\nB. TB-Profiler DR SNPs (freq >= 0.9) found as ALT in the merged VCF")
    P(f"  {found:,} of {total:,} ({found / max(1, total):.3%})")
    genes = sorted({g for g, _ in by_gene}, key=lambda g: -sum(v for (x, _), v in by_gene.items() if x == g))
    for g in genes[:12]:
        P(f"    {g:8s} " + ", ".join(f"{st} {by_gene[(g, st)]}" for st in ("ALT", "REF", "NOCALL", "noKEY") if by_gene[(g, st)]))
    P("  examples not ALT:", miss_ex)
    P(f"  RRDR, cohort draw label x merged VCF (any ALT in rpoB 761082-761162): "
      f"label1&vcf1 {rr[(1, 1)]}, label1&vcf0 {rr[(1, 0)]}, label0&vcf1 {rr[(0, 1)]}, label0&vcf0 {rr[(0, 0)]}")

    # C. old build
    os_, okeys, OG = snp_matrix(OLD)
    osset = set(os_)
    shared = [s for s in samples if s in osset]
    nj = [samples.index(s) for s in shared]; oj = [os_.index(s) for s in shared]
    common = [(i, okidx) for okidx, k in enumerate(okeys) if k in kidx and not mask[k[0] - 1] for i in [kidx[k]]]
    ni = np.array([c[0] for c in common]); oi = np.array([c[1] for c in common])
    A_ = G[np.ix_(ni, nj)]; B_ = OG[np.ix_(oi, oj)]
    both = (A_ >= 0) & (B_ >= 0)
    agree = (A_ == B_) & both
    alt_new_only = ((A_ == 1) & (B_ == 0)).sum(); alt_old_only = ((A_ == 0) & (B_ == 1)).sum()
    okset = set(okeys)
    new_only_keys = sum(1 for k in kidx if not mask[k[0] - 1] and k not in okset)
    P("\nC. Old build (gwas1000, CX333 graph) vs new build, shared isolates")
    P(f"  shared isolates {len(shared)}; unmasked SNP keys in both VCFs {len(common):,} "
      f"(only in new {new_only_keys:,}, only in old {sum(1 for k in okeys if not mask[k[0] - 1] and k not in kidx):,})")
    P(f"  cells called in both {int(both.sum()):,}: agree {int(agree.sum()) / max(1, int(both.sum())):.4%}; "
      f"ALT only in new {int(alt_new_only):,}, ALT only in old {int(alt_old_only):,}")
    open(f"{a.outdir}/summary.txt", "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
