#!/usr/bin/env python3
"""Context for the read-check cells and extra cell counts. Analysis only.
Run from the repo root after readcheck.py (main and --sample unpaired runs).
1. For each checked cell: does the build that wrote REF/ABSENT give the same
   isolate ALT at another record of the same class within 100 bp (an allele
   or representation difference rather than a missed event)?
2. IS6110 records where >= 90% of isolates are ABSENT (whole-site ABSENT).
3. NOCALL share per class and build on the shared isolates (all records).
"""
import pickle, collections
import numpy as np, pandas as pd
D = "analysis/sv_eval_1k/out"
o = pickle.load(open(f"{D}/old.pkl", "rb")); n = pickle.load(open(f"{D}/new.pkl", "rb"))
def klass(r):
    if r["chrom"] != "NC_000962.3": return None
    c, t = r["cls"], r["svtype"]
    if c == "accessory_presence": return "ACC"
    if c == "is6110" and t == "INS": return "IS"
    if t == "DEL" and c in ("sv", "small", "is6110"): return "DEL"
    if t == "INS" and c in ("sv", "small"): return "INS"
    return None
for d in (o, n):
    d["k"] = [klass(r) for r in d["recs"]]
    d["pos"] = np.array([r["pos"] for r in d["recs"]])
r = pd.read_csv(f"{D}/readcheck_results.tsv", sep="\t")
other = []
for _, c in r.iterrows():
    d = n if c.stratum == "old_alt" else o
    j = d["samples"].index(c["sample"])
    own = int(c.new_idx if c.stratum == "old_alt" else c.old_idx)
    hit = [i for i in np.where(np.abs(d["pos"] - c.chk_pos) <= 100)[0]
           if d["k"][i] == c["class"] and i != own and d["st"][i, j] == 1]
    other.append(";".join(f"{d['recs'][i]['pos']}:{d['recs'][i]['svlen']}" for i in hit))
r["other_side_alt_nearby"] = other
# 1b. clip clusters near the site: a REF verdict at the exact site is not
# trusted when the isolate has >= 10 MAPQ>=20 reads soft-clipped >= 15 bp in
# one 100 bp bin within 1 kb outside the checked span (a complex or shifted
# allele the site check cannot see); such REF verdicts become unclear.
import pysam
def clip_bins(sample, a, b):
    bam = pysam.AlignmentFile(f"runroot/refbias/work/new1k_p200_p1/{sample}.h37rv.bam")
    c = collections.Counter()
    for x in bam.fetch("NC_000962.3", max(0, a - 1000), b + 1000):
        if x.is_secondary or x.is_duplicate or x.mapping_quality < 20: continue
        ct = x.cigartuples
        for p, cl in ((x.reference_start + 1, ct[0]), (x.reference_end, ct[-1])):
            if cl[0] == 4 and cl[1] >= 15 and not (a - 40 <= p <= b + 40) and a - 1000 <= p <= b + 1000:
                c[p // 100] += 1
    return max(c.values()) if c else 0
cb, fin = [], []
for _, c in r.iterrows():
    if c["class"] in ("DEL", "INS"):
        a = int(c.chk_pos); b = int(c.chk_end) if c["class"] == "DEL" else int(c.chk_pos)
        m = clip_bins(c["sample"], a, b)
    else:
        m = -1
    cb.append(m)
    v = c.verdict
    if v == "REF" and m >= 10: v = "unclear"
    fin.append(v)
r["clip_cluster_1kb"] = cb
r["verdict_final"] = fin
alt_side = np.where(r.stratum == "old_alt", "old", "new"); oth = np.where(r.stratum == "old_alt", "new", "old")
r["right_final"] = np.where(r.verdict_final == "ALT", alt_side, np.where(r.verdict_final == "REF", oth, "unclear"))
r.to_csv(f"{D}/readcheck_results_ctx.tsv", sep="\t", index=False)
# 2
rows = []
for nm, d in (("OLD", o), ("NEW", n)):
    for i, rec in enumerate(d["recs"]):
        if d["k"][i] != "IS": continue
        s = d["st"][i]; fa = (s == 2).mean()
        rows.append((nm, rec["pos"], rec["info"].get("SITECLASS"), fa, int((s == 1).sum())))
t = pd.DataFrame(rows, columns=["build", "pos", "siteclass", "absent_frac", "alt"])
with open(f"{D}/is6110_absent_sites.tsv", "w") as f:
    t.to_csv(f, sep="\t", index=False)
print(t.assign(whole=t.absent_frac >= 0.9).groupby(["build", "siteclass", "whole"]).size())
print("cells ABSENT at whole-ABSENT sites:", {b: int((g.absent_frac[g.absent_frac >= .9] * (997 if b == "OLD" else 1000)).sum()) for b, g in t.groupby("build")})
# 3
sh = sorted(set(o["samples"]) & set(n["samples"]))
out = []
for nm, d in (("OLD", o), ("NEW", n)):
    idx = [d["samples"].index(s) for s in sh]
    st = d["st"][:, idx]
    for k in ("DEL", "INS", "IS", "ACC"):
        m = np.array([x == k for x in d["k"]])
        s = st[m]
        out.append((nm, k, int(m.sum()), s.size, int((s == 1).sum()), int((s == 0).sum()), int((s == 2).sum()), int((s == 3).sum()),
                    round((s == 3).mean(), 3)))
pd.DataFrame(out, columns=["build", "class", "records", "cells", "ALT", "REF", "ABSENT", "NOCALL", "nocall_frac"]).to_csv(
    f"{D}/state_counts_shared.tsv", sep="\t", index=False)
print(pd.read_csv(f"{D}/state_counts_shared.tsv", sep="\t").to_string())


# 4. the same clip-cluster guard on the unpaired-cell check
u = pd.read_csv(f"{D}/readcheck_results_unpaired.tsv", sep="\t")
cb = [clip_bins(c["sample"], int(c.chk_pos), int(c.chk_end) if c["class"] == "DEL" else int(c.chk_pos)) for _, c in u.iterrows()]
u["clip_cluster_1kb"] = cb
u["verdict_final"] = np.where((u.verdict == "REF") & (u.clip_cluster_1kb >= 10), "unclear", u.verdict)
u.to_csv(f"{D}/readcheck_results_unpaired_ctx.tsv", sep="\t", index=False)
print(u.groupby(["class", "stratum", "verdict_final"]).size().unstack(fill_value=0))
