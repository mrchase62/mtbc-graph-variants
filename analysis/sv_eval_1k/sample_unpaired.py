#!/usr/bin/env python3
"""Fixed-seed sample of ALT cells of records the other build has no match
for (unpaired, 'only'/'secondary' in events.tsv), DEL and INS, 15 per
build, at most 2 cells per record, among records with >= 2 ALT carriers in
the shared isolates (singleton 'called' records dominate the OLD side and
are counted separately). Strata old_only / new_only; written in the format
of readcheck_sample.tsv (the unpaired side's index in both idx columns).
"""
import pickle, random, collections
import numpy as np, pandas as pd
D = "analysis/sv_eval_1k/out"
rng = random.Random(20261017)
o = pickle.load(open(f"{D}/old.pkl", "rb")); n = pickle.load(open(f"{D}/new.pkl", "rb"))
sh = sorted(set(o["samples"]) & set(n["samples"]))
p = pd.read_csv(f"{D}/pairs.tsv", sep="\t")
def klass(r):
    if r["chrom"] != "NC_000962.3": return None
    c, t = r["cls"], r["svtype"]
    if t == "DEL" and c in ("sv", "small", "is6110"): return "DEL"
    if t == "INS" and c in ("sv", "small"): return "INS"
    return None
rows = []
for nm, d, paired in (("old", o, set(p.old_idx)), ("new", n, set(p.new_idx))):
    idx = [d["samples"].index(s) for s in sh]
    for k in ("DEL", "INS"):
        cells = []
        for i, r in enumerate(d["recs"]):
            if klass(r) != k or i in paired: continue
            a = np.where(d["st"][i, idx] == 1)[0]
            if len(a) < 2: continue
            cells += [(i, t) for t in a]
        rng.shuffle(cells); per = collections.Counter(); take = []
        for i, t in cells:
            if per[i] >= 2: continue
            per[i] += 1; take.append((i, t))
            if len(take) == 15: break
        for i, t in take:
            r = d["recs"][i]
            rows.append({"class": k, "old_idx": i, "new_idx": i, "sample": sh[t], "old_st": "ALT" if nm == "old" else "none",
                         "new_st": "ALT" if nm == "new" else "none", "new_ev": d["evnames"][int(d["ev"][i, d["samples"].index(sh[t])])],
                         "sub_old": "", "sub_new": "", "pos": r["pos"], "end": r["end"], "svlen": r["svlen"],
                         "stratum": f"{nm}_only", "seq": r["seq"] or "", "acclocus": "", "siteclass": ""})
pd.DataFrame(rows).to_csv(f"{D}/readcheck_sample_unpaired.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).groupby(["class", "stratum"]).size())
