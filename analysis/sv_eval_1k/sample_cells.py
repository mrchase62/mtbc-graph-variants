#!/usr/bin/env python3
"""Fixed-seed sample of OLD/NEW disagreement cells for the read check.
Strata = class x direction: 'old_alt' (OLD ALT, NEW REF/ABSENT) and 'new_alt'
(NEW ALT, OLD REF/ABSENT); at most 2 cells per record pair so a few
high-carrier pairs do not fill a stratum. DEL/INS/IS: 30 per stratum (all
when fewer). Accessory presence has no ALT<->REF flips, so its strata are
ALT<->NOCALL (10 each), checked as a bonus.
"""
import pickle, random, collections
import pandas as pd, numpy as np
D = "analysis/sv_eval_1k/out"
rng = random.Random(20261016)
c = pd.read_csv(f"{D}/disagree_cells.tsv.gz", sep="\t")
c["stratum"] = np.where(c.old_st == "ALT", "old_alt", "new_alt")
rows = []
for (k, s), g in c.groupby(["class", "stratum"]):
    idx = list(g.index); rng.shuffle(idx)
    per = collections.Counter(); take = []
    for i in idx:
        key = (g.at[i, "old_idx"], g.at[i, "new_idx"])
        if per[key] >= 2: continue
        per[key] += 1; take.append(i)
        if len(take) == 30: break
    rows.append(c.loc[take])
o = pickle.load(open(f"{D}/old.pkl", "rb")); n = pickle.load(open(f"{D}/new.pkl", "rb"))
p = pd.read_csv(f"{D}/pairs.tsv", sep="\t"); p = p[p["class"] == "ACC"]
shared = sorted(set(o["samples"]) & set(n["samples"]))
io = [o["samples"].index(s) for s in shared]; inn = [n["samples"].index(s) for s in shared]
acc = {"old_alt": [], "new_alt": []}
for _, r in p.iterrows():
    so = o["st"][r.old_idx, io]; sn = n["st"][r.new_idx, inn]
    for t in range(len(shared)):
        if so[t] == 1 and sn[t] == 3: acc["old_alt"].append((r, t, "ALT", "NOCALL"))
        if sn[t] == 1 and so[t] == 3: acc["new_alt"].append((r, t, "NOCALL", "ALT"))
for s, L in acc.items():
    rng.shuffle(L); per = collections.Counter(); take = []
    for r, t, a, b in L:
        if per[r.new_idx] >= 2: continue
        per[r.new_idx] += 1
        nr = n["recs"][r.new_idx]
        take.append(dict({"class": "ACC", "old_idx": r.old_idx, "new_idx": r.new_idx, "sample": shared[t], "old_st": a,
                          "new_st": b, "new_ev": ".", "sub_old": r.sub_old, "sub_new": r.sub_new, "pos": nr["pos"],
                          "end": nr["end"], "svlen": nr["svlen"], "stratum": s}))
        if len(take) == 10: break
    rows.append(pd.DataFrame(take))
out = pd.concat(rows, ignore_index=True)
# record details for the checker
out["seq"] = [n["recs"][j]["seq"] or "" for j in out.new_idx]
out["old_pos"] = [o["recs"][j]["pos"] for j in out.old_idx]
out["old_end"] = [o["recs"][j]["end"] for j in out.old_idx]
out["old_len"] = [o["recs"][j]["svlen"] for j in out.old_idx]
out["acclocus"] = [n["recs"][j]["info"].get("ACCLOCUS", "") for j in out.new_idx]
out["siteclass"] = [n["recs"][j]["info"].get("SITECLASS", "") for j in out.new_idx]
out.to_csv(f"{D}/readcheck_sample.tsv", sep="\t", index=False)
print(out.groupby(["class", "stratum"]).size())
