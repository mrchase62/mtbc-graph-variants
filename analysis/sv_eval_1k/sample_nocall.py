#!/usr/bin/env python3
"""Fixed-seed sample of matched DEL cells OLD ALT -> NEW NOCALL, 20 with
EV=no_coverage on small_allele records and 10 with EV=ambiguous_depth on
interval records (the two largest groups), at most 2 cells per record, for
the read check (stratum old_alt: checked at the OLD record's coordinates)."""
import pickle, random, collections
import numpy as np, pandas as pd
D = "analysis/sv_eval_1k/out"
rng = random.Random(20261018)
o = pickle.load(open(f"{D}/old.pkl", "rb")); n = pickle.load(open(f"{D}/new.pkl", "rb"))
sh = sorted(set(o["samples"]) & set(n["samples"]))
io = [o["samples"].index(s) for s in sh]; inn = [n["samples"].index(s) for s in sh]
p = pd.read_csv(f"{D}/pairs.tsv", sep="\t"); p = p[p["class"] == "DEL"]
want = {"no_coverage": ("small_allele", 20), "ambiguous_depth": ("interval_graph", 10)}
rows = []
for ev, (sub, N) in want.items():
    cells = []
    for _, r in p[p.sub_new == sub].iterrows():
        so = o["st"][r.old_idx, io]; sn = n["st"][r.new_idx, inn]; e = n["ev"][r.new_idx, inn]
        for t in np.where((so == 1) & (sn == 3))[0]:
            if n["evnames"][int(e[t])] == ev: cells.append((r, t))
    rng.shuffle(cells); per = collections.Counter()
    for r, t in cells:
        if per[r.new_idx] >= 2: continue
        per[r.new_idx] += 1
        rows.append({"class": "DEL", "old_idx": r.old_idx, "new_idx": r.new_idx, "sample": sh[t], "old_st": "ALT", "new_st": "NOCALL",
                     "new_ev": ev, "sub_old": r.sub_old, "sub_new": r.sub_new, "pos": r.pos_new, "end": n["recs"][r.new_idx]["end"],
                     "svlen": r.len_new, "stratum": "old_alt", "seq": "", "acclocus": "", "siteclass": "", "n_cells_group": len(cells)})
        if sum(1 for x in rows if x["new_ev"] == ev) == N: break
pd.DataFrame(rows).to_csv(f"{D}/readcheck_sample_nocall.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).groupby("new_ev").size())
