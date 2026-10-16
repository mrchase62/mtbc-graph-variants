#!/usr/bin/env python3
"""Per isolate: ALT cells of DEL interval records that its own deepest
sublineage does not share (ALT frequency among the other called isolates of
that sublineage <= 0.2, >= 3 such isolates), inside vs outside the repeat
mask. A candidate per-sample QC metric. Writes out/private_alt_per_sample.tsv.
"""
import collections
import numpy as np, pandas as pd
from common import *
d = load_records(); st = d["st"]; recs = d["recs"]; samples = d["samples"]
lin, deep = cohort(); mask = mask_array()
idx = [i for i, r in enumerate(recs) if klass(r) == "DEL" and component(r) == "interval"]
inm = np.array([mask[recs[i]["pos"] - 1:recs[i]["end"]].any() for i in idx])
S = st[idx]
grp = collections.defaultdict(list)
for j, s in enumerate(samples): grp[deep[s]].append(j)
out = np.zeros((len(samples), 2), int)
for g, js in grp.items():
    js = np.array(js)
    if len(js) < 4: continue
    sub = S[:, js]
    called = sub <= 1; alt = sub == 1
    na, nc = alt.sum(1), called.sum(1)
    for k, j in enumerate(js):
        a = alt[:, k]
        fo = (na - a) / np.maximum(1, nc - called[:, k])
        ok = (nc - called[:, k]) >= 3
        priv = a & ok & (fo <= 0.2)
        out[j, 0] += int((priv & inm).sum()); out[j, 1] += int((priv & ~inm).sum())
t = pd.DataFrame(out, columns=["private_in_mask", "private_outside"], index=samples)
t["sublineage"] = [deep[s] for s in samples]
t.to_csv(f"{OUT}/private_alt_per_sample.tsv", sep="\t")
print(t.describe(percentiles=[.5, .9, .95, .99]))
print(t.sort_values("private_in_mask", ascending=False).head(15))
odd = pd.read_csv(f"{OUT}/odd_side_samples.tsv", sep="\t", index_col=0)
print(t.loc[odd.index[:10]])
