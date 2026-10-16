#!/usr/bin/env python3
"""For each ALT/REF-discordant cell of a near-identical pair (pairs.py), which
side is the odd one out? Reference: ALT frequency of the record among the
other called isolates of the same deepest sublineage (pair excluded).
>= 0.8 -> the REF side is odd; <= 0.2 -> the ALT side is odd; else unclear.
Per isolate: number of cells where it is the odd side. Writes out/odd_side.tsv
(per class/component/mask) and out/odd_side_samples.tsv.
"""
import collections
import numpy as np, pandas as pd
from common import *
d = load_records(); st = d["st"]; samples = d["samples"]; col = {s: j for j, s in enumerate(samples)}
lin, deep = cohort()
rm = pd.read_csv(f"{R}/runroot/refbias/new1k_p200/p1/refmap.tsv", sep="\t").set_index("sample")
c = pd.read_csv(f"{OUT}/pair_discordant_cells.tsv.gz", sep="\t")
c = c[c.st_a.isin(["ALT", "REF"]) & c.st_b.isin(["ALT", "REF"])].copy()
grp = collections.defaultdict(list)
for s in samples: grp[deep[s]].append(col[s])
odd, freq = [], []
for r in c.itertuples():
    ja, jb = col[r.a], col[r.b]
    js = [j for j in grp[deep[r.a]] if j not in (ja, jb)]
    x = st[r.idx, js]; x = x[x <= 1]
    f = x.mean() if len(x) >= 3 else np.nan
    freq.append(f)
    alt_s, ref_s = (r.a, r.b) if r.st_a == "ALT" else (r.b, r.a)
    odd.append(ref_s if f >= 0.8 else alt_s if f <= 0.2 else "unclear")
c["clade_alt_freq"] = freq; c["odd_side"] = odd
c["odd_is"] = np.where(c.odd_side == "unclear", "unclear",
                       np.where(c.odd_side == np.where(c.st_a == "ALT", c.a, c.b), "ALT_side(false_ALT?)", "REF_side(false_REF?)"))
t = c.groupby(["class", "component", "in_mask", "odd_is"]).size().unstack(fill_value=0)
t.to_csv(f"{OUT}/odd_side.tsv", sep="\t")
print(t)
s = c[c.odd_side != "unclear"].groupby("odd_side").size().sort_values(ascending=False)
s = s.to_frame("odd_cells")
s["meandepth"] = rm.loc[s.index, "meandepth"].values
s["sublineage"] = [deep[x] for x in s.index]
s.to_csv(f"{OUT}/odd_side_samples.tsv", sep="\t")
print(s.head(15))
