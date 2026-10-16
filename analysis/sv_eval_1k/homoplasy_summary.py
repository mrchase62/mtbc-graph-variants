#!/usr/bin/env python3
"""Tables from out/homoplasy_records<suffix>.tsv: per set x class x unit,
records with >= 2 carriers, median changes, share with >= 2 / >= 5 changes,
and the change-count histogram. Unit 'events' = SPLIT-1 groups (DEL/INS),
'record' = each VCF record.
    homoplasy_summary.py [suffix]
"""
import sys
import numpy as np, pandas as pd
mask = np.zeros(4411600, bool)
for l in open("runroot/refbias/build/3893a3762ed1/assets/repeat_mask.bed"):
    c, a, b = l.split("\t")[:3]; mask[int(a):int(b)] = True
suf = sys.argv[1] if len(sys.argv) > 1 else ""
h = pd.read_csv(f"analysis/sv_eval_1k/out/homoplasy_records{suf}.tsv", sep="\t")
h.loc[h.unit.isin(["event_single", "event_split1"]), "unit2"] = "event"
h.loc[h.unit == "record", "unit2"] = "record"
h["masked"] = [bool(mask[p - 1:p + (L if k in ("DEL",) else 1)].any()) for p, L, k in zip(h.pos, h.svlen, h["class"])]
def by_mask(t):
    out = []
    for (s, k, u, m), g in h.groupby(["set", "class", "unit2", "masked"]):
        out.append(dict(set=s, cls=k, unit=u, masked=m, n=len(g), median=g.changes.median(), ge2=round((g.changes >= 2).mean(), 3),
                        ge5=round((g.changes >= 5).mean(), 3)))
    return pd.DataFrame(out)
rows = []
for (s, k, u), g in h.groupby(["set", "class", "unit2"]):
    c = g.changes
    rows.append(dict(set=s, cls=k, unit=u, n=len(g), median=c.median(), mean=round(c.mean(), 2),
                     ge2=round((c >= 2).mean(), 3), ge5=round((c >= 5).mean(), 3),
                     h1=(c == 1).sum(), h2=(c == 2).sum(), h3_4=c.between(3, 4).sum(), h5_9=c.between(5, 9).sum(), h10p=(c >= 10).sum(),
                     merged=int((g.unit == "event_split1").sum()) if u == "event" else "",
                     merged_changes=int(g.changes[g.unit == "event_split1"].sum()) if u == "event" else "",
                     merged_member_sum=int(g.member_change_sum[g.unit == "event_split1"].sum()) if u == "event" else "",
                     total_changes=int(c.sum()), nocall_med=g.nocall_frac.median()))
t = pd.DataFrame(rows)
t.to_csv(f"analysis/sv_eval_1k/out/homoplasy_summary{suf}.tsv", sep="\t", index=False)
print(t.to_string(index=False))
m = by_mask(h)
m.to_csv(f"analysis/sv_eval_1k/out/homoplasy_by_mask{suf}.tsv", sep="\t", index=False)
print(m.to_string(index=False))
