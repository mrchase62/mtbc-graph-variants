#!/usr/bin/env python3
"""Dump records in a window with per-lineage ALT/REF/ABSENT/NOCALL counts.
   locus.py start end [class...]"""
import sys, collections
import numpy as np
from common import *
d = load_records(); lin, deep = cohort()
L = np.array([lin[s] for s in d["samples"]])
a, b = int(sys.argv[1]), int(sys.argv[2]); ks = sys.argv[3:]
for i, r in enumerate(d["recs"]):
    if r["chrom"] != "NC_000962.3" or r["end"] < a or r["pos"] > b: continue
    k = klass(r)
    if ks and k not in ks: continue
    s = d["st"][i]
    per = " ".join(f"{g[-1]}:{(s[L==g]==1).sum()}/{(s[L==g]==0).sum()}/{(s[L==g]==2).sum()}/{(s[L==g]>=3).sum()}" for g in sorted(set(L)))
    inf = r["info"]
    print(i, k, r["pos"], r["end"], r["svlen"], component(r), inf.get("SVSOURCE","."), inf.get("AA_STATE","."), inf.get("AA_FLAG","."), "|", per)
