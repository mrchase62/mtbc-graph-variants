#!/usr/bin/env python3
"""Markdown tables for the README from the out/ files. Prints to stdout."""
import pandas as pd, numpy as np
D = "analysis/sv_eval_1k/out"
def md(df):
    cols = list(df.columns)
    s = "| " + " | ".join(map(str, cols)) + " |\n|" + "---|" * len(cols) + "\n"
    for _, r in df.iterrows(): s += "| " + " | ".join(str(r[c]) for c in cols) + " |\n"
    return s
print("### transitions\n")
t = pd.read_csv(f"{D}/transitions.tsv", sep="\t")
for k in ["DEL", "INS", "IS", "ACC"]:
    g = t[t["class"] == k].pivot_table(index="old", columns="new", values="n", aggfunc="sum", fill_value=0)
    g = g.reindex(index=[x for x in ["ALT", "REF", "ABSENT", "NOCALL"] if x in g.index],
                  columns=[x for x in ["ALT", "REF", "ABSENT", "NOCALL", "DUP"] if x in g.columns], fill_value=0)
    g.insert(0, "OLD \\ NEW", g.index)
    print(f"**{k}**\n"); print(md(g.reset_index(drop=True)))
print("### read check\n")
r = pd.read_csv(f"{D}/readcheck_results_ctx.tsv", sep="\t")
r["nearby"] = r.other_side_alt_nearby.notna()
rows = []
for (k, s), g in r.groupby(["class", "stratum"], sort=False):
    alt = "old" if s == "old_alt" else "new"
    rows.append({"class": k, "ALT side": alt.upper(), "cells": len(g),
                 "reads: ALT (ALT side right)": int((g.verdict_final == "ALT").sum()),
                 "  of which other build ALT at a nearby record": int(((g.verdict_final == "ALT") & g.nearby).sum()),
                 "reads: REF (other side right)": int((g.verdict_final == "REF").sum()),
                 "unclear": int((g.verdict_final == "unclear").sum())})
print(md(pd.DataFrame(rows)))
for f, lab in (("readcheck_results_unpaired_ctx.tsv", "unpaired"),):
    u = pd.read_csv(f"{D}/{f}", sep="\t")
    print(f"**{lab}**\n"); print(md(u.groupby(["class", "stratum", "verdict_final"]).size().unstack(fill_value=0).reset_index()))
u = pd.read_csv(f"{D}/readcheck_results_nocall.tsv", sep="\t")
print("**OLD ALT -> NEW NOCALL**\n"); print(md(u.groupby(["new_ev", "verdict"]).size().unstack(fill_value=0).reset_index()))
