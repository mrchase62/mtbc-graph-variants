#!/usr/bin/env python3
"""Group ABSENT cells into per-sample events for local reassembly.

An event is one region of one sample:
- the matched reference's deletion, where junction_check.py found one;
- otherwise the span of the sample's cells, grouping cells less than 500 bp
  apart.

Sets:
  calibration  the ctpV cases with known answers, plus core cells whose
               junction verdict is already decided: all deleted_here events,
               and a fixed random sample of present_here events and of
               supported (no depth) + deleted_here events
  unresolved   core contradicted cells with verdict unresolved,
               R_no_genotype or no_R_deletion

Output columns: event, set, sample, start, end, n_cells, expected.
`expected` is the answer the calibration checks against, or empty.
"""
import argparse
import random

import pandas as pd

CTPV = [  # sample, start, end, expected (analysis/ctpV_check.md)
    ("SAMEA2297133", 1078522, 1078756, "present_here"),
    ("SAMEA2297133", 1078811, 1080106, "present_here"),
    ("SAMEA112800746", 1078522, 1078756, "deleted_here"),
    ("SAMN07766100", 1078522, 1078756, "deleted_here"),
    ("SAMEA5542103", 1078519, 1080103, "deleted_here"),
]


def group_cells(df):
    out = []
    for s, g in df.sort_values("start").groupby("sample"):
        cur = None
        for _, r in g.iterrows():
            if pd.notna(r.r_del_start):
                key = (s, int(r.r_del_start), int(r.r_del_end))
                out.append(key)
                continue
            if cur and r.start - cur[2] < 500:
                cur[2] = max(cur[2], int(r.end))
                cur[3] += 1
            else:
                if cur:
                    out.append((cur[0], cur[1], cur[2]))
                cur = [s, int(r.start), int(r.end), 1]
        if cur:
            out.append((cur[0], cur[1], cur[2]))
    return pd.Series(out).value_counts().rename_axis("key").reset_index(name="n_cells")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junction", required=True, help="junction.masked.tsv")
    ap.add_argument("--set", choices=["calibration", "unresolved"], required=True)
    ap.add_argument("--n-present", type=int, default=50)
    ap.add_argument("--n-supported", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    j = pd.read_csv(a.junction, sep="\t")
    core = j[j.mappability == "core"]
    rows = []
    if a.set == "calibration":
        rnd = random.Random(a.seed)
        for s, st, en, exp in CTPV:
            rows.append((s, st, en, 0, exp, "ctpV"))
        for verdict, depth, n in [("deleted_here", "contradicted", None),
                                  ("present_here", "contradicted", a.n_present),
                                  ("deleted_here", "supported", a.n_supported)]:
            sub = core[(core.junction_verdict == verdict) & (core.depth_verdict == depth)]
            ev = group_cells(sub)
            keys = sorted(ev.key.tolist())
            if n is not None and len(keys) > n:
                keys = sorted(rnd.sample(keys, n))
            nc = dict(zip(ev.key, ev.n_cells))
            for k in keys:
                rows.append((k[0], k[1], k[2], nc[k], verdict,
                             f"{depth}/{verdict}"))
    else:
        sub = core[(core.depth_verdict == "contradicted")
                   & core.junction_verdict.isin(["unresolved", "R_no_genotype",
                                                 "no_R_deletion"])]
        for verdict, g in sub.groupby("junction_verdict"):
            ev = group_cells(g)
            for k, n in zip(ev.key, ev.n_cells):
                rows.append((k[0], k[1], k[2], n, "", f"contradicted/{verdict}"))
    df = pd.DataFrame(rows, columns=["sample", "start", "end", "n_cells",
                                     "expected", "source"])
    df = df.drop_duplicates(["sample", "start", "end"]).reset_index(drop=True)
    df.insert(0, "event", [f"{a.set[:3]}{i:04d}" for i in range(len(df))])
    df.to_csv(a.out, sep="\t", index=False)
    print(df.groupby("source").size().to_string())


if __name__ == "__main__":
    main()
