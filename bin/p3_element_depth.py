#!/usr/bin/env python3
"""IS6110 family copy number, from depth across H37Rv's own element spans.

Why this and not the per-locus copy-number table: reads from every copy of a
repeated element map onto each reference copy, so depth at one locus reflects the
whole family rather than that locus. That makes per-locus depth the wrong
instrument for "how many copies does this sample have" and makes the family
median the right one.

The estimate is deliberately crude and stated as such. Total copies is
approximated as (reference copies) x (median element depth / genome median
depth). It inherits every assumption of depth-based copy number -- uniform
coverage, no GC bias at the element, mapping quality unaffected by copy count --
and none of them is verified here. It is reported as a ratio and an implied count
so that the implied count can be checked against an orthogonal site count, which
is the only way to find out whether the assumptions hold.
"""
import argparse, bisect, csv, statistics, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--positions", required=True,
                    help="odgi output for each element's start and end")
    ap.add_argument("--depth", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    els = []
    for line in open(a.gff):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) > 4:
            els.append((int(f[3]), int(f[4])))
    els.sort()
    if not els:
        print("no elements in the GFF", file=sys.stderr); return 1

    # keyed by source position, never by file order: threaded odgi returns
    # results in completion order
    by_src = {}
    for line in open(a.positions):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            src = int(f[0].rsplit(",", 2)[-2]) + 1
            tgt = int(f[1].rsplit(",", 2)[-2]) + 1
            dist = int(f[2])
        except (ValueError, IndexError):
            continue
        by_src.setdefault(src, (tgt, dist))

    depth = {}
    for line in open(a.depth):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            depth[int(f[1])] = int(f[2])
    if not depth:
        print("empty depth file", file=sys.stderr); return 1
    gmed = statistics.median(depth.values())
    if gmed <= 0:
        print("genome median depth is zero", file=sys.stderr); return 1

    rows, ratios = [], []
    for s, e in els:
        ps, pe = by_src.get(s), by_src.get(e)
        if not ps or not pe or ps[1] != 0 or pe[1] != 0:
            rows.append(dict(sample=a.sample, reference=a.reference,
                             h37rv_start=s, h37rv_end=e, r_start="", r_end="",
                             median_depth="", depth_ratio="",
                             status="no_projection"))
            continue
        lo, hi = sorted((ps[0], pe[0]))
        # The termini sit on collapsed IS6110 nodes, so the two ends can land
        # in DIFFERENT copies in R; the "span" is then hundreds of kb of
        # single-copy sequence with a ratio near 1. Require the projected span
        # to be about the element's own length.
        if abs((hi - lo) - (e - s)) > max(100, (e - s) // 10):
            rows.append(dict(sample=a.sample, reference=a.reference,
                             h37rv_start=s, h37rv_end=e, r_start=lo, r_end=hi,
                             median_depth="", depth_ratio="",
                             status="span_mismatch"))
            continue
        vals = [depth.get(p, 0) for p in range(lo, hi + 1)]
        med = statistics.median(vals) if vals else 0
        ratio = med / gmed
        ratios.append(ratio)
        rows.append(dict(sample=a.sample, reference=a.reference,
                         h37rv_start=s, h37rv_end=e, r_start=lo, r_end=hi,
                         median_depth=med, depth_ratio=round(ratio, 4),
                         status="ok"))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    ok = [r for r in rows if r["status"] == "ok"]
    if ratios:
        fam = statistics.median(ratios)
        print(f"  {a.sample}: {len(ok)}/{len(els)} element spans projected; "
              f"genome median {gmed:.0f}x")
        print(f"    family depth ratio median {fam:.3f} -> implied total copies "
              f"~{fam*len(els):.1f} against {len(els)} reference copies")
    else:
        print(f"  {a.sample}: no element span projected", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
