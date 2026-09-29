#!/usr/bin/env python3
"""Collate P4b: what the SV arm actually produced, before it reaches a matrix.

Three things worth seeing, each of which would change how the records should be
used:

  caller agreement   stage 6 measured union and intersection behaving
                     differently, so the fraction both callers found is the
                     fraction a conservative consumer would keep.
  projection rate    an SV needs BOTH breakpoints on the H37Rv path to be placed
                     there. The rate says how much of the SV callset the common
                     frame can express at all.
  ubiquity           an event called in every isolate is a reference-assembly
                     artefact or a mapping hotspot, not 23 independent events --
                     the same check that exposed the 437 matrix sites.
"""
import argparse, collections, csv, os, statistics, sys


def rd(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--dir", default="refbias/p4b")
    ap.add_argument("--out", default="refbias/p4b/p4b_summary.tsv")
    a = ap.parse_args()

    rows, missing = [], []
    allev = []
    for r in rd(a.refmap):
        p = os.path.join(a.dir, f"{r['sample']}.sv_placed.tsv")
        if not os.path.exists(p):
            missing.append(r["sample"]); continue
        ev = rd(p)
        allev.extend(ev)
        c = collections.Counter(x["svtype"] for x in ev)
        rows.append(dict(sample=r["sample"], reference=r["reference"],
                         events=len(ev),
                         placed=sum(1 for x in ev if x["frame"] == "h37rv"),
                         breakend=sum(1 for x in ev if x["frame"] == "bnd"),
                         both_callers=sum(1 for x in ev if x["n_callers"] == "2"),
                         DEL=c["DEL"], DUP=c["DUP"], INS=c["INS"], INV=c["INV"]))
    if not rows:
        print("no P4b output found", file=sys.stderr); return 1
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    print(f"  {len(rows)} isolates"
          f"{'; MISSING: ' + ', '.join(missing) if missing else ''}\n")
    print(f"  {'isolate':<17s}{'events':>7s}{'placed':>7s}{'bnd':>5s}"
          f"{'both':>6s}{'DEL':>5s}{'DUP':>5s}{'INS':>5s}{'INV':>5s}")
    for r in rows:
        print(f"  {r['sample']:<17s}{r['events']:>7d}{r['placed']:>7d}"
              f"{r['breakend']:>5d}{r['both_callers']:>6d}{r['DEL']:>5d}"
              f"{r['DUP']:>5d}{r['INS']:>5d}{r['INV']:>5d}")
    n = len(allev)
    pl = sum(1 for x in allev if x["frame"] == "h37rv")
    bo = sum(1 for x in allev if x["n_callers"] == "2")
    print(f"\n  pooled {n} events over {len(rows)} isolates")
    print(f"    placed on the H37Rv path   {pl}  ({100*pl/n:.1f}%)")
    print(f"    kept as breakend pairs     {n-pl}  ({100*(n-pl)/n:.1f}%)")
    print(f"    found by BOTH callers      {bo}  ({100*bo/n:.1f}%)")
    print(f"    -- stage 6 measured union and intersection behaving differently,")
    print(f"       so this fraction is what a conservative consumer keeps")
    ln = [int(x["svlen"]) for x in allev if x["svlen"].isdigit()]
    if ln:
        ln.sort()
        print(f"    length: median {ln[len(ln)//2]}, "
              f"{sum(1 for x in ln if x >= 500)} at >= 500 bp")
    # ubiquity: the same artefact check that caught the 437 matrix sites
    loc = collections.Counter()
    for x in allev:
        if x["frame"] == "h37rv" and x["h37rv_pos"]:
            loc[(x["svtype"], int(x["h37rv_pos"]) // 500)] += 1
    ubiq = [k for k, v in loc.items() if v >= len(rows) - 1]
    print(f"\n  event loci called in >= {len(rows)-1}/{len(rows)} isolates: "
          f"{len(ubiq)} of {len(loc)} distinct loci")
    if ubiq:
        print(f"    these are reference-artefact candidates, not shared events")
        for k in ubiq[:6]:
            print(f"      {k[0]} near {k[1]*500}")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
