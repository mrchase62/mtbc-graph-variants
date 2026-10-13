#!/usr/bin/env python3
"""scale200_fix: the matched-genome call against production's H37Rv-frame
genotype for the same inherited deletion, and a summary of inherited
insertions by repeat class. Analysis only.

Deletions: an inherited row matches a catalogued interval (`svi:DEL:pos:len`
in p5/svgt_iv_states.tsv) for the same sample when the start is within 50 bp
and the lengths agree within 20%.

    compare_h37rv.py out/scale200fix.tsv svgt_iv_states.tsv repeat_mask.bed
"""
import bisect, collections, csv, sys

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
mask = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(sys.argv[3]) if l.strip())
mst = [a for a, _ in mask]


def masked(p, e):
    i = bisect.bisect_right(mst, e) - 1
    while i >= 0 and i > bisect.bisect_right(mst, p) - 50:
        if mask[i][0] < e + 1 and mask[i][1] > p:
            return True
        i -= 1
    return False


iv = collections.defaultdict(list)
for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
    _, t, p, L = r["interval"].split(":")
    ev = r["evidence"].split(",")[0]
    iv[r["sample"]].append((int(p), int(L), r["state"], ev))
for v in iv.values():
    v.sort()


def h_state(s, p, L):
    v = iv.get(s, [])
    i = bisect.bisect_left(v, (p - 50,))
    while i < len(v) and v[i][0] <= p + 50:
        q, l, st, ev = v[i]
        if min(l, L) / max(l, L) >= 0.8:
            return st, ev
        i += 1
    return "not_catalogued", ""


def group(st, ev):
    if st == "not_catalogued":
        return "not catalogued"
    if "inherited" not in ev:
        return f"{st} (non-inherited evidence)"
    return st


calls = ["carries", "lacks", "unclear", "low_depth", "unplaced"]
for part in ("all", "outside the mask", "inside the mask"):
    tab, cols = collections.Counter(), set()
    for r in rows:
        if r["svtype"] != "DEL":
            continue
        p, L = int(r["h37rv_pos"]), int(r["svlen"])
        m = masked(p, p + L)
        if part == "outside the mask" and m or part == "inside the mask" and not m:
            continue
        g = group(*h_state(r["sample"], p, L))
        tab[r["call"], g] += 1; cols.add(g)
    cols = sorted(cols)
    print(f"\nInherited DELETIONS ({part}): matched-genome call (rows) x H37Rv-frame state in the cohort VCF (columns)")
    print(f"{'':<12}" + "".join(f"{c[:26]:>28}" for c in cols))
    for c in calls:
        if any(tab[c, x] for x in cols):
            print(f"{c:<12}" + "".join(f"{tab[c, x]:>28,}" for x in cols))

print("\nInherited INSERTIONS: matched-genome call by class (today all are reported present)")
tab = collections.Counter()
for r in rows:
    if r["svtype"] != "INS":
        continue
    p = int(r["h37rv_pos"])
    cl = "masked" if masked(p, p + 1) else ("tandem" if r.get("tandem") == "yes" else "other")
    tab[cl, r["call"]] += 1
for cl in ("other", "tandem", "masked"):
    n = sum(tab[cl, c] for c in calls + ["repeat"])
    print(f"  {cl:<7} {n:>6,}: " + ", ".join(f"{c} {tab[cl, c]:,}" for c in calls + ["repeat"] if tab[cl, c]))
