#!/usr/bin/env python3
"""Score check.py's matched-genome calls on truth7 against each isolate's own
finished assembly (truth7/truth/<sample>.truth.var, paftools call against
H37Rv). Analysis only.

A row is TRUE (the isolate carries the SV) when the assembly has a variant of
the same type within 200 bp whose size is within 50% (score_sv.py's rule);
FALSE when the position lies inside an `R` region (assembly aligned once,
so the absence is observed); UNKNOWN otherwise.

    score_truth.py out/truth7.tsv truth_dir > out/truth7_score.txt
"""
import bisect, collections, csv, os, sys

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
tdir = sys.argv[2]
mask = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(sys.argv[3]) if l.strip())
mst = [a for a, _ in mask]


def masked(p, e):
    i = bisect.bisect_right(mst, e) - 1
    while i >= 0 and i > bisect.bisect_right(mst, p) - 50:
        if mask[i][0] < e + 1 and mask[i][1] > p:
            return True
        i -= 1
    return False
truth, regions = {}, {}
for s in {r["sample"] for r in rows}:
    v, R = [], []
    for line in open(os.path.join(tdir, f"{s}.truth.var")):
        f = line.rstrip("\n").split("\t")
        if f[0] == "R":
            R.append((int(f[2]), int(f[3])))
        elif f[0] == "V":
            ref, alt = f[6], f[7]
            lr, la = (0 if ref == "-" else len(ref)), (0 if alt == "-" else len(alt))
            d = la - lr
            if abs(d) >= 25:
                v.append((int(f[2]), "INS" if d > 0 else "DEL", abs(d)))
    truth[s] = sorted(v); regions[s] = sorted(R)


def state(r):
    s, p, L, t = r["sample"], int(r["h37rv_pos"]), int(r["svlen"]), r["svtype"]
    for q, tt, tl in truth[s]:
        if tt == t and abs(q - p) <= 200 and min(L, tl) / max(L, tl) >= 0.5:
            return "TRUE"
    R = regions[s]
    i = bisect.bisect_right([a for a, _ in R], p) - 1
    end = p if t == "INS" else p + L
    if i >= 0 and R[i][0] <= p and R[i][1] >= end:
        return "FALSE"
    return "UNKNOWN"


tab = collections.Counter()
part = sys.argv[4] if len(sys.argv) > 4 else "all"
for r in rows:
    p = int(r["h37rv_pos"]); e = int(r["h37rv_end"])
    m = masked(p, max(e, p + 1))
    if part == "outside" and m or part == "inside" and not m:
        continue
    if part.startswith("tandem") and (m or r.get("tandem") != ("yes" if part == "tandem" else "no")):
        continue
    tab[r["svtype"], r["call"], state(r)] += 1
calls = ["carries", "lacks", "unclear", "repeat", "low_depth", "unplaced"]
print(f"inherited events in the 7 truth isolates ({part}; 'tandem'/'tandem_no' are outside the mask): matched-genome call against the assembly")
for t in ("DEL", "INS"):
    print(f"\n{t}            truly carried  truly absent  unknown")
    for c in calls:
        n = [tab[t, c, x] for x in ("TRUE", "FALSE", "UNKNOWN")]
        if sum(n):
            print(f"  {c:<12} {n[0]:>13} {n[1]:>13} {n[2]:>8}")
    tp, fp = tab[t, "carries", "TRUE"], tab[t, "carries", "FALSE"]
    tn, fn = tab[t, "lacks", "FALSE"], tab[t, "lacks", "TRUE"]
    allT = sum(tab[t, c, "TRUE"] for c in calls)
    allF = sum(tab[t, c, "FALSE"] for c in calls)
    print(f"  present as reported today (no check): {allT} right, {allF} wrong "
          f"(precision {allT / max(1, allT + allF):.3f})")
    print(f"  'carries' precision {tp / max(1, tp + fp):.3f} ({tp}/{tp + fp}); "
          f"'lacks' correct {tn}/{tn + fn}; truly absent caught {tn}/{allF}")
