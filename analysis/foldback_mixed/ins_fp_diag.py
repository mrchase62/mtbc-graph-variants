#!/usr/bin/env python3
"""The positional SV insertion calls with no assembly insertion within 50 bp
(ctrl25_p200; user, 2026-10-14: "look into the positional SV insertion false
calls"). For each: source arm, component, length, region, nearest truth
insertion, recurrence across controls, and whether the inserted sequence
(where the record carries it) is in the control's own assembly."""
import bisect, collections, csv, gzip, sys
sys.path.insert(0, "."); import grid_score as gs
V = "../../runroot/refbias/ctrl25_p200/p5/merged.vcf.gz"
fb = {r[0]: float(r[3]) for r in csv.reader(open(gs.RATES), delimiter="\t") if r[0] != "sample"}
tf = {r["sample"]: r["file"] for r in csv.DictReader(open(gs.TRUTH), delimiter="\t")}
fh = gzip.open(V, "rt")
for l in fh:
    if l.startswith("#CHROM"): S = l.rstrip().split("\t")[9:]; break
recs = []
for l in fh:
    f = l.rstrip("\n").split("\t")
    info = dict(kv.split("=", 1) if "=" in kv else (kv, 1) for kv in f[7].split(";"))
    if not (info.get("CLASS") == "sv" and info.get("SVTYPE") == "INS"): continue
    L = abs(int(info.get("SVLEN", 0)))
    if L < 50: continue
    car = [S[j] for j, x in enumerate(f[9:]) if x.split(":")[0] == "1"]
    recs.append((int(f[1]), f[3], f[4], L, info, car))
seqs = {}
def asm(s):
    if s not in seqs:
        a = "".join(l.strip().upper() for l in gzip.open(f"{gs.ASM}/{tf[s]}", "rt") if l[0] != ">")
        seqs[s] = (a, a.translate(str.maketrans("ACGT", "TGCA"))[::-1])
    return seqs[s]
rows = []
for s in S:
    cal, tv = gs.truth(s, "GCF_000195955", f"{gs.FM}/controls_grid/truth")
    cs = sorted(cal); st = [x[0] for x in cs]
    def ok(p):
        i = bisect.bisect_right(st, p - 1) - 1
        return i >= 0 and cs[i][0] <= p - 1 < cs[i][1] and not gs.mask[p - 1]
    tins = [(p, len(a) - len(r)) for p, r, a in tv if len(a) - len(r) >= 1]
    for pos, ref, alt, L, info, car in recs:
        if s not in car or not ok(pos): continue
        if any(abs(p - pos) <= 50 and tl >= 50 for p, tl in tins): continue
        near = min((abs(p - pos), tl) for p, tl in tins) if tins else (10**9, 0)
        seq_in = ""
        if not alt.startswith("<"):
            ins = alt[len(ref):] if len(alt) > len(ref) else alt
            a, ar = asm(s)
            k = ins[:min(len(ins), 60)]
            seq_in = "yes" if (k in a or k in ar) else "no"
        rows.append(dict(sample=s, group="chimeric" if fb.get(s, 0) > 4.8 else "clean", pos=pos, len=L,
                         source=info.get("SVSOURCE", info.get("SOURCES", "")), component=info.get("COMPONENT", ""),
                         region=info.get("REGION", ""), carriers=len(car), nearest_truth_ins_bp=near[0],
                         nearest_truth_ins_len=near[1], seq_in_assembly=seq_in))
with open("out/ins_fp_diag.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print("unmatched positional SV insertion calls:", len(rows))
for k in ("group", "source", "component", "region", "seq_in_assembly"):
    print(f"  by {k}:", dict(collections.Counter(r[k] for r in rows).most_common()))
d = collections.Counter("small insertion (<50) within 50 bp" if r["nearest_truth_ins_bp"] <= 50 else
                        "insertion 51-500 bp away" if r["nearest_truth_ins_bp"] <= 500 else "nothing within 500 bp" for r in rows)
print("  nearest truth insertion:", dict(d))
loc = collections.Counter(r["pos"] for r in rows)
print("  distinct positions:", len(loc), "; in >= 5 controls:", sum(v >= 5 for v in loc.values()),
      "; calls at those:", sum(v for v in loc.values() if v >= 5))
print("  top recurrent:", [(p, v, next(r["len"] for r in rows if r["pos"] == p)) for p, v in loc.most_common(10)])
