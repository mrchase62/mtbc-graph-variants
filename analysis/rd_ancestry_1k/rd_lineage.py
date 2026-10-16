#!/usr/bin/env python3
"""Step 1. Known RDs (MtbPangenome known_RDs.tiered.tsv) against the new 1K VCF.

For each RD: H37Rv-frame deletion records (CLASS sv / is6110, SVTYPE DEL,
any COMPONENT incl. small_allele) with reciprocal overlap >= 0.5 of the
deleted bases are its matches; the best match (highest RO) is the primary
record. With no match, the record with the highest RO (or nearest) is shown.
Expected carriers: isolates whose TB-Profiler deepest sublineage is the RD's
lineage label or below it (lineage1.2.1 covers lineage1.2.1.x). Isolates
typed only above the RD's level (e.g. lineage4 for an RD of 4.3.3) are
'unresolved'. Labels naming no lineage 1-4 group (L5/6/9, animal, canettii)
expect no carriers.
Also: Fitch on the cohort tree (rooted L1 | L2-4) for the primary record.
Outputs out/rd_records.tsv, out/rd_by_lineage.tsv, out/rd_by_sublineage.tsv,
out/rd_exceptions.tsv.
"""
import collections, re
import numpy as np
from common import *

d = load_records(); recs, st, ev, evn = d["recs"], d["st"], d["ev"], d["evnames"]
samples = d["samples"]
lin, deep = cohort()
L = np.array([lin[s] for s in samples]); D = [deep[s] for s in samples]

dels = [i for i, r in enumerate(recs) if r["chrom"] == "NC_000962.3" and r["svtype"] == "DEL"
        and r["cls"] in ("sv", "is6110", "small")]
span = {i: (recs[i]["pos"], recs[i]["end"]) for i in dels}   # 0-based [pos, end) = deleted bases


def ro(a, b, c, e):
    ov = max(0, min(b, e) - max(a, c))
    return min(ov / max(1, b - a), ov / max(1, e - c))


def parse_label(lab):
    out = []
    for tok in lab.split("\\"):
        m = re.match(r"\s*(\d+(?:\.\d+)*)", tok)
        if m and m.group(1).split(".")[0] in "1234" and len(m.group(1).split(".")[0]) == 1:
            out.append(m.group(1))
    return out


def expected(lab_list, dp, excl=()):
    """1 expected carrier, 0 not, -1 unresolved."""
    x = dp.replace("lineage", "")
    if any(x == e or x.startswith(e + ".") for e in excl): return 0
    res = 0
    for lab in lab_list:
        if x == lab or x.startswith(lab + "."): return 1
        if lab.startswith(x + "."): res = -1
    return res


adj, _ = read_unrooted()
ch, root, post, rinfo = root_on_split(adj, [s for s in samples if lin[s] == "lineage1"])
F = Fitch(ch, root, post, samples)

rds = []
for l in open(f"{ANN}/known_RDs.tiered.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] == "chrom": continue
    rds.append(dict(start=int(f[1]), end=int(f[2]), rd=f[3], lab=f[4], tier=f[5]))

rows, byl, bysub, exc = [], [], [], []
prim_idx = []
for rd in rds:
    c = [(ro(rd["start"], rd["end"], *span[i]), i) for i in dels]
    c = [x for x in c if x[0] > 0]
    m = sorted([x for x in c if x[0] >= 0.5], reverse=True)
    if m:
        prim = m[0][1]; prim_ro = m[0][0]
    elif c:
        prim_ro, prim = max(c)
    else:
        # nearest by midpoint
        mid = (rd["start"] + rd["end"]) / 2
        prim = min(dels, key=lambda i: abs((span[i][0] + span[i][1]) / 2 - mid)); prim_ro = 0.0
    excl = re.findall(r"except (\d+(?:\.\d+)*)", rd["lab"])
    labs = parse_label(re.sub(r"\(except[^)]*\)", "", rd["lab"]))
    exp = np.array([expected(labs, dp, excl) for dp in D])
    inside = [i for i in range(len(recs)) if recs[i]["chrom"] == "NC_000962.3" and klass(recs[i]) is not None
              and recs[i]["pos"] >= rd["start"] and recs[i]["end"] <= rd["end"]]
    inside = [i for i in inside if (st[i] == 2).mean() < 0.5]   # skip whole-site ABSENT records
    absin = (st[inside] == 2).any(0) if inside else np.zeros(len(samples), bool)
    s = st[prim]
    r = recs[prim]; info = r["info"]
    union = np.full(len(samples), 3)
    if m:
        mm = st[[i for _, i in m]]
        union = np.where((mm == 1).any(0), 1, np.where((mm == 0).any(0), 0, np.where((mm == 2).any(0), 2, 3)))
    rows.append(dict(rd=rd["rd"], rd_start=rd["start"], rd_end=rd["end"], rd_len=rd["end"] - rd["start"],
                     label=rd["lab"], tier=rd["tier"], n_match=len(m),
                     match_pos=r["pos"], match_end=r["end"], match_len=r["end"] - r["pos"], ro=round(prim_ro, 3),
                     component=component(r), svsource=info.get("SVSOURCE", "."), svtier=info.get("SVTIER", "."),
                     is6110prox=int("IS6110PROX" in info), panelcarriers=info.get("PANELCARRIERS", "."),
                     aa_state=info.get("AA_STATE", "."), aa_flag=info.get("AA_FLAG", "."),
                     aa_changes=info.get("AA_CHANGES", "."), aa_source=info.get("AA_SOURCE", "."),
                     exp_n=int((exp == 1).sum()), unres_n=int((exp == -1).sum()),
                     exp_alt=int(((exp == 1) & (s == 1)).sum()), exp_ref=int(((exp == 1) & (s == 0)).sum()),
                     exp_abs=int(((exp == 1) & (s == 2)).sum()), exp_nc=int(((exp == 1) & (s >= 3)).sum()),
                     out_alt=int(((exp == 0) & (s == 1)).sum()), out_ref=int(((exp == 0) & (s == 0)).sum()),
                     out_abs=int(((exp == 0) & (s == 2)).sum()), out_nc=int(((exp == 0) & (s >= 3)).sum()),
                     unres_alt=int(((exp == -1) & (s == 1)).sum()),
                     absent_inside_exp=int((absin & (exp == 1)).sum()), absent_inside_out=int((absin & (exp != 1)).sum()),
                     union_exp_alt=int(((exp == 1) & (union == 1)).sum()),
                     union_out_alt=int(((exp == 0) & (union == 1)).sum()),
                     idx=prim))
    prim_idx.append(prim)
    for g in sorted(set(L)):
        k = L == g
        byl.append((rd["rd"], g, *(int((k & (s == v)).sum()) for v in (1, 0, 2)), int((k & (s >= 3)).sum())))
    for g in sorted(set(D)):
        k = np.array([x == g for x in D])
        cnt = [int((k & (s == v)).sum()) for v in (1, 0, 2)] + [int((k & (s >= 3)).sum())]
        if cnt[0] or cnt[2] or (exp[k] == 1).any():
            bysub.append((rd["rd"], g, int(exp[k][0]), *cnt))
    if m:
        for j in np.where(((exp == 1) & (s != 1)) | ((exp == 0) & (s == 1)))[0]:
            exc.append((rd["rd"], rd["lab"], samples[j], L[j], D[j], "expected_ALT" if exp[j] == 1 else "expected_REF",
                        ST[int(s[j])], evn.get(int(ev[prim, j]), "?"), ST[int(union[j])]))

# Fitch on the cohort tree for primary records
res = F.run(st[prim_idx])
for k, row in enumerate(rows):
    row.update(cohort_changes=int(res["cost"][k]), cohort_gains=int(res["gains"][k]),
               cohort_losses=int(res["losses"][k]),
               cohort_root=("TIED" if res["root_tied"][k] else ("ALT" if res["root_state"][k] else "REF")))

below = F.leaves_below()
col = {s: j for j, s in enumerate(samples)}
def change_edges(k):
    out = []
    stack = [root]
    while stack:
        p = stack.pop()
        for c in ch.get(p, []):
            a, b = F.state[p][k], F.state[c][k]
            if a != b:
                lv = below[c]
                comp = collections.Counter(deep[x].replace("lineage", "L") for x in lv).most_common(2)
                out.append(f"{'+' if b else '-'}{len(lv)}[" + ",".join(f"{g}:{n}" for g, n in comp) + "]")
            stack.append(c)
    return " ".join(sorted(out, key=lambda x: -int(re.match(r"[+-](\d+)", x).group(1))))
for k, row in enumerate(rows):
    row["change_edges"] = change_edges(k)
cols = list(rows[0].keys())
with open(f"{OUT}/rd_records.tsv", "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for r in rows: fo.write("\t".join(str(r[c]) for c in cols) + "\n")
with open(f"{OUT}/rd_by_lineage.tsv", "w") as fo:
    fo.write("rd\tlineage\tALT\tREF\tABSENT\tNOCALL\n")
    for r in byl: fo.write("\t".join(map(str, r)) + "\n")
with open(f"{OUT}/rd_by_sublineage.tsv", "w") as fo:
    fo.write("rd\tsublineage\texpected(1 yes,0 no,-1 unresolved)\tALT\tREF\tABSENT\tNOCALL\n")
    for r in bysub: fo.write("\t".join(map(str, r)) + "\n")
with open(f"{OUT}/rd_exceptions.tsv", "w") as fo:
    fo.write("rd\tlabel\tsample\tlineage\tsublineage\tdirection\tST\tEV\tunion_of_matches\n")
    for r in exc: fo.write("\t".join(map(str, r)) + "\n")
print("root:", rinfo)
print(f"{'rd':14s} {'lab':10s} tier  nm  ro   comp   expN  eA eR eAb eNC | oA oAb | chg g l root | AA flag chg")
for r in rows:
    if r["exp_n"] == 0 and r["out_alt"] == 0 and r["unres_alt"] == 0: continue
    print(f"{r['rd']:14s} {r['label'][:10]:10s} {r['tier'][:4]} {r['n_match']:2d} {r['ro']:.2f} {r['component'][:6]:6s} {r['exp_n']:4d} {r['exp_alt']:3d} {r['exp_ref']:2d} {r['exp_abs']:2d} {r['exp_nc']:3d} | {r['out_alt']:3d} {r['out_abs']:3d} | {r['cohort_changes']:3d} {r['cohort_gains']} {r['cohort_losses']} {r['cohort_root']} | {r['aa_state']} {r['aa_flag']} {r['aa_changes']} unres_alt={r['unres_alt']}")
