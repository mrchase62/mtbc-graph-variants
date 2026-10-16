#!/usr/bin/env python3
"""Step 3. Ancestral reconstruction of SV records on the rooted cohort tree.

Tree: new1k.fasttree.nwk rooted on the branch separating lineage 1 (250
isolates, monophyletic) from lineages 2-4 (see common.root_on_split).
Records: H37Rv-frame DEL and INS (CLASS sv, any COMPONENT incl. small_allele),
plus ACC (accessory presence) and IS6110 insertions for context; scored when
>= 2 ALT and >= 2 REF isolates. Leaves ALT 1, REF 0, other states unknown.
Fitch up pass gives the change count and the root set; a root set {0,1} is
TIED (the root sits between L1 and L2-4, so any L1-vs-rest pattern ties).
Down pass from the root state (REF when tied) counts changes away from the
root state ('derived') and back to it ('reversals').
The cohort root state is compared with INFO/AA_STATE (panel parsimony at the
L1-4,7 ancestor, outgroups L5/L6 used for ties).
Outputs out/ancestry_records.tsv, out/ancestry_agreement.tsv,
out/homoplastic_top.tsv.
"""
import bisect, collections
import numpy as np
from common import *

d = load_records(); recs, st = d["recs"], d["st"]; samples = d["samples"]
lin, deep = cohort(); mask = mask_array(); G = genes(); IS = is6110_h37rv()
gstarts = [g[0] for g in G]
adj, _ = read_unrooted()
ch, root, post, rinfo = root_on_split(adj, [s for s in samples if lin[s] == "lineage1"])
F = Fitch(ch, root, post, samples)

K = [klass(r) for r in recs]
sel = [i for i, r in enumerate(recs) if K[i] in ("DEL", "INS", "ACC", "IS6110")
       and (st[i] == 1).sum() >= 2 and (st[i] == 0).sum() >= 2]
S = st[sel]
res = F.run(S)
# derived vs reversal relative to the root state
der = np.where(res["root_state"] == 0, res["gains"], res["losses"])
rev = np.where(res["root_state"] == 0, res["losses"], res["gains"])
below = F.leaves_below()

# cohort IS6110 insertion sites carried by >= 1 isolate
iss = sorted(r["pos"] for i, r in enumerate(recs) if K[i] == "IS6110" and (st[i] == 1).any())


def near(sorted_pos, a, b):
    k = bisect.bisect_left(sorted_pos, a)
    best = 10 ** 9
    for j in (k - 1, k):
        if 0 <= j < len(sorted_pos):
            p = sorted_pos[j]
            best = min(best, 0 if a <= p <= b else min(abs(p - a), abs(p - b)))
    return best


def is_dist(a, b):
    best = 10 ** 9
    for s, e in IS:
        best = min(best, 0 if (s <= b and e >= a) else min(abs(s - b), abs(a - e)))
    return best


def gene_at(a, b):
    k = bisect.bisect_right(gstarts, b)
    hit = [g for g in G[max(0, k - 60):k] if g[1] >= a and g[0] <= b]
    if hit:
        names = [g[3] if g[3] and g[3] != g[2] else g[2] for g in hit]
        return ",".join(names[:4]) + (f",+{len(names)-4}" if len(names) > 4 else "")
    left = G[k - 1] if k else None
    return f"intergenic(after {left[2]})" if left else "intergenic"


def edges(k):
    out = []
    stack = [root]
    while stack:
        p = stack.pop()
        for c in ch.get(p, []):
            if F.state[p][k] != F.state[c][k]:
                lv = below[c]
                comp = collections.Counter(deep[x].replace("lineage", "L") for x in lv).most_common(1)
                out.append((len(lv), f"{'+' if F.state[c][k] else '-'}{len(lv)}[{comp[0][0]}]"))
            stack.append(c)
    return " ".join(x[1] for x in sorted(out, reverse=True)[:12])


rows = []
for k, i in enumerate(sel):
    r = recs[i]; inf = r["info"]
    a, b = r["pos"], max(r["end"], r["pos"])
    inm = bool(mask[max(0, a - 1):b].any())
    rs = "TIED" if res["root_tied"][k] else ("ALT" if res["root_state"][k] else "REF")
    rows.append(dict(idx=i, cls=K[i], component=component(r), pos=a, end=r["end"], svlen=r["svlen"],
                     in_mask=int(inm), alt=int((S[k] == 1).sum()), ref=int((S[k] == 0).sum()),
                     nocall=int((S[k] >= 2).sum()), changes=int(res["cost"][k]),
                     root=rs, derived=int(der[k]), reversals=int(rev[k]),
                     aa_state=inf.get("AA_STATE", "."), aa_flag=inf.get("AA_FLAG", "."),
                     aa_changes=inf.get("AA_CHANGES", "."), aa_source=inf.get("AA_SOURCE", "."),
                     is6110prox=int("IS6110PROX" in inf), lowconf=int("LOWCONF" in inf),
                     svtier=inf.get("SVTIER", "."), panelcarriers=inf.get("PANELCARRIERS", "."),
                     dist_h37rv_is6110=is_dist(a, b), dist_cohort_is6110=near(iss, a, b),
                     gene=gene_at(a, b)))
cols = list(rows[0])
with open(f"{OUT}/ancestry_records.tsv", "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for r in rows: fo.write("\t".join(str(r[c]) for c in cols) + "\n")

# agreement of cohort root with panel AA_STATE
agr = collections.Counter()
for r in rows:
    flag = "none" if r["aa_state"] == "." else (r["aa_flag"].split(":")[0] if r["aa_flag"] != "." else "clean")
    key = (r["cls"], r["component"], flag, "mask" if r["in_mask"] else "outside",
           "nocall>=50%" if r["nocall"] >= 500 else "nocall<50%")
    agr[key + ("n",)] += 1
    if r["aa_state"] == ".": continue
    if r["root"] == "TIED": agr[key + ("tied",)] += 1
    elif r["root"] == r["aa_state"]: agr[key + ("agree",)] += 1
    else: agr[key + ("disagree",)] += 1
keys = sorted({k[:5] for k in agr})
with open(f"{OUT}/ancestry_agreement.tsv", "w") as fo:
    fo.write("class\tcomponent\taa_flag\tmask\tnocall\trecords\tcohort_root_tied\tagree\tdisagree\n")
    for k in keys:
        fo.write("\t".join(k) + "\t" + "\t".join(str(agr[k + (x,)]) for x in ("n", "tied", "agree", "disagree")) + "\n")

top = sorted([(r, k) for k, r in enumerate(rows) if not r["in_mask"] and r["cls"] in ("DEL", "INS")],
             key=lambda x: (-x[0]["changes"], -x[0]["alt"]))[:40]
with open(f"{OUT}/homoplastic_top.tsv", "w") as fo:
    c2 = ["cls", "component", "pos", "end", "svlen", "alt", "ref", "nocall", "changes", "derived", "reversals", "root",
          "aa_state", "aa_flag", "is6110prox", "dist_h37rv_is6110", "dist_cohort_is6110", "gene"]
    fo.write("\t".join(c2) + "\tchange_edges\n")
    for r, k in top:
        fo.write("\t".join(str(r[c]) for c in c2) + "\t" + edges(k) + "\n")
print("root:", rinfo, "records scored:", len(rows))
