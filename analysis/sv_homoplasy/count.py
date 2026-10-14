#!/usr/bin/env python3
"""How many graph SVs (>= 50 bp) arose more than once on the panel tree?
Analysis only (user, 2026-10-14). P200 build 3893a3762ed1.

Fitch small parsimony on the rooted panel tree (runroot/panel/P200/trees/
p200.rooted.nwk) for every biallelic SV record of the collapsed graph VCF;
H37Rv is REF; kit genomes (L8, canettii) have no SV genotype here and are
pruned. The number of state changes is the minimum number of events; with
the ancestral state from the build's parsimony table (panel_polarity.tsv,
ancestral at the L1-4,7 MRCA), changes are split into gains and losses of
the derived allele. "Evolved multiple times" = >= 2 changes. Reported for all
SVs and for those outside the build's repeat mask.
"""
import collections, gzip, sys
import numpy as np
sys.path.insert(0, "bin")
from panel_ancestral import parse_newick
B = "runroot/refbias/build/3893a3762ed1"
T = "runroot/panel/P200/trees/p200.rooted.nwk"
children, root, leaves = parse_newick(open(T).read())
mask = np.zeros(4411600, bool)
for l in open(f"{B}/assets/repeat_mask.bed"):
    c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True
fh = gzip.open(f"{B}/assets/graph_collapsed.vcf.gz", "rt")
for l in fh:
    if l.startswith("#CHROM"): samples = l.rstrip("\n").split("\t")[9:]; break
keep = set(samples) | {"GCF_000195955"}
def prune(n):
    if n not in children: return n if n in keep else None
    k = [x for x in (prune(c) for c in children[n]) if x is not None]
    if not k: return None
    if len(k) == 1: return k[0]
    ch[n] = k; return n
ch = {}; rt = prune(root)
post = []
st = [(rt, 0)]
while st:
    n, d = st.pop()
    if d: post.append(n); continue
    st.append((n, 1))
    for c in ch.get(n, []): st.append((c, 0))
rows = []
for l in fh:
    f = l.rstrip("\n").split("\t")
    ref, alt = f[3], f[4]
    if "," in alt or max(len(ref), len(alt)) - min(len(ref), len(alt)) < 50 and max(len(ref), len(alt)) < 50: continue
    if abs(len(alt) - len(ref)) < 50 and not (len(ref) >= 50 and len(alt) >= 50): continue
    g = {s: (int(x.split(":")[0]) if x.split(":")[0] in "01" else None) for s, x in zip(samples, f[9:])}
    g["GCF_000195955"] = 0
    pos = int(f[1])
    rows.append((pos, ref, alt, g))
pol = {}
for l in open(f"{B}/assets/panel_polarity.tsv"):
    x = l.split("\t")
    if x[0] == "chrom": continue
    pol[(int(x[1]), x[2], x[3])] = x[5]
out = collections.Counter(); hist = collections.Counter(); per = []
for pos, ref, alt, g in rows:
    sets, cost = {}, 0
    for n in post:
        if n not in ch:
            v = g.get(n); sets[n] = {0, 1} if v is None else {v}; continue
        a, b = [sets[c] for c in ch[n]][:2] if len(ch[n]) == 2 else (None, None)
        ss = [sets[c] for c in ch[n]]
        inter = set.intersection(*ss)
        if inter: sets[n] = inter
        else:
            # multifurcation: count states, keep the most frequent, cost = k - max
            cnt = collections.Counter(s for x in ss for s in x)
            m = max(cnt.values()); sets[n] = {s for s, v in cnt.items() if v == m}; cost += len(ss) - m
    carriers = sum(1 for v in g.values() if v == 1)
    kind = "INS" if len(alt) > len(ref) else "DEL"
    masked = bool(mask[pos - 1:pos - 1 + max(len(ref), 1)].any())
    anc = pol.get((pos, ref, alt), "?")
    per.append((pos, kind, abs(len(alt) - len(ref)), carriers, cost, masked, anc))
    hist[min(cost, 6)] += 1
import csv
with open("analysis/sv_homoplasy/sv_changes.tsv", "w") as fo:
    fo.write("pos\ttype\tsize\tcarriers\tmin_changes\tmasked\tancestral\n")
    for r in per: fo.write("\t".join(map(str, r)) + "\n")
n = len(per)
for lab, sub in (("all", per), ("outside mask", [r for r in per if not r[5]])):
    m = [r for r in sub if r[4] >= 2]
    print(f"{lab}: {len(sub):,} SVs; >= 2 changes {len(m):,} ({len(m)/max(1,len(sub)):.0%}); >= 3 {sum(r[4]>=3 for r in sub):,}; >= 5 {sum(r[4]>=5 for r in sub):,}")
    print("   by type:", dict(collections.Counter(r[1] for r in m)), "  ancestral known (REF/ALT):", sum(r[6] in ("REF","ALT") for r in m))
print("changes histogram (6 = 6+):", dict(sorted(hist.items())))
top = sorted([r for r in per if not r[5]], key=lambda r: -r[4])[:15]
print("most homoplastic outside the mask: pos type size carriers changes ancestral")
for r in top: print("  ", r[0], r[1], r[2], r[3], r[4], r[6])
