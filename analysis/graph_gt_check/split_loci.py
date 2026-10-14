#!/usr/bin/env python3
"""What are the split SV loci? (user, 2026-10-14: "look into QC-1 and the
split loci"). Loci = SV keys (>= 50 bp, normalised and pooled as
panel_ancestral.read_panel) of one type within 50 bp of each other, with
disjoint carriers. Each locus is classed by its alleles:
  copy_number allele lengths are whole multiples of the shortest (within 3%
              or 3 bp), or the sequence is periodic: a VNTR / duplication
              copy-number series
  is6110      the allele contains >= 50% of IS6110 by 21-mers
  same_event_variants  allele lengths within 5% of each other: one event
              whose alleles differ inside (SNPs or small indels)
  other
plus masked (any allele touches the repeat mask). Writes out/split_loci.tsv.
"""
import collections, os, sys
import numpy as np
sys.path.insert(0, "bin"); import panel_ancestral as pa
B = "runroot/refbias/build/3893a3762ed1"
seq = "".join(l.strip().upper() for l in open(f"{B}/refs/GCF_000195955.fasta") if l[0] != ">")
mask = np.zeros(len(seq) + 200, bool)
for l in open(f"{B}/assets/repeat_mask.bed"):
    c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True
is_ = "".join(l.strip().upper() for l in open(os.environ["MTB_DATA"] + "/annotation/IS6110.fasta") if l[0] != ">")
RC = str.maketrans("ACGT", "TGCA")
isk = {is_[i:i + 21] for i in range(len(is_) - 20)} | {is_.translate(RC)[::-1][i:i + 21] for i in range(len(is_) - 20)}
samples, keys, copies, _ = pa.read_panel(f"{B}/assets/graph_collapsed.vcf.gz", os.environ["MTB_BCFTOOLS"], seq, "NC_000962.3")
sv = []
for k, c in zip(keys, copies):
    if pa.kind_of(k[2], k[3]) != "SV": continue
    col = pa.pooled_column(c)
    t = "INS" if len(k[3]) > len(k[2]) else "DEL"
    allele = k[3][len(k[2]):] if t == "INS" else k[2][len(k[3]):]
    if t == "INS" and not k[3].startswith(k[2][:1]): allele = k[3]
    sv.append((k[1], t, allele, {i for i, g in enumerate(col) if g == "1"}, len(k[2])))
def periodic(s):
    n = len(s)
    for p in range(1, min(200, n // 2) + 1):
        if sum(s[i] == s[i - p] for i in range(p, n)) >= 0.9 * (n - p): return p
    return 0
def isfrac(s):
    km = [s[i:i + 21] for i in range(max(0, len(s) - 20))]
    return sum(x in isk for x in km) / max(1, len(km))
out = []
for t in ("INS", "DEL"):
    x = sorted([r for r in sv if r[1] == t])
    grp, cur = [], []
    for r in x:
        if cur and r[0] - cur[-1][0] > 50: grp.append(cur); cur = []
        cur.append(r)
    if cur: grp.append(cur)
    for g in grp:
        if len(g) < 2: continue
        if any(g[i][3] & g[j][3] for i in range(len(g)) for j in range(i + 1, len(g))): continue
        lens = [len(r[2]) for r in g]
        per = [periodic(r[2]) for r in g]
        isf = max(isfrac(r[2]) for r in g)
        m = any(mask[r[0] - 1:r[0] + max(r[4], 1)].any() for r in g)
        u = min(lens)
        mult = all(abs(L / u - round(L / u)) * u <= max(3, 0.03 * L) for L in lens)
        if isf >= 0.5: cls = "is6110"
        elif max(lens) <= 1.05 * min(lens): cls = "same_event_variants"
        elif mult or all(p for p in per): cls = "copy_number"
        else: cls = "other"
        out.append((g[0][0], t, len(g), ",".join(map(str, lens)), sum(len(r[3]) for r in g), cls, int(m), max(per)))
with open("analysis/graph_gt_check/out/split_loci.tsv", "w") as fo:
    fo.write("pos\ttype\trecords\tallele_lengths\tcarriers\tclass\tmasked\tperiod\n")
    for r in out: fo.write("\t".join(map(str, r)) + "\n")
c = collections.Counter((r[5], r[6]) for r in out)
print(f"split loci: {len(out)}")
for cls in ("copy_number", "is6110", "same_event_variants", "other"):
    print(f"  {cls:12s} outside mask {c[(cls, 0)]:4d}  inside {c[(cls, 1)]:4d}")
print("records per locus:", dict(sorted(collections.Counter(min(r[2], 6) for r in out).items())))
for r in [r for r in out if r[6] == 0 and r[5] == "other"][:12]: print("  ", r)
