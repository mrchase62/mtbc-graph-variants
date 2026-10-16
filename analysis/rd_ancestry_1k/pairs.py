#!/usr/bin/env python3
"""Step 2. SV / IS6110 / accessory discordance between near-identical isolates.

SNP distance: analysis/sv_eval_1k/out/snps.fasta (build_tree.py: biallelic
H37Rv-frame SNPs of the new 1K VCF outside the build repeat mask, variable in
the cohort; GT . -> N). Pairwise distance = sites where both isolates have a
base and the bases differ (N ignored; the number of sites compared is kept).
Pairs with distance <= 5 (tier 1) and 6-12 (tier 2).

For every pair and every SV-like record (classes of common.klass), a cell pair
is 'both called' when both ST are REF or ALT, and 'discordant' when one is ALT
and the other REF. ABSENT-involving differences (ABSENT vs REF/ALT) are counted
separately. Broken down by class, COMPONENT, repeat mask, EV of the two cells,
and whether the two isolates were mapped to the same matched reference (P1
refmap). Known-RD records (RO >= 0.5 to an RD of known_RDs.tiered.tsv) are
counted on their own.
Outputs: out/pair_distances.tsv (all pairs <= 30 SNPs), out/pair_discordance.tsv
(per pair x class), out/pair_discordant_cells.tsv.gz, out/pair_summary.txt.
"""
import collections, gzip, itertools
import numpy as np
from common import *

d = load_records(); recs, st, ev, evn = d["recs"], d["st"], d["ev"], d["evnames"]
samples = d["samples"]; n = len(samples)
lin, deep = cohort(); rmap = refmap(); mask = mask_array()

# ---- SNP distances
names, seqs = [], []
for l in open(f"{R}/analysis/sv_eval_1k/out/snps.fasta"):
    if l.startswith(">"): names.append(l[1:].strip())
    else: seqs.append(np.frombuffer(l.strip().encode(), np.uint8))
assert names == samples, "snps.fasta order differs from the VCF"
X = np.vstack(seqs); del seqs
L = X.shape[1]
both = np.zeros((n, n), np.float64); same = np.zeros((n, n), np.float64)
for a in range(0, L, 4000):
    blk = X[:, a:a + 4000]
    called = (blk != ord("N")).astype(np.float32)
    both += called @ called.T
    for b in b"ACGT":
        m = (blk == b).astype(np.float32); same += m @ m.T
dist = np.rint(both - same).astype(int); comp = np.rint(both).astype(int)
iu = np.triu_indices(n, 1)
with open(f"{OUT}/pair_distances.tsv", "w") as fo:
    fo.write("a\tb\tsnps\tsites_compared\tlineage_a\tsub_a\tsub_b\tref_a\tref_b\n")
    for i, j in zip(*iu):
        if dist[i, j] <= 30:
            a, b = samples[i], samples[j]
            fo.write(f"{a}\t{b}\t{dist[i,j]}\t{comp[i,j]}\t{lin[a]}\t{deep[a]}\t{deep[b]}\t{rmap.get(a,'.')}\t{rmap.get(b,'.')}\n")
print("SNP columns", L, "median sites compared", int(np.median(comp[iu])))
pairs = [(i, j, int(dist[i, j])) for i, j in zip(*iu) if dist[i, j] <= 12]

# ---- records
K = np.array([klass(r) or "" for r in recs])
COMP = np.array([component(r) for r in recs])
INM = np.array([bool(mask[max(0, r["pos"] - 1):max(r["end"], r["pos"])].any()) if r["chrom"] == "NC_000962.3" else False
                for r in recs])
use = np.where(K != "")[0]
rds = [l.rstrip("\n").split("\t") for l in open(f"{ANN}/known_RDs.tiered.tsv")][1:]
isrd = np.zeros(len(recs), bool); rdname = {}
for f in rds:
    s0, e0 = int(f[1]), int(f[2])
    for i in use:
        r = recs[i]
        if r["chrom"] != "NC_000962.3" or r["svtype"] != "DEL": continue
        ov = max(0, min(e0, r["end"]) - max(s0, r["pos"]))
        if ov and min(ov / (e0 - s0), ov / max(1, r["end"] - r["pos"])) >= 0.5:
            isrd[i] = True; rdname.setdefault(i, []).append(f[3])
print("RD-matched records:", int(isrd.sum()))

sub = st[use]; evs = ev[use]
Ku, Cu, Mu, RDu = K[use], COMP[use], INM[use], isrd[use]
rows = []
cells = gzip.open(f"{OUT}/pair_discordant_cells.tsv.gz", "wt")
cells.write("a\tb\tsnps\tsame_ref\tidx\tclass\tcomponent\tin_mask\tknown_rd\tpos\tend\tsvlen\tst_a\tst_b\tev_a\tev_b\n")
agg = collections.Counter(); aggc = collections.Counter(); evpair = collections.Counter()
for i, j, dd in pairs:
    a, b = samples[i], samples[j]
    sr = rmap.get(a) == rmap.get(b)
    x, y = sub[:, i], sub[:, j]
    called = (x <= 1) & (y <= 1)
    disc = called & (x != y)
    absd = ((x == 2) & (y <= 1)) | ((y == 2) & (x <= 1))
    tier = "le5" if dd <= 5 else "6to12"
    for k in sorted(set(Ku)):
        m = Ku == k
        rows.append((a, b, dd, tier, int(sr), lin[a], k, int((called & m).sum()), int((disc & m).sum()),
                     int((absd & m).sum()), int((disc & m & RDu).sum()), int((called & m & RDu).sum())))
    for t in np.where(disc | absd)[0]:
        g = use[t]; r = recs[g]
        cells.write(f"{a}\t{b}\t{dd}\t{int(sr)}\t{g}\t{Ku[t]}\t{Cu[t]}\t{int(Mu[t])}\t{','.join(rdname.get(g, []))}\t{r['pos']}\t{r['end']}\t{r['svlen']}\t{ST[int(x[t])]}\t{ST[int(y[t])]}\t{evn.get(int(evs[t, i]))}\t{evn.get(int(evs[t, j]))}\n")
    for t in np.where(called)[0]:
        key = (tier, Ku[t], Cu[t], int(Mu[t]), int(sr))
        aggc[key] += 1
        if disc[t]:
            agg[key] += 1
            e1, e2 = evn.get(int(evs[t, i])), evn.get(int(evs[t, j]))
            ea, er = (e1, e2) if x[t] == 1 else (e2, e1)
            evpair[(tier, Ku[t], Cu[t], ea, er)] += 1
cells.close()
with open(f"{OUT}/pair_discordance.tsv", "w") as fo:
    fo.write("a\tb\tsnps\ttier\tsame_ref\tlineage\tclass\tboth_called\tdiscordant_alt_ref\tabsent_vs_called\trd_discordant\trd_both_called\n")
    for r in rows: fo.write("\t".join(map(str, r)) + "\n")
with open(f"{OUT}/pair_rates_by_component.tsv", "w") as fo:
    fo.write("tier\tclass\tcomponent\tin_mask\tsame_ref\tboth_called\tdiscordant\trate_per_1e4\n")
    for k in sorted(aggc):
        fo.write("\t".join(map(str, k)) + f"\t{aggc[k]}\t{agg[k]}\t{1e4*agg[k]/aggc[k]:.2f}\n")
with open(f"{OUT}/pair_discordant_ev.tsv", "w") as fo:
    fo.write("tier\tclass\tcomponent\tev_of_ALT_cell\tev_of_REF_cell\tcells\n")
    for k, v in evpair.most_common(): fo.write("\t".join(map(str, k)) + f"\t{v}\n")
print("pairs <=5:", sum(p[2] <= 5 for p in pairs), " 6-12:", sum(p[2] > 5 for p in pairs))
