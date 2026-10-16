#!/usr/bin/env python3
"""OLD (gwas1000, graph 7713a8d71d8e) vs NEW (new1k_p200, graph 3893a3762ed1)
SV-like records on the 525 shared isolates. Analysis only.

Classes (H37Rv-frame records only):
  DEL  svtype DEL of CLASS sv (any COMPONENT/SVSOURCE), CLASS is6110 DEL
       (NEW copy losses) and CLASS small deletions >= 50 bp (OLD keeps these
       as small; NEW writes them CLASS=sv SVSOURCE=small_allele).
       Match: reciprocal overlap >= 50% of deleted bases [POS+1, END].
  INS  svtype INS of CLASS sv and small insertions >= 50 bp.
       Match: start within 50 bp, length within 20%.
  IS   CLASS is6110 INS. Match: same site within 20 bp.
  ACC  CLASS accessory_presence. Match: same ACCLOCUS, else start within
       50 bp and ACCLEN within 20%.
Pairs are one-to-one, greedy by ALT-carrier Jaccard over the shared isolates
among the geometric candidates (geometry breaks ties), so that where a locus
carries several alleles each OLD record pairs with the NEW record that best
represents it. "secondary" = a record whose
candidates were all taken by a better pair (typically split records).
Outputs out/pairs.tsv (pair-level), out/events.tsv, out/transitions.tsv,
out/new_change_ev.tsv, out/disagree_cells.tsv.gz.
"""
import collections, gzip, pickle
import numpy as np
D = "analysis/sv_eval_1k/out"
o = pickle.load(open(f"{D}/old.pkl", "rb")); n = pickle.load(open(f"{D}/new.pkl", "rb"))
shared = sorted(set(o["samples"]) & set(n["samples"]))
io = np.array([o["samples"].index(s) for s in shared]); inn = np.array([n["samples"].index(s) for s in shared])
SN = ["REF", "ALT", "ABSENT", "NOCALL", "DUP", "?", "?", "?", "?", "MISSING"]

def klass(r):
    if r["chrom"] != "NC_000962.3": return None
    c, t = r["cls"], r["svtype"]
    if c == "accessory_presence": return "ACC"
    if c == "is6110" and t == "INS": return "IS"
    if t == "DEL" and c in ("sv", "small", "is6110"): return "DEL"
    if t == "INS" and c in ("sv", "small"): return "INS"
    return None

def sub(r):
    i = r["info"]
    if r["cls"] == "small": return "small"
    if r["cls"] == "is6110" and r["svtype"] == "DEL": return "is6110_loss:" + i.get("COMPONENT", i.get("SVSOURCE", "."))
    if i.get("SVSOURCE") == "small_allele": return "small_allele"
    c = i.get("COMPONENT", ".")
    if c == "interval": return "interval_" + i.get("SVSOURCE", ".")
    return c

def score(k, a, b):
    if k == "DEL":
        s1, e1, s2, e2 = a["pos"] + 1, a["end"], b["pos"] + 1, b["end"]
        ov = min(e1, e2) - max(s1, s2) + 1
        if ov <= 0: return None
        l1, l2 = e1 - s1 + 1, e2 - s2 + 1
        if ov >= 0.5 * l1 and ov >= 0.5 * l2: return ov / max(l1, l2)
        return None
    if k == "IS":
        d = abs(a["pos"] - b["pos"])
        return 1 - d / 21 if d <= 20 else None
    if k == "ACC":
        if a["info"].get("ACCLOCUS") == b["info"].get("ACCLOCUS"): return 2.0
    d = abs(a["pos"] - b["pos"]); L1, L2 = a["svlen"], b["svlen"]
    if d <= 50 and L1 and L2 and min(L1, L2) >= 0.8 * max(L1, L2):
        return 1 - d / 102 + 0.5 * min(L1, L2) / max(L1, L2) - 0.5
    return None

def recs(d, idx):
    out = collections.defaultdict(list)
    for j, r in enumerate(d["recs"]):
        k = klass(r)
        if k: out[k].append(j)
    return out
RO, RN = recs(o, io), recs(n, inn)
ost, nst = o["st"][:, io], n["st"][:, inn]
nev = n["ev"][:, inn]
pairs_out = open(f"{D}/pairs.tsv", "w")
pairs_out.write("class\told_idx\tnew_idx\tpos_old\tpos_new\tlen_old\tlen_new\tsub_old\tsub_new\tscore\told_alt\tnew_alt\tboth_alt\told_nocall\tnew_nocall\n")
ev_rows = []
trans = collections.Counter(); chg_ev = collections.Counter()
cells = gzip.open(f"{D}/disagree_cells.tsv.gz", "wt")
cells.write("class\told_idx\tnew_idx\tsample\told_st\tnew_st\tnew_ev\tsub_old\tsub_new\tpos\tend\tsvlen\n")
for k in ["DEL", "INS", "IS", "ACC"]:
    A, B = RO[k], RN[k]
    A.sort(key=lambda j: o["recs"][j]["pos"]); B.sort(key=lambda j: n["recs"][j]["pos"])
    bp = np.array([n["recs"][j]["pos"] for j in B])
    cand = []
    for ja in A:
        a = o["recs"][ja]
        lo = a["pos"] - (60 if k != "DEL" else 0); hi = (a["end"] if k == "DEL" else a["pos"]) + 60
        if k == "ACC" or k == "DEL":
            lo = a["pos"] - 200000 if k == "DEL" else 0; hi = a["end"] + 1 if k == "DEL" else 5e6
        for t in range(np.searchsorted(bp, lo), np.searchsorted(bp, hi)):
            jb = B[t]; s = score(k, a, n["recs"][jb])
            if s is not None:
                x, y = ost[ja] == 1, nst[jb] == 1
                u = int((x | y).sum()); jac = (x & y).sum() / u if u else 0.0
                cand.append((jac + 0.01 * s, ja, jb))
    cand.sort(key=lambda x: -x[0])
    pa, pb = {}, {}
    anya, anyb = set(), set()
    for s, ja, jb in cand:
        anya.add(ja); anyb.add(jb)
        if ja in pa or jb in pb: continue
        pa[ja] = jb; pb[jb] = ja
    oalt = {j: int((ost[j] == 1).sum()) for j in A}; nalt = {j: int((nst[j] == 1).sum()) for j in B}
    for ja in A:
        st = "both" if ja in pa else ("secondary" if ja in anya else "only")
        ev_rows.append((k, "old", sub(o["recs"][ja]), st, oalt[ja] > 0))
    for jb in B:
        st = "both" if jb in pb else ("secondary" if jb in anyb else "only")
        ev_rows.append((k, "new", sub(n["recs"][jb]), st, nalt[jb] > 0))
    for ja, jb in pa.items():
        a, b = o["recs"][ja], n["recs"][jb]
        so, sn = ost[ja].astype(int), nst[jb].astype(int)
        sa, sb = sub(a), sub(b)
        for x, y, e in zip(so, sn, nev[jb]):
            trans[(k, SN[x], SN[y])] += 1
            if x != y:
                chg_ev[(k, sb, SN[x], SN[y], n["evnames"][int(e)])] += 1
        pairs_out.write(f"{k}\t{ja}\t{jb}\t{a['pos']}\t{b['pos']}\t{a['svlen']}\t{b['svlen']}\t{sa}\t{sb}\t"
                        f"{[s for s,x,y in cand if x==ja and y==jb][0]:.3f}\t{oalt[ja]}\t{nalt[jb]}\t{int(((so==1)&(sn==1)).sum())}\t"
                        f"{int((so==3).sum())}\t{int((sn==3).sum())}\n")
        dis = np.where(((so == 1) & np.isin(sn, [0, 2])) | ((sn == 1) & np.isin(so, [0, 2])))[0]
        for t in dis:
            cells.write(f"{k}\t{ja}\t{jb}\t{shared[t]}\t{SN[so[t]]}\t{SN[sn[t]]}\t{n['evnames'][int(nev[jb][t])]}\t{sa}\t{sb}\t"
                        f"{b['pos']}\t{b['end']}\t{b['svlen']}\n")
cells.close(); pairs_out.close()
with open(f"{D}/events.tsv", "w") as f:
    f.write("class\tside\tsub\tstatus\thas_alt_shared\tn\n")
    for key, v in sorted(collections.Counter(ev_rows).items()): f.write("\t".join(map(str, key)) + f"\t{v}\n")
with open(f"{D}/transitions.tsv", "w") as f:
    f.write("class\told\tnew\tn\n")
    for key, v in sorted(trans.items()): f.write("\t".join(key) + f"\t{v}\n")
with open(f"{D}/new_change_ev.tsv", "w") as f:
    f.write("class\tsub_new\told\tnew\tnew_ev\tn\n")
    for key, v in sorted(chg_ev.items(), key=lambda x: -x[1]): f.write("\t".join(key) + f"\t{v}\n")
print(len(shared), "shared")
