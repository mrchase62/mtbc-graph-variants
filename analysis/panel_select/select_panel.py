#!/usr/bin/env python3
"""Choose the new panel by demand from the new1k isolates. Analysis only.

Method (as explained to the user 2026-10-13; the distance is check 3's,
analysis/panel_checks/bin/check3_coverage.py):

1. Pool: candidates.tsv rows with in_pool != "no" (qc_candidates.py).
   Forced in first: H37Rv, the two in-graph outgroups, and every pool
   lineage 7 genome (L7 stays in the graph; no L7 isolates create demand).
2. Genome states: direct alignment to H37Rv (panel_checks/states/*.npz).
3. Isolate states: P1's own H37Rv calls (work/new1k_p1/<s>.h37rv.vcf.gz).
   A single-base ALT is the isolate's base; any other record (indel, multi-
   base) makes its positions unknown; positions with depth < 5
   (lowdepth/<s>.bed.gz) are unknown; everything else is REF.
4. Distance = positions outside the build's repeat mask where both are
   called and differ.
5. Greedy k-median: add the genome that most lowers the summed distance of
   isolates to their nearest panel genome. Cross-validated: isolates are
   split in two halves (fixed seed, stratified by lineage); select on one,
   score on the other, both ways. The final order uses all isolates.

    select_panel.py --candidates candidates.tsv --outdir out
"""
import argparse, collections, csv, gzip, os, random, statistics as S
import numpy as np
import pysam

REPO = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
BUILD = f"{REPO}/runroot/refbias/build/7713a8d71d8e-fix1"
STATES = f"{REPO}/analysis/panel_checks/states"
WORK = f"{REPO}/runroot/refbias/work/new1k_p1"
L = 4411532
SIZES = [50, 100, 150, 200, 250, 300, 350, 400]


def refcode():
    seq = "".join(l.strip() for l in open(f"{BUILD}/refs/GCF_000195955.fasta") if not l.startswith(">"))
    lut = np.zeros(256, np.uint8)
    for i, b in enumerate("ACGT", 1): lut[ord(b)] = i
    return lut[np.frombuffer(seq.encode(), np.uint8)]


def isolate_state(s, rc):
    """Full-length state vector: 0 unknown, 1-4 ACGT."""
    st = rc.copy()
    code = {b: i for i, b in enumerate("ACGT", 1)}
    for r in pysam.VariantFile(f"{WORK}/{s}.h37rv.vcf.gz"):
        gt = r.samples[0]["GT"]
        if gt is None or gt[0] is None:
            st[r.pos - 1:r.pos - 1 + len(r.ref)] = 0; continue
        if gt[0] == 0:
            continue
        al = r.alleles[gt[0]]
        if len(r.ref) == 1 and len(al) == 1 and al in code:
            st[r.pos - 1] = code[al] if st[r.pos - 1] in (rc[r.pos - 1], code[al]) else 0
        else:
            st[r.pos - 1:r.pos - 1 + max(len(r.ref), 1)] = 0
    with gzip.open(f"lowdepth/{s}.bed.gz", "rt") as fh:
        for l in fh:
            _, a, b = l.split("\t"); st[int(a):int(b)] = 0
    return st


def greedy(D, forced, n_max):
    """D: isolates x genomes (int). Returns the selection order."""
    D = D.astype(np.int64)
    order = list(forced)
    cur = D[:, order].min(axis=1) if order else np.full(len(D), 10**9, np.int64)
    taken = np.zeros(D.shape[1], bool); taken[order] = True
    while len(order) < min(n_max, D.shape[1]):
        tot = np.minimum(D, cur[:, None]).sum(axis=0)
        tot[taken] = np.iinfo(np.int64).max
        j = int(np.argmin(tot)); order.append(j); taken[j] = True
        cur = np.minimum(cur, D[:, j])
    return order


def stats(d):
    return (f"{S.median(d):g}\t{int(np.percentile(d, 90))}\t{int((d > 150).sum())}\t"
            f"{int((d > 300).sum())}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--cohort", default=f"{REPO}/analysis/cohort_1k/cohort_1k.tsv")
    ap.add_argument("--refmap", default=f"{REPO}/runroot/refbias/new1k/p1/refmap.tsv")
    ap.add_argument("--seed", type=int, default=20261013)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rc = refcode()
    mask = np.zeros(L, bool)
    for l in open(f"{BUILD}/assets/repeat_mask.bed"):
        c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True

    cand = {r["id"]: r for r in csv.DictReader(open(a.candidates), delimiter="\t")}
    pool = [g for g, r in cand.items() if r["in_pool"] != "no"]
    cx_now = [g for g, r in cand.items() if r["source"] == "cx333"
              and r["lineage"].startswith(("lineage1", "lineage2", "lineage3", "lineage4", "lineage7"))]
    genomes = sorted(set(pool) | set(cx_now))
    gi = {g: j for j, g in enumerate(genomes)}
    iso = list(csv.DictReader(open(a.cohort), delimiter="\t"))
    print(f"pool {len(pool)}; current CX333 L1-4,7 {len(cx_now)}; isolates {len(iso)}", flush=True)

    # informative positions (genomes, then isolates), then isolate states there
    U = np.zeros(L, bool)
    GS = {}
    for g in genomes:
        st = np.load(f"{STATES}/{g}.npz")["state"]; GS[g] = st
        U |= (st > 0) & (st != rc)
    for r in iso:
        st = isolate_state(r["sample"], rc)
        U |= (st > 0) & (st != rc)
    U &= ~mask
    U = np.where(U)[0]
    print(f"informative positions: {len(U):,}", flush=True)
    I = np.stack([isolate_state(r["sample"], rc)[U] for r in iso])
    G = np.stack([GS[g][U] for g in genomes]); del GS
    D = np.zeros((len(iso), len(genomes)), np.int32)
    Gc = G != 0
    for i in range(len(iso)):
        D[i] = ((G != I[i]) & Gc & (I[i] != 0)).sum(axis=1)
    np.savez_compressed(f"{a.outdir}/distances.npz", D=D, genomes=np.array(genomes),
                        isolates=np.array([r["sample"] for r in iso]))

    out = []
    # validation against P1's own choice
    rm = {r["sample"]: r for r in csv.DictReader(open(a.refmap), delimiter="\t")}
    cxj = np.array([gi[g] for g in cx_now])
    agree = n = 0; pairs = []
    for i, r in enumerate(iso):
        p = rm.get(r["sample"])
        if not p: continue
        n += 1
        agree += genomes[cxj[int(np.argmin(D[i, cxj]))]] == p["reference"]
        if p["reference"] in gi:
            pairs.append((int(D[i, gi[p["reference"]]]), int(p["snp_distance"])))
    rr = np.corrcoef(np.array(pairs).T)[0, 1]
    out.append(f"VALIDATION: nearest current-CX333 genome by this distance = P1's choice for {agree}/{n}; "
               f"distance to P1's reference vs P1 snp_distance: r = {rr:.3f}, medians "
               f"{S.median(p[0] for p in pairs):g} vs {S.median(p[1] for p in pairs):g}")

    # selection
    pj = np.array([gi[g] for g in pool])
    Dp = D[:, pj]
    forced = [k for k, g in enumerate(pool) if cand[g]["in_pool"].startswith("forced")
              or cand[g]["lineage"].startswith("lineage7")]
    out.append(f"pool {len(pool)} genomes ({collections.Counter(cand[g]['source'] for g in pool)}); "
               f"forced first: {', '.join(pool[k] for k in forced)}")
    rnd = random.Random(a.seed)
    half = np.zeros(len(iso), int)
    byl = collections.defaultdict(list)
    for i, r in enumerate(iso): byl[r["lineage"]].append(i)
    for l, ii in sorted(byl.items()):
        rnd.shuffle(ii)
        for k, i in enumerate(ii): half[i] = k % 2
    out += ["", "CROSS-VALIDATED (select on one half, score on the other)",
            "select\tscore\tsize\tmedian\tp90\tover150\tover300"]
    nmax = max(SIZES)
    for h in (0, 1):
        o = greedy(Dp[half == h], forced, nmax)
        De = Dp[half != h]
        for k in SIZES:
            out.append(f"half{h}\thalf{1 - h}\t{k}\t{stats(De[:, o[:k]].min(axis=1))}")
        out.append(f"half{h}\thalf{1 - h}\tCX333 now ({len(cxj)})\t{stats(D[half != h][:, cxj].min(axis=1))}")

    o = greedy(Dp, forced, len(pool))
    out += ["", "ALL ISOLATES (selected and scored on all 1,000; optimistic)",
            "size\tmedian\tp90\tover150\tover300\tcx333\tmarin\tgenbank"]
    for k in SIZES:
        src = collections.Counter(cand[pool[j]]["source"] for j in o[:k])
        out.append(f"{k}\t{stats(Dp[:, o[:k]].min(axis=1))}\t{src['cx333']}\t{src['marin_candidate']}\t{src['genbank']}")
    out.append(f"CX333 now ({len(cxj)})\t{stats(D[:, cxj].min(axis=1))}")
    with open(f"{a.outdir}/order.tsv", "w") as fo:
        fo.write("rank\tid\tsource\tlineage\tforced\tmedian_after\tover300_after\n")
        cur = np.full(len(iso), 10**9)
        for k, j in enumerate(o, 1):
            cur = np.minimum(cur, Dp[:, j])
            fo.write(f"{k}\t{pool[j]}\t{cand[pool[j]]['source']}\t{cand[pool[j]]['lineage']}\t"
                     f"{'yes' if j in forced else ''}\t{S.median(cur):g}\t{int((cur > 300).sum())}\n")

    # gaps by sublineage at each size, and now
    def grp(lab):
        return ".".join(lab.split(".")[:3])
    by = collections.defaultdict(list)
    for i, r in enumerate(iso): by[grp(r["deepest"])].append(i)
    now = D[:, cxj].min(axis=1)
    cols = [150, 300]
    out += ["", "BY SUBLINEAGE (all isolates): median distance and isolates over 300",
            "sublineage\tisolates\tnow_median\tnow_over300\t" + "\t".join(f"k{k}_median\tk{k}_over300" for k in cols)]
    at = {k: Dp[:, o[:k]].min(axis=1) for k in cols}
    for gname in sorted(by):
        ii = by[gname]
        out.append(f"{gname}\t{len(ii)}\t{S.median(now[ii]):g}\t{int((now[ii] > 300).sum())}\t" +
                   "\t".join(f"{S.median(at[k][ii]):g}\t{int((at[k][ii] > 300).sum())}" for k in cols))
    with open(f"{a.outdir}/per_isolate.tsv", "w") as fo:
        fo.write("sample\tlineage\tdeepest\tnearest_now\tdist_now\t" + "\t".join(f"nearest_k{k}\tdist_k{k}" for k in cols) + "\n")
        for i, r in enumerate(iso):
            fo.write(f"{r['sample']}\t{r['lineage']}\t{r['deepest']}\t{genomes[cxj[int(np.argmin(D[i, cxj]))]]}\t{now[i]}\t" +
                     "\t".join(f"{pool[o[int(np.argmin(Dp[i, o[:k]]))]]}\t{at[k][i]}" for k in cols) + "\n")
    open(f"{a.outdir}/summary.txt", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
