#!/usr/bin/env python3
"""Choose the panel by demand from the whole passing pool (the population
the 10K is drawn from). Analysis only. User, 2026-10-13: "yes, start it"
(select on the fold-back-passing pool).

Pool: passes select_isolates.py (pass_standard.tsv), lineage 1-4,
fold-back (inv_chimera_per_1k, old scale) <= 2.5, as in
analysis/cohort_1k/draw_1k.py.

Isolate states come from the collection's own Mutect2 calls on H37Rv
(<s>.m2variants.pkl.gz), no new alignment:
  PASS single-base SNV with AF >= 0.9   -> that base
  any record with AF < 0.25             -> ignored (REF)
  anything else (indel, MNV, non-PASS, AF 0.25-0.9) -> unknown over its REF span
  no record                             -> REF (no per-base depth exists)
Genome states, mask, distance and greedy are as in select_panel.py. The
distance uses all positions (full-length genome states), computed sparsely:
  D = c_g - #{p in A|N: G!=0, G!=ref} + #{p in A: G!=0, G!=alt}
Isolates are weighted so each lineage carries equal weight (the 10K draws
equal numbers per lineage).

Validation (--validate): the same states for the new1k isolates, compared
with select_panel.py's P1-based distances (out/distances.npz).

    select_pool.py --candidates candidates.tsv --outdir out_pool [--threads 8]
"""
import argparse, collections, csv, gzip, os, pickle, random, statistics as S
from multiprocessing import Pool
import numpy as np

REPO = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
BUILD = f"{REPO}/runroot/refbias/build/7713a8d71d8e-fix1"
STATES = f"{REPO}/analysis/panel_checks/states"
FB = f"{REPO}/analysis/foldback_qc"
M2 = ("/n/netscratch/sfortune_lab/Lab/pculviner/notebooks/260728_sv_exploration/"
      "pipeline_outputs/{0}/results/{0}.m2variants.pkl.gz")
L = 4411532
LINEAGES = ("lineage1", "lineage2", "lineage3", "lineage4")
SIZES = [50, 100, 150, 175, 200, 250, 300, 350, 400]
CODE = {b: i for i, b in enumerate("ACGT", 1)}


def m2_calls(s):
    """(alt_pos[int32], alt_code[uint8], unknown_pos[int32]) 0-based, or None."""
    try:
        o = pickle.load(gzip.open(M2.format(s)))
    except (OSError, EOFError, pickle.UnpicklingError):
        return None
    alt_p, alt_c, unk = [], [], []
    for pos, ref, alt, ok, af in zip(o["POS"], o["REF"], o["ALT"], o["PASSF"], o["AF"]):
        if af < 0.25:
            continue
        p = int(pos) - 1
        if ok and af >= 0.9 and len(ref) == 1 and len(alt) == 1 and alt in CODE:
            alt_p.append(p); alt_c.append(CODE[alt])
        else:
            unk.extend(range(p, p + max(len(ref), 1)))
    a = np.array(alt_p, np.int32); c = np.array(alt_c, np.uint8)
    u = np.unique(np.array(unk, np.int32))
    keep = ~np.isin(a, u)
    return a[keep], c[keep], u


def refcode():
    seq = "".join(l.strip() for l in open(f"{BUILD}/refs/GCF_000195955.fasta") if not l.startswith(">"))
    lut = np.zeros(256, np.uint8)
    for i, b in enumerate("ACGT", 1): lut[ord(b)] = i
    return lut[np.frombuffer(seq.encode(), np.uint8)]


def pool_isolates():
    fb = {}
    for f in [f"{FB}/screen50k/rates_all.tsv", f"{FB}/cohorts_in_use_rates.tsv"]:
        for r in csv.DictReader(open(f), delimiter="\t"):
            fb[r["sample"]] = float(r["inv_chimera_per_1k"])
    out = []
    for r in csv.DictReader(open(f"{FB}/screen50k/pass_standard.tsv"), delimiter="\t"):
        if r["lineage"] in LINEAGES and r["sample"] in fb and fb[r["sample"]] <= 2.5:
            out.append((r["sample"], r["lineage"], r["deepest"]))
    return out


def distances(calls, G, rc, c_g):
    D = np.zeros((len(calls), G.shape[0]), np.int32)
    for i, (a, c, u) in enumerate(calls):
        au = np.concatenate([a, u])
        Gau = G[:, au]
        sub = ((Gau != 0) & (Gau != rc[au])).sum(axis=1)
        Ga = Gau[:, :len(a)]
        add = ((Ga != 0) & (Ga != c)).sum(axis=1)
        D[i] = c_g - sub + add
    return D


def greedy(D, w, forced, n_max):
    D = D.astype(np.int64)
    order = list(forced)
    cur = D[:, order].min(axis=1)
    taken = np.zeros(D.shape[1], bool); taken[order] = True
    while len(order) < min(n_max, D.shape[1]):
        tot = (np.minimum(D, cur[:, None]) * w[:, None]).sum(axis=0)
        tot[taken] = np.inf
        j = int(np.argmin(tot)); order.append(j); taken[j] = True
        cur = np.minimum(cur, D[:, j])
    return order


def stats(d):
    return (f"{S.median(d):g}\t{int(np.percentile(d, 90))}\t{(d > 150).mean():.1%}\t{(d > 300).mean():.1%}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20261013)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rc = refcode()
    mask = np.zeros(L, bool)
    for l in open(f"{BUILD}/assets/repeat_mask.bed"):
        c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True

    cand = {r["id"]: r for r in csv.DictReader(open(a.candidates), delimiter="\t")}
    pool = [g for g, r in cand.items() if r["in_pool"] != "no"]
    cx_now = [g for g, r in cand.items() if r["source"] == "cx333"
              and r["lineage"].startswith(LINEAGES + ("lineage7",))]
    genomes = sorted(set(pool) | set(cx_now))
    gi = {g: j for j, g in enumerate(genomes)}
    G = np.zeros((len(genomes), L), np.uint8)
    for j, g in enumerate(genomes):
        G[j] = np.load(f"{STATES}/{g}.npz")["state"]
    G[:, mask] = 0
    c_g = ((G != 0) & (G != rc)).sum(axis=1).astype(np.int32)
    out = []

    # validation on new1k against the P1-based distances
    z = np.load("out/distances.npz")
    v_iso = list(z["isolates"]); Dp1 = z["D"]; vg = [gi[g] for g in z["genomes"]]
    with Pool(a.threads) as P:
        vc = P.map(m2_calls, v_iso, chunksize=20)
    ok = [i for i, c in enumerate(vc) if c is not None]
    Dm = distances([vc[i] for i in ok], G, rc, c_g)[:, vg]
    Dp = Dp1[ok]
    cxv = np.array([k for k, g in enumerate(z["genomes"]) if g in set(cx_now)])
    poolv = np.array([k for k, g in enumerate(z["genomes"]) if g in set(pool)])
    same_cx = (Dm[:, cxv].argmin(1) == Dp[:, cxv].argmin(1)).mean()
    same_pool = (Dm[:, poolv].argmin(1) == Dp[:, poolv].argmin(1)).mean()
    # how much worse is the Mutect2-chosen genome, by the P1 distance
    rows = np.arange(len(ok))
    excess = Dp[:, poolv][rows, Dm[:, poolv].argmin(1)] - Dp[:, poolv].min(1)
    r = np.corrcoef(Dm.ravel(), Dp.ravel())[0, 1]
    out.append(f"VALIDATION on new1k ({len(ok)} with Mutect2 calls): same nearest genome as P1-based "
               f"{same_cx:.1%} (CX333) / {same_pool:.1%} (pool); Mutect2's pick is within 10 SNPs of the "
               f"P1-best for {(excess <= 10).mean():.1%} (median excess {np.median(excess):g}); "
               f"r = {r:.3f}, medians {np.median(Dm):g} vs {np.median(Dp):g}")
    print(out[-1], flush=True)

    # the pool
    iso = pool_isolates()
    print(f"pool isolates: {len(iso)}", flush=True)
    with Pool(a.threads) as P:
        calls = P.map(m2_calls, [s for s, _, _ in iso], chunksize=50)
    keep = [i for i, c in enumerate(calls) if c is not None]
    out.append(f"pool isolates: {len(iso)}; with Mutect2 calls: {len(keep)}")
    iso = [iso[i] for i in keep]; calls = [calls[i] for i in keep]
    D = distances(calls, G, rc, c_g)
    del calls
    lin = np.array([l for _, l, _ in iso])
    n_l = collections.Counter(lin)
    w = np.array([1.0 / n_l[l] for l in lin])
    out.append("lineages: " + ", ".join(f"{l} {n_l[l]}" for l in LINEAGES) + " (equal weight per lineage)")
    np.savez_compressed(f"{a.outdir}/distances.npz", D=D, genomes=np.array(genomes),
                        isolates=np.array([s for s, _, _ in iso]))

    pj = np.array([gi[g] for g in pool]); cxj = np.array([gi[g] for g in cx_now])
    Dq = D[:, pj]
    forced = [k for k, g in enumerate(pool) if cand[g]["in_pool"].startswith("forced")
              or cand[g]["lineage"].startswith("lineage7")]
    rnd = random.Random(a.seed)
    half = np.zeros(len(iso), int)
    byl = collections.defaultdict(list)
    for i, l in enumerate(lin): byl[l].append(i)
    for l, ii in sorted(byl.items()):
        rnd.shuffle(ii)
        for k, i in enumerate(ii): half[i] = k % 2
    out += ["", "CROSS-VALIDATED (select on one half, score on the other; unweighted stats)",
            "select\tscore\tsize\tmedian\tp90\tover150\tover300"]
    for h in (0, 1):
        o = greedy(Dq[half == h], w[half == h], forced, max(SIZES))
        De = Dq[half != h]
        for k in SIZES:
            out.append(f"half{h}\thalf{1 - h}\t{k}\t{stats(De[:, o[:k]].min(1))}")
        out.append(f"half{h}\thalf{1 - h}\tCX333 now ({len(cxj)})\t{stats(D[half != h][:, cxj].min(1))}")
        print("\n".join(out[-len(SIZES) - 1:]), flush=True)

    o = greedy(Dq, w, forced, len(pool))
    cur = np.full(len(iso), 10**9); last = 0
    with open(f"{a.outdir}/order.tsv", "w") as fo:
        fo.write("rank\tid\tsource\tlineage\tforced\tnearest_for\tmedian_after\tover300_after\n")
        nearest = collections.Counter()
        for k, j in enumerate(o, 1):
            new = np.minimum(cur, Dq[:, j])
            if (new < cur).any(): last = k
            cur = new
        # isolates whose nearest (first-chosen among ties) is each genome, at full size
        cur = np.full(len(iso), 10**9)
        for k, j in enumerate(o, 1):
            new = np.minimum(cur, Dq[:, j]); nf = int((new < cur).sum()); cur = new
            fo.write(f"{k}\t{pool[j]}\t{cand[pool[j]]['source']}\t{cand[pool[j]]['lineage']}\t"
                     f"{'yes' if j in forced else ''}\t{nf}\t{S.median(cur):g}\t{int((cur > 300).sum())}\n")
    out += ["", f"ALL POOL ISOLATES: last rank that lowers any isolate's distance: {last}",
            "size\tmedian\tp90\tover150\tover300\tcx333\tmarin\tgenbank"]
    for k in SIZES:
        src = collections.Counter(cand[pool[j]]["source"] for j in o[:k])
        out.append(f"{k}\t{stats(Dq[:, o[:k]].min(1))}\t{src['cx333']}\t{src['marin_candidate']}\t{src['genbank']}")
    out.append(f"CX333 now ({len(cxj)})\t{stats(D[:, cxj].min(1))}")

    def grp(lab):
        return ".".join(lab.split(".")[:3])
    by = collections.defaultdict(list)
    for i, (_, _, d) in enumerate(iso): by[grp(d)].append(i)
    now = D[:, cxj].min(1)
    at = {k: Dq[:, o[:k]].min(1) for k in (175, 250)}
    out += ["", "BY SUBLINEAGE: isolates, median distance, share over 300",
            "sublineage\tisolates\tnow_median\tnow_over300\tk175_median\tk175_over300\tk250_median\tk250_over300"]
    for gname in sorted(by):
        ii = by[gname]
        if len(ii) < 10: continue
        out.append(f"{gname}\t{len(ii)}\t{S.median(now[ii]):g}\t{(now[ii] > 300).mean():.0%}\t" +
                   "\t".join(f"{S.median(at[k][ii]):g}\t{(at[k][ii] > 300).mean():.0%}" for k in (175, 250)))
    open(f"{a.outdir}/summary.txt", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
