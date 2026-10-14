#!/usr/bin/env python3
"""Graph genotype check, step 2: each panel genome's graph VCF genotype
against its own direct alignment to H37Rv (user, 2026-10-14). Analysis only.

Graph side: the P200 build's collapsed graph VCF, keys normalised and pooled
exactly as bin/panel_ancestral.py does (read_panel, pooled_column).
Direct side, per genome:
  SNP keys        the genome's direct-alignment state
                  (analysis/panel_checks/states/<acc>.npz): its base = ALT ->
                  1, = REF -> 0, N or a third base -> no call
  indel/MNP/SV    route A (kit/<acc>.var, .sv.tsv) through
                  panel_ancestral.kit_event_states: the same conservative
                  matching used for the L8/canettii outgroups
A cell is scored when both sides call it. Discordant cells are graph-ALT /
direct-REF ("graph_only") or graph-REF / direct-ALT ("direct_only").

Also:
- representation splits: SV keys of one type within 50 bp of each other,
  carried by disjoint genomes (one event written as several records);
- homoplasy rerun: the SV homoplasy count (analysis/sv_homoplasy) with
  discordant cells set to no call, and restricted to SV records with no
  discordant cell.

    score.py --outdir out
"""
import argparse, collections, csv, os, sys
import numpy as np

R = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
sys.path.insert(0, f"{R}/bin")
import panel_ancestral as pa  # noqa: E402

B = f"{R}/runroot/refbias/build/3893a3762ed1"
KIT = f"{R}/analysis/graph_gt_check/kit"
STATES = f"{R}/analysis/panel_checks/states"
TREE = f"{R}/runroot/panel/P200/trees/p200.rooted.nwk"
H37 = "GCF_000195955"
CODE = {b: i for i, b in enumerate("ACGT", 1)}


def fitch_changes(ch, post, col):
    sets, cost = {}, 0
    for n in post:
        if n not in ch:
            v = col.get(n)
            sets[n] = {0, 1} if v is None else {v}
            continue
        ss = [sets[c] for c in ch[n]]
        inter = set.intersection(*ss)
        if inter:
            sets[n] = inter
        else:
            cnt = collections.Counter(s for x in ss for s in x)
            m = max(cnt.values())
            sets[n] = {s for s, v in cnt.items() if v == m}
            cost += len(ss) - m
    return cost


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    seq = "".join(l.strip().upper() for l in open(f"{B}/refs/{H37}.fasta") if not l.startswith(">"))
    contig = "NC_000962.3"
    mask = np.zeros(len(seq) + 200, bool)
    for l in open(f"{B}/assets/repeat_mask.bed"):
        c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True
    samples, keys, copies, st = pa.read_panel(f"{B}/assets/graph_collapsed.vcf.gz",
                                              os.environ["MTB_BCFTOOLS"], seq, contig)
    kinds = [pa.kind_of(k[2], k[3]) for k in keys]
    cols = [pa.pooled_column(c) for c in copies]
    print(f"keys: {len(keys):,} ({dict(collections.Counter(kinds))}); genomes {len(samples)}", flush=True)
    nk = len(keys)
    G = np.array([[1 if c[i] == "1" else 0 if c[i] == "0" else -1 for i in range(len(samples))] for c in cols], np.int8)
    D = np.full_like(G, -1)
    snp_idx = [i for i, k in enumerate(kinds) if k == "SNP"]
    for j, g in enumerate(samples):
        if not os.path.exists(f"{KIT}/{g}.var"):
            print("no kit output:", g); continue
        stt = np.load(f"{STATES}/{g}.npz")["state"]
        for i in snp_idx:
            _, pos, ref, alt = keys[i]
            b = stt[pos - 1]
            D[i, j] = 1 if b == CODE.get(alt) else 0 if b == CODE.get(ref) else -1
        ev = pa.read_kit_events(f"{KIT}/{g}", seq, contig)
        sts = pa.kit_event_states(keys, kinds, seq, contig, ev)
        for i, v in sts.items() if isinstance(sts, dict) else enumerate(sts):
            if kinds[i] != "SNP" and v is not None:
                D[i, j] = v
    both = (G >= 0) & (D >= 0)
    disc = both & (G != D)
    gonly = both & (G == 1) & (D == 0)
    donly = both & (G == 0) & (D == 1)
    masked = np.array([bool(mask[k[1] - 1:k[1] - 1 + max(len(k[2]), 1)].any()) for k in keys])
    K = np.array(kinds)
    out = ["CELL CONCORDANCE (cells called by both)", "class\tregion\tcells\tdiscordant\trate\tgraph_only\tdirect_only\trecords_with_any\trecords"]
    for cl in ("SNP", "MNP", "indel", "SV"):
        for reg, m in (("outside mask", ~masked), ("inside mask", masked)):
            sel = (K == cl) & m
            c = int(both[sel].sum()); d = int(disc[sel].sum())
            out.append(f"{cl}\t{reg}\t{c:,}\t{d:,}\t{d / max(c, 1):.3%}\t{int(gonly[sel].sum()):,}\t{int(donly[sel].sum()):,}"
                       f"\t{int((disc[sel].sum(1) > 0).sum()):,}\t{int(sel.sum()):,}")
    pg = disc.sum(0) / np.maximum(both.sum(0), 1)
    order = np.argsort(-pg)
    out += ["", "GENOMES WITH THE MOST DISCORDANCE (all classes)", "genome\tdiscordant\tcalled\trate\tSV_disc\tindel_disc"]
    for j in order[:15]:
        out.append(f"{samples[j]}\t{int(disc[:, j].sum())}\t{int(both[:, j].sum())}\t{pg[j]:.3%}\t"
                   f"{int(disc[K == 'SV', j].sum())}\t{int(disc[K == 'indel', j].sum())}")
    out.append(f"median genome discordance {np.median(pg):.3%}")

    # representation splits among SV keys
    sv = [i for i in range(nk) if kinds[i] == "SV"]
    def typ(i): return "INS" if len(keys[i][3]) > len(keys[i][2]) else "DEL"
    sv.sort(key=lambda i: keys[i][1])
    groups, cur = [], []
    for i in sv:
        if cur and (keys[i][1] - keys[cur[-1]][1] > 50):
            groups.append(cur); cur = []
        cur.append(i)
    if cur: groups.append(cur)
    n_multi = n_disjoint = 0
    for gidx in groups:
        for t in ("INS", "DEL"):
            ii = [i for i in gidx if typ(i) == t]
            if len(ii) < 2: continue
            n_multi += 1
            carr = [set(np.where(G[i] == 1)[0]) for i in ii]
            if all(not (carr[x] & carr[y]) for x in range(len(ii)) for y in range(x + 1, len(ii))):
                n_disjoint += 1
    out += ["", f"REPRESENTATION: {n_multi:,} loci with >= 2 SV records of one type within 50 bp; "
                f"{n_disjoint:,} of them with disjoint carriers (one event possibly split)"]

    # homoplasy rerun
    children, root, leaves = pa.parse_newick(open(TREE).read())
    keep = set(samples) | {H37}
    ch = {}
    def prune(n):
        if n not in children: return n if n in keep else None
        k = [x for x in (prune(c) for c in children[n]) if x is not None]
        if not k: return None
        if len(k) == 1: return k[0]
        ch[n] = k; return n
    rt = prune(root)
    post, stk = [], [(rt, 0)]
    while stk:
        n, d = stk.pop()
        if d: post.append(n); continue
        stk.append((n, 1))
        for c in ch.get(n, []): stk.append((c, 0))
    rows = []
    for i in sv:
        col = {s: (int(G[i, j]) if G[i, j] >= 0 and not disc[i, j] else None) for j, s in enumerate(samples)}
        col[H37] = 0
        raw = {s: (int(G[i, j]) if G[i, j] >= 0 else None) for j, s in enumerate(samples)}; raw[H37] = 0
        rows.append((keys[i][1], typ(i), abs(len(keys[i][3]) - len(keys[i][2])), int((G[i] == 1).sum()),
                     fitch_changes(ch, post, raw), fitch_changes(ch, post, col), int(disc[i].sum()), int(both[i].sum()),
                     bool(masked[i])))
    with open(f"{a.outdir}/sv_records.tsv", "w") as fo:
        fo.write("pos\ttype\tsize\tcarriers\tchanges_graph\tchanges_discordant_removed\tdiscordant_cells\tcells_checked\tmasked\n")
        for r in rows: fo.write("\t".join(map(str, r)) + "\n")
    out += ["", "SV HOMOPLASY (>= 2 changes on the panel tree)",
            "set\tSVs\tgraph genotypes\tdiscordant cells removed\tonly records with no discordant cell"]
    for lab, sub in (("all", rows), ("outside mask", [r for r in rows if not r[8]])):
        clean = [r for r in sub if r[6] == 0]
        out.append(f"{lab}\t{len(sub):,}\t{sum(r[4] >= 2 for r in sub):,}\t{sum(r[5] >= 2 for r in sub):,}"
                   f"\t{sum(r[4] >= 2 for r in clean):,} of {len(clean):,}")
    with open(f"{a.outdir}/records.tsv", "w") as fo:
        fo.write("chrom\tpos\tref\talt\tclass\tmasked\tcarriers\tcalled_both\tdiscordant\tgraph_only\tdirect_only\n")
        for i, k in enumerate(keys):
            fo.write(f"{k[0]}\t{k[1]}\t{k[2] if len(k[2]) <= 30 else str(len(k[2])) + 'bp'}\t"
                     f"{k[3] if len(k[3]) <= 30 else str(len(k[3])) + 'bp'}\t{kinds[i]}\t{int(masked[i])}\t"
                     f"{int((G[i] == 1).sum())}\t{int(both[i].sum())}\t{int(disc[i].sum())}\t{int(gonly[i].sum())}\t{int(donly[i].sum())}\n")
    open(f"{a.outdir}/summary.txt", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
