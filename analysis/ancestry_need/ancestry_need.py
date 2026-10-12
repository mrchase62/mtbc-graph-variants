#!/usr/bin/env python3
"""Which loci need genomes outside lineages 1-4 and 7 to know their ancestral
state? (user, 2026-10-12: "I am still on the fence about removing lineages
from the graph that may provide key ancestry for specific loci.") Report
only. Run from the repository root with $MTB_PY.

Inputs (build 7713a8d71d8e-fix1, CX333): the panel tree
(assets/ancestral_inputs/tree.nwk, the tree production's AA uses) and the
graph VCF (assets/graph_collapsed.vcf.gz: every panel genome's genotype at
every SNP, indel and SV record against H37Rv).

1. Ingroup I: the MRCA of every lineage 1-4 and 7 genome on the tree.
2. Records that vary inside I (two or more alleles among I's genomes with a
   call). Only these can give events inside lineages 1-4.
3. Fitch parsimony, upward, over I's own genomes: if I's set is one allele,
   the ingroup decides its ancestral state alone. If it is tied, outside
   genomes are needed.
4. The state at I from the whole tree (Fitch's second pass: the parent's
   final set where it meets the node's set, else the node's set), in three
   panels:
     full       all 333 CX333 genomes (today)
     proposed   I + GCF_022870225 (L5.1) + GCF_022870205 (L6.3.1), the
                lineage 5/6 pair of test graph B, + lineage 8 GCF_012923765
                and canettii ET1291 GCF_035581225 (the kit; their graph
                genotypes stand in for kit rows: the outgroup test found
                them 98-99.96% identical)
     no_L56     I + lineage 8 + ET1291 only
5. Each record is placed in an H37Rv gene (or between genes) from the H37Rv
   annotation; drug-resistance genes are flagged.
Writes out/records.tsv (every record varying inside I), out/genes.tsv
(genes with records whose ancestry the ingroup cannot decide) and prints a
summary.
"""
import bisect
import collections
import csv
import gzip
import os
import sys

import numpy as np

sys.path.insert(0, "bin")
from ancestral_alleles import parse_newick  # noqa: E402

B = "runroot/refbias/build/7713a8d71d8e-fix1"
O = "analysis/ancestry_need/out"
QC = "/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/qc/qc_master_table.tsv"
L56 = ["GCF_022870225", "GCF_022870205"]   # default: test graph B's pair
KIT = ["GCF_012923765", "GCF_035581225"]
DR = {"rpoB", "rpoC", "rpoA", "katG", "fabG1", "inhA", "ahpC", "embA", "embB", "embC", "gyrA", "gyrB",
      "rpsL", "rrs", "rrl", "gid", "pncA", "ethA", "ethR", "eis", "tlyA", "Rv0678", "atpE", "pepQ",
      "rplC", "ddn", "fgd1", "fbiA", "fbiB", "fbiC", "thyA", "folC", "ribD", "alr", "ald", "panD"}


def popcount(x):
    x = x.astype(np.uint32)
    c = np.zeros_like(x)
    while x.any():
        c += x & 1
        x >>= 1
    return c


def prune(children, root, keep):
    """Children dict of the tree restricted to `keep`, unary nodes collapsed."""
    out = {}

    def rec(n):
        if n not in children:
            return n if n in keep else None
        kids = [k for k in (rec(c) for c in children[n]) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        out[n] = kids
        return n
    r = rec(root)
    return out, r


def fitch(children, root, leaf_bits, nsite):
    """Upward sets for every node, then final sets top-down."""
    up = {}

    def post(n):
        if n not in children:
            up[n] = leaf_bits[n]
            return
        for c in children[n]:
            post(c)
        inter = np.bitwise_and.reduce([up[c] for c in children[n]])
        uni = np.bitwise_or.reduce([up[c] for c in children[n]])
        up[n] = np.where(inter != 0, inter, uni)
    sys.setrecursionlimit(10000)
    post(root)
    fin = {root: up[root]}
    stack = [root]
    while stack:
        n = stack.pop()
        for c in children.get(n, []):
            m = fin[n] & up[c]
            fin[c] = np.where(m != 0, m, up[c])
            stack.append(c)
    return up, fin


def mrca(children, root, targets):
    best = None

    def rec(n):
        nonlocal best
        if n not in children:
            return 1 if n in targets else 0
        k = sum(rec(c) for c in children[n])
        if k == len(targets) and best is None:
            best = n
        return k
    rec(root)
    return best


def leaves_under(children, n):
    if n not in children:
        return {n}
    return set().union(*(leaves_under(children, c) for c in children[n]))


def genes(gff):
    iv = []
    for line in gzip.open(gff, "rt"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
            continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        iv.append((int(f[3]), int(f[4]), a.get("Name", ""), a.get("locus_tag", "")))
    iv.sort()
    return iv


def gene_at(iv, starts, p):
    i = bisect.bisect_right(starts, p) - 1
    hits = [g for g in iv[max(0, i - 3):i + 1] if g[0] <= p <= g[1]]
    if hits:
        return hits[-1][2] or hits[-1][3], hits[-1][3]
    lo = iv[i] if i >= 0 else None
    hi = iv[i + 1] if i + 1 < len(iv) else None
    return f"{(lo[2] or lo[3]) if lo else ''}-{(hi[2] or hi[3]) if hi else ''}", "intergenic"


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--l56", default=",".join(L56), help="in-graph lineage 5/6 genomes (comma list)")
    ap.add_argument("--sweep-line", action="store_true",
                    help="print one summary line for this --l56 and write nothing")
    a = ap.parse_args()
    l56 = [x for x in a.l56.split(",") if x]
    os.makedirs(O, exist_ok=True)
    children, root, leaves = parse_newick(open(f"{B}/assets/ancestral_inputs/tree.nwk").read())
    lin = {}
    for r in csv.DictReader(open(QC), delimiter="\t"):
        if r["in_panel"] == "yes":
            lin[r["accession"] if "accession" in r else next(iter(r.values()))] = r["lineage"]
    target = {g for g in leaves if lin.get(g) in ("lineage1", "lineage2", "lineage3", "lineage4", "lineage7")}
    I = mrca(children, root, target)
    under = leaves_under(children, I)
    extra = sorted(under - target)
    outside = sorted(leaves - under)
    print(f"tree: {len(leaves)} leaves; ingroup I = MRCA of {len(target)} lineage 1-4,7 genomes; "
          f"{len(under)} leaves under I" + (f" (also {', '.join(f'{g} {lin.get(g)}' for g in extra)})" if extra else ""))
    print(f"outside I: {len(outside)} genomes: " + ", ".join(
        f"{k} {v}" for k, v in sorted(collections.Counter(lin.get(g, '?') for g in outside).items())))

    # graph VCF -> per-record allele bitmasks per genome
    fh = gzip.open(f"{B}/assets/graph_collapsed.vcf.gz", "rt")
    for line in fh:
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            break
    # H37Rv is the VCF's reference path, not a sample: REF (0) everywhere
    H37 = "GCF_000195955"
    col = {s: i for i, s in enumerate(samples)}
    col[H37] = len(samples)
    missing = [g for g in leaves if g not in col]
    if missing:
        sys.exit(f"FATAL: tree leaves without a VCF column: {missing[:5]}")
    under_idx = np.array([col[g] for g in sorted(under)])
    recs, rows = [], []
    for line in fh:
        f = line.rstrip("\n").split("\t")
        gt = f[9:] + ["0"]
        g = np.array([int(x) if x.isdigit() else -1 for x in gt], dtype=np.int16)
        ins = g[under_idx]
        if len(set(ins[ins >= 0].tolist())) < 2:
            continue
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        ref, alts = f[3], f[4].split(",")
        d = max(abs(len(a) - len(ref)) for a in alts)
        kind = "SNP" if len(ref) == 1 and all(len(a) == 1 for a in alts) else ("SV" if d >= 50 else
                                                                               ("indel" if d else "MNP"))
        recs.append(dict(pos=int(f[1]), ref=ref if len(ref) <= 20 else f"{ref[:10]}..({len(ref)})",
                         alt=f[4] if len(f[4]) <= 20 else f"..({max(len(a) for a in alts)})",
                         kind=kind, size=d, lv=info.get("LV", "")))
        rows.append(g)
    G = np.array(rows, dtype=np.int16).T          # genomes x records
    nsite = G.shape[1]
    bits = np.where(G >= 0, (1 << np.clip(G, 0, 7)).astype(np.uint8), np.uint8(0xFF))
    leaf_bits = {s: bits[col[s]] for s in leaves}
    print(f"records varying inside I: {nsite:,} ({', '.join(f'{k} {v:,}' for k, v in sorted(collections.Counter(r['kind'] for r in recs).items()))})")

    up_full, fin_full = fitch(children, root, leaf_bits, nsite)
    inside = up_full[I]
    decided = popcount(inside) == 1
    res = {"full": fin_full[I]}
    for name, keep in (("proposed", under | set(l56) | set(KIT)), ("no_L56", under | set(KIT))):
        ch, rt = prune(children, root, keep)
        _, fin = fitch(ch, rt, leaf_bits, nsite)
        res[name] = fin[I]

    def st(x):
        return int(x).bit_length() - 1 if bin(int(x)).count("1") == 1 else "tied"

    iv = genes(f"{B}/annotation/GCF_000195955.gff.gz")
    starts = [g[0] for g in iv]
    out = []
    for j, r in enumerate(recs):
        f_, p_, n_ = st(res["full"][j]), st(res["proposed"][j]), st(res["no_L56"][j])
        if decided[j]:
            cat = "ingroup_decides"
        elif f_ == "tied":
            cat = "tied_even_with_full_panel"
        elif p_ == f_:
            cat = "outgroup_needed_proposed_agrees"
        elif p_ == "tied":
            cat = "outgroup_needed_proposed_unresolved"
        else:
            cat = "outgroup_needed_proposed_differs"
        gname, tag = gene_at(iv, starts, r["pos"])
        out.append(dict(r, gene=gname, locus_tag=tag, category=cat, anc_full=f_, anc_proposed=p_,
                        anc_no_L56=n_, dr_gene="yes" if gname in DR else ""))
    if a.sweep_line:
        n = collections.Counter(r["category"] for r in out)
        el = next((r["anc_proposed"] for r in out if r["pos"] == 2352065 and r["kind"] == "SV"
                   and r["category"] != "ingroup_decides" and r["anc_full"] == 1), "")
        print(f"{a.l56}\tdiffers {n['outgroup_needed_proposed_differs']}\tunresolved "
              f"{n['outgroup_needed_proposed_unresolved']}\telement_ancestral "
              f"{'present' if el == 1 else 'absent' if el == 0 else el}")
        return
    with open(f"{O}/records.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    print("\nrecords by category and kind")
    cats = ["ingroup_decides", "outgroup_needed_proposed_agrees", "outgroup_needed_proposed_differs",
            "outgroup_needed_proposed_unresolved", "tied_even_with_full_panel"]
    kinds = ["SNP", "indel", "SV", "MNP"]
    print("category\t" + "\t".join(kinds) + "\tall")
    for c in cats:
        n = collections.Counter(r["kind"] for r in out if r["category"] == c)
        print(f"{c}\t" + "\t".join(f"{n[k]:,}" for k in kinds) + f"\t{sum(n.values()):,}")
    nl = collections.Counter(r["category"] for r in out if r["category"].startswith("outgroup_needed")
                             and r["anc_no_L56"] != r["anc_full"])
    print(f"\nwithout the two L5/6 genomes (kit only), records where the outgroup is needed and the answer "
          f"differs from the full panel or is unresolved: {sum(nl.values()):,}")

    bad = [r for r in out if r["category"] in ("outgroup_needed_proposed_differs",
                                               "outgroup_needed_proposed_unresolved")]
    by_gene = collections.defaultdict(list)
    for r in bad:
        by_gene[(r["gene"], r["locus_tag"])].append(r)
    grows = [dict(gene=g, locus_tag=t, records=len(v), kinds=",".join(sorted({x["kind"] for x in v})),
                  differs=sum(x["category"].endswith("differs") for x in v),
                  unresolved=sum(x["category"].endswith("unresolved") for x in v),
                  dr_gene="yes" if g in DR else "", first_pos=v[0]["pos"])
             for (g, t), v in by_gene.items()]
    grows.sort(key=lambda x: -x["records"])
    with open(f"{O}/genes.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(grows[0]) if grows else ["gene"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(grows)
    print(f"\nloci (genes or intergenic stretches) where the proposed panel differs from or cannot resolve "
          f"the full panel's ancestral state: {len(grows):,}")
    print("gene\tlocus_tag\trecords\tkinds\tdiffers\tunresolved\tdr_gene")
    for g in grows[:30]:
        print(f"{g['gene']}\t{g['locus_tag']}\t{g['records']}\t{g['kinds']}\t{g['differs']}\t{g['unresolved']}\t{g['dr_gene']}")
    dr = [g for g in grows if g["dr_gene"]]
    print(f"\ndrug-resistance genes among them: {', '.join(g['gene'] + ' ' + str(g['records']) for g in dr) or 'none'}")
    drall = collections.Counter(r["category"] for r in out if r["dr_gene"])
    print("all records in drug-resistance genes varying inside I, by category: " +
          ", ".join(f"{k} {v}" for k, v in sorted(drall.items())))


if __name__ == "__main__":
    main()
