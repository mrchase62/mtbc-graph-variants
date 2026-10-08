#!/usr/bin/env python3
"""Outgroup test: compare the arm G and arm K trees (canettii rows from the
graph vs from direct alignment). Analysis only.

    $MTB_PY_VT compare_trees.py --g G.contree --k K.contree \\
        --panel-lineages <W>/data/qc/lineages.all.tsv \\
        --cohort <W>/refbias/cohort.scale200.tsv

Reports, for each tree:
  * the root: the tips outside the MTBC clade (the smallest clade holding
    every tip that is not canettii) -- must be exactly the canettii tips;
  * whether each main lineage is monophyletic, whether lineages 1-4 and 7
    form one clade, and which lineage splits first inside it;
and between the two:
  * the Robinson-Foulds distance (unrooted splits);
  * every split in one tree and not the other, with its UFBoot support.
"""
import argparse, collections, csv
import dendropy

CANETTII = {"GCF_035581225", "GCF_000253375", "canettii"}


def main_lineage(s):
    s = (s or "").split(":")[0].strip()
    if not s.startswith("lineage"):
        return s or "?"
    return "lineage" + s[len("lineage"):].split(".")[0]


def load_lineages(a):
    lin = {}
    for r in csv.DictReader(open(a.panel_lineages), delimiter="\t"):
        lin[r["accession"]] = main_lineage(r.get("strain"))
    for r in csv.DictReader(open(a.cohort), delimiter="\t"):
        lin[r["sample"]] = main_lineage(r.get("lineage"))
    lin["GCF_000195955"] = "lineage4"
    return lin


def support(node):
    lab = node.label or ""
    try:
        return float(lab.split("/")[-1])
    except ValueError:
        return float("nan")


def describe(name, tree, lin):
    tips = {l.taxon.label for l in tree.leaf_node_iter()}
    ing = tips - CANETTII
    mrca = tree.mrca(taxon_labels=sorted(ing))
    inside = {l.taxon.label for l in mrca.leaf_iter()}
    out = sorted(tips - inside)
    print(f"== {name}: {len(tips)} tips; MTBC clade {len(inside)} tips "
          f"(support {support(mrca):g}); outside it: {', '.join(out)}"
          f"{'  OK' if set(out) == (CANETTII & tips) else '  ** MTBC CLADE WRONG **'}")
    by = collections.defaultdict(set)
    for t in ing:
        by[lin.get(t, "?")].add(t)
    mono = {}
    for L, ts in sorted(by.items()):
        if len(ts) < 2 or L == "?":
            continue
        m = tree.mrca(taxon_labels=sorted(ts))
        extra = {l.taxon.label for l in m.leaf_iter()} - ts
        mono[L] = not extra
        print(f"   {L:<10} {len(ts):>4} tips  "
              f"{'monophyletic' if not extra else f'NOT monophyletic ({len(extra)} others inside)'}"
              f"  support {support(m):g}")
    l147 = set().union(*(by.get(f"lineage{i}", set()) for i in (1, 2, 3, 4, 7)))
    m = tree.mrca(taxon_labels=sorted(l147))
    extra = {l.taxon.label for l in m.leaf_iter()} - l147
    print(f"   lineages 1-4,7: {len(l147)} tips; "
          f"{'one clade' if not extra else f'NOT one clade ({len(extra)} others inside: ' + ', '.join(sorted(collections.Counter(lin.get(x, '?') for x in extra))) + ')'}"
          f"  support {support(m):g}")
    kids = m.child_nodes()
    first = [sorted(collections.Counter(lin.get(l.taxon.label, "?")
                                        for l in k.leaf_iter())) for k in kids]
    print(f"   first split inside lineages 1-4,7: "
          + " | ".join("+".join(f) for f in first))
    return mono


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g", required=True)
    ap.add_argument("--k", required=True)
    ap.add_argument("--panel-lineages", required=True)
    ap.add_argument("--cohort", required=True)
    a = ap.parse_args()
    lin = load_lineages(a)
    tns = dendropy.TaxonNamespace()
    g = dendropy.Tree.get(path=a.g, schema="newick", taxon_namespace=tns,
                          preserve_underscores=True, rooting="force-rooted")
    k = dendropy.Tree.get(path=a.k, schema="newick", taxon_namespace=tns,
                          preserve_underscores=True, rooting="force-rooted")
    for t in (g, k):
        t.reroot_at_edge(t.find_node_with_taxon_label("GCF_035581225").edge,
                         update_bipartitions=False)
    describe("G (canettii from the graph)", g, lin)
    describe("K (canettii by direct alignment)", k, lin)

    from dendropy.calculate import treecompare
    for t in (g, k):
        t.is_rooted = False
        t.encode_bipartitions()
    rf = treecompare.symmetric_difference(g, k)
    n_int = len([e for e in g.internal_edges()])
    print(f"\nRobinson-Foulds (unrooted): {rf} of {2 * n_int} possible "
          f"({rf / (2 * n_int):.2%})")
    bg = {e.bipartition: e.head_node for e in g.preorder_edge_iter()
          if e.head_node.is_internal()}
    bk = {e.bipartition: e.head_node for e in k.preorder_edge_iter()
          if e.head_node.is_internal()}
    for name, x, y in (("only in G", bg, bk), ("only in K", bk, bg)):
        sup = sorted(support(x[b]) for b in set(x) - set(y))
        hi = [s for s in sup if s >= 95]
        print(f"  splits {name}: {len(sup)}; UFBoot median "
              f"{sup[len(sup) // 2] if sup else float('nan'):g}; >= 95: {len(hi)}")
        for b in set(x) - set(y):
            if support(x[b]) >= 95:
                ts = [l.taxon.label for l in x[b].leaf_iter()]
                print(f"    support {support(x[b]):g}: {len(ts)} tips, lineages "
                      f"{dict(collections.Counter(lin.get(t, '?') for t in ts))}")


if __name__ == "__main__":
    main()
