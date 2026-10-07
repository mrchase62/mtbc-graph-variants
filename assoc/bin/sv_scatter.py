#!/usr/bin/env python3
"""Tests A, B and E of assoc/SV_SCATTER_PLAN.md: why SV carriers are scattered.

Reads tables that already exist and recomputes the reconstruction directly from
them, rather than rebuilding the merged VCF three times.

  A  split the ALT cells by how they were earned -- `called` by a caller,
     `inherited` from the matched reference, `depth_absent` from coverage --
     and reconstruct with each set in turn. Which arm the scatter appears in
     says whether p5svgt caused it or revealed it.
  B  reference assignment as a character on the tree. An `inherited` ALT is a
     fact about which panel genome the sample was matched to, so if that
     assignment does not itself track the tree, anything inherited from it is
     scattered by construction.
  E  the fair null. Gains per carrier depends strongly on carrier count -- at
     two carriers it can only be 0.5 or 1.0 -- so SV rows are compared against
     small variants MATCHED on carrier count, and the consistency index is
     reported alongside because it is the standard homoplasy measure.

H37Rv is included as a leaf in the REFERENCE state for every row. That is not
an assumption: these rows are keyed on H37Rv coordinates, so a deletion means
sequence H37Rv has and the sample does not, and H37Rv by definition does not
carry it. It anchors the root, which is otherwise the weak point.
"""
import argparse, collections, csv, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from write_event_matrix import (load_tree, fitch, UNK, A_BIT, D_BIT)


def tree_frame(path):
    tree = load_tree(path)
    nodes = list(tree.preorder_node_iter())
    order = {n.label: i for i, n in enumerate(nodes)}
    # EVERY NODE NEEDS ITS OWN LABEL. The index is keyed on it, so two nodes
    # sharing a label -- or several with none, as a raw IQ-TREE tree has --
    # collapse to one entry and silently rewire children and parents.
    # labelled.nwk is unique by construction; anything else must be refused.
    if None in order or len(order) != len(nodes):
        sys.exit(f"FATAL: {path}: node labels are not unique "
                 f"({len(order)} distinct among {len(nodes)} nodes, "
                 f"{sum(1 for n in nodes if n.label is None)} unlabelled); "
                 f"pass the event matrix's labelled.nwk")
    children = [[order[c.label] for c in n.child_nodes()] for n in nodes]
    parent = [None if n.parent_node is None else order[n.parent_node.label]
              for n in nodes]
    is_leaf = [n.is_leaf() for n in nodes]
    labels = [n.label for n in nodes]
    leaves = [labels[i] for i in range(len(nodes)) if is_leaf[i]]
    leaf_index = [None] * len(nodes)
    for j, i in enumerate([i for i in range(len(nodes)) if is_leaf[i]]):
        leaf_index[i] = j
    root_i = next(i for i in range(len(nodes)) if parent[i] is None)
    post = list(range(len(nodes)))[::-1]
    pre = list(range(len(nodes)))
    return dict(nodes=nodes, labels=labels, children=children, parent=parent,
                is_leaf=is_leaf, leaves=leaves, leaf_index=leaf_index,
                root=root_i, post=post, pre=pre, n=len(nodes))


def events(T, bits):
    """Gains, losses and derived-leaf counts per column."""
    down = fitch(T["post"], T["pre"], T["children"], T["parent"],
                 T["leaf_index"], bits)
    n_n, n_v = down.shape
    gains = np.zeros(n_v, dtype=np.int32)
    losses = np.zeros(n_v, dtype=np.int32)
    for i in range(n_n):
        p = T["parent"][i]
        if p is None:
            continue
        pv, cv = down[p], down[i]
        det = ((pv == A_BIT) | (pv == D_BIT)) & ((cv == A_BIT) | (cv == D_BIT))
        gains += (det & (pv == A_BIT) & (cv == D_BIT))
        losses += (det & (pv == D_BIT) & (cv == A_BIT))
    leaf_rows = [i for i in range(n_n) if T["is_leaf"][i]]
    carriers = (down[leaf_rows, :] == D_BIT).sum(axis=0)
    unresolved = (down[T["root"], :] == UNK) | (down[T["root"], :] == 3)
    return gains, losses, carriers, unresolved


def ratio_stats(g, c, label, out):
    keep = (g >= 2)
    if not keep.any():
        out.append(f"  {label:<34} no testable rows")
        return
    r = (g[keep] / np.maximum(c[keep], 1))
    ci = 1.0 / np.maximum(g[keep] + 0, 1)      # min changes for a binary char = 1
    out.append(f"  {label:<34}{keep.sum():>7,}{np.median(r):>10.2f}"
               f"{(r < 0.5).mean():>9.0%}{np.median(ci):>10.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="assoc/scale200/events")
    ap.add_argument("--sv-states", default="refbias/scale200/p5/sv_states.tsv")
    ap.add_argument("--refmap", default="refbias/scale200/p1/refmap.tsv")
    ap.add_argument("--ref-sample", default="H37Rv")
    ap.add_argument("--out", default="assoc/SV_SCATTER_RESULTS.md")
    a = ap.parse_args()

    T = tree_frame(os.path.join(a.events, "labelled.nwk"))
    li = {l: j for j, l in enumerate(T["leaves"])}
    print(f"  tree: {T['n']} nodes, {len(T['leaves'])} leaves")

    rows = collections.defaultdict(dict)     # key -> sample -> (state, evidence)
    svtype = {}
    for r in csv.DictReader(open(a.sv_states), delimiter="\t"):
        rows[r["key"]][r["sample"]] = (r["state"], r["evidence"].split(":")[0])
        svtype[r["key"]] = r["svtype"]
    keys = sorted(rows)
    print(f"  sv_states: {len(keys):,} rows x "
          f"{len({s for v in rows.values() for s in v}):,} samples")

    def build(accept):
        """Leaf bit masks, taking ALT only from the accepted evidence classes."""
        b = np.zeros((len(T["leaves"]), len(keys)), dtype=np.uint8)
        if a.ref_sample in li:
            b[li[a.ref_sample], :] = A_BIT
        for j, k in enumerate(keys):
            for s, (st, ev) in rows[k].items():
                i = li.get(s)
                if i is None:
                    continue
                if st == "ALT":
                    b[i, j] = D_BIT if ev in accept else UNK
                elif st == "REF":
                    b[i, j] = A_BIT
                # ABSENT and NOCALL stay UNK
        return b

    out = []
    out.append("# Why SV carriers do not form clades — tests A, B and E\n")
    out.append(f"Run on `{a.sv_states}` and `{a.events}/labelled.nwk`, "
               f"{len(keys):,} SV rows over {len(T['leaves'])} tips "
               f"(the cohort, the H37Rv reference tip read as REF, and the "
               f"outgroup, which has no SV data).\n")

    # ---------- Test A
    out.append("## Test A — is it the genotyping, or was the callset "
               "always like this?\n")
    out.append("ALT cells are taken only from the evidence classes named, "
               "and an ALT that is withheld becomes no-data rather than REF, "
               "so the arms differ only in which ALTs they are allowed to "
               "see.\n")
    out.append(f"  {'arm':<34}{'testable':>7}{'med g/c':>10}"
               f"{'<0.5':>9}{'med CI':>10}")
    arms = [("called only", {"called"}),
            ("called + inherited", {"called", "inherited"}),
            ("called + depth_absent", {"called", "depth_absent"}),
            ("all three (the shipped matrix)",
             {"called", "inherited", "depth_absent"})]
    per_arm = {}
    for name, acc in arms:
        g, l, c, unres = events(T, build(acc))
        per_arm[name] = (g, l, c, unres)
        ratio_stats(g, c, name, out)
    out.append("")
    for name, acc in arms:
        g, l, c, unres = per_arm[name]
        out.append(f"  {name:<34} rows with >=1 carrier {int((c>0).sum()):>6,}"
                   f"   root unresolved {int(unres.sum()):>6,}"
                   f"   total gains {int(g.sum()):>7,}")
    out.append("")

    # ---------- Test B
    out.append("## Test B — does the pattern follow reference assignment?\n")
    ref = {}
    for r in csv.DictReader(open(a.refmap), delimiter="\t"):
        ref[r["sample"]] = r.get("reference") or list(r.values())[4]
    refs = sorted({v for k, v in ref.items() if k in li})
    out.append(f"{len(refs)} distinct matched references across "
               f"{sum(1 for k in ref if k in li)} cohort tips.\n")

    # parsimony score of reference assignment as a multistate character
    idx = {r: i for i, r in enumerate(refs)}
    sets = [None] * T["n"]
    for i in range(T["n"]):
        if T["is_leaf"][i]:
            r = ref.get(T["labels"][i])
            sets[i] = {idx[r]} if r in idx else set()
    score = 0
    for i in T["post"]:
        if T["is_leaf"][i]:
            continue
        kid = [sets[c] for c in T["children"][i] if sets[c]]
        if not kid:
            sets[i] = set(); continue
        inter = set.intersection(*kid)
        if inter:
            sets[i] = inter
        else:
            sets[i] = set.union(*kid); score += 1
    out.append(f"**Reference assignment needs {score} changes on this tree** "
               f"for {len(refs)} states. A perfectly clade-consistent "
               f"assignment would need {len(refs) - 1}; the excess is "
               f"{score - (len(refs) - 1)}.\n")

    g, l, c, unres = per_arm["all three (the shipped matrix)"]
    nref, ninh = [], []
    for j, k in enumerate(keys):
        car = [s for s, (st, ev) in rows[k].items() if st == "ALT"]
        nref.append(len({ref.get(s) for s in car if s in ref}))
        ninh.append(sum(1 for s in car if rows[k][s][1] == "inherited"))
    nref, ninh = np.asarray(nref), np.asarray(ninh)
    keep = g >= 2
    out.append(f"  {'':<34}{'rows':>7}{'med gains':>11}"
               f"{'med distinct refs':>19}{'corr':>8}")
    for nm, m in (("all testable rows", keep),
                  ("no inherited carriers", keep & (ninh == 0)),
                  ("some inherited carriers", keep & (ninh > 0))):
        if not m.any():
            out.append(f"  {nm:<34} none"); continue
        cc = np.corrcoef(g[m], nref[m])[0, 1] if m.sum() > 2 else float("nan")
        out.append(f"  {nm:<34}{int(m.sum()):>7,}{np.median(g[m]):>11.1f}"
                   f"{np.median(nref[m]):>19.1f}{cc:>8.2f}")
    out.append("")

    # ---------- Test E
    out.append("## Test E — the null, matched on carrier count\n")
    sm = [r for r in csv.DictReader(
        open(os.path.join(a.events, "variants.tsv")), delimiter="\t")
        if r["class"] == "small"]
    by_k = collections.defaultdict(list)
    for r in sm:
        k = int(r["n_derived_leaves"]); gg = int(r["n_gain"])
        if gg >= 2 and k > 0:
            by_k[k].append(gg / k)
    sv_by_k = collections.defaultdict(list)
    for j in range(len(keys)):
        if g[j] >= 2 and c[j] > 0:
            sv_by_k[int(c[j])].append(g[j] / c[j])
    bins = [(2, 2), (3, 4), (5, 9), (10, 19), (20, 49), (50, 999)]
    out.append(f"  {'carriers':<12}{'SV n':>7}{'SV med':>9}"
               f"{'small n':>9}{'small med':>11}{'diff':>8}")
    for lo, hi in bins:
        sv = [x for k, v in sv_by_k.items() if lo <= k <= hi for x in v]
        s2 = [x for k, v in by_k.items() if lo <= k <= hi for x in v]
        if not sv and not s2:
            continue
        a1 = np.median(sv) if sv else float("nan")
        a2 = np.median(s2) if s2 else float("nan")
        out.append(f"  {f'{lo}-{hi}':<12}{len(sv):>7,}{a1:>9.2f}"
                   f"{len(s2):>9,}{a2:>11.2f}{a1 - a2:>8.2f}")
    out.append("")
    open(a.out, "w").write("\n".join(out) + "\n")
    print("\n".join(out))
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
