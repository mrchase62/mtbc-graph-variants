#!/usr/bin/env python3
"""Tree consistency of SV-like records, OLD vs NEW. Analysis only.

Tree: $SV_EVAL_TREE, default out/new1k.fasttree.nwk (FastTree -nt -gtr -gamma on the NEW VCF's
biallelic H37Rv-frame SNPs outside the repeat mask, all 1,000 isolates);
pruned to the 525 shared isolates for OLD vs NEW, and used whole for NEW.
Fitch small parsimony (unrooted count; a multifurcation costs k - max like
analysis/sv_homoplasy/count.py). Leaf state: ALT 1, REF 0, anything else
(NOCALL, ABSENT, DUP) unknown {0,1}.

Records: the classes of compare.py (DEL, INS, IS, ACC; H37Rv frame), counted
when informative: >= 2 ALT and >= 2 REF leaves (NOCALL/ABSENT are unknown,
so a record that is ALT or NOCALL everywhere scores 0 and is left out). Reported: change-count distribution and the share with >= 2
changes, per class and build.

SPLIT-1 effect: records grouped with the event writer's own rule
(assoc/bin/write_event_matrix.py split_sv_groups: same type, consecutive
starts within 50 bp, net lengths within 5%, not a copy-number step; the
truncated-small-allele -> graph-allele extension is not applied) over
DEL + INS records (OLD's CLASS=small >= 50 bp records are given CLASS=sv for
this, as NEW writes them); merged genotype = ALT if any member ALT, else REF
if any member REF, else unknown.
"""
import collections, os, pickle, sys
import numpy as np
R = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
sys.path.insert(0, f"{R}/bin"); sys.path.insert(0, f"{R}/assoc/bin")
sys.path.insert(0, f"{R}/analysis/sv_eval_1k")
from panel_ancestral import parse_newick  # noqa: E402
import write_event_matrix as wem  # noqa: E402
import pysam  # noqa: E402
D = f"{R}/analysis/sv_eval_1k/out"
G = pysam.FastaFile(f"{R}/runroot/refbias/build/3893a3762ed1/refs/GCF_000195955.fasta").fetch("NC_000962.3")


def klass(r):
    if r["chrom"] != "NC_000962.3": return None
    c, t = r["cls"], r["svtype"]
    if c == "accessory_presence": return "ACC"
    if c == "is6110" and t == "INS": return "IS"
    if t == "DEL" and c in ("sv", "small", "is6110"): return "DEL"
    if t == "INS" and c in ("sv", "small"): return "INS"
    return None


TREE = os.environ.get("SV_EVAL_TREE", f"{D}/new1k.fasttree.nwk")
SUFFIX = os.environ.get("SV_EVAL_SUFFIX", "")


def tree_for(keep):
    children, root, leaves = parse_newick(open(TREE).read())
    ch = {}
    def prune(nd):
        if nd not in children: return nd if nd in keep else None
        k = [x for x in (prune(c) for c in children[nd]) if x is not None]
        if not k: return None
        if len(k) == 1: return k[0]
        ch[nd] = k; return nd
    sys.setrecursionlimit(100000)
    rt = prune(root)
    post, st = [], [(rt, 0)]
    while st:
        nd, d = st.pop()
        if d: post.append(nd); continue
        st.append((nd, 1))
        for c in ch.get(nd, []): st.append((c, 0))
    return ch, post


def fitch(ch, post, leafstate):
    """leafstate: dict leaf -> 0/1/None. Returns change count."""
    sets, cost = {}, 0
    for nd in post:
        if nd not in ch:
            v = leafstate.get(nd); sets[nd] = 3 if v is None else (1 << v); continue
        ss = [sets[c] for c in ch[nd]]
        inter = 3
        for x in ss: inter &= x
        if inter: sets[nd] = inter
        else:
            cnt = [sum(1 for x in ss if x & 1), sum(1 for x in ss if x & 2)]
            m = max(cnt); sets[nd] = (1 if cnt[0] == m else 0) | (2 if cnt[1] == m else 0); cost += len(ss) - m
    return cost


def as_variant(r, i, cls_override=None):
    """The dict split_sv_groups / copy_number_step read."""
    if r["svtype"] == "INS":
        ref, alt = ("N", "N" + r["seq"]) if r["seq"] else ("N", "<INS>")
    else:
        resolved = r["cls"] == "small" or r["info"].get("SVSOURCE") == "small_allele"
        ref, alt = (G[r["pos"] - 1:r["end"]], G[r["pos"] - 1]) if resolved else ("N", "<DEL>")
    return dict(cls=cls_override or r["cls"], svtype=r["svtype"], key=("NC_000962.3", r["pos"]), chrom=r["chrom"],
                acclocus=r["info"].get("ACCLOCUS"), pos=r["pos"], svlen=r["svlen"], ref=ref, alt=alt, idx=i)


def run(d, samples, ch, post, label, out):
    """Rows (set, class, unit, idx, pos, svlen, carriers, refs, nocall_frac,
    changes, member_change_sum). A record or event is scored when it has
    >= 2 ALT and >= 2 REF leaves (informative); NOCALL / ABSENT leaves are
    unknown, so heavy no-calling lowers counts and is reported alongside."""
    leaves = set(ch_leaves(ch, post))
    keep = np.array([s in leaves for s in samples])
    samples = [s for s in samples if s in leaves]
    idx = np.array([d["samples"].index(s) for s in samples])
    st = d["st"][:, idx]
    recs = d["recs"]
    res = []
    by = collections.defaultdict(list)
    for i, r in enumerate(recs):
        k = klass(r)
        if not k: continue
        by[k].append(i)
    def leafst(row):
        return {s: (1 if x == 1 else 0 if x == 0 else None) for s, x in zip(samples, row)}
    rec_c = {}
    for k, L in by.items():
        for i in L:
            na, nr = int((st[i] == 1).sum()), int((st[i] == 0).sum())
            if na < 2 or nr < 2: continue
            c = fitch(ch, post, leafst(st[i])); rec_c[i] = c
            res.append((label, k, "record", i, recs[i]["pos"], recs[i]["svlen"], na, nr,
                        round(float((st[i] == 3).mean()), 3), c, c))
    # SPLIT-1 groups over DEL + INS (CLASS sv; OLD small given class sv)
    vs = []
    for k in ("DEL", "INS"):
        for i in by[k]:
            r = recs[i]
            if r["cls"] not in ("sv", "small"): continue
            vs.append(as_variant(r, i, "sv"))
    groups = wem.split_sv_groups(vs)
    for gidx in groups:
        members = [vs[j]["idx"] for j in gidx]
        m = st[members]
        row = np.where((m == 1).any(0), 1, np.where((m == 0).any(0), 0, 3))
        na, nr = int((row == 1).sum()), int((row == 0).sum())
        k = vs[gidx[0]]["svtype"]
        res_c = fitch(ch, post, leafst(row)) if na >= 2 and nr >= 2 else None
        if res_c is not None and len(members) > 1:
            for x in members:
                if x not in rec_c:
                    rec_c[x] = fitch(ch, post, leafst(st[x])) if (st[x] == 1).any() else 0
        msum = sum(rec_c.get(x, 0) for x in members)
        if res_c is not None:
            res.append((label, k, "event_split1" if len(members) > 1 else "event_single", members[0],
                        recs[members[0]]["pos"], recs[members[0]]["svlen"], na, nr,
                        round(float((row == 3).mean()), 3), res_c, msum))
        out.write(f"{label}\t{k}\t{len(members)}\t{','.join(str(recs[x]['pos']) + ':' + str(recs[x]['svlen']) for x in members)}\t{na}\t{res_c}\t{msum}\n")
    return res


def ch_leaves(ch, post):
    return [x for x in post if x not in ch]


def main():
    o = pickle.load(open(f"{D}/old.pkl", "rb")); n = pickle.load(open(f"{D}/new.pkl", "rb"))
    shared = sorted(set(o["samples"]) & set(n["samples"]))
    rows = []
    gf = open(f"{D}/split1_groups{SUFFIX}.tsv", "w")
    gf.write("set\ttype\tn_members\tmembers(pos:len)\tcarriers\tchanges\tmember_change_sum\n")
    ch, post = tree_for(set(shared))
    rows += run(o, shared, ch, post, "OLD_525", gf)
    rows += run(n, shared, ch, post, "NEW_525", gf)
    if not SUFFIX:
        ch, post = tree_for(set(n["samples"]))
        rows += run(n, n["samples"], ch, post, "NEW_1000", gf)
    gf.close()
    with open(f"{D}/homoplasy_records{SUFFIX}.tsv", "w") as f:
        f.write("set\tclass\tunit\tidx\tpos\tsvlen\tcarriers\trefs\tnocall_frac\tchanges\tmember_change_sum\n")
        for r in rows: f.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
