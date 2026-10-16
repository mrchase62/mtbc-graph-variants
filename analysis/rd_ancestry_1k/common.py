"""Shared inputs for analysis/rd_ancestry_1k (analysis only, 2026-10-16).

- records: analysis/sv_eval_1k/out/new.pkl (extract.py over
  runroot/refbias/new1k_p200/p5/merged.vcf.gz: every non-small record plus
  CLASS=small indels >= 50 bp; per cell ST code 0 REF 1 ALT 2 ABSENT 3 NOCALL
  4 DUP, and an EV code).
- tree: analysis/sv_eval_1k/out/new1k.fasttree.nwk (unrooted FastTree),
  rooted here on the branch that separates lineage 1 from lineages 2-4.
- vectorised Fitch parsimony (binary states, unknown leaves = {0,1}),
  up pass with the k - max rule at multifurcations, then a down pass from a
  chosen root state that keeps the parent's state where the child's set
  allows it; gains (0 -> 1) and losses (1 -> 0) are counted on that path.
"""
import collections, pickle, sys
import numpy as np

R = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
B = f"{R}/runroot/refbias/build/3893a3762ed1"
ANN = "/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/annotation"
OUT = f"{R}/analysis/rd_ancestry_1k/out"
TREE = f"{R}/analysis/sv_eval_1k/out/new1k.fasttree.nwk"
ST = {0: "REF", 1: "ALT", 2: "ABSENT", 3: "NOCALL", 4: "DUP", 9: "MISSING"}


def load_records():
    return pickle.load(open(f"{R}/analysis/sv_eval_1k/out/new.pkl", "rb"))


def klass(r):
    """Analysis class of a record."""
    c, t = r["cls"], r["svtype"]
    if r["chrom"] != "NC_000962.3":
        return "IS6110_node" if c == "is6110" else None
    if c == "accessory_presence": return "ACC"
    if c == "is6110" and t == "INS": return "IS6110"
    if c == "is6110" and t == "DEL": return "IS6110_loss"
    if t == "DEL" and c in ("sv", "small"): return "DEL"
    if t == "INS" and c in ("sv", "small"): return "INS"
    if t == "REPL": return "REPL"
    return None


def component(r):
    i = r["info"]
    if i.get("COMPONENT"): return i["COMPONENT"]
    if i.get("SVSOURCE") == "small_allele" or r["cls"] == "small": return "small_allele"
    return "."


def cohort():
    lin, deep = {}, {}
    for l in open(f"{R}/analysis/cohort_1k/cohort_1k.tsv"):
        f = l.rstrip("\n").split("\t")
        if f[0] == "sample": continue
        lin[f[0]], deep[f[0]] = f[1], f[2]
    return lin, deep


def refmap():
    m = {}
    hdr = None
    for l in open(f"{R}/runroot/refbias/new1k_p200/p1/refmap.tsv"):
        f = l.rstrip("\n").split("\t")
        if hdr is None: hdr = f; continue
        m[f[0]] = f[hdr.index("reference")]
    return m


def mask_array():
    m = np.zeros(4411600, bool)
    for l in open(f"{B}/assets/repeat_mask.bed"):
        c, s, e = l.split("\t")[:3]; m[int(s):int(e)] = True
    return m


def genes():
    g = []
    for l in open(f"{ANN}/H37Rv_snpeff_dump.txt"):
        f = l.rstrip("\n").split("\t")
        if len(f) > 6 and f[4] == "Gene":
            g.append((int(f[1]), int(f[2]), f[5], f[6]))
    g.sort()
    return g


def is6110_h37rv():
    out = []
    for l in open(f"{B}/assets/is6110_elements.gff"):
        if l.startswith("#"): continue
        f = l.split("\t")
        if len(f) > 4: out.append((int(f[3]), int(f[4])))
    return sorted(out)


# ---------------------------------------------------------------- tree
def read_unrooted(path=TREE):
    """Undirected adjacency over leaves and internal nodes."""
    s = open(path).read().strip().rstrip(";")
    adj = collections.defaultdict(list)
    stack, tok, n, i, top = [], "", 0, 0, None
    while i < len(s):
        c = s[i]
        if c == "(":
            n += 1; nd = f"_n{n}"
            if stack: adj[stack[-1]].append(nd); adj[nd].append(stack[-1])
            else: top = nd
            stack.append(nd); tok = ""
        elif c in ",)":
            name = tok.split(":")[0].strip()
            if name:
                adj[stack[-1]].append(name); adj[name].append(stack[-1])
            tok = ""
            if c == ")":
                stack.pop(); i += 1
                while i < len(s) and s[i] not in "(),": i += 1
                continue
        else:
            tok += c
        i += 1
    return dict(adj), top


def root_on_split(adj, groupA):
    """Root on the edge whose one side best matches leaf set groupA
    (fewest leaves misplaced). Returns (children, root, post, info)."""
    sys.setrecursionlimit(100000)
    leaves = [x for x in adj if not x.startswith("_n")]
    A = set(groupA) & set(leaves)
    start = next(x for x in adj if x.startswith("_n"))
    # orient from start, compute per-node subtree leaf counts of A and total
    parent = {start: None}; order = [start]
    for nd in order:
        for c in adj[nd]:
            if c not in parent: parent[c] = nd; order.append(c)
    na, nt = {}, {}
    for nd in reversed(order):
        kids = [c for c in adj[nd] if parent.get(c) == nd]
        if not kids:
            na[nd] = int(nd in A); nt[nd] = 1
        else:
            na[nd] = sum(na[c] for c in kids); nt[nd] = sum(nt[c] for c in kids)
    TA, TT = len(A), len(leaves)
    best = None
    for nd in order[1:]:
        # edge (parent[nd], nd): side below nd vs rest
        a, t = na[nd], nt[nd]
        err1 = (TA - a) + (t - a)            # below == A
        err2 = a + ((TT - t) - (TA - a))     # rest == A
        e = min(err1, err2)
        if best is None or e < best[0]: best = (e, parent[nd], nd)
    e, u, v = best
    ch = {"ROOT": [u, v]}
    seen = {"ROOT", u, v}
    q = [(u, v), (v, u)]
    while q:
        nd, frm = q.pop()
        kids = [c for c in adj[nd] if c != frm]
        if kids: ch[nd] = kids
        for c in kids: q.append((c, nd))
    post, st = [], [("ROOT", 0)]
    while st:
        nd, d = st.pop()
        if d: post.append(nd); continue
        st.append((nd, 1))
        for c in ch.get(nd, []): st.append((c, 0))
    return ch, "ROOT", post, dict(misplaced=e, groupA=TA, leaves=TT)


class Fitch:
    """Vectorised Fitch over many binary characters on one rooted tree."""

    def __init__(self, ch, root, post, samples):
        self.ch, self.root, self.post = ch, root, post
        self.leaf_col = {s: j for j, s in enumerate(samples)}
        self.leaves = [x for x in post if x not in ch]
        missing = [x for x in self.leaves if x not in self.leaf_col]
        if missing: raise SystemExit(f"tree leaves not in VCF: {missing[:5]}")

    def run(self, st, root_pref=None):
        """st: (nrec, nsample) ST codes. root_pref: per record 0/1 used when
        the root set is tied (default 0). Returns dict of arrays."""
        nrec = st.shape[0]
        S = {}
        cost = np.zeros(nrec, np.int32)
        for nd in self.post:
            if nd not in self.ch:
                x = st[:, self.leaf_col[nd]]
                S[nd] = np.where(x == 1, 2, np.where(x == 0, 1, 3)).astype(np.uint8)
                continue
            ss = [S[c] for c in self.ch[nd]]
            inter = ss[0].copy()
            for x in ss[1:]: inter &= x
            c0 = sum((x & 1) > 0 for x in ss).astype(np.int32)
            c1 = sum((x & 2) > 0 for x in ss).astype(np.int32)
            m = np.maximum(c0, c1)
            uni = ((c0 == m) * 1 | (c1 == m) * 2).astype(np.uint8)
            empty = inter == 0
            S[nd] = np.where(empty, uni, inter)
            cost += np.where(empty, len(ss) - m, 0)
        rs = S[self.root]
        tied = rs == 3
        pref = np.zeros(nrec, np.uint8) if root_pref is None else np.asarray(root_pref, np.uint8)
        rstate = np.where(rs == 1, 0, np.where(rs == 2, 1, pref)).astype(np.uint8)
        state = {self.root: rstate}
        gains = np.zeros(nrec, np.int32); losses = np.zeros(nrec, np.int32)
        stack = [self.root]   # pre-order traversal
        while stack:
            p = stack.pop()
            for c in self.ch.get(p, []):
                ps = state[p]
                bit = (S[c] >> ps) & 1
                other = np.where(S[c] & 1, 0, 1).astype(np.uint8)
                cs = np.where(bit == 1, ps, other).astype(np.uint8)
                state[c] = cs
                gains += (ps == 0) & (cs == 1)
                losses += (ps == 1) & (cs == 0)
                stack.append(c)
        self.state = state
        return dict(cost=cost, root_set=rs, root_tied=tied, root_state=rstate,
                    gains=gains, losses=losses)

    def leaves_below(self):
        memo = {}
        for nd in self.post:
            memo[nd] = [nd] if nd not in self.ch else [x for c in self.ch[nd] for x in memo[c]]
        return memo
