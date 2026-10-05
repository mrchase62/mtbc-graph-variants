#!/usr/bin/env python3
"""Reconstruct the ancestral allele at every panel SNP site.

Item 14b. The merged VCF's REF is whichever base the reference carries, and
H37Rv is lineage 4.9 -- not ancestral and not central -- so at any site where
H37Rv holds the derived allele the polarity is inverted for every other
lineage. `ALT` currently means "not H37Rv", not "derived". VCF forbids changing
REF, but reserves the INFO key `AA` for exactly this, so the fix is an
annotation rather than a re-polarisation.

WHY THIS IS A PANEL PROPERTY. The ancestral state at a position is a fact about
the 333 genomes and their tree, not about whichever isolates are in a cohort.
It is therefore reconstructed once here and joined onto every merged VCF as a
lookup -- no per-cohort reconstruction, and no cohort needs re-running to gain
the annotation.

THE METHOD, and its limits. Fitch's algorithm on the rooted tree: one upward
pass computing state sets, then a downward step to the INGROUP node -- the
MTBC ancestor, the most recent common ancestor of every leaf that is not an
outgroup -- whose state is reported. AA is the state of the common ancestor of
the genomes a cohort is drawn from, which is what "ancestral" means for them.

WHY THE MRCA AND NOT THE ROOT'S OTHER CHILD (audit Fault B). The panel tree
has TWO canettii genomes. It is rooted on GCF_035581225, and the root's other
child is GCF_000253375 plus the MTBC, so taking "the root's child that is not
the outgroup" reported the canettii-plus-MTBC node: 132 of the sites variable
within lineages 1-4 came out differently from the MTBC ancestor's state. The
outgroups are therefore a list, the ingroup is the MRCA of everything else,
and the script stops if that MRCA has an outgroup beneath it.

WHY THE INGROUP AND NOT THE ROOT (review 3.7). The root is the canettii/MTBC
split and has one outgroup genome on one side, so its set ties whenever that
genome differs from the MTBC ancestor: 5,130 sites were reported `.`, and 4,759
of them have one unambiguous state at the ingroup node. The downward step is
Fitch's own: the ingroup takes its upward set if that is one state; otherwise
the outgroups polarise it, NEAREST FIRST -- the first outgroup with a call
whose state leaves exactly one of the ingroup's tied states decides. What is
left is genuinely tied inside the MTBC and is reported `.` with TIED, not
resolved by a coin toss.

This is point-estimate parsimony, deliberately. samarray's `AncestorArray`
computes the fraction of tied optimal reconstructions in which each node takes
each state, which is strictly more informative; adopting it is a larger
decision than this annotation, and the tied flag here carries the same warning
in the cases that matter.
"""
import argparse, collections, csv, sys

def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

def parse_newick(s):
    """Return (children dict, root, leaves). Enough for Fitch: topology only."""
    s = s.strip().rstrip(";")
    children, stack, node = collections.defaultdict(list), [], None
    counter = [0]
    def newnode():
        counter[0] += 1
        return f"_n{counter[0]}"
    i, cur = 0, None
    parent_stack = []
    token = ""
    root = None
    while i < len(s):
        c = s[i]
        if c == "(":
            n = newnode()
            if parent_stack:
                children[parent_stack[-1]].append(n)
            else:
                root = n
            parent_stack.append(n); token = ""
        elif c in ",)":
            name = token.split(":")[0].strip()
            if name:
                children[parent_stack[-1]].append(name)
            token = ""
            if c == ")":
                closed = parent_stack.pop()
                if parent_stack:
                    pass
                i += 1
                # consume the internal label / branch length
                while i < len(s) and s[i] not in "(),":
                    i += 1
                continue
        else:
            token += c
        i += 1
    if token.split(":")[0].strip() and parent_stack:
        children[parent_stack[-1]].append(token.split(":")[0].strip())
    leaves = {n for n in
              {x for v in children.values() for x in v} | ({root} if root else set())
              if n not in children}
    return children, root, leaves

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--outgroups", "--outgroup", dest="outgroups",
                    default="GCF_035581225,GCF_000253375",
                    help="comma-separated outgroup leaves (default the two "
                         "canettii). The ingroup whose state is reported is "
                         "the MRCA of every other leaf; the outgroups break a "
                         "tie there, nearest to it first")
    a = ap.parse_args()

    children, root, leaves = parse_newick(open(a.tree).read())
    outgroups = [x for x in a.outgroups.split(",") if x]
    absent = [x for x in outgroups if x not in leaves]
    if not outgroups or absent:
        sys.exit(f"FATAL: outgroups must be tree leaves; not found: {absent}")
    parent = {c: p for p, kids in children.items() for c in kids}

    def path_up(n):
        out = [n]
        while n in parent:
            n = parent[n]; out.append(n)
        return out

    def under(n):
        if n not in children:
            return {n}
        return set().union(*(under(c) for c in children[n]))

    # the ingroup: the deepest node on every non-outgroup leaf's path to root
    ing_leaves = sorted(leaves - set(outgroups))
    common = None
    for l in ing_leaves:
        p = set(path_up(l))
        common = path_up(l) if common is None else [n for n in common if n in p]
    if not common:
        sys.exit("FATAL: every leaf is an outgroup")
    ingroup = common[0]
    inside = under(ingroup) & set(outgroups)
    if inside:
        sys.exit(f"FATAL: the MRCA of the ingroup ({len(ing_leaves)} leaves) "
                 f"also contains the outgroup(s) {sorted(inside)}, so its state "
                 f"would not be the ingroup ancestor's. Is the tree rooted on "
                 f"an outgroup?")
    # nearest first: the outgroup that joins the ingroup's lineage lowest
    anc = path_up(ingroup)

    def join_depth(o):
        return next(i for i, n in enumerate(anc) if o in under(n))
    outgroups.sort(key=join_depth)          # stable: ties keep the list order
    print(f"  ingroup: MRCA of {len(ing_leaves)} leaves; outgroups, nearest "
          f"first: {', '.join(outgroups)}")
    seqs = read_fasta(a.alignment)
    missing = leaves - set(seqs)
    extra = set(seqs) - leaves
    print(f"  tree: {len(leaves)} leaves, root {root}")
    print(f"  alignment: {len(seqs)} sequences, {len(next(iter(seqs.values())))} sites")
    if missing:
        sys.exit(f"  {len(missing)} tree leaves absent from the alignment: "
                 f"{sorted(missing)[:5]}")
    if extra:
        print(f"  note: {len(extra)} alignment sequences are not tree leaves; ignored")

    sites = list(csv.DictReader(open(a.sites), delimiter="\t"))
    n = len(sites)
    order = sorted(leaves)
    cols = [seqs[l] for l in order]

    # post-order once, reused for every site
    post = []
    seen = set()
    stack = [(root, False)]
    while stack:
        node, done = stack.pop()
        if done:
            post.append(node); continue
        if node in seen: continue
        seen.add(node)
        stack.append((node, True))
        for ch in children.get(node, []):
            stack.append((ch, False))
    idx = {l: i for i, l in enumerate(order)}

    out = []
    tied = unamb = nocall = 0
    for s_i in range(n):
        state = {}
        for node in post:
            kids = children.get(node, [])
            if not kids:
                b = cols[idx[node]][s_i] if node in idx else "N"
                state[node] = set() if b == "N" else {b}
                continue
            sets = [state[k] for k in kids if state[k]]
            if not sets:
                state[node] = set(); continue
            inter = set.intersection(*sets)
            state[node] = inter if inter else set.union(*sets)
        ing = state.get(ingroup, set())
        if len(ing) > 1:
            for o in outgroups:           # the outgroups polarise a tie
                pol_ = ing & state.get(o, set())
                if len(pol_) == 1:
                    ing = pol_
                    break
        if len(ing) == 1:
            aa = next(iter(ing)); flag = ""; unamb += 1
        elif len(ing) > 1:
            aa = "."; flag = "TIED:" + "".join(sorted(ing)); tied += 1
        else:
            aa = "."; flag = "NODATA"; nocall += 1
        out.append((sites[s_i], aa, flag))

    with open(a.out, "w") as fh:
        fh.write("chrom\tpos\tref\talt\tAA\tflag\tpolarity\n")
        for site, aa, flag in out:
            if aa == site["ref"]:
                pol = "ref_ancestral"
            elif aa == site["alt"]:
                pol = "alt_ancestral"
            else:
                pol = "unknown"
            fh.write(f"{site['chrom']}\t{site['pos']}\t{site['ref']}\t"
                     f"{site['alt']}\t{aa}\t{flag}\t{pol}\n")

    pol = collections.Counter(
        "ref_ancestral" if aa == s["ref"] else
        "alt_ancestral" if aa == s["alt"] else "unknown"
        for s, aa, _ in out)
    print(f"\n  sites               {n:,}")
    print(f"    unambiguous       {unamb:,}")
    print(f"    tied in the MTBC  {tied:,}")
    print(f"    no data           {nocall:,}")
    print(f"  polarity:")
    for k in ("ref_ancestral", "alt_ancestral", "unknown"):
        print(f"    {k:16s} {pol[k]:,}")
    if pol["alt_ancestral"]:
        print(f"\n  {pol['alt_ancestral']:,} sites where the REFERENCE carries the "
              f"DERIVED allele -- {pol['alt_ancestral']/n:.1%} of sites, and the "
              f"reason AA is worth having")
    print(f"  -> {a.out}")

if __name__ == "__main__":
    main()
