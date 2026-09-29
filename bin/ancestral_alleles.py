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
pass computing state sets, then the root's set read as the ancestral state.
Where the root set holds one state the call is unambiguous; where it holds two
the site is genuinely tied under parsimony and `AA` is reported as `.` with
`AA_TIED` rather than being resolved by a coin toss. Tied sites are the ones
where a reader most needs to know the reconstruction is uncertain, which is
why they are not silently broken.

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
    a = ap.parse_args()

    children, root, leaves = parse_newick(open(a.tree).read())
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
        r = state.get(root, set())
        if len(r) == 1:
            aa = next(iter(r)); flag = ""; unamb += 1
        elif len(r) > 1:
            aa = "."; flag = "TIED:" + "".join(sorted(r)); tied += 1
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
    print(f"    tied at the root  {tied:,}")
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
