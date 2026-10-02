#!/usr/bin/env python3
"""Which panel genomes traverse each graph node? A property of the graph alone.

WHY THIS IS THE FIX. `bin/p5_states.py` writes NOCALL for every node-frame key a
sample produced no record for, and its own comment says why: "its own reference
may or may not traverse that node, and this pass has not asked". The census in
accessory/ACCESSORY_CALL_CENSUS.md measures what that costs -- 91 REF calls in
6.9 million accessory cells, and 40 of 63 loci with zero REF anywhere.

This asks. Node-to-path membership depends on the graph and nothing else, so it
is cohort-independent and computed once. With it, a sample whose reference does
NOT traverse a node is ABSENT -- a measured statement, the same one the
H37Rv-frame keys already get when a position does not project -- instead of
no-call. That alone turns most of those cells into observations.

It does not settle the samples whose reference DOES traverse the node; those
need depth at the corresponding position, which is the second half of the fix.

ONLY THE NODES THAT MATTER. The graph has 84,171 segments; 20,299 of them carry
a node-frame VCF record. Restricting to those keeps the table small enough to
load per task, which is what lets p5_states.py use it without a projection.
"""
import argparse, collections, gzip, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", default="")
    ap.add_argument("--nodes", default="",
                    help="one node id per line: the nodes to record. Default is "
                         "every node, which is larger but complete.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-positions", default="",
                    help="also write node,path,start -- the node's offset along "
                         "each path that traverses it. That is what lets a "
                         "node-frame key be looked up in a sample's own "
                         "reference frame, where its gVCF depth decides REF.")
    a = ap.parse_args()

    gfa = a.gfa
    if not gfa:
        import glob
        g = sorted(glob.glob("graphs/CX333.s10k.k23.K15/*.smooth.final.gfa"))
        gfa = g[0] if g else ""
    if not gfa or not os.path.exists(gfa):
        sys.exit("FATAL: no GFA; pass --gfa")

    want = None
    if a.nodes:
        want = {l.strip() for l in open(a.nodes) if l.strip()}
        print(f"  restricting to {len(want):,} nodes")

    # node -> list of path names. Built from P lines, whose walk is a comma
    # separated list of <node><orientation>.
    # Node lengths, needed only if positions are asked for. A path's walk is an
    # ordered node list, so a node's offset along the path is the sum of the
    # lengths before it -- no projection and no aligner involved.
    nlen = {}
    if a.out_positions:
        with open(gfa) as fh:
            for line in fh:
                if line[0] != "S":
                    continue
                f = line.split("\t", 3)
                nlen[f[1]] = len(f[2].strip())
        print(f"  {len(nlen):,} segment lengths read")

    mem = collections.defaultdict(list)
    pos_rows = []
    npath = 0
    with open(gfa) as fh:
        for line in fh:
            if line[0] != "P":
                continue
            f = line.split("\t", 3)
            name = f[1]
            npath += 1
            off = 0
            for tok in f[2].split(","):
                nid = tok[:-1] if tok and tok[-1] in "+-" else tok
                if want is None or nid in want:
                    mem[nid].append(name)
                    if a.out_positions:
                        # 1-based start, and the strand, because a node walked
                        # in reverse has its offsets counted from the other end
                        # and a caller that ignores that lands on the wrong base
                        pos_rows.append((nid, name.split("#")[0], off + 1,
                                         tok[-1] if tok and tok[-1] in "+-" else "+"))
                if a.out_positions:
                    off += nlen.get(nid, 0)
    print(f"  {npath} paths walked, {len(mem):,} nodes recorded")

    counts = collections.Counter(len(v) for v in mem.values())
    tot = sum(counts.values()) or 1
    print("\n  how many panel genomes traverse a node:")
    for lab, lo, hi in (("1 only", 1, 1), ("2-9", 2, 9), ("10-99", 10, 99),
                        ("100-300", 100, 300), ("301+", 301, 10**9)):
        n = sum(v for k, v in counts.items() if lo <= k <= hi)
        print(f"    {lab:<10}{n:>8,}  {n/tot:6.1%}")

    with open(a.out, "w") as fh:
        fh.write("node\tn_paths\tpaths\n")
        for nid, paths in sorted(mem.items(), key=lambda kv: int(kv[0])
                                 if kv[0].isdigit() else 0):
            # accessions only, dropping the PanSN suffix, so a caller can test
            # its refmap value directly
            accs = sorted({p.split("#")[0] for p in paths})
            fh.write(f"{nid}\t{len(accs)}\t{','.join(accs)}\n")
    print(f"\n  -> {a.out}")
    if a.out_positions:
        # One row per (node, accession) pair. A node a path visits more than
        # once gets its FIRST occurrence: a repeat cannot be genotyped from one
        # position anyway, and silently keeping the last would be arbitrary.
        # COUNT THE OCCURRENCES, do not just keep the first. A node a path
        # visits more than once is a repeat in that genome, and depth at one
        # copy says nothing about the state of the locus -- the reads are shared
        # between copies. So the count is carried and the caller is expected to
        # refuse REF where it exceeds one, the same reasoning that makes
        # MAPQ-filtered depth necessary for element-proximal intervals.
        occ = collections.Counter((nid, acc) for nid, acc, _, _ in pos_rows)
        first = {}
        for nid, acc, start, strand in pos_rows:
            first.setdefault((nid, acc), (start, strand))
        multi = sum(1 for v in occ.values() if v > 1)
        with open(a.out_positions, "w") as fh:
            # `length` is the node's own length, recorded so a reader can
            # check that an offset falls inside the node. (P4's node offsets
            # already run along the path, so the base is start + offset on
            # either strand; see p5_states.py.)
            fh.write("node\taccession\tstart\tstrand\tn_occurrences\tlength\n")
            for (nid, acc), (start, strand) in first.items():
                fh.write(f"{nid}\t{acc}\t{start}\t{strand}\t{occ[(nid, acc)]}"
                         f"\t{nlen.get(nid, '')}\n")
        print(f"  {len(first):,} (node, accession) pairs; {multi:,} "
              f"({multi/max(len(first),1):.1%}) are repeated within their own "
              f"path and cannot be genotyped from one position")
        print(f"  -> {a.out_positions}")


if __name__ == "__main__":
    main()
