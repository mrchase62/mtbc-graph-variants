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

A BUILD ASSET NOW (audit P4P5-7). P0 (`p0_prepare.sh --step nodes`) runs this
on the build's own graph with --exclude-path <H37Rv path>: every node-frame key
lies on a node H37Rv does not traverse, so that restriction is complete for
any cohort, where the hand-made node list covered only the cohorts it was made
from. On CX333: 95,828 nodes and 1,379,748 (node, accession) pairs, covering
every node-frame key of scale200 (11,533) and gwas1000 (21,743), and equal to
the hand-made tables on all 21,400 nodes and 363,898 pairs they share.
"""
import argparse, collections, gzip, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True,
                    help="the build's own graph as GFA (P0 writes it with "
                         "`odgi view -g` from the build's .og). There is no "
                         "default: a glob over graphs/CX333... handed every "
                         "build the CX333 node ids (audit P4P5-7)")
    ap.add_argument("--nodes", default="",
                    help="one node id per line: the nodes to record. Default is "
                         "every node, which is larger but complete.")
    ap.add_argument("--exclude-path", default="",
                    help="leave out every node this path traverses. P0 passes "
                         "the H37Rv path: a node-frame key is by construction "
                         "on a node H37Rv does not traverse (p4_place.py "
                         "writes one only where the H37Rv distance is not "
                         "0), so this keeps the table complete for P5 while "
                         "dropping the shared backbone")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-positions", default="",
                    help="also write node,path,start -- the node's offset along "
                         "each path that traverses it. That is what lets a "
                         "node-frame key be looked up in a sample's own "
                         "reference frame, where its gVCF depth decides REF.")
    a = ap.parse_args()

    gfa = a.gfa
    if not os.path.exists(gfa):
        sys.exit(f"FATAL: no GFA at {gfa}")

    want = None
    if a.nodes:
        want = {l.strip() for l in open(a.nodes) if l.strip()}
        print(f"  restricting to {len(want):,} nodes")

    # Node lengths, needed only if positions are asked for. A path's walk is an
    # ordered node list, so a node's offset along the path is the sum of the
    # lengths before it -- no projection and no aligner involved. The excluded
    # path's nodes are collected in the same first pass.
    nlen = {}
    skip = set()
    if a.out_positions or a.exclude_path:
        found = False
        with open(gfa) as fh:
            for line in fh:
                if line[0] == "S" and a.out_positions:
                    f = line.split("\t", 3)
                    nlen[f[1]] = len(f[2].strip())
                elif line[0] == "P" and a.exclude_path:
                    f = line.split("\t", 3)
                    if f[1] == a.exclude_path:
                        found = True
                        skip = {t[:-1] if t and t[-1] in "+-" else t
                                for t in f[2].split(",")}
        if a.out_positions:
            print(f"  {len(nlen):,} segment lengths read")
        if a.exclude_path:
            if not found:
                sys.exit(f"FATAL: path {a.exclude_path} is not in {gfa}")
            print(f"  excluding {len(skip):,} nodes traversed by "
                  f"{a.exclude_path}")

    def keep(nid):
        return (want is None or nid in want) and nid not in skip

    # node -> list of path names. Built from P lines, whose walk is a comma
    # separated list of <node><orientation>.
    mem = collections.defaultdict(list)
    npath = 0
    # POSITIONS ARE WRITTEN ONE ACCESSION AT A TIME. Holding every (node,
    # accession) row of a whole graph at once is tens of millions of tuples;
    # one accession's rows are all that a first occurrence and an occurrence
    # count need, provided an accession's paths are adjacent, which is checked.
    pos_fh = open(a.out_positions + ".tmp", "w") if a.out_positions else None
    if pos_fh:
        # `length` is the node's own length, recorded so a reader can
        # check that an offset falls inside the node, and so it can turn a
        # key's offset into a walk offset. (A node key's offset is the
        # node's FORWARD offset -- P4 and the IS6110 arm restate odgi's
        # walking offset -- so the base along this path is start + off where
        # the path walks the node '+', and start + (length-1-off) where it
        # walks it '-': mtb_norm.forward_offset, as p5_states.py reads it.)
        pos_fh.write("node\taccession\tstart\tstrand\tn_occurrences\tlength\n")
    npairs = multi = 0
    cur, first, occ, done_accs = None, {}, collections.Counter(), set()

    def flush():
        # One row per (node, accession) pair. A node a path visits more than
        # once gets its FIRST occurrence: a repeat cannot be genotyped from one
        # position anyway, and silently keeping the last would be arbitrary.
        # COUNT THE OCCURRENCES, do not just keep the first. A node a path
        # visits more than once is a repeat in that genome, and depth at one
        # copy says nothing about the state of the locus -- the reads are shared
        # between copies. So the count is carried and the caller is expected to
        # refuse REF where it exceeds one, the same reasoning that makes
        # MAPQ-filtered depth necessary for element-proximal intervals.
        nonlocal npairs, multi
        for nid, (start, strand) in first.items():
            pos_fh.write(f"{nid}\t{cur}\t{start}\t{strand}\t{occ[nid]}"
                         f"\t{nlen.get(nid, '')}\n")
            npairs += 1
            multi += occ[nid] > 1
        first.clear(); occ.clear()

    with open(gfa) as fh:
        for line in fh:
            if line[0] != "P":
                continue
            f = line.split("\t", 3)
            name = f[1]
            acc = name.split("#")[0]
            npath += 1
            if pos_fh and acc != cur:
                if cur is not None:
                    flush()
                    done_accs.add(cur)
                if acc in done_accs:
                    sys.exit(f"FATAL: the paths of {acc} are not adjacent in "
                             f"{gfa}; occurrence counts would be split")
                cur = acc
            off = 0
            for tok in f[2].split(","):
                nid = tok[:-1] if tok and tok[-1] in "+-" else tok
                if keep(nid):
                    mem[nid].append(name)
                    if pos_fh:
                        # 1-based start, and the strand, because a node walked
                        # in reverse has its offsets counted from the other end
                        # and a caller that ignores that lands on the wrong base
                        occ[nid] += 1
                        if nid not in first:
                            first[nid] = (off + 1, tok[-1] if tok and tok[-1]
                                          in "+-" else "+")
                if pos_fh:
                    off += nlen.get(nid, 0)
    if pos_fh:
        if cur is not None:
            flush()
        pos_fh.close()
    print(f"  {npath} paths walked, {len(mem):,} nodes recorded")

    counts = collections.Counter(len(v) for v in mem.values())
    tot = sum(counts.values()) or 1
    print("\n  how many panel genomes traverse a node:")
    for lab, lo, hi in (("1 only", 1, 1), ("2-9", 2, 9), ("10-99", 10, 99),
                        ("100-300", 100, 300), ("301+", 301, 10**9)):
        n = sum(v for k, v in counts.items() if lo <= k <= hi)
        print(f"    {lab:<10}{n:>8,}  {n/tot:6.1%}")

    with open(a.out + ".tmp", "w") as fh:
        fh.write("node\tn_paths\tpaths\n")
        for nid, paths in sorted(mem.items(), key=lambda kv: int(kv[0])
                                 if kv[0].isdigit() else 0):
            # accessions only, dropping the PanSN suffix, so a caller can test
            # its refmap value directly
            accs = sorted({p.split("#")[0] for p in paths})
            fh.write(f"{nid}\t{len(accs)}\t{','.join(accs)}\n")
    os.replace(a.out + ".tmp", a.out)
    print(f"\n  -> {a.out}")
    if a.out_positions:
        os.replace(a.out_positions + ".tmp", a.out_positions)
        print(f"  {npairs:,} (node, accession) pairs; {multi:,} "
              f"({multi/max(npairs,1):.1%}) are repeated within their own "
              f"path and cannot be genotyped from one position")
        print(f"  -> {a.out_positions}")


if __name__ == "__main__":
    main()
