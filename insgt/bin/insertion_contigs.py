#!/usr/bin/env python3
"""Turn the graph's insertions into an alignment target, one contig each.

WHY DEPTH WORKS HERE. Depth across an insertion POINT says nothing -- the
inserted sequence is not in the reference, so a carrier and a non-carrier show
the same coverage there. Depth ON the inserted sequence is a different question
with a clean answer: a carrier maps reads to it at roughly genome-median depth,
a non-carrier maps essentially none. The graph supplies the sequence, so the
question becomes answerable.

The IS6110 arm already proves the mechanism: H37Rv.isclean.fasta carries the
canonical element as its own contig so copy number reads as depth on that
contig over depth on the chromosome. This generalises it from one element to
every insertion the panel carries.

AND IT IS EASIER THAN THE ELEMENT ARM. The hard part there is that H37Rv holds
sixteen copies competing for the same reads, which is why they had to be
excised. An insertion absent from H37Rv has no competing copy by construction.

WHAT THIS SCRIPT DOES NOT DO. It writes the target; it does not genotype. And a
read mapping to a contig proves the SEQUENCE is present in the sample, not that
it sits at the coordinate the graph puts it at. Placement needs the junction,
which is the clip channel, not this.

THE FOUR THINGS THAT HAD TO BE SETTLED, and how they are settled here:

  length floor      a 50 bp contig cannot be placed uniquely, so --min-len
                    drops the short ones rather than pretending they are
                    genotypable. The count kept at each floor is printed so the
                    choice is made on the distribution rather than by taste.
  shared sequence   IS6110 family members, PE/PPE fragments and other mobile
                    sequence recur across many insertions. Left alone, reads
                    multimap and one event is counted many times -- the same
                    failure as the allele ladder. Identical and contained
                    sequences are collapsed here; near-identical ones are
                    flagged for the aligner's own mapping quality to handle,
                    and the flag is carried into the catalogue so a
                    multi-mapping contig is never read as an independent locus.
  normalisation     not this script's job, but the catalogue records each
                    contig's length so depth can be normalised per base.
  shared flanks     the sequence written is the inserted interval ONLY, taken
                    as ALT minus the anchoring REF base, never with flanking
                    reference attached.
"""
import argparse, collections, gzip, hashlib, sys


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph-vcf",
                    default="graphs/CX333.s10k.k23.K15/all_variants.collapsed.vcf.gz",
                    help="the collapsed graph VCF: one record per trimmed "
                         "allele, carriers unioned, so an insertion's "
                         "sequence is ALT minus its single anchor base")
    ap.add_argument("--min-len", type=int, default=300,
                    help="shorter contigs cannot be mapped uniquely by a short "
                         "read, so they are counted and dropped")
    ap.add_argument("--max-len", type=int, default=100000)
    ap.add_argument("--min-len-ratio", type=float, default=0.70,
                    help="two sequences only merge if the shorter is at least "
                         "this fraction of the longer. Without it a short "
                         "element sitting inside a long insertion merges with "
                         "it and the carrier set is inflated.")
    ap.add_argument("--min-containment", type=float, default=0.80,
                    help="31-mer containment at which two sequences are the "
                         "same insertion differing between carriers")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-table", required=True)
    a = ap.parse_args()

    op = gzip.open if a.graph_vcf.endswith(".gz") else open
    panel, raw = [], []
    with op(a.graph_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                panel = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            ref, alt = f[3], f[4]
            if "," in alt or alt.startswith("<"):
                continue
            d = len(alt) - len(ref)
            if d < 50:
                continue
            # The inserted interval only: a VCF insertion is left-anchored, so
            # the novel sequence is ALT minus the REF prefix. Attaching flanking
            # reference would pull reads from non-carriers at the edges.
            seq = alt[len(ref):].upper()
            if not seq or len(seq) != d:
                continue
            carriers = [panel[i] for i, g in enumerate(f[9:])
                        if g.split(":")[0] == "1"]
            raw.append(dict(pos=int(f[1]), length=d, seq=seq,
                            carriers=carriers))
    print(f"  {len(raw):,} insertion records >= 50 bp over {len(panel)} panel genomes")

    band = collections.Counter()
    for r in raw:
        L = r["length"]
        band["50-149" if L < 150 else "150-299" if L < 300 else
             "300-999" if L < 1000 else "1k-5k" if L < 5000 else ">5k"] += 1
    print("  by length: " + "  ".join(
        f"{k}={band[k]:,}" for k in ("50-149", "150-299", "300-999", "1k-5k", ">5k")))

    kept = [r for r in raw if a.min_len <= r["length"] <= a.max_len]
    print(f"  {len(kept):,} clear the {a.min_len} bp floor")

    # ---- cluster by SEQUENCE SIMILARITY, not exact containment ------------
    # The first version of this collapsed only exact duplicates and exact
    # substrings, and that was wrong in a way worth recording. TbD1 came out as
    # NINETEEN separate contigs, all named ins:1761791:2153, with carrier sets
    # of 32, 20, 13, 6, 3, 2, 2, 2, 2, 1, 1 and so on summing to the 94 panel
    # genomes that carry it. Two carriers' copies of a 2,153 bp insertion differ
    # by a SNP or two, so neither contains the other and exact matching keeps
    # them apart -- the allele-ladder failure from section 39 one level up.
    #
    # It also produced duplicate FASTA names, which would have made an invalid
    # alignment target.
    #
    # Clustering on shared k-mers fixes both. Two sequences whose 31-mer sets
    # overlap by --min-containment are the same insertion with sequence-level
    # variation between carriers, and the union of their carriers is the answer
    # the catalogue wants.
    K = 31
    def kmers(s):
        return {s[i:i + K] for i in range(len(s) - K + 1)}
    kept.sort(key=lambda r: -r["length"])
    prof = [kmers(r["seq"]) for r in kept]
    # index k-mers to candidate members, so each sequence is only compared with
    # sequences it actually shares a k-mer with
    idx = collections.defaultdict(set)
    groups, gk = [], []
    for i, r in enumerate(kept):
        cand = collections.Counter()
        for km in prof[i]:
            for j in idx.get(km, ()):
                cand[j] += 1
        # CONTAINMENT ALONE OVER-MERGES. At 80% containment with min() as the
        # denominator, a 1,358 bp IS6110 insertion merged into TbD1's 2,153 bp
        # cluster at the same position, because a short sequence sits inside a
        # long one and the metric cannot tell "the same event with a SNP" from
        # "a different, shorter element here too". The TbD1 cluster came out at
        # 133 panel carriers against the 94 that actually carry it.
        #
        # Requiring similar LENGTH as well separates the two: carriers' copies
        # of one insertion differ by SNPs and indels, not by 40% of their
        # length. Two events at one position stay two.
        best, best_c = None, 0.0
        for j, shared_n in cand.items():
            c = shared_n / max(1, min(len(prof[i]), len(gk[j])))
            lo, hi = sorted((kept[i]["length"], groups[j]["length"]))
            if hi and lo / hi < a.min_len_ratio:
                continue
            if c > best_c:
                best, best_c = j, c
        if best is not None and best_c >= a.min_containment:
            g = groups[best]
            g["carriers"] |= set(r["carriers"])
            g["members"] += 1
            g["lengths"].append(r["length"])
            continue
        groups.append(dict(seq=r["seq"], pos=r["pos"], length=r["length"],
                           carriers=set(r["carriers"]), members=1,
                           lengths=[r["length"]]))
        gk.append(prof[i])
        gi = len(groups) - 1
        for km in prof[i]:
            idx[km].add(gi)
    print(f"  {len(groups):,} clusters at >= {a.min_containment:.0%} 31-mer "
          f"containment ({len(kept) - len(groups):,} merged)")
    big = [g for g in groups if g["members"] >= 5]
    if big:
        print(f"    {len(big)} clusters merged 5 or more records; largest "
              f"{max(g['members'] for g in groups)}")

    # ---- flag clusters that still share sequence with another cluster -----
    # These are different events that a read could still map to either of --
    # element families, PE/PPE fragments. Not merged, but never to be read as
    # independent loci without mapping quality behind them.
    shared = collections.Counter()
    seen = collections.defaultdict(set)
    for i, g in enumerate(groups):
        s = g["seq"]
        for j in range(0, max(1, len(s) - 200 + 1), 50):
            seen[hashlib.blake2b(s[j:j + 200].encode(),
                                 digest_size=8).digest()].add(i)
    for members in seen.values():
        if len(members) > 1:
            for i in members:
                shared[i] += 1
    print(f"  {len(shared):,} clusters share a 200-mer with another and are "
          f"flagged multimap")

    with open(a.out_fasta, "w") as fa, open(a.out_table, "w") as tb:
        tb.write("contig\th37rv_pos\tlength\tn_panel_carriers\tmembers_merged\t"
                 "multimap\tlen_min\tlen_max\tpanel_carriers\n")
        # Names must be unique or the FASTA is not a usable target; the first
        # version emitted ins:1761791:2153 nineteen times.
        used = collections.Counter()
        order = sorted(range(len(groups)), key=lambda i: groups[i]["pos"])
        for i in order:
            g = groups[i]
            base = f"ins:{g['pos']}:{g['length']}"
            used[base] += 1
            name = base if used[base] == 1 else f"{base}#{used[base]}"
            mm = 1 if shared.get(i) else 0
            fa.write(f">{name}\n")
            for k in range(0, len(g["seq"]), 60):
                fa.write(g["seq"][k:k + 60] + "\n")
            tb.write(f"{name}\t{g['pos']}\t{g['length']}\t{len(g['carriers'])}\t"
                     f"{g['members']}\t{mm}\t{min(g['lengths'])}\t"
                     f"{max(g['lengths'])}\t{','.join(sorted(g['carriers']))}\n")
    print(f"  -> {a.out_fasta}\n  -> {a.out_table}")


if __name__ == "__main__":
    main()
