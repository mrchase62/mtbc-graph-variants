#!/usr/bin/env python3
"""Cut the cohort + panel tree down to the tips the event writer reads.

The combined tree (assoc/bin/combined_alignment.py, then bin/build_snp_tree.sh)
takes its topology and branch lengths from every panel genome's SNPs as well as
the cohort's. The event reconstruction (write_event_matrix.py) reads genotypes
from the cohort's merged VCF, which the panel genomes are not in, so the tree it
gets holds exactly the cohort's isolates, the reference tip and the outgroup.

Pruning removes the other panel tips and joins the branches they leave, summing
lengths, so the remaining tips keep their patristic distances. The root stays on
the outgroup: the root must have two children, one of them the outgroup.

    $MTB_PY assoc/bin/prune_for_cohort.py --tree data/trees/<c>.combined.rooted.nwk \\
        --vcf refbias/<c>/p5/merged.vcf.gz --outgroup <outgroup> \\
        --out data/trees/<c>.rooted.nwk

Brought into the chain from analysis/combined_tree/ (review 2, R2-TREES-1). The
reference tip is the cohort alignment's own reference row (--ref-sample), so
no renaming is needed.
"""
import argparse, sys
import dendropy, pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--ref-sample", default="H37Rv")
    ap.add_argument("--outgroup", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    samples = set(pysam.VariantFile(a.vcf).header.samples)
    t = dendropy.Tree.get(path=a.tree, schema="newick", preserve_underscores=True,
                          rooting="force-rooted")
    tips = {n.taxon.label for n in t.leaf_node_iter()}
    missing = samples - tips
    if missing:
        sys.exit(f"FATAL: {len(missing)} VCF samples are not tips of {a.tree}: "
                 f"{sorted(missing)[:5]}")
    for need in (a.ref_sample, a.outgroup):
        if need not in tips:
            sys.exit(f"FATAL: {need} is not a tip of {a.tree}")
    keep = samples | {a.ref_sample, a.outgroup}
    t.prune_taxa([n.taxon for n in t.leaf_node_iter() if n.taxon.label not in keep])
    t.suppress_unifurcations()
    kids = t.seed_node.child_nodes()
    if len(kids) != 2 or not any(k.is_leaf() and k.taxon.label == a.outgroup
                                 for k in kids):
        sys.exit(f"FATAL: after pruning the root's children are not the "
                 f"outgroup {a.outgroup} and the ingroup")
    t.write(path=a.out, schema="newick", suppress_rooting=True,
            unquoted_underscores=True)
    print(f"  {len(tips)} tips -> {len(list(t.leaf_node_iter()))} "
          f"({len(samples)} isolates + {a.ref_sample} + {a.outgroup}); "
          f"root on {a.outgroup}")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    sys.exit(main())
