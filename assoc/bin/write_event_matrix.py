#!/usr/bin/env python3
"""Write a merged cohort VCF and a tree out as mtbvartools CallBytestreams.

WHY THIS EXISTS. phyoverlap2 asks a convergence question -- are the leaves
descending from a set of independent event branches enriched for a second
label -- and it reads those events from a `mtbvartools.CallBytestream`
directory laid out variant x node. Our reconstruction already exists but is
published as an `AA` INFO key inside the merged VCF, which records the state at
the ROOT and nothing about the branches. This script is the missing writer: it
reconstructs a state for every node, turns the state changes into per-branch
events, and writes both matrices in the format the library reads.

    <out>/ancestor/   states, per node   0 ancestral, 1 derived, 2 unknown
    <out>/event/      events, per branch 0 none, 1 gain, 2 undetermined, 3 loss
    <out>/labelled.nwk  the tree with unique internal labels -- USE THIS ONE
    <out>/variants.tsv  row key -> the pipeline's own key, class, polarity
    <out>/nodes.tsv     label, parent, edge length, discarded support value
    <out>/summary.txt

The two encodings are not arbitrary. `mtbvartools.trees.getVariantDistance`
counts positions where two nodes differ and neither is 2, so 2 must be the
unknown sentinel in the STATE stream. `phyoverlap2.variant`'s event lookup
selects nodes where the value equals 1, so a gain must be 1 and a loss must be
something else in the EVENT stream -- hence 3 for a loss, which keeps losses
recorded rather than silently merged into "no event".

POLARITY IS READ, NOT ASSUMED. A gain means a gain of the DERIVED allele, and
which allele that is comes from the VCF: `AA_INVERTED` marks the 8.9% of panel
sites where the reference carries the derived allele, so there GT=0 is the
derived state and GT=1 the ancestral one. Getting this backwards would invert
the direction of every event at those sites. Records with no `AA` -- every
is6110 and sv record, since those are not panel SNP sites -- are polarised
ALT-as-derived, which for an insertion is a definition rather than a guess:
presence is the derived state. For a CLASS=small record with a tied or absent
AA it IS an assumption, and summary.txt counts those separately.

THE RECONSTRUCTION is Fitch: an upward pass for state sets, then a downward
pass assigning each node the parent's state where that is consistent with its
own set, and its own set otherwise. That is one most-parsimonious assignment,
not a distribution over them; samarray's `AncestorArray` computes the fraction
of tied optimal reconstructions taking each state and is the better instrument
if the event set turns out to hinge on ambiguous branches. Every branch whose
two endpoints are not both resolved is written as 2, undetermined, rather than
being resolved by a coin toss -- the same rule the `AA` annotation already
follows for tied root states.

THE ROOT IS THE WEAK POINT, and `--root` is where the choice is made. Fitch
infers the root state from the cohort's own leaves, and where the cohort splits
evenly the root is genuinely tied -- every branch below an unresolved node is
then written undetermined, so the variant contributes no events at all. It is
tempting to pin the root to the panel-ancestral allele that `AA` already gives,
and `--root ancestral` does exactly that, but it is NOT the default and should
not be used casually: a cohort is a subtree of the panel, so a derived allele
that arose ABOVE the cohort's most recent common ancestor would be forced to
arise again inside it, manufacturing parallel gains that never happened. On the
synthetic test this turns one true loss into two spurious gains plus a loss.
The sound fix is not a flag: include an outgroup in the cohort tree, or graft
the cohort onto the rooted panel tree, so the root state is resolved by data.

WHY THE INTERNAL LABELS ARE REWRITTEN. Our trees come from IQ-TREE, whose
internal labels are support values like `100/100` and are therefore neither
unique nor node identifiers. phyoverlap2 keys every one of its dictionaries on
`node.label`, so those labels would collide and quietly merge unrelated
clades. This script assigns `n00001`-style labels, keeps the support value in
nodes.tsv, and writes the labelled tree back out; downstream work must load
THAT tree, or the labels in the bytestream will not resolve.

INTERPRETER. Needs `mtbvartools` (for the bytestream format) and `dendropy`
(so the tree is parsed exactly as the library will parse it), which the
pipeline environment does not carry. Use MTB_PY_VT from config/project_env.sh:

    $MTB_PY_VT assoc/bin/write_event_matrix.py \
        --vcf refbias/<cohort>/p5/merged.vcf.gz \
        --tree data/trees/<cohort>.rooted.nwk \
        --out assoc/<cohort>/events
"""
import argparse, collections, csv, glob, gzip, os, sys
import numpy as np

# State sets during the passes are two-bit masks, not the written encoding.
UNK, A_BIT, D_BIT, AMB = 0, 1, 2, 3        # ancestral bit, derived bit
ST_ANC, ST_DER, ST_UNK = 0, 1, 2           # written STATE encoding
EV_NONE, EV_GAIN, EV_UNDET, EV_LOSS = 0, 1, 2, 3   # written EVENT encoding


def eprint(*a):
    print(*a, file=sys.stdout, flush=True)


def read_fasta_seq(path, want):
    """One named sequence from a FASTA, without holding the rest."""
    seq, on = [], False
    for line in open(path):
        if line.startswith(">"):
            if on:
                break
            on = line[1:].split()[0] == want
        elif on:
            seq.append(line.strip())
    if not seq:
        sys.exit(f"FATAL: {want} is not in {path}")
    return "".join(seq)


def load_tree(path):
    import dendropy
    tree = dendropy.Tree.get(path=path, schema="newick",
                             preserve_underscores=True)
    for leaf in tree.leaf_nodes():          # what mtbvartools.loadTree does
        leaf.label = leaf.taxon.label
    return tree


def label_internal_nodes(tree, prefix="n"):
    """Give every internal node a unique label, keeping the old one as support.

    Returns {label: support_string}. Leaf labels are left alone: dendropy
    writes them from the taxon, and phyoverlap2 reads them from the taxon too.
    """
    support = {}
    i = 0
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        i += 1
        support[f"{prefix}{i:05d}"] = node.label or ""
        node.label = f"{prefix}{i:05d}"
    return support


def read_vcf(path, want_class=None, want_frame=None):
    """One pass over the merged VCF.

    Returns (samples, variants, gt_rows) where a variant is a dict and a
    gt_row is a uint8 array of the raw GT characters in sample order.
    """
    op = gzip.open if path.endswith(".gz") else open
    samples, variants, gt_rows = None, [], []
    n_skipped = collections.Counter()
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            info = {}
            for kv in f[7].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info[k] = v
                elif kv:
                    info[kv] = True
            cls = info.get("CLASS", "")
            frame = info.get("FRAME", "")
            if want_class and cls not in want_class:
                n_skipped[f"class={cls}"] += 1
                continue
            if want_frame and frame not in want_frame:
                n_skipped[f"frame={frame}"] += 1
                continue
            chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
            # THE ROW KEY, and why it is not simply (POS, REF, ALT). It must
            # be a 3-tuple so mtbvartools.maskVCB can unpack it, but POS alone
            # does not identify a record in this file, in two ways that both
            # bit on the first real cohort:
            #
            #   Node-frame records all sit at POS=1 on a contig of length 1,
            #   with the real offset inside the graph node carried in the ID
            #   (node:<node>:<offset>:REF>ALT). Twelve records on one node
            #   therefore share (1, ...) until the offset is used.
            #
            #   SV records use a symbolic ALT, so two different insertions at
            #   one position are both (pos, N, <INS>). SVLEN separates them.
            #
            # So: the locus is an int POS on the H37Rv contig, which keeps a
            # genome mask usable there, and the ID's node:<node>:<offset>
            # prefix otherwise -- which also makes an H37Rv mask fail loudly
            # on a record it cannot describe rather than mask the wrong base.
            # Verified unique: 85,970 records -> 85,970 keys on scale200.
            if vid.startswith("node:"):
                locus = ":".join(vid.split(":")[:3])
            else:
                locus = pos
            akey = alt
            svlen = info.get("SVLEN", "")
            if akey.startswith("<") and akey.endswith(">") and svlen:
                akey = f"{akey[:-1]}:{svlen}>"
            # A level-1 presence record sits at the locus's H37Rv anchor with a
            # symbolic <INS> and no SVLEN, so it would key as (pos, N, <INS>)
            # and collide with any caller insertion anchored at the same base.
            # The locus id separates them while POS stays an int, so a genome
            # mask remains usable at that coordinate.
            acc = info.get("ACCLOCUS", "")
            if acc and akey.startswith("<"):
                akey = f"{akey[:-1]}:acc:{acc}>"
            key = (locus, ref, akey)
            if "AA_INVERTED" in info:
                polarity, derived = "alt_ancestral", "REF"
            elif info.get("AA", ".") not in (".", ""):
                polarity, derived = "ref_ancestral", "ALT"
            elif cls in ("is6110", "sv", "accessory_presence"):
                # For accessory_presence the polarity is not an assumption: the
                # locus is sequence absent from H37Rv by construction, so
                # carrying it is the derived state by definition of the frame.
                polarity, derived = "presence_is_derived", "ALT"
            else:
                polarity, derived = "unpolarised", "ALT"
            variants.append(dict(
                key=key, chrom=chrom, pos=pos, id=vid, ref=ref, alt=alt,
                cls=cls, frame=frame, aa=info.get("AA", ""), panel_af="",
                aa_flag=info.get("AA_FLAG", ""), polarity=polarity,
                derived=derived, region=info.get("REGION", ""),
                evidence=info.get("EVIDENCE", ""),
                siteclass=info.get("SITECLASS", "")))
            gt_rows.append(np.frombuffer(
                "".join([s[0] for s in f[9:]]).encode(), dtype=np.uint8))
    return samples, variants, gt_rows, n_skipped


def leaf_bits(gt_chars, derived_is_ref, absent):
    """Raw GT characters -> two-bit state masks, in the derived-allele frame."""
    b = np.zeros(gt_chars.shape, dtype=np.uint8)
    is0 = gt_chars == ord("0")
    is1 = gt_chars == ord("1")
    is2 = gt_chars == ord("2")
    if derived_is_ref:
        b[is0] = D_BIT
        b[is1] = A_BIT
    else:
        b[is0] = A_BIT
        b[is1] = D_BIT
    if absent == "ref":            # GT=2 read as the REF allele's state
        b[is2] = D_BIT if derived_is_ref else A_BIT
    elif absent == "ancestral":
        b[is2] = A_BIT
    # absent == "unknown" leaves GT=2 at UNK, which is the default: a deleted
    # region is an absence of observation about the allele, not an observation
    # of the reference allele.
    return b


def fitch(post, pre, children, parent, leaf_index, bits, pin_root=None):
    """Vectorised Fitch over a block of variants.

    bits is (n_leaves, n_block); returns the resolved per-node state masks as
    (n_nodes, n_block).
    """
    n_nodes = len(post)
    n_block = bits.shape[1]
    up = np.zeros((n_nodes, n_block), dtype=np.uint8)
    for node in post:                      # upward pass
        kids = children[node]
        if not kids:
            up[node] = bits[leaf_index[node]]
            continue
        acc_and = np.full(n_block, AMB, dtype=np.uint8)
        acc_or = np.zeros(n_block, dtype=np.uint8)
        for ch in kids:
            b = up[ch]
            known = b != UNK
            acc_and = np.where(known, acc_and & b, acc_and)
            acc_or |= b
        # An empty union means no descendant had data, which is unknown --
        # not the all-ones mask acc_and would otherwise carry.
        up[node] = np.where(acc_or == UNK, UNK,
                            np.where(acc_and != UNK, acc_and, acc_or))
    down = np.zeros((n_nodes, n_block), dtype=np.uint8)
    for node in pre:                       # downward pass
        p = parent[node]
        u = up[node]
        if p is None:
            if pin_root is None:
                down[node] = u
            else:
                # Pin only where the caller could resolve polarity AND the
                # variant has some data; an unpolarised variant has no
                # panel-ancestral allele to pin to, and a variant with no data
                # anywhere should not be reported as resolved at the root.
                down[node] = np.where(pin_root & (u != UNK), A_BIT, u)
            continue
        cand = down[p] & u
        down[node] = np.where(u == UNK, UNK, np.where(cand != UNK, cand, u))
    return down


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="merged cohort VCF")
    ap.add_argument("--tree", required=True,
                    help="rooted newick; its leaves must be VCF samples")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--class", dest="cls", default="",
                    help="comma-separated CLASS values to keep (default all)")
    ap.add_argument("--frame", default="",
                    help="comma-separated FRAME values to keep (default all). "
                         "Pass h37rv if a downstream genome mask will be used")
    ap.add_argument("--absent", choices=("unknown", "ref", "ancestral"),
                    default="unknown",
                    help="how to read GT=2, a deleted region (default unknown)")
    ap.add_argument("--root", choices=("fitch", "ancestral"), default="fitch",
                    help="state at the tree root. `fitch` (the default) infers "
                         "it from this cohort's own data and leaves it tied "
                         "where parsimony cannot resolve it. `ancestral` pins "
                         "it to the panel-ancestral allele wherever AA "
                         "resolved -- READ THE WARNING in this script's "
                         "header before using it")
    ap.add_argument("--chunk", type=int, default=20000,
                    help="variants per Fitch block")
    ap.add_argument("--panel-polarity",
                    default="refbias/assets/panel_polarity.tsv",
                    help="ancestral allele read off the outgroup's own "
                         "genotype in the panel VCF, from bin/panel_polarity.py. "
                         "Consumed ONLY where AA did not resolve, because Fitch "
                         "over 333 genomes and a tree is the stronger "
                         "inference. It is what polarises indels at all: the AA "
                         "table comes from a SNP alignment, so every indel is "
                         "unpolarised without it and defaults to ALT-as-derived "
                         "-- backwards wherever the reference is the odd genome "
                         "out.")
    ap.add_argument("--min-panel-af", type=float, default=0.05,
                    help="the outgroup's allele is only taken as ANCESTRAL "
                         "against the reference when the panel also carries it "
                         "at least this often. Below it, calling the allele "
                         "ancestral implies it was lost in 95%+ of the panel, "
                         "which is far less parsimonious than the outgroup "
                         "having acquired it independently -- 216 of 288 "
                         "candidate inversions on scale200 are in that "
                         "position. Those stay unpolarised rather than being "
                         "flipped on one genome's say-so.")
    ap.add_argument("--no-collapse-ladders", action="store_true",
                    help="keep every rung of an insertion ladder as its own "
                         "variant. Off by default: N nested ALT alleles at one "
                         "position are one insertion read to N depths, and "
                         "treating them as independent manufactures convergence."
                         " See assoc/IDSB_ON_GWAS1000.md.")
    ap.add_argument("--dedupe", choices=("fail", "first", "drop", "suffix"),
                    default="fail",
                    help="what to do when two records share a row key. A "
                         "merged cohort VCF cannot, since its keys are unique "
                         "by construction, but a graph deconstruct VCF emits "
                         "one record per bubble and nested bubbles repeat a "
                         "site. `fail` (default) stops, `first` keeps the "
                         "first copy, `drop` removes every copy")
    ap.add_argument("--allow-extra-samples", action="store_true",
                    help="tolerate VCF samples that are not tree leaves")
    ap.add_argument("--accessory-presence", default="",
                    help="LEVEL 2. Directory of <sample>.presence.tsv from "
                         "accessory/bin/locus_presence.py. A node-frame variant "
                         "that lies inside an accessory locus is reconstructed "
                         "ONLY over the samples that carry that locus; for a "
                         "sample that does not carry the insert the variant is "
                         "INAPPLICABLE -- not reference, not unknown -- so its "
                         "leaf is made uninformative and the reconstruction is "
                         "confined to the carrier clades. Needs --node-locus.")
    ap.add_argument("--node-locus",
                    default="accessory/assets/node_locus.tsv",
                    help="node -> accessory locus, from "
                         "accessory/bin/node_locus_from_p4.py. This is p4's own "
                         "coordinate anchor, not sequence containment.")
    ap.add_argument("--polarise-by-root", action="store_true",
                    help="A SECOND PASS that fixes polarity where the "
                         "reconstruction contradicts the assumption. 70% of "
                         "variants have no resolved ancestral allele and are "
                         "polarised alternate-as-derived by default. That "
                         "assumption is usually right, because a rare minor "
                         "allele is almost always the derived one, but where it "
                         "is wrong it does specific damage: it moves the "
                         "convergence signal from the gain column to the loss "
                         "column, so a variant with repeated independent "
                         "ORIGINS reads as one with repeated losses and scores "
                         "zero gains, which makes it untestable rather than "
                         "wrong. Measured on gwas1000: 2,439 variants "
                         "reconstruct with the derived allele AT THE ROOT, and "
                         "1,255 have two or more losses and fewer than two "
                         "gains. This flips those and reconstructs them again.")
    ap.add_argument("--prune-tips", default="",
                    help="comma-separated tree leaves to remove before "
                         "anything else. For an outgroup whose only job was to "
                         "place the root and which has no genotype in the "
                         "cohort VCF; the rooting it established survives its "
                         "removal")
    ap.add_argument("--outgroup-fasta", default="",
                    help="the alignment the tree was built from, including the "
                         "outgroup sequence assoc/bin/add_outgroup.py "
                         "appended. Requires --outgroup-sites and "
                         "--outgroup-name. WITHOUT THIS the outgroup tip has "
                         "no data, and pruning it or leaving it unknown are "
                         "equally useless: the root state falls back to what "
                         "the cohort alone can say, which is what the outgroup "
                         "was added to avoid")
    ap.add_argument("--outgroup-sites", default="",
                    help="the sites table written beside that alignment, one "
                         "row per column in order")
    ap.add_argument("--outgroup-name", default="",
                    help="the outgroup's label, a tree leaf that is not a VCF "
                         "sample")
    ap.add_argument("--ref-sample", default="",
                    help="a tree leaf that is the VCF's REFERENCE rather than "
                         "one of its samples, and so carries no genotype "
                         "column. It is given the reference allele at every "
                         "record. The panel tree's H37Rv tip needs this, and "
                         "bin/vcf_to_alignment.py takes the same option for "
                         "the same reason")
    a = ap.parse_args()

    try:
        from mtbvartools.KeyedByteArray import KeyedByteArray
    except ImportError:
        sys.exit("FATAL: mtbvartools is not importable. Use $MTB_PY_VT "
                 "(see config/project_env.sh), not the pipeline python.")

    os.makedirs(a.out, exist_ok=True)
    want_class = {x for x in a.cls.split(",") if x}
    want_frame = {x for x in a.frame.split(",") if x}

    # ---- tree
    tree = load_tree(a.tree)
    drop = [x for x in a.prune_tips.split(",") if x]
    if drop:
        have = {l.taxon.label for l in tree.leaf_nodes()}
        absent = [d for d in drop if d not in have]
        if absent:
            sys.exit(f"FATAL: --prune-tips names leaves that are not in the "
                     f"tree: {absent}")
        tree.prune_taxa_with_labels(drop)
        eprint(f"  pruned {len(drop)} tip(s): {', '.join(drop)}")
    support = label_internal_nodes(tree)
    labelled = os.path.join(a.out, "labelled.nwk")
    tree.write(path=labelled, schema="newick",
               suppress_internal_node_labels=False, suppress_rooting=False,
               unquoted_underscores=True)
    # Round trip through the library's own loader: if the labels do not come
    # back, nothing downstream will resolve and it is better to fail here.
    check = load_tree(labelled)
    want = {n.label for n in tree.preorder_node_iter()}
    got = {n.label for n in check.preorder_node_iter()}
    if want != got:
        sys.exit(f"FATAL: {len(want - got)} labels did not survive the newick "
                 f"round trip, e.g. {sorted(want - got)[:5]}")

    nodes = list(tree.preorder_node_iter())
    order = {n.label: i for i, n in enumerate(nodes)}
    if len(order) != len(nodes):
        sys.exit("FATAL: node labels are not unique after relabelling")
    children = [[order[c.label] for c in n.child_nodes()] for n in nodes]
    parent = [None if n.parent_node is None else order[n.parent_node.label]
              for n in nodes]
    is_leaf = [n.is_leaf() for n in nodes]
    labels = [n.label for n in nodes]
    leaves = [labels[i] for i in range(len(nodes)) if is_leaf[i]]
    root_i = next(i for i in range(len(nodes)) if parent[i] is None)
    # post-order = reverse of pre-order for a tree, which is all the upward
    # pass needs (children before parents).
    post = list(range(len(nodes)))[::-1]
    pre = list(range(len(nodes)))
    eprint(f"  tree: {len(nodes)} nodes, {len(leaves)} leaves, "
           f"root {labels[root_i]}")

    # ---- VCF
    eprint(f"  reading {a.vcf} ...")
    samples, variants, gt_rows, skipped = read_vcf(
        a.vcf, want_class or None, want_frame or None)

    # ---- polarity from the outgroup, where AA left the record unpolarised
    if a.panel_polarity and os.path.exists(a.panel_polarity):
        pol_tab = {}
        for r in csv.DictReader(open(a.panel_polarity), delimiter="\t"):
            pol_tab[(int(r["pos"]), r["ref"].upper(), r["alt"].upper())] = (
                r["ancestral"], r["panel_af"])
        n_fix = collections.Counter()
        for v in variants:
            if v["polarity"] != "unpolarised":
                continue
            hit = pol_tab.get((v["pos"], v["ref"].upper(), v["alt"].upper()))
            if hit is None:
                continue
            anc, af = hit
            v["panel_af"] = af
            if anc == "ALT":
                try:
                    ok = float(af) >= a.min_panel_af
                except ValueError:
                    ok = False
                if not ok:
                    n_fix["ALT ancestral but too rare in the panel to trust"] += 1
                    continue
                v["polarity"], v["derived"] = "alt_ancestral_outgroup", "REF"
                n_fix["ALT is ancestral -- the polarity was backwards"] += 1
            else:
                v["polarity"], v["derived"] = "ref_ancestral_outgroup", "ALT"
                n_fix["REF is ancestral -- the polarity is confirmed"] += 1
            n_fix["  of which indels" if len(v["ref"]) != len(v["alt"])
                  else "  of which SNPs"] += 1
        if n_fix:
            eprint("  polarity from the outgroup, for records AA left "
                   "unresolved:")
            for k in sorted(n_fix):
                eprint(f"    {k:<48}{n_fix[k]:>8,}")
    if samples is None:
        sys.exit("FATAL: no #CHROM line in the VCF")
    # ---- COLLAPSE INSERTION LADDERS ---------------------------------------
    # One mobile-element insertion read to N different depths by the
    # small-variant caller becomes N VCF records whose ALT alleles are nested
    # prefixes of each other. Each then gets its own carrier set and scores its
    # own independent gains, which is manufactured convergence. idsB was the
    # case that exposed it: twelve alleles at 3,797,827, all beginning
    # CTGAACCGCCCCGGCATGTCCGGAGACTCC, whose reverse complement matches the
    # canonical element at offset 1,283 of 1,375 -- the IS6110 3' terminus. Its
    # carriers had the same phenotype rate as everyone else, 0.45 against 0.46.
    #
    # Measured on gwas1000: 68 positions carry five or more distinct insertion
    # alleles and 40 of them reach the scan, the largest ladder being 42 alleles
    # at one position.
    #
    # A sample carrying ANY rung carries the event, so the rungs are merged into
    # the longest allele and the genotypes unioned. That is the biologically
    # correct reading -- the insertion is present and the caller resolved
    # different amounts of it -- and it is the conservative one, because it
    # removes origins rather than adding them.
    if not a.no_collapse_ladders and variants:
        by_pos = collections.defaultdict(list)
        for i, v in enumerate(variants):
            if v["cls"] == "small" and len(v["alt"]) > len(v["ref"]) + 10:
                by_pos[(v["chrom"], v["pos"], v["ref"])].append(i)
        # NOT `drop`: that name already holds the pruned tip labels and is
        # read again at the end to write the summary. Shadowing it made the
        # run write every bytestream and then die on
        # `', '.join(drop)` with ints in it.
        rung_drop, merged, ladders = set(), 0, 0
        for idxs in by_pos.values():
            if len(idxs) < 2:
                continue
            # longest first; a shorter allele that is a prefix of a longer one
            # is the same insertion, read less far
            idxs.sort(key=lambda i: -len(variants[i]["alt"]))
            claimed = set()
            for a_i in idxs:
                if a_i in claimed:
                    continue
                grp = [a_i]
                for b_i in idxs:
                    if b_i == a_i or b_i in claimed:
                        continue
                    if variants[a_i]["alt"].startswith(variants[b_i]["alt"]):
                        grp.append(b_i); claimed.add(b_i)
                if len(grp) < 2:
                    continue
                ladders += 1
                keep = grp[0]
                row = gt_rows[keep].copy()
                for b_i in grp[1:]:
                    other = gt_rows[b_i]
                    # ALT anywhere wins; otherwise a call beats a no-call
                    row = np.where(other == ord("1"), other,
                                   np.where(row == ord("."), other, row))
                    rung_drop.add(b_i); merged += 1
                gt_rows[keep] = row
                variants[keep]["ladder_rungs"] = len(grp)
        if rung_drop:
            keep_i = [i for i in range(len(variants))
                      if i not in rung_drop]
            variants = [variants[i] for i in keep_i]
            gt_rows = [gt_rows[i] for i in keep_i]
            eprint(f"  collapsed {ladders:,} insertion ladders, merging "
                   f"{merged:,} rungs into their longest allele")

    eprint(f"  vcf: {len(samples)} samples, {len(variants)} records kept")
    for k, v in sorted(skipped.items()):
        eprint(f"    skipped {v:,} records with {k}")
    if not variants:
        sys.exit("FATAL: no records kept -- check --class and --frame")

    col = {s: i for i, s in enumerate(samples)}
    if a.ref_sample and a.ref_sample in col:
        eprint(f"  note: --ref-sample {a.ref_sample} IS a genotyped sample "
               f"here; its own column is used and the option ignored")
        a.ref_sample = ""
    if a.ref_sample and a.ref_sample not in set(leaves):
        sys.exit(f"FATAL: --ref-sample {a.ref_sample} is not a tree leaf")
    missing = [l for l in leaves if l not in col and l != a.ref_sample
               and l != a.outgroup_name]
    extra = [s for s in samples if s not in set(leaves)]
    if missing:
        sys.exit(f"FATAL: {len(missing)} tree leaves are not VCF samples, "
                 f"e.g. {missing[:5]}")
    if extra and not a.allow_extra_samples:
        sys.exit(f"FATAL: {len(extra)} VCF samples are not tree leaves, e.g. "
                 f"{extra[:5]}. Pass --allow-extra-samples to ignore them")
    if extra:
        eprint(f"  note: ignoring {len(extra)} VCF samples absent from the tree")
    # -1 marks the reference tip, which has no genotype column and is read as
    # the reference allele at every record; -2 marks the outgroup, whose
    # alleles come from the alignment rather than the VCF.
    take = np.asarray([-2 if (a.outgroup_name and l == a.outgroup_name)
                       else -1 if l == a.ref_sample else col[l]
                       for l in leaves])
    take_safe = np.maximum(take, 0)
    is_ref_tip = take == -1
    is_og_tip = take == -2
    if a.ref_sample:
        eprint(f"  note: tree leaf {a.ref_sample} is the VCF reference; "
               f"reading it as REF at every record")
    leaf_index = [None] * len(nodes)
    for j, i in enumerate([i for i in range(len(nodes)) if is_leaf[i]]):
        leaf_index[i] = j

    # ---- the outgroup's alleles, from the alignment the tree was built from
    og_gt = None
    if a.outgroup_name:
        if not (a.outgroup_fasta and a.outgroup_sites):
            sys.exit("FATAL: --outgroup-name needs --outgroup-fasta and "
                     "--outgroup-sites")
        og_seq = read_fasta_seq(a.outgroup_fasta, a.outgroup_name)
        og_sites = list(csv.DictReader(open(a.outgroup_sites), delimiter="\t"))
        if len(og_seq) != len(og_sites):
            sys.exit(f"FATAL: outgroup sequence is {len(og_seq)} long but the "
                     f"sites table has {len(og_sites)} rows")
        at = {(r["chrom"], int(r["pos"]), r["ref"].upper(), r["alt"].upper()): i
              for i, r in enumerate(og_sites)}
        og_gt = np.full(len(variants), ord("."), dtype=np.uint8)
        hit = collections.Counter()
        for j, v in enumerate(variants):
            i = at.get((v["chrom"], v["pos"], v["ref"].upper(),
                        v["alt"].upper()))
            if i is None:
                hit["site not in the alignment"] += 1
                continue
            b = og_seq[i].upper()
            if b == v["ref"].upper():
                og_gt[j] = ord("0"); hit["reference allele"] += 1
            elif b == v["alt"].upper():
                og_gt[j] = ord("1"); hit["alternate allele"] += 1
            else:
                hit["N or a third base"] += 1
        eprint(f"  outgroup {a.outgroup_name}:")
        for k, n in sorted(hit.items(), key=lambda x: -x[1]):
            eprint(f"    {k:<26s}{n:>8,}  {n / len(variants):6.2%}")

    keys = [v["key"] for v in variants]
    counts = collections.Counter(keys)
    dup = [k for k, c in counts.items() if c > 1]
    if dup and a.dedupe == "fail":
        sys.exit(f"FATAL: {len(dup)} duplicate row keys, e.g. {dup[:3]}. "
                 f"Keys must be unique -- the bytestream index is a dict. "
                 f"Pass --dedupe suffix to keep both, or first or drop")
    if dup and a.dedupe == "suffix":
        # KEEP BOTH RECORDS. The row key is (locus, ref, alt-with-SVLEN), which
        # is deliberately independent of the VCF ID so that the same event
        # reached in two frames or two cohorts maps to one key. The cost is
        # that two structural records at the same position, with the same REF
        # and the same type and length, collide however unique their IDs are.
        # gwas1000 has one such pair: two INS clusters at 802491, both of
        # length 108, carried by 35 and by 567 isolates.
        #
        # `first` would keep the 35-carrier row and silently discard the
        # 567-carrier one; `drop` would discard both. Neither is right when the
        # carrier sets differ by a factor of 16, so this appends an occurrence
        # number to the SECOND and later keys and keeps every record. The first
        # occurrence keeps its key unchanged, so nothing that already matches a
        # key elsewhere stops matching.
        n_seen = collections.Counter()
        renamed = []
        for v in variants:
            k = v["key"]
            n_seen[k] += 1
            if n_seen[k] > 1:
                v["key"] = k + (f"#{n_seen[k]}",) if isinstance(k, tuple) \
                    else f"{k}#{n_seen[k]}"
                renamed.append((k, v["key"], v.get("id", "")))
        eprint(f"  --dedupe suffix: {len(renamed):,} record(s) of "
               f"{len(variants):,} given an occurrence suffix; none removed")
        for k, nk, vid in renamed[:10]:
            eprint(f"    {k} -> {nk}   (VCF ID {vid})")
        keys = [v["key"] for v in variants]
        counts = collections.Counter(keys)
        dup = [k for k, c in counts.items() if c > 1]
    if dup:
        seen, keep = set(), []
        for i, k in enumerate(keys):
            if a.dedupe == "drop" and counts[k] > 1:
                continue
            if k in seen:
                continue
            seen.add(k)
            keep.append(i)
        eprint(f"  --dedupe {a.dedupe}: {len(dup):,} repeated keys, "
               f"{len(variants) - len(keep):,} of {len(variants):,} records "
               f"removed")
        variants = [variants[i] for i in keep]
        gt_rows = [gt_rows[i] for i in keep]
        keys = [v["key"] for v in variants]

    # ---- LEVEL 2: which samples is an inside-the-insert variant even about? -
    # The two levels, in Michael's framing: level 1 is whether the insert is
    # there, level 2 is the variation inside it among the samples that have it.
    # Conflating them is what put reference calls on samples that do not carry
    # the sequence. Here level 2 is made conditional.
    #
    # The mechanism is the Fitch pass itself and needs no separate tree. A leaf
    # left at UNK contributes nothing: the upward pass returns UNK for a node
    # whose descendants are all UNK, and the event pass only calls a gain where
    # parent and child are BOTH resolved. So marking every non-carrier leaf UNK
    # confines the reconstruction, and every event it can report, to the carrier
    # clades -- which is what "inapplicable" means operationally.
    #
    # What is NOT free is the bookkeeping. An inapplicable cell and a missing
    # cell both read as UNK downstream, so the count of each is recorded per
    # variant; a variant applicable to 40 samples is a different object from one
    # applicable to 997 and its event branches must be pooled against its own
    # kind by the region null, not against the whole tree.
    acc_carriers = {}
    node_locus = {}
    n_l2 = 0
    if a.accessory_presence:
        if not os.path.exists(a.node_locus):
            sys.exit(f"FATAL: --accessory-presence needs --node-locus; "
                     f"{a.node_locus} does not exist")
        with open(a.node_locus, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("locus"):
                    node_locus[r["node"]] = r["locus"]
        per = collections.defaultdict(set)
        pf = sorted(glob.glob(os.path.join(a.accessory_presence,
                                           "*.presence.tsv")))
        for f in pf:
            with open(f, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("state") == "PRESENT":
                        per[r["locus"]].add(r["sample"])
        acc_carriers = dict(per)
        eprint(f"  level 2: {len(node_locus):,} nodes placed in "
               f"{len({v for v in node_locus.values()}):,} accessory loci; "
               f"{len(pf):,} presence tables")
        # attach the locus to every variant that sits on a placed node
        for v in variants:
            nd = ""
            if v["id"].startswith("node:"):
                nd = v["id"].split(":")[1]
            lid = node_locus.get(nd, "")
            v["acc_locus"] = lid
            if lid:
                n_l2 += 1
        car_leaf = {}
        for lid, ss in acc_carriers.items():
            car_leaf[lid] = np.asarray([l in ss for l in leaves])
        eprint(f"  {n_l2:,} variants are inside a placed accessory locus and "
               f"become conditional")
    else:
        for v in variants:
            v["acc_locus"] = ""
        car_leaf = {}

    n_v, n_n = len(variants), len(nodes)
    eprint(f"  matrices: {n_n} x {n_v} uint8 x2 = "
           f"{2 * n_n * n_v / 1e9:.2f} GB in memory")
    states = np.empty((n_n, n_v), dtype=np.uint8)
    events = np.empty((n_n, n_v), dtype=np.uint8)

    # ---- reconstruct, in blocks
    for s in range(0, n_v, a.chunk):
        e = min(s + a.chunk, n_v)
        block = np.empty((len(leaves), e - s), dtype=np.uint8)
        for j in range(s, e):
            v = variants[j]
            g = gt_rows[j][take_safe]
            if is_ref_tip.any():
                g = np.where(is_ref_tip, ord("0"), g)
            if is_og_tip.any():
                g = np.where(is_og_tip, og_gt[j], g)
            b = leaf_bits(g, v["derived"] == "REF", a.absent)
            # LEVEL 2: a non-carrier has no sequence for this variant to be in,
            # so its leaf is made uninformative rather than called reference.
            lid = v.get("acc_locus") or ""
            if lid:
                keep_l = car_leaf.get(lid)
                if keep_l is None:
                    b[:] = UNK            # locus never level-1 genotyped
                    v["n_applicable"] = 0
                else:
                    v["n_inapplicable"] = int((~keep_l & (b != UNK)).sum())
                    b = np.where(keep_l, b, UNK)
                    v["n_applicable"] = int(keep_l.sum())
            block[:, j - s] = b
        pin = None
        if a.root == "ancestral":
            pin = np.asarray(
                [variants[j]["polarity"] != "unpolarised" for j in range(s, e)])
        down = fitch(post, pre, children, parent, leaf_index, block, pin)
        # states: written encoding
        st = np.full(down.shape, ST_UNK, dtype=np.uint8)
        st[down == A_BIT] = ST_ANC
        st[down == D_BIT] = ST_DER
        states[:, s:e] = st
        # events: compare each node with its parent
        ev = np.full(down.shape, EV_UNDET, dtype=np.uint8)
        for i in range(n_n):
            p = parent[i]
            if p is None:
                continue                    # no branch above the root
            pv, cv = down[p], down[i]
            det = ((pv == A_BIT) | (pv == D_BIT)) & \
                  ((cv == A_BIT) | (cv == D_BIT))
            row = np.full(pv.shape, EV_UNDET, dtype=np.uint8)
            row[det & (pv == cv)] = EV_NONE
            row[det & (pv == A_BIT) & (cv == D_BIT)] = EV_GAIN
            row[det & (pv == D_BIT) & (cv == A_BIT)] = EV_LOSS
            ev[i] = row
        events[:, s:e] = ev
        eprint(f"    reconstructed {e:,}/{n_v:,}")

    # ---- SECOND PASS: polarity the reconstruction disagrees with -----------
    # Only variants whose ancestral allele was ASSUMED are eligible. A variant
    # with a resolved ancestral allele from the panel or the outgroup keeps it;
    # external evidence outranks a parsimony inference, and silently overriding
    # it would make the `polarity` column a lie.
    #
    # The criterion is the reconstructed root. If the allele called derived sits
    # at the root, then on this tree it is the ancestral one and every event
    # counted for this variant runs backwards. Flipping is not circular in a way
    # that manufactures signal: the root state is inferred from the tip
    # distribution and the topology, both independent of the phenotype, and the
    # flip changes which of gain and loss a transition is called, not how many
    # transitions there are.
    n_flip = 0
    if a.polarise_by_root:
        # WHICH CLASSES ARE ELIGIBLE, and why it is not simply "all of them".
        # A variant reconstructing derived-at-root can mean three different
        # things and only one is a polarity error. Measured on gwas1000:
        #
        #   caller INS records      1,503 of 7,313 (21%). NOT polarity. That arm
        #                           has zero reference calls in 7.29M cells, so
        #                           Fitch sees only ALT tips and unknown tips and
        #                           the root can only come out ALT. Flipping
        #                           would dress a missing-genotype problem as a
        #                           resolved one.
        #   is6110 sites            30 of 2,250 (1.3%). An element inserts; it
        #                           does not un-insert and leave the site. So
        #                           presence really is derived there and this
        #                           rate is artefact-level. Not eligible.
        #   catalogued deletions    0 of 354 in the coherent tier, 39 of 854 in
        #                           the scattered-repeat tier. The assumption
        #                           holds where the evidence is good, so the few
        #                           exceptions are genotype noise. Not eligible.
        #   accessory presence      24 of 802 (3%). GENUINE polarity error. An
        #                           insert can be ancestral and lost, which is
        #                           the regions-of-difference pattern Behruznia
        #                           et al. report as dominating MTBC. TbD1 is
        #                           the type case: 415 carriers, root derived,
        #                           ZERO gains and 4 losses, so under
        #                           presence-is-derived it is untestable while
        #                           its 4 independent losses are the whole
        #                           signal. Eligible.
        _elig = {"unpolarised", "presence_is_derived"}
        cand = [j for j, v in enumerate(variants)
                if v["polarity"] in _elig
                and states[root_i, j] == ST_DER
                and (v["polarity"] == "unpolarised"
                     or v["cls"] == "accessory_presence")]
        _bycls = collections.Counter(
            f'{variants[j]["polarity"]}/{variants[j]["cls"]}' for j in cand)
        eprint(f"  polarise-by-root: {len(cand):,} eligible variants reconstruct "
               f"with the derived allele at the root: "
               + ", ".join(f"{k} {n:,}" for k, n in _bycls.most_common()))
        for s2 in range(0, len(cand), a.chunk):
            idx = cand[s2:s2 + a.chunk]
            block = np.empty((len(leaves), len(idx)), dtype=np.uint8)
            for c, j in enumerate(idx):
                v = variants[j]
                g = gt_rows[j][take_safe]
                if is_ref_tip.any():
                    g = np.where(is_ref_tip, ord("0"), g)
                if is_og_tip.any():
                    g = np.where(is_og_tip, og_gt[j], g)
                # the flip: derived becomes REF instead of ALT
                b = leaf_bits(g, True, a.absent)
                lid = v.get("acc_locus") or ""
                if lid:
                    keep_l = car_leaf.get(lid)
                    if keep_l is None:
                        b[:] = UNK
                    else:
                        b = np.where(keep_l, b, UNK)
                block[:, c] = b
            down = fitch(post, pre, children, parent, leaf_index, block, None)
            st = np.full(down.shape, ST_UNK, dtype=np.uint8)
            st[down == A_BIT] = ST_ANC
            st[down == D_BIT] = ST_DER
            ev = np.full(down.shape, EV_UNDET, dtype=np.uint8)
            for i in range(n_n):
                pp = parent[i]
                if pp is None:
                    continue
                pv, cv = down[pp], down[i]
                det = ((pv == A_BIT) | (pv == D_BIT)) & \
                      ((cv == A_BIT) | (cv == D_BIT))
                row = np.full(pv.shape, EV_UNDET, dtype=np.uint8)
                row[det & (pv == cv)] = EV_NONE
                row[det & (pv == A_BIT) & (cv == D_BIT)] = EV_GAIN
                row[det & (pv == D_BIT) & (cv == A_BIT)] = EV_LOSS
                ev[i] = row
            for c, j in enumerate(idx):
                states[:, j] = st[:, c]
                events[:, j] = ev[:, c]
                variants[j]["polarity"] = "root_inferred"
                variants[j]["derived"] = "REF"
                n_flip += 1
            eprint(f"    repolarised {min(s2 + a.chunk, len(cand)):,}"
                   f"/{len(cand):,}")
        if n_flip:
            br = np.asarray([i for i in range(n_n) if i != root_i])
            still = sum(1 for j in cand if states[root_i, j] == ST_DER)
            gains = sum(int((events[br, j] == EV_GAIN).sum()) for j in cand)
            eprint(f"  flipped {n_flip:,} variants to root_inferred; "
                   f"{gains:,} gains now counted in them, and {still:,} still "
                   f"reconstruct derived-at-root (a flip cannot fix those -- "
                   f"they are genuinely ambiguous)")

    # ---- write the two bytestreams, both orientations
    def write_pair(name, mat):
        d = os.path.join(a.out, name)
        os.makedirs(d, exist_ok=True)
        kba = KeyedByteArray(f"{d}/by_variant.kba", mode="w",
                             columns=labels, dtype="uint8")
        for j, k in enumerate(keys):
            kba.write(k, mat[:, j])
        kba.close()
        kba = KeyedByteArray(f"{d}/by_node.kba", mode="w",
                             columns=keys, dtype="uint8")
        for i, lab in enumerate(labels):
            kba.write(lab, mat[i, :])
        kba.close()
        eprint(f"  -> {d}/by_variant.kba  ({len(keys)} rows x {len(labels)} cols)")
        eprint(f"  -> {d}/by_node.kba     ({len(labels)} rows x {len(keys)} cols)")

    write_pair("ancestor", states)
    write_pair("event", events)

    # ---- side tables
    ROOT_STATE = {ST_ANC: "ancestral", ST_DER: "derived", ST_UNK: "unresolved"}
    with open(os.path.join(a.out, "variants.tsv"), "w") as fh:
        fh.write("row_key\tid\tchrom\tpos\tref\talt\tclass\tframe\tregion\t"
                 "aa\taa_flag\tpolarity\tderived\tpanel_af\tevidence\t"
                 "siteclass\tacc_locus\tn_applicable\tn_inapplicable\t"
                 "root_state\tn_gain\tn_loss\tn_undet\tn_derived_leaves\n")
        leaf_rows = np.asarray([i for i in range(n_n) if is_leaf[i]])
        branch_rows = np.asarray([i for i in range(n_n) if i != root_i])
        for j, v in enumerate(variants):
            # the root row carries no branch, so it is not an undetermined
            # event -- counting it inflated every variant's n_undet by one
            ev = events[branch_rows, j]
            fh.write("\t".join(str(x) for x in (
                "|".join(str(p) for p in v["key"]), v["id"], v["chrom"],
                v["pos"], v["ref"], v["alt"], v["cls"], v["frame"],
                v["region"], v["aa"], v["aa_flag"], v["polarity"],
                v["derived"], v.get("panel_af", ""), v["evidence"],
                v["siteclass"], v.get("acc_locus", ""),
                v.get("n_applicable", ""), v.get("n_inapplicable", ""),
                ROOT_STATE[int(states[root_i, j])],
                int((ev == EV_GAIN).sum()), int((ev == EV_LOSS).sum()),
                int((ev == EV_UNDET).sum()),
                int((states[leaf_rows, j] == ST_DER).sum()))) + "\n")
    with open(os.path.join(a.out, "nodes.tsv"), "w") as fh:
        fh.write("label\tparent\tedge_length\tis_leaf\tsupport\tn_leaves\n")
        n_leaves_under = [0] * n_n
        for i in post:
            n_leaves_under[i] = 1 if is_leaf[i] else sum(
                n_leaves_under[c] for c in children[i])
        for i, n in enumerate(nodes):
            el = n.edge.length if n.edge is not None and n.edge.length else 0
            fh.write(f"{labels[i]}\t"
                     f"{'' if parent[i] is None else labels[parent[i]]}\t"
                     f"{el}\t{int(is_leaf[i])}\t{support.get(labels[i], '')}\t"
                     f"{n_leaves_under[i]}\n")

    # ---- summary
    branch = np.ones(n_n, dtype=bool)
    branch[root_i] = False
    ev_b = events[branch, :]
    pol = collections.Counter(v["polarity"] for v in variants)
    cls = collections.Counter(v["cls"] for v in variants)
    gains = (ev_b == EV_GAIN).sum(axis=0)
    lines = []
    lines.append(f"cohort vcf        {a.vcf}")
    lines.append(f"tree              {a.tree}")
    lines.append(f"labelled tree     {labelled}")
    lines.append(f"absent (GT=2) as  {a.absent}")
    lines.append(f"root policy       {a.root}")
    if drop:
        lines.append(f"pruned tips       {', '.join(drop)}")
    if a.outgroup_name:
        lines.append(f"outgroup          {a.outgroup_name}, alleles from "
                     f"{a.outgroup_fasta}")
    lines.append(f"nodes             {n_n}  ({len(leaves)} leaves, "
                 f"{n_n - len(leaves)} internal)")
    lines.append(f"variants          {n_v}")
    for k, v in sorted(cls.items()):
        lines.append(f"  class {k:<10s}{v}")
    lines.append("polarity source")
    # The list is explicit so the order is stable, which means a NEW polarity
    # label silently vanishes from the tally -- `root_inferred` did exactly
    # that on its first run, leaving 551 variants unaccounted between this
    # table and the total. Anything not named above is appended rather than
    # dropped.
    _named = ("ref_ancestral", "alt_ancestral", "ref_ancestral_outgroup",
              "alt_ancestral_outgroup", "presence_is_derived",
              "root_inferred", "unpolarised")
    for k in _named:
        if pol[k]:
            lines.append(f"  {k:<22s}{pol[k]:>9,}")
    for k in sorted(set(pol) - set(_named)):
        if pol[k]:
            lines.append(f"  {k:<22s}{pol[k]:>9,}  (unlisted label)")
    _shown = sum(pol[k] for k in pol)
    if _shown != len(variants):
        lines.append(f"  WARNING tally {_shown:,} != {len(variants):,} variants")
    if pol["unpolarised"]:
        lines.append(f"  NOTE {pol['unpolarised']:,} records have no resolved AA "
                     f"and were polarised ALT-as-derived by assumption")
    if pol["root_inferred"]:
        lines.append(f"  NOTE {pol['root_inferred']:,} had that assumption "
                     f"CONTRADICTED by the reconstruction -- the assumed-derived "
                     f"allele sat at the root -- and were flipped and rebuilt")
    rs = collections.Counter(int(x) for x in states[root_i, :])
    lines.append("state at the tree root")
    for code, name in ((ST_ANC, "ancestral"), (ST_DER, "derived"),
                       (ST_UNK, "unresolved")):
        if rs[code]:
            lines.append(f"  {name:<22s}{rs[code]:>9,}")
    if rs[ST_UNK]:
        lines.append(f"  NOTE every branch below an unresolved node is written "
                     f"undetermined, so those {rs[ST_UNK]:,} variants "
                     f"contribute no events")
    lines.append("branch events (root branch excluded)")
    for name, code in (("none", EV_NONE), ("gain", EV_GAIN),
                       ("loss", EV_LOSS), ("undetermined", EV_UNDET)):
        n = int((ev_b == code).sum())
        lines.append(f"  {name:<22s}{n:>12,}  "
                     f"{n / ev_b.size:.4%} of {ev_b.size:,} branch x variant")
    lines.append("independent gains per variant")
    lines.append(f"  variants with 0 gains {int((gains == 0).sum()):>9,}"
                 f"   (untestable: no independent origin)")
    lines.append(f"  variants with 1 gain  {int((gains == 1).sum()):>9,}"
                 f"   (single origin: no convergence to test)")
    lines.append(f"  variants with >=2     {int((gains >= 2).sum()):>9,}"
                 f"   (testable by phyoverlap2)")
    if (gains >= 2).any():
        lines.append(f"  max gains in one variant {int(gains.max()):>6,}")
    txt = "\n".join(lines) + "\n"
    open(os.path.join(a.out, "summary.txt"), "w").write(txt)
    eprint()
    eprint(txt)
    eprint(f"  -> {a.out}/summary.txt")


if __name__ == "__main__":
    main()
