#!/usr/bin/env python3
"""Ancestral allele of every panel variant, by Fitch parsimony over the panel
tree, for every variant class. Writes the build's two polarity tables.

THE DECISION (user, 2026-10-13: "yes, use parsimony over all outgroups for
every variant class"). One reconstruction over the panel tree, from each panel
genome's graph-VCF genotype, for SNPs, indels, SVs, IS6110 and accessory
records alike. It replaces two scripts, and the two ways they could disagree:

  - bin/ancestral_alleles.py: Fitch over a SNP ALIGNMENT, so SNPs only, and a
    tie at the ancestral node broken by the NEAREST outgroup that decides it;
  - bin/panel_polarity.py: one outgroup genome's own genotype, read off the
    graph VCF, for every class, used wherever the first left a gap.

On a reduced panel both rules fail at the same kind of site. The tatC-helY
element (68 bp, H37Rv 2,352,065): lineage 5.1 GCF_022870265 carries it and
lineage 6.1 GCF_000253355 does not, and the two join the lineage 1-4,7 clade
at the same node. Nearest-first then lets the order of the outgroup list
decide the ancestral state, and a single outgroup read gives whichever answer
that genome has. Parsimony over the whole tree uses every outgroup at once
(analysis/ancestry_need/RESULTS.md, the reference implementation:
analysis/ancestry_need/ancestry_need.py).

THE METHOD, step by step.
  1. Keys. Every biallelic record of the build's collapsed graph VCF, keyed
     (chrom, pos, ref, alt) after left-alignment against the build's H37Rv
     with bin/mtb_norm.normalise, the function the cohort VCF's keys go
     through (review 2, R2-INT-1). Copies of one normalised key are pooled.
  2. Each genome's state at each key: 1 (carries the ALT) if its GT is 1 in
     any copy, 0 if it is 0 in some copy and 1 in none, otherwise no call.
     At a SNP key, a genome at 0 that carries another allele changing the
     same base (a sibling record at that position) is no call, not REF
     (audit TP-7, R2-BUILD-3; vcf_to_alignment.py's rule).
       - H37Rv, the VCF's reference path, has no sample column: it is 0 at
         every key on its own contig (--ref-sample).
       - Outgroups outside the graph ("the outgroup kit": lineage 8
         GCF_012923765, canettii GCF_035581225) have rows only for SNPs, from
         their direct-alignment states (analysis/panel_checks/bin/
         genome_states.py: one uint8 per H37Rv base, 0 = uncalled, 1-4 =
         ACGT). Their base is 0 if it is the REF, 1 if the ALT, and no call if
         uncalled or a third base. At indel, MNP and SV keys they take route
         A's event states when --kit-events gives them (user, 2026-10-14,
         "option A": bin/outgroup_kit.py, the assembly aligned straight to
         H37Rv; OUT-1's pseudo-isolates may replace them), by the rules of
         kit_event_states(): the same normalised allele -> 1; the span in
         one callable alignment with no kit variant there -> 0; SVs matched
         by type, start within 50 bp and net length within 20%; anything
         ambiguous -> no call. Without --kit-events they are no call there.
         A kit genome that is also a VCF sample takes the kit's states (for a
         check run on a panel that still contains it).
  3. Fitch's upward pass over the rooted tree, one bitmask per key (1 = REF,
     2 = ALT, 3 = no call, i.e. either): a node's set is the intersection of
     its children's sets if that is not empty, else their union.
  4. The ancestral node is the MRCA of every leaf not named in --outgroups
     (for the P200 panel: the MRCA of the lineage 1-4 and 7 genomes). If that
     MRCA has an outgroup beneath it, the run stops (audit fault B: the
     canettii-plus-MTBC node was reported as the MTBC ancestor).
  5. Fitch's downward pass, from the root to the ancestral node, over the
     whole tree: a node's final set is its parent's final set where that
     meets the node's upward set, else its upward set. All outgroups act
     together and their list order is never read; there is no nearest-first
     step. Where the ancestral node's final set is one allele, that is the
     ancestral allele. Where it is both, the key is TIED and stays
     unpolarised: no coin toss.
  6. Output. --out-ancestral: ancestral.tsv, one row per SNP key on the H37Rv
     contig with ACGT alleles (the merged VCF's AA tag; columns chrom pos ref
     alt AA flag polarity, flag TIED:<bases> or NODATA as before). It now
     covers every SNP key of the graph, not only the tree alignment's columns.
     --out-polarity: panel_polarity.tsv, one row per key of every class whose
     ancestral allele was resolved (the event writer's --panel-polarity;
     columns chrom pos ref alt outgroup_gt ancestral panel_af n_records
     method). `outgroup_gt` is kept for the reader and now holds the
     ancestral node's state, 0 or 1. `panel_af` is unchanged: carriers over
     genomes called, over the VCF's sample columns, each genome once.
     `method` is `ingroup` where the ancestral node's upward set alone
     decided it, `outgroups` where only the downward pass did. TIED and
     NODATA keys are left out, so they stay unpolarised in the event writer.
  7. How far to trust it (fix list AA-1, user 2026-10-08). Both tables also
     carry, per key:
       changes  Fitch's change count (the fewest gains plus losses that
                explain the leaves) over the panel's GRAPH genomes: H37Rv,
                every VCF sample, the in-graph outgroups. Kit genomes are left
                out of the count (set to "either", which is the same as
                pruning them): their rows are direct-alignment calls, not
                graph genotypes, and the canettii branch's changes are not
                MTBC homoplasy. They still take part in the reconstruction.
       flag     a comma-separated list (the merged VCF's AA_FLAG):
                TIED:<bases> (ancestral.tsv; TIED is implied in
                panel_polarity.tsv, which leaves tied keys out) or NODATA as
                before; OUTGROUP where method is `outgroups`, i.e. the
                ingroup's own upward set was a tie and the outgroups decided
                it; HOMOPLASTIC where changes >= --homoplastic-min (default 2,
                confirmed by the user 2026-10-14). HOMOPLASTIC combines with
                any of the others.
     ancestral.tsv also gains `method` (ingroup, outgroups, tied, nodata).
     New columns are appended, so readers by column name are unaffected.
     The probabilistic reconstruction (AA-1 point 4) is not done.

WHAT IS NOT CHANGED. The event writer still consults panel_polarity.tsv only
where AA did not resolve (they agree at SNPs now, being one reconstruction),
still refuses an ALT-ancestral call below --min-panel-af (5%), and accessory
sequence absent from H37Rv is still presence-is-derived there.

DROPPED, deliberately. ancestral_alleles.py's nearest-first tie break. And
the alignment's site filter (--max-missing 0.10 in vcf_to_alignment.py): the
tree is still built from that alignment, but the states here are the VCF's
own, so a SNP key with many missing genomes gets an AA where the data allow.

    $MTB_PY bin/panel_ancestral.py --tree <panel>.rooted.nwk \\
        --panel-vcf <build>/assets/graph_collapsed.vcf.gz \\
        --h37rv <build>/refs/GCF_000195955.fasta --ref-sample GCF_000195955 \\
        --outgroups GCF_022870265,GCF_000253355,GCF_012923765,GCF_035581225 \\
        --kit GCF_012923765=<states>/GCF_012923765.npz \\
        --kit GCF_035581225=<states>/GCF_035581225.npz \\
        --kit-events GCF_012923765=analysis/outgroup_kit/out/GCF_012923765 \\
        --kit-events GCF_035581225=analysis/outgroup_kit/out/GCF_035581225 \\
        --out-ancestral <build>/assets/ancestral.tsv \\
        --out-polarity <build>/assets/panel_polarity.tsv

P0 step `ancestral` runs it (docs/PANEL_TREE.md section 1).
"""
import argparse
import collections
import csv
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mtb_norm  # noqa: E402

MISSING = {".", "./.", ".|."}
REF_BIT, ALT_BIT, ANY = 1, 2, 3          # Fitch state sets, per key
ACGT = "ACGT"


# --------------------------------------------------------------------------
# tree
# --------------------------------------------------------------------------
def parse_newick(s):
    """(children dict, root, leaves). Topology only: branch lengths and
    internal labels (IQ-TREE's support values) are skipped."""
    s = s.strip().rstrip(";")
    children = collections.defaultdict(list)
    stack, token, root, n = [], "", None, 0
    i = 0
    while i < len(s):
        c = s[i]
        if c == "(":
            n += 1
            node = f"_n{n}"
            if stack:
                children[stack[-1]].append(node)
            else:
                if root is not None:
                    sys.exit("FATAL: the tree has more than one top-level clade")
                root = node
            stack.append(node)
            token = ""
        elif c in ",)":
            name = token.split(":")[0].strip().strip("'\"")
            if name:
                if not stack:
                    sys.exit("FATAL: malformed newick (a leaf outside any clade)")
                children[stack[-1]].append(name)
            token = ""
            if c == ")":
                if not stack:
                    sys.exit("FATAL: malformed newick (unbalanced ')')")
                stack.pop()
                i += 1
                while i < len(s) and s[i] not in "(),":   # label / length
                    i += 1
                continue
        else:
            token += c
        i += 1
    if stack or root is None:
        sys.exit("FATAL: malformed newick (unbalanced '(' or no clade)")
    children = dict(children)
    leaves = {x for v in children.values() for x in v if x not in children}
    return children, root, leaves


def postorder(children, root):
    out, stack = [], [(root, False)]
    while stack:
        node, done = stack.pop()
        if done:
            out.append(node)
            continue
        stack.append((node, True))
        for ch in children.get(node, []):
            stack.append((ch, False))
    return out


def ancestral_node(children, root, leaves, outgroups):
    """The MRCA of every leaf not in `outgroups`, and the path root -> it.
    Exits if the MRCA has an outgroup beneath it (audit fault B)."""
    parent = {c: p for p, kids in children.items() for c in kids}
    ingroup = sorted(leaves - set(outgroups))
    if not ingroup:
        sys.exit("FATAL: every leaf is an outgroup")

    def path_up(x):
        out = [x]
        while x in parent:
            x = parent[x]
            out.append(x)
        return out
    common = path_up(ingroup[0])
    for leaf in ingroup[1:]:
        p = set(path_up(leaf))
        common = [x for x in common if x in p]
    node = common[0]
    below = set()
    stack = [node]
    while stack:
        x = stack.pop()
        if x in children:
            stack.extend(children[x])
        else:
            below.add(x)
    inside = sorted(below & set(outgroups))
    if inside:
        sys.exit(f"FATAL: the MRCA of the {len(ingroup)} non-outgroup leaves "
                 f"also contains the outgroup(s) {inside}, so its state would "
                 f"not be the ingroup ancestor's (audit fault B). Is the tree "
                 f"rooted inside the ingroup, or is an outgroup missing from "
                 f"--outgroups?")
    return node, list(reversed(path_up(node))), len(ingroup)


def fitch_at(children, root, path, leaf_bits, nkey):
    """Upward pass over the whole tree, then the downward pass along `path`
    (root ... ancestral node). Returns (upward set, final set) at the last
    node of the path, as uint8 arrays over keys."""
    up = {}
    for node in postorder(children, root):
        kids = children.get(node)
        if not kids:
            up[node] = leaf_bits.get(node, np.full(nkey, ANY, np.uint8))
            continue
        inter = np.bitwise_and.reduce([up[k] for k in kids])
        uni = np.bitwise_or.reduce([up[k] for k in kids])
        up[node] = np.where(inter != 0, inter, uni).astype(np.uint8)
    fin = up[path[0]]
    for node in path[1:]:
        m = fin & up[node]
        fin = np.where(m != 0, m, up[node]).astype(np.uint8)
    return up[path[-1]], fin


def fitch_changes(children, root, leaf_bits, nkey, skip=()):
    """Fitch's change count per key: the fewest state changes on the tree that
    explain the leaves (one per node whose children's sets do not meet).
    Leaves in `skip`, and leaves with no data, are the set {REF, ALT}, which
    never meets a child set emptily, so they add nothing: the count is the
    count over the remaining leaves, as if the skipped ones were pruned."""
    up, cost = {}, np.zeros(nkey, dtype=np.int32)
    for node in postorder(children, root):
        kids = children.get(node)
        if not kids:
            up[node] = (np.full(nkey, ANY, np.uint8) if node in skip else
                        leaf_bits.get(node, np.full(nkey, ANY, np.uint8)))
            continue
        cur = up[kids[0]]
        for k in kids[1:]:          # a multifurcation, resolved pairwise
            inter = cur & up[k]
            cost += inter == 0
            cur = np.where(inter != 0, inter, cur | up[k]).astype(np.uint8)
        up[node] = cur
    return cost


def aa_flags(how, changes, homoplastic_min, tied_bases=""):
    """AA-1's flag list for one key, comma-joined (a VCF Number=. list):
      TIED:<bases> / TIED  the ancestral node's final set is both alleles
      NODATA               no genome had a call
      OUTGROUP             resolved, but only by the downward pass: the
                           ingroup's own upward set was a tie and the
                           outgroups decided it (method `outgroups`)
      HOMOPLASTIC          the panel needs `homoplastic_min` or more changes
    HOMOPLASTIC combines with any of the others."""
    f = []
    if how == "tied":
        f.append("TIED:" + tied_bases if tied_bases else "TIED")
    elif how == "nodata":
        f.append("NODATA")
    elif how == "outgroups":
        f.append("OUTGROUP")
    if changes >= homoplastic_min:
        f.append("HOMOPLASTIC")
    return ",".join(f)


# --------------------------------------------------------------------------
# panel genotypes
# --------------------------------------------------------------------------
def read_panel(vcf, bcftools, seq, h37rv_contig):
    """Pooled normalised keys of the panel VCF.

    Returns (samples, keys, copies, stats): keys is a list of (chrom, pos,
    ref, alt) in VCF order of first appearance, copies[i] the list of GT
    lists of every record normalising to keys[i]."""
    r = subprocess.run([bcftools, "query", "-l", vcf], capture_output=True,
                       text=True)
    if r.returncode != 0:
        sys.exit(f"FATAL: bcftools query -l {vcf}: {r.stderr.strip()}")
    samples = r.stdout.split()
    proc = subprocess.Popen(
        [bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", vcf],
        stdout=subprocess.PIPE, text=True)
    pooled = collections.OrderedDict()
    st = collections.Counter()
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        chrom, pos, ref, alt, gts = f[0], int(f[1]), f[2].upper(), f[3].upper(), f[4:]
        st["records"] += 1
        if "," in alt:
            st["multiallelic"] += 1
            continue
        if len(gts) != len(samples):
            sys.exit(f"FATAL: {vcf} {chrom}:{pos} has {len(gts)} genotypes "
                     f"for {len(samples)} samples")
        if chrom.split("#")[-1] == h37rv_contig:
            k = mtb_norm.normalise(seq, pos, ref, alt)
            st["moved"] += k != (pos, ref, alt)
            pos, ref, alt = k
        else:
            st["offpath"] += 1
        pooled.setdefault((chrom, pos, ref, alt), []).append(gts)
    if proc.wait() != 0:
        sys.exit(f"FATAL: bcftools query failed on {vcf}")
    return samples, list(pooled), list(pooled.values()), st


def pooled_column(copies):
    """Per sample over every copy of a key: '1' carrier, '0' called, '.'."""
    if len(copies) == 1:
        return copies[0]
    return ["1" if "1" in c else "0" if c - MISSING else "."
            for c in map(set, zip(*copies))]


def kind_of(ref, alt):
    if alt.startswith("<"):
        return "SV"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    d = abs(len(ref) - len(alt))
    return "SV" if d >= 50 else "indel" if d else "MNP"


# --------------------------------------------------------------------------
# kit events: outgroup kit route A (bin/outgroup_kit.py)
# --------------------------------------------------------------------------
SV_POS_TOL = 50       # bp between starts (compare_h37rv.py, score.py)
SV_LEN_FRAC = 0.8     # smaller / larger net length (compare_h37rv.py)
SV_GUARD_MIN = 10     # a kit indel this long near an SV key makes it missing


def key_span(pos, ref, alt):
    """H37Rv bases a normalised key touches, its anchor and one flanking base
    included: an insertion sits between pos and pos+1."""
    return pos, pos + max(len(ref), 1)


def read_kit_events(prefix, seq, contig):
    """Route A outputs of one kit genome: <prefix>.var (paftools call: R =
    H37Rv stretches covered by exactly one alignment, V = SNPs and indels,
    0-based half-open) and <prefix>.sv.tsv (outgroup_kit.py's >= 50 bp rows,
    1-based). Returns (callable, vkeys, small_mask, sv_mask, svs):
      callable    int32 per base (1-based index): R interval id, 0 = not
      vkeys       the set of the kit's V lines as normalised (pos, ref, alt)
      small_mask  bool per base: touched by any kit V line or SV row
      sv_mask     bool per base: within SV_POS_TOL of a kit SV row or of a
                  kit V indel of SV_GUARD_MIN bp or more
      svs         candidate SVs: (type, start, net_len, key or None); key is
                  the normalised V line for an indel of >= 50 bp (it has its
                  sequence), None for a between-alignment row (no sequence)"""
    n = len(seq)
    callable_ = np.zeros(n + 2, dtype=np.int32)
    small = np.zeros(n + 2, dtype=bool)
    svm = np.zeros(n + 2, dtype=bool)
    vkeys, svs = set(), []
    nR = 0

    def mark(arr, a, b, pad=0):
        arr[max(1, a - pad):min(n, b + pad) + 1] = True

    for line in open(prefix + ".var"):
        c = line.rstrip("\n").split("\t")
        if c[0] not in ("R", "V") or c[1].split("#")[-1] != contig:
            continue
        s, e = int(c[2]), int(c[3])
        if c[0] == "R":
            nR += 1
            if callable_[s + 1:e + 1].any():
                sys.exit(f"FATAL: {prefix}.var: R intervals overlap at {s + 1}")
            callable_[s + 1:e + 1] = nR
            continue
        r = "" if c[6] == "-" else c[6].upper()
        x = "" if c[7] == "-" else c[7].upper()
        if r and x:                              # SNP (paftools: one base)
            pos, ref, alt = s + 1, r, x
        else:                                    # indel, anchored before it
            if s < 1:
                continue
            anchor = seq[s - 1]
            pos, ref, alt = s, anchor + r, anchor + x
        k = mtb_norm.normalise(seq, pos, ref, alt)
        vkeys.add(k)
        a_, b_ = key_span(*k)
        mark(small, a_, b_)
        d = len(k[2]) - len(k[1])
        if abs(d) >= SV_GUARD_MIN:
            mark(svm, a_, b_, SV_POS_TOL)
        if abs(d) >= 50:
            svs.append(("INS" if d > 0 else "DEL", k[0] + 1, abs(d), k))
    for row in csv.DictReader(open(prefix + ".sv.tsv"), delimiter="\t"):
        st, en = int(row["h37rv_start"]), int(row["h37rv_end"])
        mark(small, st - 1, en + 1)
        mark(svm, st - 1, en + 1, SV_POS_TOL)
        if row["source"] == "cs":
            continue                             # the V line, with sequence
        rl, al = int(row["ref_len"]), int(row["alt_len"])
        if row["type"] == "DEL" and rl - al >= 50:
            svs.append(("DEL", st, rl - al, None))
        elif row["type"] == "INS" and al - rl >= 50:
            svs.append(("INS", st, al - rl, None))
    if not nR:
        sys.exit(f"FATAL: {prefix}.var has no R (callable) lines on {contig}")
    return callable_, vkeys, small, svm, svs


def kit_event_states(keys, kinds, seq, contig, ev):
    """The kit genome's state (1, 0 or -1 no call) at every non-SNP key, by
    route A. Rules (any ambiguity is no call):
      indel / MNP (< 50 bp): 1 if the kit has the same normalised allele; 0
        if the key's span (anchor and flank included) lies in ONE callable
        R interval and no kit variant touches it; else no call.
      SV (>= 50 bp): each kit SV is matched to the graph SV keys of its type
        whose start is within SV_POS_TOL and whose net length is within
        SV_LEN_FRAC. If one candidate is the kit's exact allele, it gets 1
        and the others no call; else a single candidate gets 1 and several
        get no call. An unmatched SV key is 0 if its span lies in one
        callable interval and no kit SV or indel of >= SV_GUARD_MIN bp lies
        within SV_POS_TOL of it; else no call."""
    callable_, vkeys, small, svm, svs = ev
    cs_small = np.concatenate([[0], np.cumsum(small)])
    cs_sv = np.concatenate([[0], np.cumsum(svm)])

    def any_in(cs, a, b):
        return cs[b + 1] - cs[a] > 0

    def one_interval(a, b):
        return callable_[a] != 0 and callable_[a] == callable_[b]

    n = len(seq)
    out = np.full(len(keys), -1, dtype=np.int8)
    sv_ix = collections.defaultdict(list)          # type -> [(start, len, j)]
    for j, (chrom, pos, ref, alt) in enumerate(keys):
        if kinds[j] == "SNP" or chrom.split("#")[-1] != contig \
                or alt.startswith("<") or "*" in alt:
            continue
        a, b = key_span(pos, ref, alt)
        a, b = max(1, a - 1), min(n, b + 1)
        if kinds[j] != "SV":
            if (pos, ref, alt) in vkeys:
                out[j] = 1
            elif one_interval(a, b) and not any_in(cs_small, a, b):
                out[j] = 0
            continue
        d = len(alt) - len(ref)
        sv_ix["INS" if d > 0 else "DEL"].append((pos + 1, abs(d), j))
        if one_interval(a, b) and not any_in(cs_sv, a, b):
            out[j] = 0
    for typ, start, ln, k in svs:
        cand = [j for s, l, j in sv_ix[typ]
                if abs(s - start) <= SV_POS_TOL and min(l, ln) / max(l, ln) >= SV_LEN_FRAC]
        exact = [j for j in cand if k is not None and keys[j][1:] == k]
        if exact:
            for j in cand:
                if out[j] != 1:
                    out[j] = 1 if j in exact else -1
        elif len(cand) == 1:
            out[cand[0]] = 1
        else:
            for j in cand:
                if out[j] != 1:
                    out[j] = -1
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # no defaults: they were CX333's tree, VCF and canettii; P0 step
    # ancestral passes the build's
    ap.add_argument("--tree", required=True,
                    help="rooted newick; leaves = the VCF's samples, the "
                         "reference (--ref-sample) and the kit genomes")
    ap.add_argument("--panel-vcf", required=True,
                    help="<build>/assets/graph_collapsed.vcf.gz")
    ap.add_argument("--h37rv", required=True,
                    help="H37Rv FASTA (<build>/refs/<acc>.fasta); keys are "
                         "left-aligned against it as the cohort VCF's are")
    ap.add_argument("--h37rv-contig", default="NC_000962.3",
                    help="H37Rv contig name; panel contigs are matched on the "
                         "part after the last '#' (PanSN)")
    ap.add_argument("--ref-sample", required=True,
                    help="the tree leaf that is the VCF's reference path "
                         "(H37Rv's accession): REF at every key on its contig")
    ap.add_argument("--outgroups", required=True,
                    help="comma-separated outgroup leaves, in-graph and kit. "
                         "The ancestral node is the MRCA of every other leaf; "
                         "their order is not used")
    ap.add_argument("--kit", action="append", default=[],
                    metavar="ACC=STATES.npz",
                    help="an outgroup outside the graph and its direct-"
                         "alignment states (genome_states.py's npz, array "
                         "'state', one uint8 per H37Rv base). SNP keys only. "
                         "Repeatable, or comma-separated")
    ap.add_argument("--kit-events", action="append", default=[],
                    metavar="ACC=PREFIX",
                    help="route A event states of a kit genome "
                         "(bin/outgroup_kit.py): PREFIX.var (paftools call, "
                         "-L 10000 -l 1000) and PREFIX.sv.tsv. Gives it states "
                         "at indel, MNP and SV keys (kit_event_states). The "
                         "genome must also be a --kit. Repeatable, or "
                         "comma-separated")
    ap.add_argument("--accessions", default="",
                    help="the build's accessions.txt; if given, the tree's "
                         "leaves other than kit genomes must be exactly these")
    ap.add_argument("--bcftools",
                    default=os.environ.get("MTB_BCFTOOLS", "bcftools"))
    ap.add_argument("--out-ancestral", required=True,
                    help="<build>/assets/ancestral.tsv (SNP AA table)")
    ap.add_argument("--out-polarity", required=True,
                    help="<build>/assets/panel_polarity.tsv (all classes)")
    # AA-1. The user confirmed 2 (2026-10-14: "use 2 changes as the
    # homoplastic threshold"); a parameter so a check run can vary it.
    ap.add_argument("--homoplastic-min", type=int, default=2,
                    help="flag a key HOMOPLASTIC when its Fitch change count "
                         "over the panel's graph genomes is at least this "
                         "(default 2: the allele arose, or was lost, more "
                         "than once on the panel tree)")
    a = ap.parse_args()
    if a.homoplastic_min < 1:
        sys.exit(f"FATAL: --homoplastic-min must be >= 1, not "
                 f"{a.homoplastic_min} (0 would flag every key)")

    # ---- tree, outgroups, ancestral node
    children, root, leaves = parse_newick(open(a.tree).read())
    outgroups = [x for x in a.outgroups.split(",") if x]
    if len(set(outgroups)) != len(outgroups):
        sys.exit(f"FATAL: --outgroups repeats a name: {outgroups}")
    absent = [x for x in outgroups if x not in leaves]
    if not outgroups or absent:
        sys.exit(f"FATAL: outgroups must be leaves of {a.tree}; not found: "
                 f"{absent or '(none given)'}")
    anc, path, n_in = ancestral_node(children, root, leaves, outgroups)

    # ---- kit
    kit = {}
    for item in (x for arg in a.kit for x in arg.split(",") if x):
        if "=" not in item:
            sys.exit(f"FATAL: --kit takes ACC=STATES.npz, not {item!r}")
        acc, path_ = item.split("=", 1)
        if acc in kit:
            sys.exit(f"FATAL: --kit names {acc} twice")
        if acc not in outgroups:
            sys.exit(f"FATAL: kit genome {acc} is not in --outgroups; the kit "
                     f"holds outgroups only")
        if not os.path.exists(path_):
            sys.exit(f"FATAL: no kit states at {path_}")
        kit[acc] = path_
    kit_ev = {}
    for item in (x for arg in a.kit_events for x in arg.split(",") if x):
        if "=" not in item:
            sys.exit(f"FATAL: --kit-events takes ACC=PREFIX, not {item!r}")
        acc, pre = item.split("=", 1)
        if acc in kit_ev:
            sys.exit(f"FATAL: --kit-events names {acc} twice")
        if acc not in kit:
            sys.exit(f"FATAL: --kit-events genome {acc} is not a --kit genome")
        for suf in (".var", ".sv.tsv"):
            if not os.path.exists(pre + suf):
                sys.exit(f"FATAL: no kit events at {pre + suf}")
        kit_ev[acc] = pre

    seq = "".join(l.strip() for l in open(a.h37rv) if not l.startswith(">"))
    if not seq:
        sys.exit(f"FATAL: no sequence in {a.h37rv}")
    seq = seq.upper()
    kit_state = {}
    for acc, p in kit.items():
        d = np.load(p)
        if "state" not in d.files:
            sys.exit(f"FATAL: {p} has no 'state' array (genome_states.py's)")
        s = d["state"]
        if s.dtype != np.uint8 or s.ndim != 1 or len(s) != len(seq):
            sys.exit(f"FATAL: {p}: 'state' is {s.dtype} {s.shape}, not uint8 "
                     f"of H37Rv's length {len(seq):,}")
        if s.max() > 4:
            sys.exit(f"FATAL: {p}: state codes above 4 (0 = N, 1-4 = ACGT)")
        kit_state[acc] = s

    # ---- panel
    samples, keys, copies, st = read_panel(a.panel_vcf, a.bcftools, seq,
                                           a.h37rv_contig)
    if not keys:
        sys.exit(f"FATAL: no biallelic records in {a.panel_vcf}")
    graph = set(samples) | {a.ref_sample}
    have = graph | set(kit)
    extra, miss = sorted(leaves - have), sorted(have - leaves)
    if extra or miss:
        sys.exit(f"FATAL: the tree's leaves must be exactly the VCF's samples, "
                 f"the reference {a.ref_sample} and the kit genomes. Leaves "
                 f"with no data: {extra[:5]}{' ...' if len(extra) > 5 else ''}"
                 f" ({len(extra)}); genomes not in the tree: {miss[:5]}"
                 f"{' ...' if len(miss) > 5 else ''} ({len(miss)})")
    if a.accessions:
        want = {l.strip() for l in open(a.accessions) if l.strip()}
        got = leaves - (set(kit) - set(samples) - {a.ref_sample})
        if got != want:
            sys.exit(f"FATAL: the tree's graph genomes are not the build's "
                     f"accessions ({a.accessions}): {len(got - want)} extra, "
                     f"{len(want - got)} missing")
    both = sorted(set(kit) & set(samples))

    in_graph = [o for o in outgroups if o not in kit]
    print(f"  tree: {len(leaves)} leaves; ancestral node = MRCA of the {n_in} "
          f"non-outgroup leaves, {len(path) - 1} nodes below the root")
    print(f"  outgroups in the graph: {', '.join(in_graph) or 'none'}; "
          f"kit (SNP states only): {', '.join(kit) or 'none'}")
    if both:
        print(f"  note: {', '.join(both)} also in the VCF; the kit's states "
              f"are used and the VCF column ignored for Fitch")
    print(f"  panel: {st['records']:,} records ({st['multiallelic']:,} "
          f"multiallelic skipped), {len(keys):,} keys; {st['moved']:,} "
          f"changed by left-alignment"
          + (f"; {st['offpath']:,} off the H37Rv path, not normalised"
             if st["offpath"] else ""))

    # ---- leaf state sets
    nkey = len(keys)
    on_path = np.array([k[0].split("#")[-1] == a.h37rv_contig for k in keys])
    kinds = [kind_of(k[2], k[3]) for k in keys]
    G = np.empty((len(samples), nkey), dtype=np.int8)      # 1, 0, -1 no call
    panel_af = np.empty(nkey)
    n_called = np.zeros(nkey, dtype=np.int64)
    for j, cp in enumerate(copies):
        col = pooled_column(cp)
        g = np.array([1 if x == "1" else 0 if x not in MISSING else -1
                      for x in col], dtype=np.int8)
        G[:, j] = g
        called = int((g >= 0).sum())
        n_called[j] = called
        panel_af[j] = (g == 1).sum() / called if called else float("nan")

    # sibling rule at SNP keys: GT 0 while carrying another allele that
    # changes the same base is no call
    by_pos = collections.defaultdict(list)
    for j, (c, p, r, x) in enumerate(keys):
        by_pos[(c, p)].append(j)
    n_sib = 0
    for js in by_pos.values():
        if len(js) < 2:
            continue
        change = [j for j in js if keys[j][3][:1] != keys[j][2][:1]
                  and not keys[j][3].startswith("<")]
        for j in js:
            if kinds[j] != "SNP":
                continue
            others = [o for o in change if o != j]
            if not others:
                continue
            carries = np.any(G[:, others] == 1, axis=1) & (G[:, j] == 0)
            n_sib += int(carries.sum())
            G[carries, j] = -1

    def bits(g):
        return np.where(g == 1, ALT_BIT, np.where(g == 0, REF_BIT, ANY)
                        ).astype(np.uint8)
    leaf_bits = {}
    for i, s in enumerate(samples):
        if s not in kit:
            leaf_bits[s] = bits(G[i])
    if a.ref_sample not in leaf_bits:
        leaf_bits[a.ref_sample] = np.where(on_path, REF_BIT, ANY).astype(np.uint8)
    snp_ix = np.array([j for j in range(nkey) if kinds[j] == "SNP" and on_path[j]
                       and keys[j][2] in ACGT and keys[j][3] in ACGT], dtype=np.int64)
    snp_pos = np.array([keys[j][1] for j in snp_ix], dtype=np.int64)
    snp_ref = np.array([ACGT.index(keys[j][2]) + 1 for j in snp_ix], dtype=np.uint8)
    snp_alt = np.array([ACGT.index(keys[j][3]) + 1 for j in snp_ix], dtype=np.uint8)
    for acc, s in kit_state.items():
        b = np.full(nkey, ANY, dtype=np.uint8)
        if len(snp_ix):
            base = s[snp_pos - 1]
            b[snp_ix] = np.where(base == snp_ref, REF_BIT,
                                 np.where(base == snp_alt, ALT_BIT, ANY))
        print(f"  kit {acc}: called at {int((b != ANY).sum()):,} of "
              f"{len(snp_ix):,} SNP keys")
        if acc in kit_ev:
            ev = kit_event_states(keys, kinds, seq, a.h37rv_contig,
                                  read_kit_events(kit_ev[acc], seq, a.h37rv_contig))
            nonsnp = np.array([k != "SNP" for k in kinds])
            b = np.where(nonsnp, bits(ev), b).astype(np.uint8)
            c = collections.Counter()
            for j in np.flatnonzero(nonsnp):
                c[(kinds[j], int(ev[j]))] += 1
            print(f"  kit {acc} events: " + "; ".join(
                f"{k} {c[(k, 1)]:,} ALT / {c[(k, 0)]:,} REF / {c[(k, -1)]:,} no call"
                for k in ("MNP", "indel", "SV")))
        leaf_bits[acc] = b

    # ---- Fitch
    up_anc, fin = fitch_at(children, root, path, leaf_bits, nkey)
    # AA-1: the change count over the PANEL. Kit genomes (outside the graph)
    # are left out: their states are direct-alignment SNP calls and route A
    # event calls, not graph genotypes, and the canettii branch's changes lie
    # outside the MTBC (AA-1's fix1 measurement counted them, overstating
    # homoplasy). H37Rv and the in-graph outgroups (L5.1, L6.1) are panel
    # genomes and count.
    changes = fitch_changes(children, root, leaf_bits, nkey, skip=set(kit))
    anydata = np.zeros(nkey, dtype=bool)
    for b in leaf_bits.values():
        anydata |= b != ANY

    def state(j):
        if not anydata[j]:
            return None, "nodata"
        f = int(fin[j])
        if f == ANY:
            return None, "tied"
        return (1 if f == ALT_BIT else 0), \
            ("ingroup" if int(up_anc[j]) == f else "outgroups")

    # ---- write
    order = sorted(range(nkey), key=lambda j: (keys[j][0], keys[j][1],
                                               keys[j][2], keys[j][3]))
    tally = collections.Counter()
    snp_set = set(snp_ix.tolist())
    tmp_a, tmp_p = a.out_ancestral + ".tmp", a.out_polarity + ".tmp"
    with open(tmp_a, "w") as fa, open(tmp_p, "w") as fp:
        fa.write("chrom\tpos\tref\talt\tAA\tflag\tpolarity\tchanges\tmethod\n")
        fp.write("chrom\tpos\tref\talt\toutgroup_gt\tancestral\tpanel_af"
                 "\tn_records\tmethod\tchanges\tflag\n")
        for j in order:
            chrom, pos, ref, alt = keys[j]
            g, how = state(j)
            nch = int(changes[j])
            homo = nch >= a.homoplastic_min
            tally[(kinds[j], how)] += 1
            tally[(kinds[j], "homoplastic")] += homo
            if g == 1:
                tally[(kinds[j], "alt_ancestral")] += 1
            if j in snp_set:
                if g is None:
                    aa, pol = ".", "unknown"
                else:
                    aa = alt if g else ref
                    pol = "alt_ancestral" if g else "ref_ancestral"
                flag = aa_flags(how, nch, a.homoplastic_min,
                                "".join(sorted({ref, alt})))
                fa.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{aa}\t{flag}\t{pol}"
                         f"\t{nch}\t{how}\n")
            if g is None:
                continue
            af = panel_af[j]
            fp.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{g}\t"
                     f"{'ALT' if g else 'REF'}\t{af:.4f}\t{len(copies[j])}"
                     f"\t{how}\t{nch}\t{aa_flags(how, nch, a.homoplastic_min)}\n")
    os.replace(tmp_a, a.out_ancestral)
    os.replace(tmp_p, a.out_polarity)

    print(f"  {n_sib:,} SNP cells set to no call: a sibling allele changes "
          f"the same base")
    print("\n  ancestral state at the ancestral node, by class")
    print(f"    (outgroups: flagged OUTGROUP; homoplastic: flagged HOMOPLASTIC, "
          f">= {a.homoplastic_min} Fitch changes over the graph genomes, "
          f"kit genomes left out; any state)")
    hdr = ["ingroup", "outgroups", "tied", "nodata", "alt_ancestral",
           "homoplastic"]
    print("    class  " + "".join(f"{h:>14}" for h in hdr))
    for k in ("SNP", "MNP", "indel", "SV"):
        if any(tally[(k, h)] for h in hdr):
            print(f"    {k:<6} " + "".join(f"{tally[(k, h)]:>14,}" for h in hdr))
    print(f"  -> {a.out_ancestral}  ({len(snp_ix):,} SNP keys)")
    print(f"  -> {a.out_polarity}  (resolved keys of every class; tied and "
          f"no-data keys are left out and stay unpolarised)")


if __name__ == "__main__":
    main()
