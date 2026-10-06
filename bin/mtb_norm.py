#!/usr/bin/env python3
"""Allele normalisation and cohort-invariant record keys.

The matrix keyed rows on `h37rv:<pos>`, which collapses distinct alleles at one
position into one row. Such a row needs a single canonical REF chosen from among
them, and `p5_keys.py` chose the first-seen record's -- so the choice depended on
sample iteration order. Adding 77 samples to a 23-sample cohort changed the
canonical REF at 17 positions and flipped 68 existing calls between ALT and REF.
A deterministic rule such as longest-REF does not fix it either, because adding a
sample can introduce a longer REF.

**Cohort-invariance requires the key and the canonical REF to be properties of the
allele alone.** So each record is normalised independently -- right-trimmed,
left-aligned, left-trimmed, the same operations `bcftools norm` performs -- and
the key carries the normalised REF and ALT. Adding samples then only ever adds
keys; it can never alter an existing one.

The cost is that two callers' different representations of the same event stay
separate rows. Left-alignment removes most of that, which is why it is done here
rather than trimming alone: without it, `TA->T` at 100 and `AT->A` at 101 are the
same deletion in two rows.
"""


def normalise(seq, pos, ref, alt):
    """Left-align and trim one allele. `pos` is 1-based, `seq` 0-based H37Rv.

    Returns (pos, ref, alt). Mirrors bcftools norm: trim the shared suffix, shift
    a pure indel as far left as the reference allows, then trim the shared prefix
    down to a single anchor base.
    """
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    if not ref or not alt or ref == alt:
        return pos, ref, alt

    # 1. trim the shared suffix, keeping at least one base on each side
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]

    # 2. left-shift a pure indel while the base before it matches the base that
    #    would be lost from the end. Only applies where one allele is a single
    #    base, which is what "pure indel" means after step 1.
    if seq and (len(ref) == 1) != (len(alt) == 1):
        while pos > 1 and len(ref) > 0 and len(alt) > 0 \
                and ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
            prev = seq[pos - 2].upper()
            ref = prev + ref[:-1]
            alt = prev + alt[:-1]
            pos -= 1

    # 3. trim the shared prefix, keeping one anchor base
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def h37rv_key(pos, ref, alt):
    return f"h37rv:{pos}:{ref}>{alt}"


def node_key(node, offset, ref, alt):
    """`offset` must be the node's FORWARD offset (forward_offset below)."""
    return f"node:{node}:{offset}:{ref}>{alt}"


def forward_offset(offset, strand, length):
    """A base's offset in the node's own (forward, GFA S-line) orientation.

    `odgi position -v` counts the offset along the WALKING direction of the
    path it was asked about, so where two references walk a node in opposite
    directions one base gets two offsets, off and L-1-off: R 201016,83,+ and
    H37Rv 201016,46,- are the same base of a 130 bp node. Every node key is
    therefore written with the forward offset, and every reader that needs a
    path position converts back with this same function (it is its own
    inverse). A `-` walk with no node length cannot be converted: None.
    """
    off = int(offset)
    if strand != "-":
        return off
    if length is None or length == "":
        return None
    return int(length) - 1 - off


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def node_forward_alleles(off, strand, ref, alt, rseq="", r_pos=0):
    """A node-frame record restated on the node's FORWARD strand (D41).

    `off` is the forward offset (forward_offset above) of the record's first
    base as the sample's reference R reads it; `strand` is `-` when R's
    sequence reads the node reverse complemented (odgi's walk flag combined
    with R being stored flipped in the panel: graph_frame.complement_needed);
    `ref`/`alt` are R's alleles; `rseq` is R's sequence (refs frame) and
    `r_pos` the record's 1-based position in it.

    Without this, one event read by references walking the node in opposite
    directions became two keys, C>T and G>A (scale200: 11 bases).

    Returns (offset, ref, alt, status). On `+` nothing changes. On `-`:
    - a SNP or MNP is reverse complemented and starts at off-(len-1);
    - a left-anchored indel, D the deleted and I the inserted bases, keeps
      its anchor on the left of the FORWARD strand: that is the base after
      R's span, complemented, at off-len(D)-1;
    - any other shape is reverse complemented as a block at off-(len(ref)-1).
    Status `off_node` (the restated record would start before the node) or
    `no_anchor_sequence` returns the record unchanged, for the caller to
    count; the key then stays in R's orientation.
    """
    if strand != "-":
        return off, ref, alt, "forward"
    rc = lambda x: x.translate(_COMP)[::-1]   # noqa: E731
    if len(ref) == len(alt):
        new = off - (len(ref) - 1)
        if new < 0:
            return off, ref, alt, "off_node"
        return new, rc(ref), rc(alt), "reversed"
    if ref[:1] == alt[:1] and (len(ref) == 1 or len(alt) == 1):
        D, I = ref[1:], alt[1:]
        new = off - len(D) - 1
        if new < 0:
            return off, ref, alt, "off_node"
        i = r_pos + len(D)            # 0-based index of the base after R's span
        if not rseq or not (0 <= i < len(rseq)):
            return off, ref, alt, "no_anchor_sequence"
        anc = rseq[i].translate(_COMP).upper()
        return new, anc + rc(D), anc + rc(I), "reanchored"
    new = off - (len(ref) - 1)
    if new < 0:
        return off, ref, alt, "off_node"
    return new, rc(ref), rc(alt), "reversed_block"


def node_forward_restate(node, off, strand, ref, alt, rseq="", r_pos=0,
                         locate=None):
    """node_forward_alleles, plus the record whose forward start is on
    ANOTHER node. Returns (node, offset, ref, alt, status).

    Read in reverse, an indel's forward anchor (or an MNP's forward first
    base) is R's base AFTER the span, and where the span reaches the node's
    forward start that base is on the neighbouring node -- the node a
    reference reading forward keys the same event on. 64% of CX333's nodes
    are 1 bp, so this is 1.1% of scale200's node-frame records, not a corner.

    `locate(q)`, for R's 1-based position q, returns (node, forward offset,
    reading strand) of the base there, or None where R's node table cannot
    say (a node R visits more than once). The record moves to that node
    when R reads it reverse complemented too; otherwise (an inversion
    junction, or no answer) it is left as R reads it, status `off_node`.
    """
    o, r, a, st = node_forward_alleles(off, strand, ref, alt, rseq, r_pos)
    if st != "off_node" or locate is None:
        return node, o, r, a, st
    rc = lambda x: x.translate(_COMP)[::-1]   # noqa: E731
    if len(ref) != len(alt) and ref[:1] == alt[:1] and (
            len(ref) == 1 or len(alt) == 1):
        D, I = ref[1:], alt[1:]
        q = r_pos + len(D) + 1                 # the forward anchor, in R
        if not rseq or not (1 <= q <= len(rseq)):
            return node, o, r, a, "off_node"
        anc = rseq[q - 1].translate(_COMP).upper()
        new_ref, new_alt = anc + rc(D), anc + rc(I)
    else:
        q = r_pos + len(ref) - 1               # the forward first base, in R
        new_ref, new_alt = rc(ref), rc(alt)
    hit = locate(q)
    if hit is None or hit[2] != "-":
        return node, o, r, a, "off_node"
    return hit[0], hit[1], new_ref, new_alt, "moved_to_start_node"


def load_path_nodes(path, accession):
    """R's walk from the build's node table: a sorted list of
    (start, length, node, strand, n_occurrences), `start` 1-based along the
    PANEL path. For node_forward_restate's `locate`."""
    import csv
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["accession"] == accession and r.get("length"):
                out.append((int(r["start"]), int(r["length"]), r["node"],
                            r["strand"], int(r.get("n_occurrences") or 1)))
    out.sort()
    return out


def load_node_lengths(path, want=None):
    """{node id (str): length} from the build's node table
    (assets/node_positions.tsv, columns node and length; P0 step nodes).
    `want`, a set of node ids, restricts what is kept."""
    import csv
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n = r["node"]
            if (want is None or n in want) and r.get("length"):
                out[n] = int(r["length"])
    return out


def _selftest():
    # H37Rv-like context; index 0 is position 1
    seq = "ACGTTTTGCA" * 10
    cases = [
        # a shared-suffix pair collapses to the minimal allele
        ((5, "TTTG", "TTG"), "one base deleted from the T-run"),
        # the same deletion written at a different offset must normalise alike
        ((6, "TTG", "TG"), "same deletion, shifted representation"),
        # a SNP is untouched
        ((3, "G", "A"), "SNP"),
        # a long REF with a shared prefix trims to an anchor
        ((1, "ACGT", "A"), "3 bp deletion"),
    ]
    out = []
    for (p, r, a), label in cases:
        out.append((label, normalise(seq, p, r, a)))
    a = normalise(seq, 5, "TTTG", "TTG")
    b = normalise(seq, 6, "TTG", "TG")
    same = (a == b)
    for label, v in out:
        print(f"    {label:<38s} -> {v}")
    print(f"    two representations of one deletion agree: {same}")
    # keys must be stable and allele-specific
    k1 = h37rv_key(*normalise(seq, 5, "TTTG", "TTG"))
    k2 = h37rv_key(*normalise(seq, 5, "T", "TA"))
    print(f"    distinct alleles at one position get distinct keys: {k1 != k2}")
    print(f"      {k1}\n      {k2}")
    return same and k1 != k2


if __name__ == "__main__":
    import sys
    print("  mtb_norm self-test")
    sys.exit(0 if _selftest() else 1)
