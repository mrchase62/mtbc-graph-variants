#!/usr/bin/env python3
"""Stage 1 of P1h: enumerate junctions from the ELEMENT side of a P1g alignment.

WHAT THIS IS, AND HOW IT DIFFERS FROM is6110_junctions.py
Both read exactly the same evidence out of exactly the same BAM: a read that
crosses an insertion boundary is split by bwa, and the SA tag records the contig,
coordinate, strand and CIGAR of the other piece. Nothing here assembles a
consensus, matches a seed, or calls an aligner. No external detector is involved.

The difference is only which end you start from.

  is6110_junctions.py  starts at clipped positions ON THE CHROMOSOME contig and
                       asks whether their SA lands on the element contig.
  this script          starts at clipped reads ON THE ELEMENT contig and reads
                       the chromosomal coordinate out of their SA.

The two enumerate the same set of junctions, so either is a check on the other.
Starting from the element side matters because the chromosome-side path reaches
its answer through two filters that P1H_PLAN.md section 1 shows are dropping real
junctions: a 1 kb single-linkage clustering radius that chains four separate
junction signatures in the 3,105-3,109 kb interval of SAMN13208090 into one site,
and a fixed 0.10 read-through fraction that discards that isolate's largest
junction (888,833/888,836, a 3 bp pair carrying 68 and 60 element SA records).
Neither filter exists on this path: a read on the element contig is element
sequence by where it aligned, so attribution needs no threshold at all.

THE JUNCTION COORDINATE IS COMPUTED, NOT APPROXIMATED
The SA tag's POS is the leftmost base of the far-side block, which for a junction
read is the boundary only when the far side is clipped on its left. Taking POS
unconditionally spreads one junction over a read length. Here the SA CIGAR is
parsed and the boundary facing the clip is used, so a junction returns one
coordinate and the two flanks of an insertion sit a target-site duplication
apart, not a read length apart.

WHAT IS COUNTED, AND WHAT IS DELIBERATELY NOT FILTERED
Secondary and duplicate alignments are excluded. SUPPLEMENTARY alignments are
KEPT: for a junction read one piece is primary and the other supplementary, and
which piece lands on the element contig is arbitrary, so excluding them would
discard about half the evidence. No read-through, depth or clustering criterion
is applied at this stage -- stage 2 of P1H_PLAN.md is the accounting that decides
which criteria are justified, and applying them here first would prejudge it.

ELEMENT-ELEMENT JUNCTIONS are reported separately, as the direct test of
P1G_RESULTS.md section 4 explanation 1: copies arranged next to each other put a
boundary inside element sequence, and a read crossing it clips on the element
contig with its SA returning to the element contig.
"""
import argparse, collections, csv, os, re, sys

import pysam

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_REF_OPS = set("MDN=X")
_EXCLUDE = 0x100 | 0x400            # secondary, duplicate. 0x800 is kept.


def sa_blocks(sa_field):
    """Parse an SA:Z: field into (rname, pos, strand, cigar) tuples."""
    out = []
    for rec in sa_field.rstrip(";").split(";"):
        if not rec:
            continue
        f = rec.split(",")
        if len(f) >= 5:
            out.append((f[0], int(f[1]), f[2], f[3], int(f[4])))
    return out


def cigar_stats(cigar):
    """Return (left_clip, right_clip, reference_length) for a CIGAR string."""
    ops = _CIGAR_RE.findall(cigar)
    if not ops:
        return 0, 0, 0
    reflen = sum(int(n) for n, op in ops if op in _REF_OPS)
    left = int(ops[0][0]) if ops[0][1] in "SH" else 0
    right = int(ops[-1][0]) if ops[-1][1] in "SH" else 0
    return left, right, reflen


def junction_of(pos, cigar):
    """The 1-based boundary of an SA block on the side facing its clip.

    START  the block is clipped on its left, so the junction is its leftmost base
    END    the block is clipped on its right, so it is its rightmost base
    """
    left, right, reflen = cigar_stats(cigar)
    if left >= right:
        return pos, "START"
    return pos + reflen - 1, "END"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--element-contig", default="IS6110")
    ap.add_argument("--chrom-contig", default="NC_000962.3_isclean")
    ap.add_argument("--crossmap", default="is6110/assets/H37Rv.isclean.crossmap.tsv",
                    help="clean -> original coordinates, for the reported column only")
    ap.add_argument("--terminus-slop", type=int, default=20,
                    help="an element-side block reaching this close to a contig "
                         "end is at that terminus; further in is INTERNAL and is "
                         "reported but flagged, not silently dropped")
    ap.add_argument("--merge", type=int, default=10,
                    help="chromosomal junction coordinates within this distance "
                         "are one stack. Sized for a target-site duplication "
                         "plus a base or two of alignment wobble -- NOT for "
                         "repeat-locus scatter, which is what the 1 kb radius in "
                         "is6110_isclean_summary.py was absorbing")
    ap.add_argument("--min-sa-mapq", type=int, default=20,
                    help="threshold for the reads_q column only; nothing is "
                         "dropped at this stage")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    xmap = []
    if a.crossmap and os.path.exists(a.crossmap):
        xmap = [{k: int(v) for k, v in r.items()}
                for r in csv.DictReader(open(a.crossmap), delimiter="\t")]

    def to_orig(p):
        shift = 0
        for r in xmap:
            if p > r["clean_junction"]:
                shift = r["cum_deleted"]
            else:
                break
        return p + shift

    af = pysam.AlignmentFile(a.alignment)
    try:
        ellen = af.get_reference_length(a.element_contig)
    except KeyError:
        sys.exit(f"{a.alignment}: no contig {a.element_contig}")

    # chromosomal junction coordinate -> evidence
    hits = collections.defaultdict(lambda: collections.Counter())
    tally = collections.Counter()

    for rd in af.fetch(a.element_contig):
        if rd.flag & _EXCLUDE or rd.is_unmapped:
            continue
        tally["reads_on_element"] += 1
        ct = rd.cigartuples or []
        left = ct[0][1] if ct and ct[0][0] in (4, 5) else 0
        right = ct[-1][1] if ct and ct[-1][0] in (4, 5) else 0
        if not (left or right):
            continue
        tally["clipped"] += 1
        sa = rd.get_tag("SA") if rd.has_tag("SA") else None
        if not sa:
            tally["clipped_no_sa"] += 1
            continue

        # which terminus of the element contig does this block reach through?
        start0, end0 = rd.reference_start, rd.reference_end        # 0-based, half open
        if start0 <= a.terminus_slop and left >= right:
            el_end = "START"
        elif end0 >= ellen - a.terminus_slop and right >= left:
            el_end = "END"
        else:
            el_end = "INTERNAL"

        for rname, pos, strand, cigar, mapq in sa_blocks(sa):
            if rname == a.element_contig:
                tally["sa_to_element"] += 1
                continue
            if rname != a.chrom_contig:
                tally["sa_offtarget"] += 1
                continue
            jpos, jside = junction_of(pos, cigar)
            c = hits[jpos]
            c["reads"] += 1
            c["el_" + el_end] += 1
            c["chr_" + jside] += 1
            c["strand_" + strand] += 1
            # THE JOINT COUNT, AND WHY THE MARGINALS WERE NOT ENOUGH.
            # Orientation is in none of the three counters above. A
            # non-stranded library hits every junction from both strands, so
            # fwd/rev measures the library rather than the element: over
            # scale200, 97.3% of stacks sit between 0.3 and 0.7 forward and
            # only 2 of 2,478 are unambiguous. Both termini are reached at
            # 2,175 of those 2,478, so el_start/el_end alone say nothing
            # either. What carries orientation is WHICH TERMINUS SITS ON WHICH
            # SIDE -- the element's 5' end against the left flank is one
            # orientation and against the right flank is the other -- and that
            # pairing is available right here and was being discarded.
            c[f"pair_{el_end}_{jside}"] += 1
            # SA MAPPING QUALITY, recorded rather than applied. In the direct-
            # repeat interval around 3,106-3,109 kb a single copy scatters its
            # chromosomal side over several placements: SAMEA104061150 carries
            # one copy by every other measure and yields five stacks there, of
            # which the true junction has mean SA MAPQ 58.8 and the other four
            # have 0 to 13. That is ambiguity of placement, not extra copies,
            # and the SA record already carries the number that says so. It is
            # kept as a column because stage 2 of P1H_PLAN.md is the accounting
            # that decides which criteria are justified; filtering here would
            # prejudge it.
            c["mapq_sum"] += mapq
            if mapq > c["mapq_max"]:
                c["mapq_max"] = mapq
            if mapq >= a.min_sa_mapq:
                c["reads_q"] += 1
            tally["sa_to_chrom"] += 1

    # merge neighbouring junction coordinates into stacks
    stacks, cur = [], []
    for p in sorted(hits):
        if cur and p - cur[-1] > a.merge:
            stacks.append(cur); cur = []
        cur.append(p)
    if cur:
        stacks.append(cur)

    rows = []
    for st in stacks:
        tot = sum(hits[p]["reads"] for p in st)
        peak = max(st, key=lambda p: hits[p]["reads"])
        agg = collections.Counter()
        for p in st:
            agg.update(hits[p])
        agg["mapq_max"] = max(hits[p]["mapq_max"] for p in st)
        rows.append(dict(
            sample=a.sample, clean_pos=peak, orig_pos=to_orig(peak),
            reads=tot, reads_q=agg["reads_q"], positions=len(st),
            span=st[-1] - st[0],
            sa_mapq_max=agg["mapq_max"],
            sa_mapq_mean=round(agg["mapq_sum"] / tot, 1) if tot else 0,
            el_start=agg["el_START"], el_end=agg["el_END"],
            el_internal=agg["el_INTERNAL"],
            chr_start=agg["chr_START"], chr_end=agg["chr_END"],
            both_el_termini=int(agg["el_START"] > 0 and agg["el_END"] > 0),
            both_chr_sides=int(agg["chr_START"] > 0 and agg["chr_END"] > 0),
            fwd=agg["strand_+"], rev=agg["strand_-"],
            p_start_left=agg["pair_START_START"],
            p_start_right=agg["pair_START_END"],
            p_end_left=agg["pair_END_START"],
            p_end_right=agg["pair_END_END"]))
    rows.sort(key=lambda r: r["clean_pos"])

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else
                           ["sample", "clean_pos", "orig_pos", "reads", "reads_q",
                            "positions", "span", "sa_mapq_max", "sa_mapq_mean",
                            "el_start", "el_end", "el_internal", "chr_start",
                            "chr_end", "both_el_termini", "both_chr_sides",
                            "fwd", "rev", "p_start_left", "p_start_right",
                            "p_end_left", "p_end_right"],
                           delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"[elside] {a.sample}: {tally['reads_on_element']} reads on "
          f"{a.element_contig}, {tally['clipped']} clipped, "
          f"{tally['sa_to_chrom']} SA to chromosome, "
          f"{tally['sa_to_element']} SA back to the element contig, "
          f"{tally['clipped_no_sa']} clipped with no SA, "
          f"{tally['sa_offtarget']} SA elsewhere; "
          f"{len(rows)} stacks -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
