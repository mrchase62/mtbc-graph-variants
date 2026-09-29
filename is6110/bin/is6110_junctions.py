#!/usr/bin/env python3
"""Resolve the far side of a clipped-read junction from the SA tag, and attribute
it to IS6110 by where that far side lands.

Ported from Peter Culviner's `260728_sv_exploration/260825_check_ends.ipynb`
(`fetchClipEndSA` / `fetchJunctions`). The notebook's version is written against
`samarray.CRAMArray`. samarray IS available -- Michael's clone is at
/n/netscratch/sfortune_lab/Lab/mchase/samarray and its CLAUDE.md names the
interpreter that carries the compiled extension
(/n/boslfs02/LABS/sfortune_lab/Lab/culviner/envs/samarray_env/bin/python) -- but
this port deliberately uses plain `pysam` so the junction reader runs on any BAM
or CRAM under the pipeline's own interpreter, with no TileDB collection needed. The port is
validated by reproducing the notebook's two worked examples exactly -- see
`--selftest`.

WHY THIS EXISTS
`is6110/bin/panisa_score.py` attributes a panISa site to IS6110 by matching the
clipped consensus sequence against the element's termini: a 20 bp seed with a
mismatch budget, anchored at one offset. That is a sequence-similarity test with
tunable parameters and its own failure modes, and it had been blamed for an
attribution defect it did not cause (`is6110/docs/COPY_NUMBER_FIX.md` section 1).

bwa has already solved the same problem and written the answer into the BAM. When
a read spans a junction, the part that does not fit the local alignment is
realigned and reported in the SA tag, giving the contig, position, strand and
CIGAR of the far side. So the far side can be READ rather than inferred, and a
junction is IS6110 when its far side lands at the terminus of a known IS6110
copy. No seed length, no mismatch budget, no consensus assembly.

Worked example, from the notebook, isolate SAMEA1118335 at H37Rv position 1593:

    clips END 266 -> (3796415,'END') x116  (3552708,'START') x114  (2634039,'START') x2

H37Rv's IS6110 copies span 3,795,058-3,796,412 and 3,552,713-3,554,067, so the
far side sits 3 bp past one element's end and 5 bp before another's start. That
is the signature: the clipped reads at the insertion point continue INTO an
element, and the SA tag says which one.

WHAT THIS DOES AND DOES NOT MEASURE
It detects NON-REFERENCE junctions. A copy the sample shares with H37Rv produces
no junction at all, so this cannot count total copies and must not be summed with
the family-depth estimate (`is6110_family_depth.py`). The two answer different
questions; see `is6110/docs/COPY_NUMBER_FIX.md` section 4.2.

TERMINOLOGY, kept identical to the notebook so the two can be compared
  START   the clipped end is the read's leftmost aligned base
  END     the clipped end is its rightmost aligned base
The same labelling is applied to the far side, so the label says which direction
the alignment leaves the reference through.
"""
import argparse, collections, csv, re, sys

import pysam

JunctionSummary = collections.namedtuple(
    "JunctionSummary", "clips junctions readthrough")

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_CLIP_OPS = {4, 5}                          # CIGAR S, H
_EXCLUDE_FLAGS = 0x4 | 0x100 | 0x400        # unmapped, secondary, duplicate


def _ref_len(cigar):
    """Reference bases consumed by a CIGAR string."""
    return sum(int(n) for n, op in _CIGAR_RE.findall(cigar) if op in "MDN=X")


def fetch_junctions(aln, contig, position, min_flank=20, min_mapq=0,
                    exclude_flags=_EXCLUDE_FLAGS):
    """Predicted junction sites for clipped reads ending at *position* (0-based).

    For each SA record the far side falls on whichever of the supplementary
    segment's two ends the read enters or leaves through, accounting for a strand
    flip between the two segments.

    Read-through flanks are measured in reference span, so a read carrying a
    deletion near *position* counts the deleted bases toward its flank.

    Returns
    -------
    JunctionSummary
        ``clips``       Counter keyed 'START'/'END', clipped ends at *position*
        ``junctions``   Counter keyed (contig, pos, 'START'/'END'), one count per
                        SA record
        ``readthrough`` int, reads whose aligned block extends more than
                        *min_flank* reference bases beyond *position* both sides
    """
    clips = collections.Counter()
    junctions = collections.Counter()
    readthrough = 0

    for read in aln.fetch(contig=contig, start=position, stop=position + 1):
        if read.flag & exclude_flags or read.mapping_quality < min_mapq:
            continue
        cigar = read.cigartuples
        if not cigar:
            continue

        # a CIGAR spanning no reference base is treated as spanning one (bam_endpos)
        span = max(read.reference_length or 0, 1)
        left, right = read.reference_start, read.reference_start + span - 1

        if left + min_flank < position < right - min_flank:
            readthrough += 1

        for at_left, (op, _), end_pos in ((True, cigar[0], left),
                                          (False, cigar[-1], right)):
            if end_pos != position or op not in _CLIP_OPS:
                continue
            clips["START" if at_left else "END"] += 1

            # the clipped end is the read's 3' end -- so the SA segment continues
            # from it -- whenever the clip sits on the strand-facing far side
            sa_follows = at_left == read.is_reverse
            entries = read.get_tag("SA").split(";") if read.has_tag("SA") else []
            for entry in filter(None, entries):
                sa_contig, sa_pos, sa_strand, sa_cigar = entry.split(",")[:4]
                sa_pos = int(sa_pos) - 1
                if sa_follows != (sa_strand == "-"):
                    junctions[sa_contig, sa_pos, "START"] += 1
                else:
                    junctions[sa_contig, sa_pos + _ref_len(sa_cigar) - 1, "END"] += 1

    return JunctionSummary(clips, junctions, readthrough)


def load_elements(gff):
    """(start, end) of each IS6110 copy, converted to 0-based inclusive."""
    els = []
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) > 4:
            els.append((int(f[3]) - 1, int(f[4]) - 1))
    return sorted(els)


def attribute_by_contig(junctions, element_contig):
    """Attribution on an IS-clean reference: the SA target IS the element contig.

    On is6110/assets/H37Rv.isclean.fasta the chromosome carries no element
    sequence at all, so a junction read's clipped part can only realign to the
    separate IS6110 contig. That makes attribution a contig-name test rather than
    a "within tol bp of one of sixteen termini" test -- no tolerance, no ambiguity
    about which copy, and no dependence on the element annotation being right.
    """
    hit = other = 0
    where = collections.Counter()
    for (contig, pos, _side), n in junctions.items():
        if contig == element_contig:
            hit += n
            where[pos // 100 * 100] += n
        else:
            other += n
    return hit, other, where


def attribute(junctions, els, tol):
    """Split junction counts into those landing at an IS6110 terminus and the rest.

    A far side counts as IS6110 when it sits within *tol* bp of either terminus of
    any reference copy. Both termini qualify because a read may enter an element
    from either end, and the element is present in both orientations in H37Rv.
    """
    hit = other = 0
    elements = collections.Counter()
    for (contig, pos, _side), n in junctions.items():
        matched = None
        for s, e in els:
            if abs(pos - s) <= tol or abs(pos - e) <= tol:
                matched = s
                break
        if matched is None:
            other += n
        else:
            hit += n
            elements[matched] += n
    return hit, other, elements


SELFTEST = [
    # (sample, position, expected clips, expected junctions) from the notebook
    ("SAMEA1118335", 1593,
     {"END": 266},
     {("NC_000962.3", 3796415, "END"): 116,
      ("NC_000962.3", 3552708, "START"): 114,
      ("NC_000962.3", 2634039, "START"): 2}),
    ("SAMEA1118335", 1591,
     {"START": 291, "END": 4},
     {("NC_000962.3", 3554069, "END"): 146,
      ("NC_000962.3", 3795054, "START"): 106,
      ("NC_000962.3", 3122551, "END"): 2,
      ("NC_000962.3", 3552711, "START"): 2,
      ("NC_000962.3", 3796412, "END"): 1,
      ("NC_000962.3", 1543864, "START"): 1,
      ("NC_000962.3", 3797939, "END"): 1,
      ("NC_000962.3", 889062, "START"): 1,
      ("NC_000962.3", 492731, "END"): 1,
      ("NC_000962.3", 1543810, "START"): 1}),
]


def selftest(a):
    """Reproduce the notebook's two worked examples exactly.

    This is the whole validation of the port: same isolate, same positions, same
    counts. If it does not match, nothing downstream should be believed.
    """
    aln = pysam.AlignmentFile(a.alignment, reference_filename=a.reference)
    ok = True
    for sample, pos, want_clips, want_junc in SELFTEST:
        got = fetch_junctions(aln, a.contig, pos, min_flank=10)
        c_ok = dict(got.clips) == want_clips
        j_ok = dict(got.junctions) == want_junc
        ok = ok and c_ok and j_ok
        print(f"  {sample} at {pos}:")
        print(f"    clips     {'OK ' if c_ok else 'MISMATCH'} {dict(got.clips)}")
        if not c_ok:
            print(f"      expected {want_clips}")
        print(f"    junctions {'OK ' if j_ok else 'MISMATCH'} "
              f"{len(got.junctions)} distinct, {sum(got.junctions.values())} records")
        if not j_ok:
            for k in sorted(set(want_junc) | set(got.junctions)):
                if want_junc.get(k) != got.junctions.get(k):
                    print(f"      {k}: expected {want_junc.get(k)}, got {got.junctions.get(k)}")
    print("\n  selftest: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--alignment", required=True, help="BAM or CRAM")
    ap.add_argument("--reference", help="reference FASTA, required for CRAM")
    ap.add_argument("--contig", default="NC_000962.3")
    ap.add_argument("--sample", default="")
    ap.add_argument("--sites", help="TSV with a 'pos' column, 1-based")
    ap.add_argument("--gff", help="H37Rv IS6110 element GFF, enables attribution")
    ap.add_argument("--min-element-sa", type=int, default=5,
                    help="with --element-contig, SA records on the element contig "
                         "required to call a junction IS6110")
    ap.add_argument("--element-contig", default=None,
                    help="attribute by SA target contig instead of by proximity to "
                         "an element terminus; use on the IS-clean reference, where "
                         "the element is a separate contig (e.g. IS6110)")
    ap.add_argument("--tol", type=int, default=20,
                    help="bp from an element terminus that still counts as that "
                         "element; the notebook's examples land within 5 bp")
    ap.add_argument("--min-flank", type=int, default=10)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--min-clips", type=int, default=3,
                    help="clipped ends required at a position before it is reported")
    ap.add_argument("--out")
    ap.add_argument("--selftest", action="store_true",
                    help="reproduce the notebook's worked examples and exit")
    a = ap.parse_args()

    if a.selftest:
        return selftest(a)
    if not a.sites:
        ap.error("--sites is required unless --selftest is given")

    aln = pysam.AlignmentFile(a.alignment, reference_filename=a.reference)
    els = load_elements(a.gff) if (a.gff and not a.element_contig) else []

    rows = []
    for r in csv.DictReader(open(a.sites), delimiter="\t"):
        pos1 = int(r["pos"])
        j = fetch_junctions(aln, a.contig, pos1 - 1,
                            min_flank=a.min_flank, min_mapq=a.min_mapq)
        n_clip = sum(j.clips.values())
        if n_clip < a.min_clips:
            continue
        if a.element_contig:
            hit, other, elements = attribute_by_contig(j.junctions, a.element_contig)
        elif els:
            hit, other, elements = attribute(j.junctions, els, a.tol)
        else:
            hit, other, elements = (0, 0, collections.Counter())
        top = elements.most_common(2)
        rows.append(dict(
            sample=a.sample, pos=pos1,
            clips_start=j.clips.get("START", 0), clips_end=j.clips.get("END", 0),
            readthrough=j.readthrough,
            two_sided=int(bool(j.clips.get("START")) and bool(j.clips.get("END"))),
            sa_records=sum(j.junctions.values()),
            sa_at_is6110=hit, sa_elsewhere=other,
            # On the IS-clean reference the element contig is reachable ONLY by
            # genuine element sequence, so any solid stack of SA records there is
            # decisive and a majority rule is wrong: a real junction routinely has
            # more records pointing elsewhere on the chromosome for unrelated
            # reasons. Measured on SAMEA104061150, the true site had 74 element
            # records against 125 others and the majority rule scored it 0.
            # On normal H37Rv the majority rule stands, because "near a terminus"
            # is a much weaker test that noise can satisfy.
            is6110=int(hit >= a.min_element_sa) if a.element_contig
                   else int(hit > 0 and hit >= other),
            top_elements=";".join(f"{p+1}:{n}" for p, n in top)))

    if not rows:
        print(f"  {a.sample}: no site reached --min-clips {a.min_clips}", file=sys.stderr)
        return 0
    out = a.out or "/dev/stdout"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    n_is = sum(r["is6110"] for r in rows)
    n_two = sum(r["two_sided"] for r in rows)
    print(f"  {a.sample}: {len(rows)} sites with >= {a.min_clips} clipped ends; "
          f"{n_is} attributed to IS6110 by SA tag, {n_two} two-sided")
    return 0


if __name__ == "__main__":
    sys.exit(main())
