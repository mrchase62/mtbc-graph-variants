#!/usr/bin/env python3
"""Genotype a catalogued interval from BOTH reference frames, with depth and clips.

WHY TWO FRAMES. A structural variant is invisible in the frame that already
carries it and visible in the frame that does not. A deletion absorbed into a
sample's matched panel reference presents there as contiguous sequence -- no
depth drop, no clipped reads, nothing a caller can find -- and that class is 62%
of every carrier call the interval arm makes. The same event against H37Rv is a
depth cliff with clipped reads at both ends. Symmetrically, a deletion H37Rv
carries and the matched reference lacks is measured directly in the matched
frame and not at all against H37Rv.

Worked example, `svi:DEL:453367:2606` in SAMEA1016021 (matched reference
GCF_013010385): the matched frame projects 0 of 9 probes and the existing arm
calls ALT by `inherited_reference`, which is an inference from the projection
failing. In the H37Rv frame the flanks sit at depth 80 and 87, 2,564 of the
2,606 interval bases are at zero, and 30 and 44 soft-clipped reads cluster
within 60 bp of the two breakpoints. Same call, direct measurement.

THREE CHANNELS, BECAUSE THEY FAIL DIFFERENTLY.

  depth              presence or absence across an interval. Read-length blur
                     at the edges means it cannot place a breakpoint, and it is
                     blind to insertions, where coverage is identical either
                     way.
  clips, H37Rv       base-resolution breakpoints for anything H37Rv lacks --
                     the whole inherited class, and insertions.
  clips, matched     base-resolution breakpoints for anything the matched
                     reference lacks.

THE ASYMMETRY THIS SCRIPT HAS TO RESPECT. H37Rv depth is confounded by reference
bias and matched-reference depth is not: an isolate far from H37Rv loses depth
in divergent regions for mapping reasons, not because sequence is missing, which
is the bias this project exists to characterise. So H37Rv depth alone over-calls
deletions in exactly the lineages furthest from the reference.

Clips are what separate the two. A real deletion gives a SHARP depth cliff with
clips CLUSTERED at one coordinate; divergence gives a gradual loss with clips
scattered. This script therefore reports the clip count in a tight window at
each breakpoint next to the depth, and the promoted `strong` call requires both.

WHAT IS AND IS NOT COMPUTED HERE. The H37Rv frame needs no projection at all,
because the catalogue is already in H37Rv coordinates -- no odgi, no probe
projection, nothing but the CRAM. The matched frame's call is not recomputed; it
is read from the existing `svgt_iv` output so the two are directly comparable
and this script adds a frame rather than replacing an arm.

Insertions are carried through with depth columns left empty and clip columns
filled, because clips are the only instrument that reaches them. No insertion is
promoted on depth.

DIAGNOSING A CONTRADICTION. When the two frames disagree, three columns say why,
and they are recorded for every interval rather than only the disagreeing ones:

  h_depth_q0     depth with NO mapping-quality filter. A region whose Q30 depth
                 is zero and whose Q0 depth is normal is not missing sequence --
                 it is sequence that maps ambiguously, a repeat. This separates
                 "the reads are not there" from "the reads are there and cannot
                 be placed", which depth at one threshold cannot.
  h_edge_left/right   median depth in the 100 bp just INSIDE each breakpoint
                 over the 100 bp just OUTSIDE it. A real deletion is a cliff and
                 gives a ratio near zero at a precise coordinate; divergence
                 gives a gradual slope and a middling ratio.
  h_clip_left/right   clustered clips, which only a breakpoint produces.

Expected signatures. Reference bias -- an isolate far from H37Rv losing depth
because its sequence is too divergent to map -- is low Q30 depth, low Q0 depth,
no clustered clips, soft edges. A repeat is Q30 depth zero with Q0 depth normal.
A real deletion is both depths zero, clips at both ends, hard edges.
"""
import argparse, collections, csv, os, re, subprocess, sys, tempfile

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")


def _finish(p, cmd, errf):
    """Wait for p and stop on a nonzero exit, showing samtools' own message.

    Both readers below discarded stderr and ignored the exit status. A CRAM
    decode error partway through the contig (a reference MD5 mismatch, say)
    left zero depth over the rest of it -- false depth_absent ALT calls -- and
    a failed `view` left no clips, so no call could ever be marked strong.
    """
    rc = p.wait()
    if rc != 0:
        errf.seek(0)
        sys.exit(f"FATAL: {' '.join(cmd)} exited {rc}:\n{errf.read()[-2000:]}")


def depth_array(samtools, aln, ref, contig, length, min_mapq):
    """Per-base depth over the whole contig, as a list indexed by position."""
    cmd = [samtools, "depth", "-a", "-r", contig, "--reference", ref]
    if min_mapq:
        cmd += ["-Q", str(min_mapq)]
    cmd.append(aln)
    d = [0] * (length + 2)
    errf = tempfile.TemporaryFile("w+")
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errf,
                         text=True, bufsize=1 << 20)
    n = 0
    for line in p.stdout:
        f = line.split("\t")
        if len(f) < 3:
            continue
        d[int(f[1])] = int(f[2])
        n += 1
    _finish(p, cmd, errf)
    if n == 0:
        sys.exit(f"FATAL: samtools depth returned nothing for {aln}. "
                 f"Check the contig name and the reference.")
    return d


def clip_positions(samtools, aln, ref, min_mapq, min_clip):
    """Junction coordinates of soft-clipped reads, split by which end clipped.

    A read that crosses the LEFT breakpoint of a deletion aligns up to the base
    before it and then clips, so its junction is its LAST aligned base and it
    lands in `right_end`. A read crossing the RIGHT breakpoint aligns from the
    base after it, so its junction is its FIRST aligned base, in `left_end`.
    Both are counted separately; a deletion shows one of each, at the two ends.
    """
    cmd = [samtools, "view", "--reference", ref]
    if min_mapq:
        cmd += ["-q", str(min_mapq)]
    cmd.append(aln)
    lead, trail = collections.Counter(), collections.Counter()
    errf = tempfile.TemporaryFile("w+")
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errf,
                         text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t", 6)
        if len(f) < 6:
            continue
        cig = f[5]
        if "S" not in cig:
            continue
        parts = CIGAR.findall(cig)
        if not parts:
            continue
        pos = int(f[3])
        if parts[0][1] == "S" and int(parts[0][0]) >= min_clip:
            lead[pos] += 1
        if parts[-1][1] == "S" and int(parts[-1][0]) >= min_clip:
            span = sum(int(n) for n, o in parts if o in REF_CONSUMING)
            trail[pos + span - 1] += 1
    _finish(p, cmd, errf)
    return lead, trail


def load_crossmap(path):
    """Original H37Rv position -> IS-clean position, from the excision table.

    Excising sixteen copies shifts everything downstream, so a coordinate means
    nothing in the IS-clean frame until it is converted. The table gives each
    excised span and the bases removed. A position INSIDE an excised span has no
    IS-clean equivalent at all -- it WAS the element -- and the converter says so
    rather than returning a neighbouring base.
    """
    spans = []
    if not os.path.exists(path):
        return spans
    for r in csv.DictReader(open(path, newline=""), delimiter="\t"):
        spans.append((int(r["orig_start"]), int(r["orig_end"]),
                      int(r["deleted_len"])))
    spans.sort()
    return spans


def to_clean(spans, pos):
    cum = 0
    for s, e, d in spans:
        if s <= pos <= e:
            return None
        if e < pos:
            cum += d
        else:
            break
    return pos - cum


def overlaps_excision(spans, s, e):
    """Does [s, e] touch any excised copy?

    TESTING THE ENDPOINTS IS NOT ENOUGH, and getting this wrong produced a
    clean-looking artefact. svi:DEL:2784613:1358 spans 2784613-2785974 and the
    excised copy is 2784614-2785970 -- strictly inside it. Both endpoints
    convert, so the endpoint test passed, but the 1,357 excised bases between
    them are gone in the IS-clean frame and the interval collapses to about
    five bases of ordinary depth. Every one of 198 isolates then read REF with
    consistency 1.00, which looks like a finding and is a coordinate error.
    An interval that overlaps an excised copy has no IS-clean representation at
    all and stays on plain H37Rv, where depth over the copy is exactly the
    measurement wanted: does this sample carry the element H37Rv carries?
    """
    return any(not (e < a or s > b) for a, b, _ in spans)


def _med(x):
    return sorted(x)[len(x) // 2] if x else 0


def window_sum(counter, centre, half):
    return sum(counter.get(p, 0) for p in range(centre - half, centre + half + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--intervals", required=True,
                    help="the catalogue, H37Rv coordinates, from bin/sv_intervals.py")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cram", required=True, help="the H37Rv alignment")
    ap.add_argument("--h37rv-fasta", required=True)
    ap.add_argument("--contig", default="NC_000962.3")
    ap.add_argument("--contig-length", type=int, default=4411532)
    ap.add_argument("--matched-states", default="",
                    help="this sample's svgt_iv table, so the matched frame's "
                         "call is carried alongside rather than recomputed")
    ap.add_argument("--min-dp", type=int, default=5)
    ap.add_argument("--present-frac", type=float, default=0.8)
    ap.add_argument("--deleted-frac", type=float, default=0.2)
    ap.add_argument("--flank", type=int, default=200)
    ap.add_argument("--clip-window", type=int, default=25,
                    help="half-width for counting clips at a breakpoint. Tight "
                         "on purpose: a clustered clip is the signal, a "
                         "scattered one is divergence")
    ap.add_argument("--min-clip", type=int, default=20,
                    help="a clip shorter than this cannot be placed and is "
                         "background, so it is not counted")
    ap.add_argument("--min-clips", type=int, default=3,
                    help="clipped reads needed at a breakpoint to call it "
                         "supported")
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--isclean-bam", default="",
                    help="this sample's alignment to the IS-clean H37Rv, from "
                         "pass p1g, used for element-proximal intervals ONLY. "
                         "On plain H37Rv a read carrying element sequence has "
                         "sixteen places to align, so the aligner can hide a "
                         "junction read inside a reference copy and depth there "
                         "counts mismapped copies rather than the locus. The "
                         "IS-clean backbone carries no element sequence, so such "
                         "a read must clip at its true junction.")
    ap.add_argument("--isclean-crossmap",
                    default="is6110/assets/H37Rv.isclean.crossmap.tsv")
    ap.add_argument("--isclean-contig", default="NC_000962.3_isclean")
    ap.add_argument("--isclean-length", type=int, default=4389830)
    ap.add_argument("--isclean-fasta", default="is6110/assets/H37Rv.isclean.fasta")
    ap.add_argument("--samtools", default=os.environ.get("MTB_SAMTOOLS", "samtools"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ivs = list(csv.DictReader(open(a.intervals, newline=""), delimiter="\t"))
    if not ivs:
        sys.exit(f"FATAL: no intervals in {a.intervals}")

    matched = {}
    if a.matched_states and os.path.exists(a.matched_states):
        for r in csv.DictReader(open(a.matched_states, newline=""), delimiter="\t"):
            matched[r["interval"]] = (r["state"], r["evidence"])
        print(f"  matched-frame calls read for {len(matched):,} intervals")

    print(f"  {a.sample}: reading H37Rv depth at Q>={a.min_mapq} ...")
    dep = depth_array(a.samtools, a.cram, a.h37rv_fasta, a.contig,
                      a.contig_length, a.min_mapq)
    # The unfiltered pass is what tells a repeat from an absence, and it costs
    # one more sweep of a 4.4 Mb contig, so it is always taken.
    print(f"  {a.sample}: reading H37Rv depth unfiltered ...")
    dep0 = depth_array(a.samtools, a.cram, a.h37rv_fasta, a.contig,
                       a.contig_length, 0)
    gmed = sorted(x for x in dep[1:] if x)
    gmed = gmed[len(gmed) // 2] if gmed else 0
    print(f"  genome median depth {gmed}")
    print(f"  {a.sample}: collecting clip junctions ...")
    lead, trail = clip_positions(a.samtools, a.cram, a.h37rv_fasta,
                                 a.min_mapq, a.min_clip)
    print(f"  clip junctions: {sum(lead.values()):,} leading, "
          f"{sum(trail.values()):,} trailing")

    # ---- the IS-clean frame, for element-proximal intervals only ----------
    icd = icd0 = iclead = ictrail = None
    spans = []
    if a.isclean_bam and os.path.exists(a.isclean_bam):
        spans = load_crossmap(a.isclean_crossmap)
        if not spans:
            sys.exit(f"FATAL: --isclean-bam given but no crossmap at "
                     f"{a.isclean_crossmap}; a coordinate in the IS-clean frame "
                     f"is meaningless without it")
        print(f"  {a.sample}: reading the IS-clean frame, {len(spans)} excised spans ...")
        icd = depth_array(a.samtools, a.isclean_bam, a.isclean_fasta,
                          a.isclean_contig, a.isclean_length, a.min_mapq)
        icd0 = depth_array(a.samtools, a.isclean_bam, a.isclean_fasta,
                           a.isclean_contig, a.isclean_length, 0)
        iclead, ictrail = clip_positions(a.samtools, a.isclean_bam,
                                         a.isclean_fasta, a.min_mapq, a.min_clip)
        print(f"  IS-clean clip junctions: {sum(iclead.values()):,} leading, "
              f"{sum(ictrail.values()):,} trailing")

    counts, conc = collections.Counter(), collections.Counter()
    rows = []
    n_clean = 0
    for r in ivs:
        s, e = int(r["start"]), int(r["end"])
        if s < 1 or e > a.contig_length or e < s:
            continue
        # Element-proximal intervals are measured on the IS-clean backbone
        # when one is given, in ITS coordinates. An interval whose breakpoints
        # fall inside an excised copy has no IS-clean equivalent -- the interval
        # IS the element -- and stays on plain H37Rv, flagged in frame_used.
        D, D0, LEAD, TRAIL = dep, dep0, lead, trail
        S, E = s, e
        frame_used = "h37rv"
        if icd is not None and str(r.get("is6110_prox", "0")) in ("1", "True"):
            cs, ce = to_clean(spans, s), to_clean(spans, e)
            if (not overlaps_excision(spans, s, e)
                    and cs is not None and ce is not None and ce > cs):
                D, D0, LEAD, TRAIL = icd, icd0, iclead, ictrail
                S, E = cs, ce
                frame_used = "isclean"
                n_clean += 1
            else:
                frame_used = "h37rv_inside_excision"
        inner = D[S:E + 1]
        inner0 = D0[S:E + 1]
        covered = sum(1 for x in inner if x >= a.min_dp)
        frac = covered / max(len(inner), 1)
        lf = D[max(1, S - a.flank):S]
        rf = D[E + 1:E + 1 + a.flank]
        flank_ok = (sum(1 for x in lf if x >= a.min_dp) >= 0.5 * max(len(lf), 1)
                    or sum(1 for x in rf if x >= a.min_dp) >= 0.5 * max(len(rf), 1))
        # A deletion's left breakpoint is a TRAILING clip just before `start`;
        # its right breakpoint is a LEADING clip just after `end`.
        cl = window_sum(TRAIL, S - 1, a.clip_window)
        cr = window_sum(LEAD, E + 1, a.clip_window)
        clipped = (cl >= a.min_clips) + (cr >= a.min_clips)
        # Edge shape: median depth just inside each breakpoint over the median
        # just outside it. Near zero at a cliff, middling on a slope, -1 when
        # the outside window has no depth to divide by.
        w = 100
        out_l, in_l = _med(D[max(1, S - w):S]), _med(D[S:S + w])
        in_r, out_r = _med(D[max(S, E - w):E + 1]), _med(D[E + 1:E + 1 + w])
        edge_l = (in_l / out_l) if out_l else -1.0
        edge_r = (in_r / out_r) if out_r else -1.0

        if not flank_ok:
            st, why = "NOCALL", "flanks not covered in the H37Rv frame"
        elif frac >= a.present_frac:
            st, why = "REF", "depth_present"
        elif frac <= a.deleted_frac:
            st = "ALT"
            why = ("depth_absent+clips_both" if clipped == 2 else
                   "depth_absent+clips_one" if clipped == 1 else
                   "depth_absent_no_clips")
        elif clipped == 2:
            st, why = "ALT", "clips_both_ambiguous_depth"
        else:
            st, why = "NOCALL", "ambiguous depth, no clustered clips"
        counts[st] += 1

        mst, mwhy = matched.get(r["interval"], ("", ""))
        if mst:
            if mst == st:
                agree = "agree"
            elif "NOCALL" in (mst, st) or "ABSENT" in (mst, st):
                agree = "one_frame_only"
            else:
                agree = "conflict"
            conc[f"{mst}->{st}:{agree}"] += 1
        else:
            agree = ""
        # STRONG means two independent instruments, not two restatements of one.
        # Depth and clips in the H37Rv frame are independent of each other; the
        # matched frame agreeing is a third. A depth-only ALT in the H37Rv frame
        # is exactly the call reference bias can fake, so it is never strong.
        strong = int(st == "ALT" and (clipped >= 1 or
                                      (mst == "ALT" and frac <= a.deleted_frac)))
        rows.append(dict(
            sample=a.sample, interval=r["interval"], svtype=r["svtype"],
            source=r.get("source", ""), support_tier=r.get("support_tier", ""),
            start=s, end=e, svlen=r.get("svlen", ""),
            h_state=st, h_evidence=why,
            h_covered_frac=f"{frac:.3f}",
            h_depth_median=_med(inner), h_depth_q0=_med(inner0),
            h_edge_left=f"{edge_l:.2f}", h_edge_right=f"{edge_r:.2f}",
            h_genome_median=gmed, frame_used=frame_used,
            h_clip_left=cl, h_clip_right=cr, h_clip_sides=clipped,
            m_state=mst, m_evidence=mwhy, concordance=agree, strong=strong))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    if icd is not None:
        print(f"  {n_clean:,} element-proximal intervals measured on the "
              f"IS-clean backbone")
    tot = sum(counts.values())
    print(f"  H37Rv-frame states over {tot:,} intervals: " +
          "  ".join(f"{k} {v:,} ({v/tot:.1%})" for k, v in sorted(counts.items())))
    print(f"  strong ALT (depth and clips, or both frames): "
          f"{sum(r['strong'] for r in rows):,}")
    if conc:
        print("  matched -> H37Rv, most common transitions:")
        for k, v in conc.most_common(8):
            print(f"    {k:<34} {v:,}")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
