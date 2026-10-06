#!/usr/bin/env python3
"""P4: place one isolate's calls into the common frame, routed by region.

The routing is what the tests measured, not a preference:

  core sequence      DIRECT arm, calls against H37Rv     PPV 0.9546 vs 0.9135
  PE/PPE and masked  COMPOSED arm, calls against R       PPV 0.8664 vs 0.7217

T15's hybrid of those two rules scored 0.9366 with 34% fewer false positives than
either arm alone. T16 tried to make the routing unnecessary by moving to the graph
frame and failed: composition is uniformly ~13% wrong on inherited differences and
no region rule separates those errors, so the repeat mask -- which does separate
them -- stays the boundary.

The two arms partition, they do not overlap: a direct-arm record in masked
sequence is dropped because the composed arm covers it, and a composed-arm record
in core sequence is dropped because the direct arm covers it. Every kept record
carries the arm that produced it, so the routing is auditable rather than baked in.

"Composed" means what it meant in T14 and T15: the sample's own calls against R,
plus the differences R itself carries from H37Rv, which the sample inherits by
not calling against them. Dropping the inherited half would lose recall, which is
why composition is used in PE/PPE at all despite being the less precise mechanism.
Each record says which half it came from, since the inherited half is the one
T16 showed carries 89.2% of composition's false positives.

Off-path records are keyed on the graph NODE. Two samples matched to different
references carrying the same accessory locus land on the same node, which is
exact; the H37Rv anchor agrees only to within tens of bp and the accessory
catalogue has no locus within 50 bp for 29% of accessory-scale calls. Node ids are
build-scoped, which is what the build id in every row is for.
"""
import argparse, bisect, collections, csv, gzip, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mtb_norm  # noqa: E402
import importlib.util as _ilu  # noqa: E402
_gfs = _ilu.spec_from_file_location(
    "graph_frame", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "graphframe", "bin", "graph_frame.py"))
graph_frame = _ilu.module_from_spec(_gfs); _gfs.loader.exec_module(graph_frame)


def op(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def strip_gap(x):
    return x.replace("-", "").replace("*", "")


def classify(ref, alt):
    r, a = strip_gap(ref), strip_gap(alt)
    return "SNP" if (len(r) == 1 and len(a) == 1) else "INDEL"


def load_vcf(path):
    out = []
    for line in op(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        if f[6] not in (".", "PASS"):
            continue
        ref = f[3].upper()
        for alt in f[4].upper().split(","):
            if set(ref + alt) - set("ACGTN"):
                continue
            out.append((int(f[1]), ref, alt, f[5]))
    return out


def load_intervals(bed):
    """PE/PPE and other masked intervals, each MERGED and sorted.

    The class is the mask's own first name field (`pe_ppe|PPE1|named`), which
    build_repeat_mask.py writes for every row; the old test, "PE" anywhere in
    the upper-cased name, is kept only for a BED without that field, since it
    also matches unrelated names that merely contain the letters.
    """
    pe, other = [], []
    for line in open(bed):
        f = line.split()
        if len(f) < 3:
            continue
        s, e = int(f[1]), int(f[2])
        name = f[3] if len(f) > 3 else ""
        klass = name.split("|", 1)[0].lower() if "|" in name else ""
        is_pe = (klass == "pe_ppe") if klass else ("PE" in name.upper())
        (pe if is_pe else other).append((s, e))
    return merge_intervals(pe), merge_intervals(other)


def merge_intervals(iv):
    """Union of half-open intervals, sorted and non-overlapping.

    The mask is heavily nested -- PE_PGRS4 336359-339273 contains
    336559-339142 -- and a lookup that tests only the last interval starting
    at or before a position misses every base of the outer interval past the
    inner one's end: 13% of masked bases were routed as core that way.
    Merging first makes the single-interval test exact.
    """
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def make_hit(iv):
    """Membership test for a 1-based position against 0-based half-open BED.

    BED [s, e) covers 1-based positions s+1 .. e, so the test is s < p <= e.
    The old `s <= p < e` compared a VCF POS as if it were 0-based, which
    misassigned the base on each side of every interval. `iv` must be merged.
    """
    starts = [x[0] for x in iv]
    def hit(p):
        i = bisect.bisect_left(starts, p) - 1      # last start strictly < p
        return i >= 0 and iv[i][0] < p <= iv[i][1]
    return hit


def load_seq(path):
    op2 = gzip.open if path.endswith(".gz") else open
    with op2(path, "rt") as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">")).upper()


def parse_pos_file(path, want_node=False):
    """odgi position output, keyed by SOURCE position.

    NOT by file order. `odgi position -t 4` returns results in thread-completion
    order, not input order -- verified directly: feeding positions 1000..8000
    returns 4000 5000 6000 7000 8000 1000 2000 3000. Pairing by order therefore
    mispairs every coordinate whenever more than one thread is used, and a
    count-equality check cannot detect a reordering. The source position is in
    column 1, so there is no reason to infer it from position in the file.
    """
    out = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        try:
            src = int(f[0].rsplit(",", 2)[-2]) + 1
            if want_node:
                # odgi's node offset runs along R's WALK of the node; the
                # strand is kept so main() can restate it as the node's
                # forward offset (mtb_norm.forward_offset)
                nid, off, nst = f[1].split(",")
                val = (int(nid), int(off), nst.strip())
            else:
                tgt = f[1].rsplit(",", 2)
                # f[3] is strand.vs.ref as REPLACED by graphframe/bin/
                # frame_convert.py: the relation between the two refs
                # sequences, which is `-` only when the panel stores one of the
                # two accessions reverse complemented. odgi's own flag is
                # carried in f[4].
                t = int(tgt[-2]) + 1
                st = f[3] if len(f) > 3 and f[3] in ("+", "-") else "+"
                # WHERE f[4] IS `-` (R and H37Rv walk the node in opposite
                # directions) THE TARGET IS ONE BASE HIGH AND R READS
                # COMPLEMENTED to what f[3] says. The earlier note here, that
                # odgi already accounts for the inverted step, rested on 4
                # pilot SNPs. Measured instead by reading R's 31-mer around the
                # source against H37Rv's around the target:
                #   f[3] f[4]  homolog            GCF_000193185  scale200 P4
                #   +    -     t-1, complemented  8,684 / 8,694   167 / 228
                #   -    -     t-1, same strand        -              7 / 7
                #   +    +     t, same strand    34,298 / 34,465
                # (bin/p5_states.py corrects the H37Rv -> R direction the same
                # way.) The shift is -1 because the target here is always
                # H37Rv, which the panel stores forward; for a target stored
                # reverse complemented it would be +1.
                if len(f) > 4 and f[4].strip() == "-":
                    t -= 1
                    st = "-" if st == "+" else "+"
                val = (t, int(f[2]) if len(f) > 2 else 0, st)
        except (ValueError, IndexError):
            continue
        out.setdefault(src, val)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--build-id", required=True)
    ap.add_argument("--direct", required=True)
    ap.add_argument("--matched", required=True)
    ap.add_argument("--h37rv-pos", required=True)
    ap.add_argument("--node-pos", required=True)
    ap.add_argument("--mask", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--graph-vcf", required=True,
                    help="source of the inherited half of the composed arm")
    ap.add_argument("--acc-tol", type=int, default=50)
    ap.add_argument("--near-tol", type=int, default=50,
                    help="dist.to.ref at or below this is a small-indel interior, "
                         "not accessory sequence")
    ap.add_argument("--h37rv", default="",
                    help="H37Rv FASTA. Needed to re-anchor an indel that "
                         "projects onto the other strand; without it those "
                         "records are left as called and reported.")
    ap.add_argument("--ref-fasta", default="",
                    help="the matched reference's FASTA, for the anchor base of "
                         "an OFF-path record, where H37Rv is not the reference")
    ap.add_argument("--node-lengths", default="",
                    help="the build's node table (assets/node_positions.tsv, "
                         "P0 step nodes), for node lengths. A node key's "
                         "offset is the node's FORWARD offset, and where R "
                         "walks the node in reverse that is L-1-offset; "
                         "without the length such a record cannot be keyed "
                         "and is dropped, counted")
    ap.add_argument("--frames", default="",
                    help="the build's graph frame table (default: "
                         "MTB_GRAPH_FRAMES or MTB_BUILD_DIR's), to tell whether "
                         "R is stored reverse complemented in the panel; needed "
                         "to write node-frame alleles on the node's forward "
                         "strand")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")

    def rc(x):
        return x.translate(comp)[::-1]

    hseq = load_seq(a.h37rv) if a.h37rv else ""
    rseq = load_seq(a.ref_fasta) if a.ref_fasta else ""

    def on_strand(pos_h, r_pos, ref, alt, strand, on_path):
        """Alleles and anchor position as they read on H37Rv's strand.

        Returns (pos_h, ref, alt, status).

        A SNP complements in place. An indel cannot: the VCF anchor base is the
        base BEFORE the change, and "before" is the other side once the strand
        is reversed, so the record has to move as well as complement. Writing
        the called alleles at the called position, which is what this did
        before, leaves the anchor matching H37Rv 0.0% of the time against
        100.0% for the same-strand records.

        For a left-anchored record at R position `r_pos`, with D the deleted
        bases REF[1:] and I the inserted bases ALT[1:], exactly one of which is
        empty: the R interval the record touches maps to the H37Rv interval
        ending at `pos_h`, reversed, so the H37Rv-strand anchor sits at
        pos_h - len(D) - 1 and the changed content is reverse complemented.

        THE ANCHOR BASE COMES FROM DIFFERENT PLACES ON AND OFF THE PATH. On the
        H37Rv path the record will be read against H37Rv, so the anchor must be
        H37Rv's own base or the record is malformed. Off the path H37Rv is not
        the reference the record describes -- its H37Rv column is an anchor, not
        a position -- so the anchor comes from R's own sequence, complemented.
        The two are compared whenever both are available and the disagreement
        rate is reported."""
        if strand != "-":
            return pos_h, ref, alt, "same_strand"
        if len(ref) == 1 and len(alt) == 1:
            return pos_h, ref.translate(comp), alt.translate(comp), "snp"
        if not (ref[0] == alt[0] and (len(ref) == 1 or len(alt) == 1)):
            return pos_h, ref, alt, "unhandled_shape"
        D, I = ref[1:], alt[1:]
        newpos = pos_h - len(D) - 1
        # the upper bound needs the H37Rv sequence; without --h37rv every
        # reverse-strand indel used to fail it (len("") == 0) and was marked
        # off_the_end rather than handled from R's sequence as documented
        if newpos < 1 or (hseq and newpos + len(D) > len(hseq)):
            return pos_h, ref, alt, "off_the_end"
        h_anchor = hseq[newpos - 1] if hseq else ""
        r_anchor = (rseq[r_pos + len(D)].translate(comp)
                    if rseq and r_pos + len(D) < len(rseq) else "")
        anchor = h_anchor if on_path else r_anchor
        if not anchor:
            return pos_h, ref, alt, "no_anchor_sequence"
        st = "reanchored"
        if h_anchor and r_anchor and h_anchor != r_anchor:
            st = "reanchored_anchor_disagrees"
        return newpos, anchor + rc(D), anchor + rc(I), st

    pe_iv, other_iv = load_intervals(a.mask)
    in_pe, in_other = make_hit(pe_iv), make_hit(other_iv)

    def region(p):
        if in_pe(p):
            return "pe_ppe"
        if in_other(p):
            return "masked"
        return "core"

    acc = sorted((int(r["pos"]), r["locus_id"])
                 for r in csv.DictReader(open(a.loci, newline=""), delimiter="\t"))
    acc_pos = [x[0] for x in acc]

    def acc_name(anchor):
        if not acc_pos:
            return ""
        i = bisect.bisect_left(acc_pos, anchor)
        best = ("", 10 ** 9)
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(acc):
                d = abs(acc[j][0] - anchor)
                if d < best[1]:
                    best = (acc[j][1], d)
        return best[0] if best[1] <= a.acc_tol else ""

    matched = load_vcf(a.matched)
    hmap = parse_pos_file(a.h37rv_pos, want_node=False)
    nmap = parse_pos_file(a.node_pos, want_node=True)
    src_positions = {pos for pos, _, _, _ in matched}
    missing = sum(1 for p in src_positions if p not in hmap)
    if src_positions and missing > len(src_positions) * 0.05:
        print(f"  FATAL: {missing} of {len(src_positions)} source positions have "
              f"no projection -- odgi output does not cover the input",
              file=sys.stderr)
        return 1
    hpos = [hmap.get(pos) for pos, _, _, _ in matched]
    npos = [nmap.get(pos) for pos, _, _, _ in matched]
    # ONE KEY PER BASE. odgi counts the node offset along the walking
    # direction, so two references walking a node in opposite directions gave
    # one base two keys (scale200: 13 of 11,163 keyed nodes had records from
    # both directions). Offsets are restated in the node's forward
    # orientation, which needs the node's length where the walk is `-`.
    want_len = {str(v[0]) for v in npos if v is not None and v[2] == "-"}
    nlen = (mtb_norm.load_node_lengths(a.node_lengths, want_len)
            if a.node_lengths and want_len else {})

    rows = []
    n_near = 0
    n_anchor = collections.Counter()
    n_nodestrand = collections.Counter()
    # D41: R's alleles read the node reverse complemented when R walks it `-`
    # in the panel, or when R itself is stored flipped in the panel (22 of 333
    # accessions), but not both. Asked once, and only if R has node rows.
    r_flipped = None
    frames = None
    walk = None

    def locate(q):
        """(node, forward offset, reading strand) of R's 1-based refs
        position q, from R's walk in the build's node table; None where the
        table cannot say. Loaded the first time it is needed."""
        nonlocal walk
        if not a.node_lengths:
            return None
        if walk is None:
            walk = mtb_norm.load_path_nodes(a.node_lengths, a.reference)
        pp = frames.to_panel(a.reference, q - 1) + 1
        i = bisect.bisect_right(walk, (pp, float("inf"))) - 1
        if i < 0:
            return None
        start, ln, nd, st, occ = walk[i]
        if occ != 1 or not (start <= pp < start + ln):
            return None
        fo = mtb_norm.forward_offset(pp - start, st, ln)
        return nd, fo, "-" if (st == "-") != r_flipped else "+"

    # --- direct arm: core sequence only ---------------------------------------
    for pos, ref, alt, qual in load_vcf(a.direct):
        reg = region(pos)
        if reg != "core":
            continue          # the composed arm owns masked sequence
        rows.append(dict(
            sample=a.sample, reference=a.reference, build_id=a.build_id,
            arm="direct", component="called", region=reg, frame="h37rv",
            key=f"h37rv:{pos}", r_pos="", h37rv_pos=pos, dist_to_ref=0, node="",
            node_offset="", frame_strand="+", ref=ref, alt=alt,
            kind=classify(ref, alt),
            size=abs(len(strip_gap(alt)) - len(strip_gap(ref))),
            acc_locus="", qual=qual))

    # --- composed arm, called half: masked sequence and all off-path ----------
    n_off = 0
    n_noproj = collections.Counter()
    for (pos, ref, alt, qual), hp, np_ in zip(matched, hpos, npos):
        r_pos = pos      # the matched reference's own coordinate
        if hp is None or np_ is None:
            # counted, not silently skipped: the 5% guard above covers only
            # the H37Rv projection, and a record with no NODE projection was
            # dropped here without appearing in any tally
            n_noproj["no H37Rv projection" if hp is None else "no node projection"] += 1
            continue
        h, dist, strand = hp[0], hp[1], hp[2] or "+"
        node = np_[0]
        off = mtb_norm.forward_offset(np_[1], np_[2], nlen.get(str(node)))
        if off is None and dist != 0:
            # a node key in the walking direction would be a second key for
            # the same base, so it is not written
            n_noproj["no node length for a reverse-walked node"] += 1
            continue
        if off is None:
            off = ""          # on-path: the column is informational; unknown
        r_ref, r_alt = ref, alt          # as R reads them, for the node frame
        h, ref, alt, anchor_status = on_strand(h, r_pos, ref, alt, strand,
                                               dist == 0)
        n_anchor[anchor_status] += 1
        if dist != 0:
            # "Off the H37Rv path" is not the same as "accessory". dist.to.ref
            # fires for the interior of any insertion in R relative to H37Rv,
            # including a 1 bp one, and section 11 measured 96.3% of off-path
            # calls within 50 bp of the path -- small-indel interiors, which sit
            # inside the same H37Rv ORF as their anchor and are annotatable from
            # H37Rv. Only the remainder is accessory sequence needing the matched
            # reference's own GFF. One label for both would send 96% of these
            # records down the wrong annotation route.
            n_off += 1
            if dist <= a.near_tol:
                reg_off = "off_path_near"
                n_near += 1
                nm = ""
            else:
                reg_off = "off_path_accessory"
                nm = acc_name(h)
            # ONE KEY PER EVENT, NOT ONE PER WALK DIRECTION (D41). The
            # alleles were H37Rv-strand (on_strand above), which says nothing
            # about the node: references walking the node in opposite
            # directions wrote one event as C>T and G>A, two keys. They are
            # restated from R's own alleles on the node's forward strand.
            if r_flipped is None:
                frames = graph_frame.Frames(a.frames or None)
                r_flipped = frames.flipped(a.reference)
            rs = "-" if (np_[2] == "-") != r_flipped else "+"
            node, off, ref, alt, nst = mtb_norm.node_forward_restate(
                node, off, rs, r_ref, r_alt, rseq, r_pos, locate)
            n_nodestrand[nst] += 1
            rows.append(dict(
                sample=a.sample, reference=a.reference, build_id=a.build_id,
                arm="composed", component="called", region=reg_off,
                frame="node", key=f"node:{node}:{off}", r_pos=r_pos, h37rv_pos=h,
                dist_to_ref=dist, node=node, node_offset=off,
                frame_strand=strand, ref=ref, alt=alt, kind=classify(ref, alt),
                size=abs(len(strip_gap(alt)) - len(strip_gap(ref))),
                acc_locus=nm, qual=qual))
            continue
        reg = region(h)
        if reg == "core":
            continue          # the direct arm owns core sequence
        rows.append(dict(
            sample=a.sample, reference=a.reference, build_id=a.build_id,
            arm="composed", component="called", region=reg, frame="h37rv",
            key=f"h37rv:{h}", r_pos=r_pos, h37rv_pos=h, dist_to_ref=0, node=node,
            node_offset=off, frame_strand=strand, ref=ref, alt=alt,
            kind=classify(ref, alt),
            size=abs(len(strip_gap(alt)) - len(strip_gap(ref))),
            acc_locus="", qual=qual))

    # --- composed arm, inherited half: R's own differences, masked only ------
    # This is the half T16 showed carries 89.2% of composition's false positives.
    # It is kept because dropping it loses the recall that makes composition worth
    # using in PE/PPE at all, and it is labelled so downstream can weigh it.
    # SUPPRESSED WHERE A CALLED RECORD OVERLAPS IT, not only where one sits at
    # the same position (audit P4P5-8). The sample's own call replaced R's
    # bases over its REF span; an inherited allele inside that span is R's
    # sequence the sample no longer has. Equality of positions kept 782 such
    # records in scale200, e.g. a called TGGG>T at 976,895 beside an
    # inherited TTG>GGG at 976,896.
    called_spans = sorted(
        (int(r["h37rv_pos"]), int(r["h37rv_pos"]) + max(len(r["ref"]), 1) - 1)
        for r in rows if r["frame"] == "h37rv")
    c_starts = [s for s, _ in called_spans]
    c_reach, m = [], 0
    for _, e in called_spans:
        m = max(m, e); c_reach.append(m)

    def overlaps_called(s, e):
        j = bisect.bisect_right(c_starts, e) - 1
        while j >= 0 and c_reach[j] >= s:
            if called_spans[j][1] >= s:
                return True
            j -= 1
        return False
    n_inh = n_inh_over = 0
    if not os.path.exists(a.graph_vcf):
        sys.exit(f"FATAL: --graph-vcf {a.graph_vcf} does not exist; the "
                 f"inherited half would be silently empty")
    hdr = None
    for line in op(a.graph_vcf):
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t"); break
    if hdr and a.reference in hdr[9:]:
        col = 9 + hdr[9:].index(a.reference)
        for line in op(a.graph_vcf):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= col:
                continue
            g = f[col].split(":", 1)[0]
            if not g.isdigit() or g == "0":
                continue
            pos, ref = int(f[1]), f[3].upper()
            alts = f[4].upper().split(",")
            ai = int(g)
            if ai > len(alts):
                continue
            alt = alts[ai - 1]
            if set(ref + alt) - set("ACGTN"):
                continue
            if region(pos) == "core":
                continue
            if overlaps_called(pos, pos + max(len(ref), 1) - 1):
                n_inh_over += 1
                continue
            n_inh += 1
            rows.append(dict(
                sample=a.sample, reference=a.reference, build_id=a.build_id,
                arm="composed", component="inherited", region=region(pos),
                frame="h37rv", key=f"h37rv:{pos}", r_pos="", h37rv_pos=pos,
                dist_to_ref=0, node="", node_offset="",
                frame_strand="+", ref=ref, alt=alt,
                kind=classify(ref, alt),
                size=abs(len(strip_gap(alt)) - len(strip_gap(ref))),
                acc_locus="", qual="."))
    else:
        print(f"  WARNING: {a.reference} has no column in {a.graph_vcf}; "
              f"the inherited half of the composed arm is EMPTY for this "
              f"sample, which understates its recall", file=sys.stderr)

    if not rows:
        print("  no records placed -- refusing to emit an empty table",
              file=sys.stderr)
        return 1
    # to a temporary name and renamed: p4_place.sh skips any sample whose
    # table is non-empty, so a task killed mid-write must not leave one
    with open(a.out + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    os.replace(a.out + ".tmp", a.out)

    c = collections.Counter((r["arm"], r["component"], r["region"]) for r in rows)
    print(f"  {len(rows)} records placed")
    for why, n in sorted(n_noproj.items()):
        print(f"  matched-arm records dropped, {why}: {n}")
    for k in sorted(c):
        print(f"    {k[0]:<9s}{k[1]:<10s}{k[2]:<9s}{c[k]:>7d}")
    if n_nodestrand:
        print("  node-frame alleles on the node's forward strand (D41): "
              + ", ".join(f"{k} {v}" for k, v in sorted(n_nodestrand.items())))
    n_acc = n_off - n_near
    # count names among accessory-scale records only: counting them across all
    # off-path rows produced "11 of 5 carry a name", which is the kind of ratio
    # that should never reach a report
    named = sum(1 for r in rows
                if r["region"] == "off_path_accessory" and r["acc_locus"])
    print(f"  off-path: {n_off} records -- {n_near} small-indel interiors "
          f"(<= {a.near_tol} bp from the path), {n_acc} accessory-scale")
    print(f"    of the {n_acc} accessory-scale, {named} carry an accessory "
          f"locus name within {a.acc_tol} bp")
    rev = [r for r in rows if r["frame_strand"] == "-"]
    if rev:
        nre = sum(1 for r in rev if r["kind"] != "SNP")
        print(f"  opposite strand: {len(rev)} records project onto H37Rv's other "
              f"strand, because the panel stores {a.reference} reverse "
              f"complemented or the step is locally inverted; {nre} of them "
              f"are indels")
        for k in sorted(n_anchor):
            if k != "same_strand":
                print(f"    {k}: {n_anchor[k]}")
    print(f"  inherited: {n_inh} records from {a.reference}'s own differences; "
          f"{n_inh_over} dropped where a called record overlaps them")
    print(f"  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
