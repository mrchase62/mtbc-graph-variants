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
    pe, other = [], []
    for line in open(bed):
        f = line.split()
        if len(f) < 3:
            continue
        s, e = int(f[1]), int(f[2])
        name = f[3].upper() if len(f) > 3 else ""
        (pe if ("PE" in name or "PPE" in name) else other).append((s, e))
    pe.sort(); other.sort()
    return pe, other


def make_hit(iv):
    starts = [x[0] for x in iv]
    def hit(p):
        i = bisect.bisect_right(starts, p) - 1
        return i >= 0 and iv[i][0] <= p < iv[i][1]
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
                nid, off, _ = f[1].split(",")
                val = (int(nid), int(off), None)
            else:
                tgt = f[1].rsplit(",", 2)
                # f[3] is strand.vs.ref as REPLACED by graphframe/bin/
                # frame_convert.py: the relation between the two refs
                # sequences, which is `-` only when the panel stores one of the
                # two accessions reverse complemented. odgi's own flag is
                # carried in f[4] and is deliberately not used here -- it
                # already accounts for locally inverted steps, and folding it
                # in again complements the allele at exactly those sites.
                val = (int(tgt[-2]) + 1, int(f[2]) if len(f) > 2 else 0,
                       f[3] if len(f) > 3 and f[3] in ("+", "-") else "+")
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
    ap.add_argument("--graph-vcf",
                    default="graphs/CX333.s10k.k23.K15/all_variants.nolab.vcf.gz",
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
        if newpos < 1 or newpos + len(D) > len(hseq):
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

    rows = []
    n_near = 0
    n_anchor = collections.Counter()

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
    for (pos, ref, alt, qual), hp, np_ in zip(matched, hpos, npos):
        r_pos = pos      # the matched reference's own coordinate
        if hp is None or np_ is None:
            continue
        h, dist, strand = hp[0], hp[1], hp[2] or "+"
        node, off = np_[0], np_[1]
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
    called_at = {r["h37rv_pos"] for r in rows if r["frame"] == "h37rv"}
    n_inh = 0
    if os.path.exists(a.graph_vcf):
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
                if region(pos) == "core" or pos in called_at:
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
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    c = collections.Counter((r["arm"], r["component"], r["region"]) for r in rows)
    print(f"  {len(rows)} records placed")
    for k in sorted(c):
        print(f"    {k[0]:<9s}{k[1]:<10s}{k[2]:<9s}{c[k]:>7d}")
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
              f"strand because the panel stores {a.reference} reverse "
              f"complemented; {nre} of them are indels")
        for k in sorted(n_anchor):
            if k != "same_strand":
                print(f"    {k}: {n_anchor[k]}")
    print(f"  inherited: {n_inh} records from {a.reference}'s own differences")
    print(f"  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
