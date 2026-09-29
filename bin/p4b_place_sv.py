#!/usr/bin/env python3
"""P4b: merge one isolate's SV calls across callers and place them in the frame.

P2 has been calling structural variants since it was written -- 5,966 delly and
5,479 dysgu records over the pilot -- and nothing consumed them: P4 reads only the
small-variant VCF, so the matrix contains zero SV records. This closes that.

Two things make SVs unlike small variants here.

**An SV has two breakpoints and they can behave differently.** A deletion whose
start projects cleanly onto H37Rv and whose end lands off-path is not one
liftable record; placing only the start would put a 1.3 kb deletion at a
coordinate whose partner is fictional. Stage 9 settled the principle for
reachability -- an SV counts only if BOTH breakpoints are usable -- and the same
rule applies with projection substituted for callability. A record failing it is
kept as a node-keyed breakend pair rather than discarded, which is what VCF
already provides for the case.

**Two callers are not two observations of equal weight.** Stage 6 measured union
and intersection behaving differently, so which callers found a record is part of
the record. Calls are merged within a sample on the project's own tolerance --
position within max(200 bp, 20% of length), length within 50%, type must agree --
and `SRC` records the callers that agreed.
"""
import argparse, bisect, collections, csv, gzip, os, subprocess, sys


def op(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_sv_vcf(path, caller):
    out = []
    if not path or not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        info = dict((kv.split("=", 1) + [""])[:2] for kv in f[7].split(";") if kv)
        svtype = info.get("SVTYPE", "")
        if not svtype:
            continue
        try:
            pos = int(f[1])
        except ValueError:
            continue
        try:
            end = int(info.get("END", pos))
        except ValueError:
            end = pos
        svlen = info.get("SVLEN", "")
        try:
            svlen = abs(int(svlen.split(",")[0])) if svlen else abs(end - pos)
        except ValueError:
            svlen = abs(end - pos)
        out.append(dict(caller=caller, chrom=f[0], pos=pos, end=end,
                        svtype=svtype, svlen=svlen, qual=f[5], filt=f[6],
                        pe=info.get("PE", ""), sr=info.get("SR", ""),
                        cipos=info.get("CIPOS", ""), ciend=info.get("CIEND", ""),
                        alt=f[4]))
    return out


def tol(length):
    """The project's own SV matching tolerance."""
    return max(200, int(0.2 * max(length, 1)))


def same_event(a, b):
    if a["svtype"] != b["svtype"]:
        return False
    if abs(a["pos"] - b["pos"]) > tol(max(a["svlen"], b["svlen"])):
        return False
    lo, hi = sorted((max(a["svlen"], 1), max(b["svlen"], 1)))
    return hi / lo <= 1.5


def merge_callers(records):
    """Cluster records from all callers into events, recording which agreed."""
    records = sorted(records, key=lambda r: (r["svtype"], r["pos"]))
    events = []
    for r in records:
        for e in events:
            if same_event(e["rep"], r):
                e["members"].append(r)
                # keep the call with the most split-read support as the
                # representative, since that is the one whose breakpoints are
                # best resolved
                try:
                    if int(r["sr"] or 0) > int(e["rep"]["sr"] or 0):
                        e["rep"] = r
                except ValueError:
                    pass
                break
        else:
            events.append(dict(rep=r, members=[r]))
    return events


FIELDS = ["sample", "reference", "build_id", "svtype", "svlen", "frame",
          "key", "h37rv_pos", "h37rv_end", "r_pos", "r_end", "node_pos",
          "node_end", "src", "n_callers", "sr", "pe", "qual", "filter",
          "component"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--build-id", required=True)
    ap.add_argument("--delly")
    ap.add_argument("--dysgu")
    ap.add_argument("--positions", required=True,
                    help="odgi output for every breakpoint, R -> H37Rv")
    ap.add_argument("--mask", required=True)
    ap.add_argument("--graph-vcf", default="",
                    help="source of the INHERITED half, as p4_place.py uses "
                         "it. Required unless --no-inherited: a relative "
                         "default that resolved only from one working "
                         "directory used to skip the half without a word")
    ap.add_argument("--min-sv", type=int, default=50)
    ap.add_argument("--no-inherited", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if not a.no_inherited and not (a.graph_vcf and os.path.exists(a.graph_vcf)):
        sys.exit(f"FATAL: --graph-vcf {a.graph_vcf!r} does not exist; pass it, "
                 f"or --no-inherited to skip the inherited half deliberately")

    recs = parse_sv_vcf(a.delly, "delly") + parse_sv_vcf(a.dysgu, "dysgu")
    # No early return on zero calls. It used to write a header-only table here,
    # BEFORE the inherited half ran, so R's own large differences from H37Rv --
    # which need no calls at all -- were dropped for exactly the samples
    # closest to their reference. Zero calls is zero events; carry on.
    if not recs:
        print(f"  {a.sample}: no SV records to place", file=sys.stderr)
    events = merge_callers(recs)

    # projections keyed by SOURCE position: threaded odgi returns results in
    # completion order, never in input order
    proj = {}
    for line in open(a.positions):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            src = int(f[0].rsplit(",", 2)[-2]) + 1
            tgt = int(f[1].rsplit(",", 2)[-2]) + 1
            dist = int(f[2])
        except (ValueError, IndexError):
            continue
        # column 4 is the relation between R's and H37Rv's refs sequences,
        # restated by frame_convert.py from-panel; `-` means R runs reverse to
        # H37Rv here, so an event's breakpoints swap sides (see below)
        strand = f[3].strip() if len(f) > 3 and f[3].strip() in "+-" else "+"
        proj.setdefault(src, (tgt, dist, strand))

    pe_iv = []
    for line in open(a.mask):
        f = line.split()
        if len(f) >= 3:
            pe_iv.append((int(f[1]), int(f[2])))
    pe_iv.sort()

    rows = []
    counts = collections.Counter()
    for e in events:
        r = e["rep"]
        callers = sorted({m["caller"] for m in e["members"]})
        ps, pe_ = proj.get(r["pos"]), proj.get(r["end"])
        both = ps is not None and pe_ is not None and ps[1] == 0 and pe_[1] == 0
        h1 = h2 = None
        if both:
            # DIRECTION FROM GEOMETRY, not from the strand label. Column 4 is
            # the relation between the two whole refs sequences; locally the
            # graph can run either way (measured: a DEL on the flipped
            # GCF_965124535 labelled `-` whose breakpoints project increasing).
            # So when an event has two distinct breakpoints, their projected
            # order says whether R runs reverse to H37Rv here. The label is
            # used only for a point insertion, which has one breakpoint.
            if r["end"] > r["pos"]:
                reverse = pe_[0] < ps[0]
            else:
                reverse = ps[2] == "-"
            if not reverse:
                h1, h2 = ps[0], pe_[0]
            elif r["svtype"] == "INS":
                # An insertion sits between R's POS and POS+1. Reversed, that
                # is between f(POS)-1 and f(POS), so its anchor is f(POS)-1.
                h1 = ps[0] - 1
                h2 = h1 + (r["end"] - r["pos"])
                counts["reverse_strand"] += 1
            else:
                # Symbolic SVs cover R bases POS+1 .. END (POS is the padding
                # base). Reversed onto H37Rv they cover f(END) .. f(POS)-1, so
                # the padding base is f(END)-1 and the last base is f(POS)-1.
                # Using ps/pe as they came gave start > end, shifted by the
                # whole event length, and every depth probe landed outside it.
                h1, h2 = pe_[0] - 1, ps[0] - 1
                counts["reverse_strand"] += 1
            # A placed event must keep roughly its called span. Endpoints that
            # project to distant, unrelated places are not one H37Rv interval.
            if r["svtype"] != "INS":
                span_r, span_h = abs(r["end"] - r["pos"]), abs(h2 - h1)
                if abs(span_h - span_r) > tol(max(span_r, r["svlen"])):
                    both = False
                    counts["span_mismatch"] += 1
        if both:
            frame = "h37rv"
            key = f"sv:{r['svtype']}:{h1}:{h2}"
            counts["placed"] += 1
        else:
            # one or both breakpoints have no H37Rv equivalent: a breakend pair.
            # The key is built from R coordinates and the reference; h1/h2 keep
            # whatever single breakpoint did project, for inspection only.
            frame = "bnd"
            key = f"bnd:{r['svtype']}:{r['pos']}:{r['end']}:{a.reference}"
            h1, h2 = (ps[0] if ps and ps[1] == 0 else ""), \
                     (pe_[0] if pe_ and pe_[1] == 0 else "")
            counts["breakend"] += 1
        counts[f"type_{r['svtype']}"] += 1
        counts[f"src_{'+'.join(callers)}"] += 1
        rows.append(dict(sample=a.sample, reference=a.reference,
                         build_id=a.build_id, svtype=r["svtype"],
                         svlen=r["svlen"], frame=frame, key=key,
                         h37rv_pos=h1, h37rv_end=h2, r_pos=r["pos"],
                         r_end=r["end"], node_pos="", node_end="",
                         src=",".join(callers), n_callers=len(callers),
                         sr=r["sr"], pe=r["pe"], qual=r["qual"],
                         filter=r["filt"], component="called"))

    # --- the inherited half ---------------------------------------------------
    # P4 composes a sample's small variants from TWO halves: its own calls
    # against R, and R's own differences from H37Rv out of the graph VCF. This
    # script had only the first, so a placed SV record stated sample-versus-R
    # at H37Rv coordinates while H37Rv-framed truth states sample-versus-H37Rv.
    # An SV that R and the sample SHARE was in truth and could not be in the
    # output at any caller quality -- measured against the seven truth
    # assemblies at deletion sensitivity 0.04, which is a design gap and not a
    # caller failure.
    #
    # The nesting rule is p4_place.py's: vg deconstruct emits LV=0 parents and
    # LV=1 children for one underlying difference, and summing both counts it
    # twice.
    #
    # LABELLED, NOT BLENDED. P4's own inherited half measured PPV 0.631 to
    # 0.823 against assembly truth depending on how far R sits from the
    # isolate, so it buys recall at a cost a consumer must be able to weigh.
    # `component` says which half a row came from, and `filter` is INHERITED
    # because these carry no caller FILTER of their own.
    n_inh = 0
    if not a.no_inherited:
        called_at = sorted(int(r["h37rv_pos"]) for r in rows
                           if r["frame"] == "h37rv" and r["h37rv_pos"] != "")
        child, hdr = [], None
        for line in op(a.graph_vcf):
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    hdr = line.rstrip("\n").split("\t")
                continue
            f = line.rstrip("\n").split("\t", 8)
            info = dict((kv.split("=", 1) + [""])[:2]
                        for kv in f[7].split(";") if kv)
            if info.get("LV") != "0":
                pos = int(f[1])
                child.append((pos, pos + max(len(f[3]), 1)))
        child.sort()
        cstart = [c[0] for c in child]

        def has_child(pos, end):
            i = bisect.bisect_left(cstart, end) - 1
            while i >= 0 and child[i][1] > pos:
                if child[i][0] < end:
                    return True
                i -= 1
            return False

        col = None
        if hdr and a.reference in hdr[9:]:
            col = 9 + hdr[9:].index(a.reference)
        if col is None:
            print(f"    WARNING: {a.reference} has no column in "
                  f"{a.graph_vcf}; the inherited half is EMPTY for this "
                  f"sample, which understates its recall", file=sys.stderr)
        else:
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
                d = len(alt) - len(ref)
                if abs(d) < a.min_sv:
                    continue
                info = dict((kv.split("=", 1) + [""])[:2]
                            for kv in f[7].split(";") if kv)
                if info.get("LV") == "0" and has_child(pos, pos + max(len(ref), 1)):
                    continue
                i = bisect.bisect_left(called_at, pos - 200)
                if i < len(called_at) and called_at[i] <= pos + 200:
                    continue          # the called half already states this
                svt = "INS" if d > 0 else "DEL"
                end = pos + (0 if d > 0 else abs(d))
                n_inh += 1
                rows.append(dict(
                    sample=a.sample, reference=a.reference,
                    build_id=a.build_id, svtype=svt, svlen=abs(d),
                    frame="h37rv", key=f"sv:{svt}:{pos}:{end}",
                    h37rv_pos=pos, h37rv_end=end, r_pos="", r_end="",
                    node_pos="", node_end="", src="graph_vcf", n_callers=0,
                    sr="", pe="", qual="", filter="INHERITED",
                    component="inherited"))

    # header-only when there is nothing, never an IndexError on rows[0]; an
    # empty table and a missing file must not look the same downstream
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"    inherited: {n_inh} SV records from R's own differences "
          f"from H37Rv")
    print(f"  {a.sample}: {len(recs)} caller records -> {len(events)} events; "
          f"{counts['placed']} placed on the H37Rv path, "
          f"{counts['breakend']} as breakend pairs")
    for k in ("reverse_strand", "span_mismatch"):
        if counts[k]:
            print(f"    {k.replace('_', ' ')}: {counts[k]}")
    agree = sum(1 for r in rows if r["n_callers"] > 1)
    print(f"    both callers agreed on {agree} of {len(rows)} "
          f"({100*agree/max(len(rows), 1):.1f}%)")
    for k in sorted(counts):
        if k.startswith("type_"):
            print(f"    {k[5:]:<8s}{counts[k]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
