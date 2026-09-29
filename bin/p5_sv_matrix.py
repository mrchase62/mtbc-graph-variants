#!/usr/bin/env python3
"""Cluster SV events across samples into a matrix, carrying their confidence.

Kept as a SEPARATE table from the small-variant matrix, joinable by sample,
because three things about SV rows are genuinely different and merging them into
one table would force each to carry the other's empty columns:

  keys        small variants key on an exact normalised allele. SVs cannot: two
              samples' breakpoints for one event differ by tens of bp, so rows
              come from tolerance clustering -- position within
              max(200 bp, 20% of length), length within 50%, type agreeing,
              which is this project's own convention from its scoring work.
  states      **REF is deliberately not asserted.** For a small variant, covered
              and uncalled means reference. For an SV it does not: stage 6
              measured delly's sensitivity at 0.29 for 50-499 bp deletions and
              0.01-0.05 for insertions, so "no call" is mostly "not detected".
              Cells are ALT, NOCALL or ABSENT only.
  confidence  a raw union callset is mostly false here. Stage 6 on the matched
              arm: deletions 50-499 bp PPV 0.889, 500-4999 bp 0.706, 5 kb+ 0.099;
              insertions 0.02-0.45. Stage 11: delly's bottom three QUAL quartiles
              are PPV <= 0.17 against 0.814 for the top. Every row therefore
              carries type, length band, caller agreement and a QUAL band, and
              the recommended filter is stated rather than applied.

Only 8.9% of the pilot's 7,837 events were found by both callers, so caller
agreement is the single most discriminating column here, not a footnote.

WHY DUP IS DROPPED BY DEFAULT, AND ON A DIFFERENT ARGUMENT FROM INV
Dropped 2026-09-21. An earlier note here called DUP "the same argument as INV
one step weaker", and measuring it showed that was wrong. By the internal
measures DUP is not worse than the classes that are kept:

    class   n     singletons   two-caller   QUAL high   recurrent >=5
    DUP     320      71.6%        19.4%          3            10
    DEL    2622      73.5%        24.0%         30           200
    INS    1461      69.0%         2.3%         30           205

Its singleton rate is LOWER than DEL's and its caller agreement HIGHER than
INS's. So the INV argument -- that the class looks like artefact on its own
statistics -- does not transfer.

The actual case against DUP is that it is entirely unmeasured:

  * truth from `paftools call` reports DEL and INS, so DUP has never been
    scored against a finished assembly and has no sensitivity or PPV at all;
  * it is the ONLY class with no composed rows -- 320 of 320 are `called` --
    so it gets none of the sensitivity the composed half carries, which is
    what makes the rest of this table usable;
  * 10 of 320 recur in five or more isolates, against 200 and 205.

That is a reason to keep it out of a table whose recommended filters are all
backed by a measurement, not a claim that the calls are false. `--drop-svtype
INV` keeps them if wanted.

WHY INV IS DROPPED BY DEFAULT
Dropped 2026-09-21 on the evidence in `refbias/SV_ARM_LIMIT.md`. Inversions
were 6,708 of 11,111 rows at cohort scale -- 65% of the table -- and every
internal measure said they were artefact: 88% singletons, only 10 recurring in
five or more isolates, 11.2% with two-caller support, 0.2% reaching QUAL band
high, and 70% of them 50-499 bp.

The decisive test is that the callers do not describe the one inversion this
project has independently confirmed. The lineage-2.2.1 / W148 3 Mb nested
rearrangement has four verified breakpoints, and all four have an INV call
within 5 kb -- which is not detection, because the INV callset is one event per
658 bp and a random 5 kb window holds 15.3 of them. The calls near the
breakpoints are 75 to 512 bp against a real event of 3 Mb.

So nothing real is lost: the one confirmed inversion was not represented in
the table before this change either. INV is also unscorable here, since truth
from `paftools call` reports DEL and INS only, and the composed half
contributes no inversions at all.

`--drop-svtype ''` puts them back. The capability is kept rather than deleted
because a caller that can describe a 3 Mb rearrangement would change this, and
DUP -- 320 rows, also unscored, also 72% singletons -- is the same argument
one step weaker and is deliberately NOT dropped, so the choice stays visible.

WHY --keep-filter DEFAULTS TO PASS,INHERITED AND NOT PASS
Changed 2026-09-21. P4b has two halves: rows a caller produced
(`component=called`, carrying a caller FILTER and QUAL) and rows composed from
the isolate's reference through the graph (`component=inherited`,
`FILTER=INHERITED`, no caller QUAL at all). Scoring against the seven finished
assemblies settles which matters -- see
`truth7/docs/SV_SCORED_AGAINST_ASSEMBLIES.md`:

    called half only     DEL sensitivity 0.041   INS sensitivity 0.000
    with composed half   DEL sensitivity 0.622   INS sensitivity 0.557

A default of PASS therefore returned a table that found 9 of 222 true
deletions and none of 79 true insertions. The composed half is not an extra;
it is the arm's sensitivity.

BUT IT IS LABELLED, NOT BLENDED, WHICH IS WHY THIS IS NOT A ONE-LINE CHANGE
A composed row has no caller QUAL, so banding it on a caller's calibration
would attribute a confidence to it that was never measured -- and with a blank
QUAL it would land in the lowest band, making the rows that carry the
sensitivity look like the least reliable. So:

  * `component` is a column: called, inherited, or mixed for a cluster with
    both;
  * an inherited-only cluster gets `qual_band=composed` and a PPV from P4's
    own measurement of its inherited half, 0.631 to 0.823 against assembly
    truth, taken at the conservative end and settable with --composed-ppv;
  * the summary and the recommended filter are reported per component,
    because their precision differs: DEL PPV 0.573 for the union against
    0.118 for the called half alone.
"""
import argparse, collections, csv, os, statistics, sys

# Per-caller QUAL bands, read from stage 11's calibration table rather than
# hardcoded.
#
# An earlier version hardcoded DELLY's quartile boundaries (242 / 719 / 2760) and
# applied them to every row. dysgu's QUAL runs 5 to 27, so every dysgu row fell
# below the lowest delly boundary and was banded `very_low` -- and the conclusion
# drawn from that was that dysgu was uncalibrated and calibrating it was the
# highest-value change available. It was already calibrated: stage 11 scored
# dysgu, lumpy, manta and wham alongside delly. The defect was mine, applying one
# caller's scale to another.
#
# It matters in the direction that flatters dysgu: its calibrated PPV runs
# 0.53 at the bottom quartile to 0.97 at the top, against delly's 0.03 to 0.81.


def load_bands(path):
    """caller -> [(lower bound, band name, measured PPV)], highest first."""
    import csv as _csv
    by = {}
    try:
        for r in _csv.DictReader(open(path, newline=""), delimiter="\t"):
            if r.get("bin_lo") in (None, "", "ALL"):
                continue
            try:
                lo = float(r["bin_lo"]); ppv = float(r["ppv"])
            except (ValueError, KeyError):
                continue
            by.setdefault(r["caller"], []).append((lo, ppv))
    except OSError:
        return {}
    out = {}
    for caller, rows in by.items():
        rows.sort(reverse=True)
        names = ["high", "medium", "low", "very_low"]
        out[caller] = [(lo, names[min(i, 3)], ppv)
                       for i, (lo, ppv) in enumerate(rows)]
    return out


def qual_band(q, caller, bands):
    tab = bands.get(caller)
    if not tab:
        return "unknown", ""
    try:
        v = float(q)
    except (TypeError, ValueError):
        return "unknown", ""
    for lo, name, ppv in tab:
        if v >= lo:
            return name, ppv
    return "very_low", tab[-1][2] if tab else ""


def size_band(n):
    if n < 50:
        return "<50"
    if n < 500:
        return "50-499"
    if n < 5000:
        return "500-4999"
    return "5k+"


# The tolerance is CAPPED. This project's SV convention is
# max(200 bp, 20% of length), which is scale-free and was written for scoring
# events of plausible size. Applied to this callset it fails: the 90th percentile
# of placed events is 339 kb and the largest is 4,289,795 bp in a 4.4 Mb genome,
# so one such "event" opens a half-megabase merge window and absorbs every real
# event whose length happens to fall within 50% of its own. Uncapped, 7,569
# events clustered into 22 rows.
TOL_CAP = 5000


def tol(length):
    return min(max(200, int(0.2 * max(length, 1))), TOL_CAP)


def same_event(a, b):
    if a["svtype"] != b["svtype"]:
        return False
    if abs(a["pos"] - b["pos"]) > tol(max(a["svlen"], b["svlen"])):
        return False
    lo, hi = sorted((max(a["svlen"], 1), max(b["svlen"], 1)))
    return hi / lo <= 1.5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--dir", default="refbias/p4b")
    ap.add_argument("--max-svlen", type=int, default=100000,
                    help="events longer than this are mapping artefacts in a "
                         "4.4 Mb clonal genome, not structural variants")
    ap.add_argument("--keep-filter", default="PASS,INHERITED",
                    help="comma-separated FILTER values to keep, or 'any'")
    ap.add_argument("--drop-svtype", default="INV,DUP",
                    help="comma-separated SVTYPEs to exclude. INV and DUP by "
                         "default, on different grounds: INV because it looks "
                         "like artefact on its own statistics and the callers "
                         "do not describe the one confirmed inversion; DUP "
                         "because it is unmeasured -- no truth class, no "
                         "composed rows, and 10 of 320 recurrent. Pass an "
                         "empty string to keep both.")
    ap.add_argument("--composed-ppv", type=float, default=0.631,
                    help="PPV assigned to an inherited-only cluster. P4 "
                         "measured its inherited half at 0.631-0.823 against "
                         "assembly truth depending on how far the matched "
                         "reference sits from the isolate; the conservative "
                         "end is the default. It is NOT a caller calibration "
                         "and is deliberately not read from --calibration.")
    ap.add_argument("--calibration", default="refbias/STAGE11.calibration.tsv")
    ap.add_argument("--out", default="refbias/p5/sv_matrix.tsv")
    a = ap.parse_args()

    bands = load_bands(a.calibration)
    if not bands:
        print(f"  WARNING: no calibration at {a.calibration}; every band will "
              f"read 'unknown'", file=sys.stderr)
    else:
        print("  QUAL bands per caller, from stage 11:")
        for c, tab in sorted(bands.items()):
            print(f"    {c:<7s}" + "  ".join(f"{n}>={lo:g}(ppv {p:.2f})"
                                             for lo, n, p in tab))
    keep = None if a.keep_filter == "any" else set(a.keep_filter.split(","))
    drop_ty = {x.strip().upper() for x in a.drop_svtype.split(",") if x.strip()}
    n_filt_by = collections.Counter()
    n_dropped_ty = collections.Counter()
    n_filt = 0
    samples, events = [], []
    for r in csv.DictReader(open(a.refmap, newline=""), delimiter="\t"):
        p = os.path.join(a.dir, f"{r['sample']}.sv_placed.tsv")
        if not os.path.exists(p):
            continue
        samples.append(r["sample"])
        for x in csv.DictReader(open(p, newline=""), delimiter="\t"):
            if x["frame"] != "h37rv" or not x["h37rv_pos"]:
                continue          # breakend pairs are not clusterable by position
            if x["svtype"].upper() in drop_ty:
                n_dropped_ty[x["svtype"]] += 1
                continue
            if keep and x.get("filter") not in keep:
                n_filt += 1
                # x is the SV row; r is the refmap row for the sample. An
                # earlier version of this line read r and so reported every
                # excluded row as blank.
                n_filt_by[x.get("filter", "") or "(blank)"] += 1
                continue
            try:
                events.append(dict(sample=r["sample"], svtype=x["svtype"],
                                   pos=int(x["h37rv_pos"]),
                                   end=int(x["h37rv_end"] or x["h37rv_pos"]),
                                   svlen=int(x["svlen"]), src=x["src"],
                                   n_callers=int(x["n_callers"]),
                                   qual=x["qual"], sr=x["sr"], pe=x["pe"],
                                   component=x.get("component", "called")))
            except ValueError:
                continue
    n_big = sum(1 for e in events if e["svlen"] > a.max_svlen)
    events = [e for e in events if e["svlen"] <= a.max_svlen]
    if not events:
        print("no placeable SV events found", file=sys.stderr); return 1
    print(f"  {len(events)} events kept over {len(samples)} samples")
    if n_dropped_ty:
        tot_d = sum(n_dropped_ty.values())
        print(f"    dropped {tot_d} on SVTYPE (--drop-svtype "
              f"{a.drop_svtype}): "
              + ", ".join(f"{k} {v}" for k, v in n_dropped_ty.most_common()))
        print(f"      not a quality filter and not scored away -- see the note "
              f"at the top of this script and refbias/SV_ARM_LIMIT.md.")
    print(f"    excluded {n_filt} on FILTER (kept: {a.keep_filter}) and "
          f"{n_big} longer than {a.max_svlen} bp")
    # Break the FILTER exclusion down, because two very different things hide
    # in one number. A row failing a caller's quality filter is a weak call.
    # A row with FILTER=INHERITED is P4b's COMPOSED half, added 2026-09-20:
    # it carries no caller FILTER of its own and is excluded here by the
    # default --keep-filter PASS, which is a design decision rather than a
    # quality judgement. P4b's own header says "labelled, not blended", and
    # nothing in this script reads `component` yet, so leaving the exclusion
    # invisible would hide an open question. See refbias/archive/README.md.
    if n_filt_by:
        for k, v in sorted(n_filt_by.items(), key=lambda kv: -kv[1]):
            tag = ("  <- P4b's composed half, excluded by design, not by "
                   "quality" if k == "INHERITED" else "")
            print(f"      {k:<12s} {v:6d}{tag}")

    # cluster across samples on the project's own tolerance
    # A for/else here is a trap and an earlier version fell into it: `break` was
    # used both to stop scanning when a cluster matched AND to stop scanning
    # when the remaining clusters were too far away, so the `else` never ran on
    # the second case and every event that needed a NEW cluster was silently
    # dropped. 7,569 events became 22 rows. An explicit flag says what is meant.
    events.sort(key=lambda e: (e["svtype"], e["pos"]))
    by_type = collections.defaultdict(list)
    for e in events:
        by_type[e["svtype"]].append(e)
    clusters = []
    for svtype, evs in by_type.items():
        open_clusters = []
        for e in evs:
            placed = False
            # only clusters whose representative is still within reach need
            # checking, and the list is walked newest-first because events
            # arrive in position order
            still_open = []
            for c in open_clusters:
                if e["pos"] - c["rep"]["pos"] <= TOL_CAP + 1000:
                    still_open.append(c)
            open_clusters = still_open
            for c in reversed(open_clusters):
                if same_event(c["rep"], e):
                    c["members"].append(e)
                    placed = True
                    break
            if not placed:
                c = dict(rep=e, members=[e])
                clusters.append(c)
                open_clusters.append(c)
    print(f"  clustered into {len(clusters)} distinct events")

    rows = []
    for i, c in enumerate(clusters):
        by_s = collections.defaultdict(list)
        for m in c["members"]:
            by_s[m["sample"]].append(m)
        pos = int(statistics.median(m["pos"] for m in c["members"]))
        ln = int(statistics.median(m["svlen"] for m in c["members"]))
        scored = [m for m in c["members"]
                  if m["qual"].replace(".", "", 1).isdigit()]
        best = max(scored, key=lambda m: float(m["qual"])) if scored else None
        best_q = float(best["qual"]) if best else 0.0
        qual_src = best["src"] if best else ""
        # band against the caller that produced this QUAL, not against delly
        primary = "delly" if "delly" in qual_src else (
            "dysgu" if "dysgu" in qual_src else qual_src.split(",")[0])
        comps = {m.get("component", "called") for m in c["members"]}
        component = ("mixed" if len(comps) > 1 else
                     (comps.pop() if comps else "called"))
        if best is None and component == "inherited":
            # no caller QUAL exists for a composed row; banding it on a
            # caller's scale would invent a confidence, and a blank QUAL would
            # put the rows carrying the sensitivity in the lowest band
            band, ppv = "composed", a.composed_ppv
            qual_src = "graph_vcf"
        else:
            band, ppv = qual_band(best_q, primary, bands)
        both = sum(1 for s in by_s if any(m["n_callers"] > 1 for m in by_s[s]))
        cells = []
        for s in samples:
            cells.append("ALT" if s in by_s else "NOCALL")
        rows.append(dict(
            key=f"sv:{c['rep']['svtype']}:{pos}:{ln}", svtype=c["rep"]["svtype"],
            h37rv_pos=pos, svlen=ln, size_band=size_band(ln),
            n_alt=len(by_s), n_both_callers=both,
            src=",".join(sorted({x for m in c["members"] for x in m["src"].split(",")})),
            max_qual=round(best_q, 1) if best else "",
            qual_band=band, stage11_ppv=ppv,
            qual_from=qual_src, component=component,
            cells=cells))

    # KEYS MUST BE UNIQUE, and clustering does not guarantee it. The key is
    # built from a cluster's MEDIAN position and MEDIAN length, so two
    # clusters the tolerance kept apart can still land on the same pair and
    # produce the same string. gwas1000 had exactly one such collision out of
    # 23,910 rows -- two INS clusters at 802491, both of median length 108,
    # carried by 35 and by 567 isolates. The merged VCF then held two records
    # with the same ID, and write_event_matrix.py refused it, correctly, since
    # its index is a dict.
    #
    # The colliding rows are NOT merged. Their carrier sets differ by a factor
    # of 16, so folding them together would assert that one event is the other
    # and discard whichever came second. They are made distinct instead, by
    # appending an occurrence number, and the collision is printed so that a
    # clustering which splits one event stays visible rather than being
    # papered over here.
    seen, dupes = {}, []
    for r in rows:
        k = r["key"]
        n = seen.get(k, 0) + 1
        seen[k] = n
        if n > 1:
            r["key"] = f"{k}#{n}"
            dupes.append((r["key"], r["n_alt"]))
    if dupes:
        print(f"\n  {len(dupes)} key collision(s) made distinct with a #n "
              f"suffix; two clusters shared a median position and length:")
        for k, na in dupes[:10]:
            print(f"    {k}  n_alt={na}")

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["key", "svtype", "h37rv_pos", "svlen", "size_band", "n_alt",
                    "n_both_callers", "src", "max_qual", "qual_band",
                    "stage11_ppv", "qual_from", "component"] + samples)
        for r in rows:
            w.writerow([r["key"], r["svtype"], r["h37rv_pos"], r["svlen"],
                        r["size_band"], r["n_alt"], r["n_both_callers"],
                        r["src"], r["max_qual"], r["qual_band"],
                        r["stage11_ppv"], r["qual_from"],
                        r["component"]] + r["cells"])

    n = len(rows)
    print(f"\n  {n} SV rows x {len(samples)} samples")
    for k, lab in (("svtype", "by type"), ("size_band", "by size"),
                   ("qual_band", "by QUAL band"), ("component", "by component")):
        c = collections.Counter(r[k] for r in rows)
        print(f"    {lab}: " + "  ".join(f"{a_}={b_}" for a_, b_ in c.most_common()))
    sing = sum(1 for r in rows if r["n_alt"] == 1)
    ubiq = sum(1 for r in rows if r["n_alt"] == len(samples))
    both = sum(1 for r in rows if r["n_both_callers"] > 0)
    print(f"    {sing} singletons, {ubiq} carried by all {len(samples)}, "
          f"{both} with two-caller support in at least one sample")
    c = collections.Counter(r["qual_from"] for r in rows)
    print(f"    QUAL source: " + "  ".join(f"{k or 'none'}={v}"
                                           for k, v in c.most_common()))
    # TWO filters, not one, because the halves have different calibrations and
    # a single QUAL-band rule silently drops the whole composed half -- it has
    # no caller QUAL, so it can never be band 'high'. Stating one filter here
    # was what made the composed rows invisible in the first place.
    sized = [r for r in rows if r["size_band"] in ("50-499", "500-4999")]
    prec = [r for r in sized if r["svtype"] == "DEL"
            and r["component"] in ("called", "mixed")
            and r["qual_band"] == "high"]
    sens = [r for r in sized if r["svtype"] in ("DEL", "INS")
            and (r["component"] in ("inherited", "mixed")
                 or r["qual_band"] == "high")]
    print(f"\n  recommended filters, from this project's own measurements.")
    print(f"  Two, because the halves are calibrated differently and one")
    print(f"  QUAL-band rule would discard every composed row:")
    print(f"\n    precision-first:  DEL, 50-4999 bp, called half, QUAL 'high'")
    print(f"                      ->  {len(prec)} of {n} rows")
    print(f"      stage 6 put that class at PPV 0.706-0.889 and stage 11 at")
    print(f"      0.814. It is the most trustworthy set and it is very small.")
    print(f"\n    sensitivity-first: DEL or INS, 50-4999 bp, composed half OR")
    print(f"                      called half at QUAL 'high'")
    print(f"                      ->  {len(sens)} of {n} rows")
    print(f"      measured against the seven finished assemblies at DEL")
    print(f"      sensitivity 0.622 / PPV 0.573 and INS 0.557 / 0.288.")
    kept_ty = sorted({r["svtype"] for r in rows})
    print(f"\n    Types in the table: {', '.join(kept_ty)}. Both filters above")
    print(f"    are backed by a measurement against the seven finished")
    print(f"    assemblies, which is why the unmeasured classes are dropped")
    print(f"    rather than kept and caveated -- see --drop-svtype.")
    called = [r for r in rows if r["component"] == "called"]
    inh = [r for r in rows if r["component"] == "inherited"]
    mix = [r for r in rows if r["component"] == "mixed"]
    print(f"\n  the two halves, and why the default keeps both:")
    print(f"    called    {len(called):6d} rows   caller FILTER and QUAL, "
          f"banded on stage 11")
    print(f"    inherited {len(inh):6d} rows   composed through the graph, no "
          f"caller QUAL, PPV {a.composed_ppv} from P4's own measurement")
    print(f"    mixed     {len(mix):6d} rows   a cluster containing both")
    print(f"    Against the seven finished assemblies the called half alone "
          f"gives DEL\n    sensitivity 0.041 and INS 0.000; with the composed "
          f"half, 0.622 and 0.557.\n    See truth7/docs/"
          f"SV_SCORED_AGAINST_ASSEMBLIES.md. Precision is the trade:\n"
          f"    DEL PPV 0.573 for the union against 0.118 for the called half.")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
