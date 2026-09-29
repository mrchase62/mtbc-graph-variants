#!/usr/bin/env python3
"""Site calling on two independent axes: strength of evidence, and site class.

WHY THE LABELS CHANGED.  The old tiers were A for span 2-6 with both element
termini, B for any other two-sided span, and C for one-sided.  That conflates
two things that have nothing to do with each other, and reading A as better
than B was wrong.

  EVIDENCE is a real quality ordering.  A two-sided call sees both flanks; a
  one-sided call sees one.  The one-sided population carries half the
  depth-normalised read support of the two-sided one -- 0.50 and 0.51 in the
  two arms, derived independently -- and is contradicted by an assembly more
  often, 8.8% against 1.7%.  So two_sided ranks above one_sided, and
  is6110/docs/SATURATION_CHECK.md is where that comes from.

  SITE CLASS is not a quality ordering at all.  It records whether the
  reference happens to carry the same insertion, and it was being smuggled in
  as A against B.  Where the reference lacks the insertion, both flanks align
  onto a single target site from opposite directions and the span is the
  duplication length, 2-6 bp.  Where the reference carries it, the build
  removed the element alone and left the target site in duplicate, so the
  flanks land on adjacent bases and the span collapses to 0-1.  A matched
  reference shares most of the sample's insertions and a distant one shares
  few, which is the whole of the TSD difference between the arms.
  is6110/docs/TSD_LOSS_EXPLAINED.md has the measurement.

SITE CLASS IS READ FROM THE CROSSMAP, NOT INFERRED FROM SPAN.  Span is a good
proxy at one end -- 233 of 234 stacks at span 0-1 are at a shared site -- but
it is not diagnostic in the middle, since span 7 and above occurs in both
classes at roughly 20-26%.  So the class is determined by asking whether the
call lies within --class-slop of a coordinate where that isolate's own
reference had an interval removed, and span stays a reported column rather
than becoming a label.

PROMOTION, AND WHY IT IS NOT A PLAIN SUM.  Half-unit support is also what one
insertion split across two one-sided stacks looks like, and 29% of them are
exactly that: within 1 kb of another one-sided stack, reaching through
OPPOSITE termini of the element contig 50 times against 5, with combined
support near one unit.  Counting each stack separately would count those
insertions twice, and unequally between arms.  So flank pairs are collapsed
first: opposite termini within --pair-window, or the same terminus within
--scatter-window.  Binning the gaps puts opposite-terminus enrichment at
essentially 1.00 from 20 bp out to 2 kb with no decay to take an edge from,
while the 0-20 bp bin is 0.00 opposite, which is the same flank scattered over
nearby placements.  The window is a judgement, so --sweep prints it.

Arm-agnostic: point --dir at any p1g or p1i output.
"""
import argparse, collections, csv, glob, os, statistics as st, sys

TWO, ONE = "two_sided", "one_sided"
SHARED, LACKING, UNKNOWN = "ref_shared", "ref_lacking", "unknown"


def junctions(path):
    """Clean coordinates where a build removed one of its intervals."""
    if not path or not os.path.exists(path):
        return None
    return sorted(int(r["clean_junction"])
                  for r in csv.DictReader(open(path), delimiter="\t"))


def resolve_crossmaps(dir_, samples, crossmap, crossmap_dir, refmap):
    """sample -> removed-interval coordinates, for whichever arm this is.

    A single --crossmap covers the fixed-reference arm.  The matched arm needs
    one per sample, reached through the reference assignment in --refmap.
    """
    if crossmap:
        j = junctions(crossmap)
        return {s: j for s in samples}
    out = {}
    assign = {}
    if refmap and os.path.exists(refmap):
        assign = {r["sample"]: r["reference"]
                  for r in csv.DictReader(open(refmap), delimiter="\t")}
    cache = {}
    for s in samples:
        ref = assign.get(s)
        if not ref:
            out[s] = None
            continue
        if ref not in cache:
            cache[ref] = junctions(os.path.join(crossmap_dir,
                                                f"{ref}.crossmap.tsv"))
        out[s] = cache[ref]
    return out


def classify_site(pos, junc, slop):
    if junc is None:
        return UNKNOWN
    return SHARED if any(abs(pos - j) <= slop for j in junc) else LACKING


def load(path, min_reads_q):
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        if int(r["reads_q"]) < min_reads_q:
            continue
        both = int(r["both_el_termini"])
        out.append(dict(
            pos=int(r["clean_pos"]), orig=int(r["orig_pos"]),
            reads=int(r["reads"]), reads_q=int(r["reads_q"]),
            span=int(r["span"]), evidence=TWO if both else ONE,
            term=("S" if int(r["el_start"]) >= int(r["el_end"]) else "E")))
    out.sort(key=lambda x: x["pos"])
    return out


def collapse(stacks, pair_window, scatter_window):
    """Merge one-sided stacks that are two views of a single insertion.

    Two-sided stacks are left alone: they already carry both termini, so they
    are not a half view of anything.
    """
    ones = [x for x in stacks if x["evidence"] == ONE]
    sites, used = [], set()
    for i, x in enumerate(ones):
        if i in used:
            continue
        group = [x]
        used.add(i)
        for j in range(i + 1, len(ones)):
            if j in used:
                continue
            gap = ones[j]["pos"] - group[-1]["pos"]
            opp = ones[j]["term"] != group[-1]["term"]
            if gap > (pair_window if opp else scatter_window):
                if gap > max(pair_window, scatter_window):
                    break
                continue
            group.append(ones[j])
            used.add(j)
        sites.append(dict(
            pos=group[0]["pos"], orig=group[0]["orig"],
            reads=sum(g["reads"] for g in group),
            span=group[-1]["pos"] - group[0]["pos"],
            evidence=ONE, merged=len(group),
            paired=int(len(group) > 1 and
                       len({g["term"] for g in group}) > 1)))
    for o in stacks:
        if o["evidence"] == TWO:
            o = dict(o); o["merged"] = 1; o["paired"] = 0
            sites.append(o)
    sites.sort(key=lambda s: s["pos"])
    return sites


def sites_for(dir_, sample, min_reads_q, pw, sw, junc, slop):
    stacks = load(os.path.join(dir_, f"{sample}.elstacks.tsv"), min_reads_q)
    sites = collapse(stacks, pw, sw)
    for s in sites:
        s["site_class"] = classify_site(s["pos"], junc, slop)
    return stacks, sites


def count(dir_, min_reads_q, pw, sw, cmaps, slop):
    rows = []
    for p in sorted(glob.glob(os.path.join(dir_, "*.elstacks.tsv"))):
        s = os.path.basename(p).split(".")[0]
        stacks, sites = sites_for(dir_, s, min_reads_q, pw, sw,
                                  cmaps.get(s), slop)
        ev = collections.Counter(x["evidence"] for x in sites)
        cl = collections.Counter(x["site_class"] for x in sites)
        dp = os.path.join(dir_, f"{s}.elementdepth.tsv")
        cn = ""
        if os.path.exists(dp):
            v = float(next(csv.DictReader(open(dp), delimiter="\t"))["element_ratio"])
            cn = round(v, 2) if v > 0 else ""
        rows.append(dict(
            sample=s, depth_cn=cn,
            sites=ev[TWO] + ev[ONE], two_sided=ev[TWO], one_sided=ev[ONE],
            ref_shared=cl[SHARED], ref_lacking=cl[LACKING],
            unknown_class=cl[UNKNOWN],
            one_sided_stacks=sum(1 for x in stacks if x["evidence"] == ONE),
            collapsed=sum(x["merged"] - 1 for x in sites
                          if x["evidence"] == ONE),
            flank_pairs=sum(x["paired"] for x in sites)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="refbias/scale/p1g")
    ap.add_argument("--min-reads-q", type=int, default=10)
    ap.add_argument("--pair-window", type=int, default=2000)
    ap.add_argument("--scatter-window", type=int, default=50)
    ap.add_argument("--crossmap", default="",
                    help="one crossmap for every sample, for the "
                         "fixed-reference arm")
    ap.add_argument("--crossmap-dir", default="is6110/assets/isclean_matched",
                    help="per-reference crossmaps, for the matched arm")
    ap.add_argument("--refmap", default="refbias/scale/p1/refmap.tsv",
                    help="sample to reference assignment, for the matched arm")
    ap.add_argument("--class-slop", type=int, default=25,
                    help="how close a call must be to a removed-interval "
                         "coordinate to count as ref_shared")
    ap.add_argument("--sweep", action="store_true",
                    help="print site totals across a range of pair windows")
    ap.add_argument("--out", default="")
    ap.add_argument("--calls-out", default="",
                    help="per-call rows, for scoring against contig breaks")
    a = ap.parse_args()

    samples = sorted(os.path.basename(p).split(".")[0]
                     for p in glob.glob(os.path.join(a.dir, "*.elstacks.tsv")))
    if not samples:
        sys.exit(f"no *.elstacks.tsv under {a.dir}")
    cmaps = resolve_crossmaps(a.dir, samples, a.crossmap, a.crossmap_dir,
                              a.refmap)
    nc = sum(1 for s in samples if cmaps.get(s) is None)

    rows = count(a.dir, a.min_reads_q, a.pair_window, a.scatter_window,
                 cmaps, a.class_slop)
    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
    if a.calls_out:
        os.makedirs(os.path.dirname(a.calls_out) or ".", exist_ok=True)
        with open(a.calls_out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "r_pos", "evidence", "site_class", "span",
                        "reads", "merged", "placement"])
            for s in samples:
                for x in sites_for(a.dir, s, a.min_reads_q, a.pair_window,
                                   a.scatter_window, cmaps.get(s),
                                   a.class_slop)[1]:
                    w.writerow([s, x["orig"], x["evidence"], x["site_class"],
                                x["span"], x["reads"], x["merged"],
                                "on_path"])

    tot = lambda k: sum(r[k] for r in rows)
    print(f"{a.dir}: {len(rows)} isolates, pair window {a.pair_window} bp, "
          f"scatter window {a.scatter_window} bp"
          + (f", {nc} without a crossmap" if nc else ""))
    print("\n  evidence -- this IS a quality ordering")
    print(f"    two_sided                {tot('two_sided'):6d}")
    print(f"    one_sided (promoted)     {tot('one_sided'):6d}   "
          f"from {tot('one_sided_stacks')} stacks, {tot('collapsed')} collapsed "
          f"into {tot('flank_pairs')} flank pairs")
    print(f"    SITES                    {tot('sites'):6d}")
    print("\n  site class -- this is NOT a quality ordering")
    print(f"    ref_lacking              {tot('ref_lacking'):6d}   "
          f"the reference does not carry this insertion")
    print(f"    ref_shared               {tot('ref_shared'):6d}   "
          f"the reference carries it too, and the build removed it")
    if tot("unknown_class"):
        print(f"    unknown                  {tot('unknown_class'):6d}   "
              f"no crossmap for the isolate's reference")
    s = tot("sites")
    if s:
        print(f"    shared fraction          {tot('ref_shared')/s:6.3f}")

    cn = [r["depth_cn"] for r in rows if r["depth_cn"] != ""]
    if cn:
        sub = [r for r in rows if r["depth_cn"] != ""]
        m = st.mean(r["sites"] for r in sub)
        print(f"\n  mean depth-ratio estimate  {st.mean(cn):6.2f}")
        print(f"  mean sites                 {m:6.2f}   "
              f"recovery {m/st.mean(cn):.3f}")

    if a.sweep:
        print(f"\n  --- pair-window sweep ---\n  {'window':>8} {'sites':>7} "
              f"{'collapsed':>10}")
        for w in (0, 10, 50, 100, 250, 500, 1000, 1500, 2000, 3000, 5000):
            r2 = count(a.dir, a.min_reads_q, w, a.scatter_window, cmaps,
                       a.class_slop)
            print(f"  {w:8d} {sum(x['sites'] for x in r2):7d} "
                  f"{sum(x['collapsed'] for x in r2):10d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
