#!/usr/bin/env python3
"""Collate P1's per-isolate reference choices into one refmap, and report.

Written in Python rather than shell because the shell version silently corrupted
a row: `IFS=$'\t' read -r a b c` collapses CONSECUTIVE tabs, since tab is an IFS
whitespace character, so an empty field is lost and every later column shifts
left. The isolate with no main-lineage call (canettii) had an empty lineage
field, and its `deepest` column came out holding the mean depth. csv.DictReader
preserves empty fields.

WHY THE ELEMENT-CONTENT COLUMNS ARE HERE, AND WHY `margin` ALONE IS NOT ENOUGH
`margin` is the SNP gap to the runner-up, and on its own it does not say
whether the choice matters. Measured 2026-09-22 in
`refbias/REFERENCE_TIE_BREAK.md`: all 34 lineage 7 isolates were run through
the matched arm twice, forced to each of two candidates 350 SNPs apart at a
median margin of 1 SNP, and the call sets were indistinguishable -- same site
count in 34 of 34 isolates, same (evidence, site_class) multiset in 34 of 34.

The reason is that an element-free arm's calls depend on the reference's
element CONTENT, not on its SNP distance: substitutions do not move an element
boundary. Both lineage 7 candidates carried exactly one interval, at the
homologous locus, so identical content gave identical calls.

That makes a small margin harmless when content agrees and dangerous when it
does not. Across the 100-isolate cohort, 44 of the 49 isolates with both
crossmaps available have a chosen and runner-up pair whose interval counts
DIFFER, by 1 to 8, and 9 of those have a margin of 5 SNPs or less. Those 9 are
the at-risk set, and `margin` cannot distinguish them from the 30 isolates
with a small margin overall. `ref_intervals`, `rank2_intervals` and
`interval_delta` are what does.

A blank in those columns means no crossmap exists for that reference, which is
not zero intervals -- 51 of 100 cohort isolates are missing a crossmap for one
of their two candidates.

THE TIE-BREAK, AND WHAT IT IS AND IS NOT JUSTIFIED BY
When the top two candidates are within --tie-margin SNPs and one carries more
element intervals, the richer reference is chosen. Four experiments, written up
in `refbias/REFERENCE_TIE_BREAK.md`:

    interval gap 0, 34 isolates   no difference in the calls at all
    interval gap 2,  5 isolates   no difference; read support identical 5 of 5
    interval gap 4,  1 isolate    same count, 1 of 19 calls reclassified
                                  ref_shared -> ref_lacking and 1 downgraded
                                  two_sided -> one_sided on the poorer reference
    interval gap 6,  1 isolate    counts DIFFER, 16 against 15, and the richer
                                  reference is 15/15 two_sided against 14/16

THIS IS A QUALITY PREFERENCE, NOT A DOMINANCE ARGUMENT.
An earlier version of this note claimed weak dominance: that a reference
interval the sample does not share is silent -- true, and measured, 8 of 14
intervals produced no call for one isolate -- so extra intervals could only add
correct ref_shared calls and never remove anything. The gap-6 experiment
falsifies the conclusion. Swapping GCF_014901095 for GCF_009730235 upgraded one
one_sided/ref_shared call to two_sided, exactly as the mechanism predicts, but
also made a one_sided/ref_lacking call disappear, so the richer reference
returned ONE FEWER site.

The reason the argument fails is that it describes adding intervals in
isolation, and a reference swap is not that experiment: the two genomes are 315
and 319 SNPs from the sample and differ from each other throughout, so the
alignments differ for reasons unrelated to element content. The read-support
vectors show it plainly.

What survives is the quality comparison, and it favours the richer reference in
all four experiments: it never produced worse evidence, and at gap 6 it gave
15 of 15 two_sided against 14 of 16 with two one_sided calls. The call it lost
was one_sided/ref_lacking, the cell measured at 0.106 contradicted by an
assembly against 0.000 for one_sided/ref_shared, so it is more likely to have
dropped a false positive than a real site -- likely, not demonstrated.

So the rule is applied because it improves evidence quality on the evidence
available, not because a richer reference cannot hurt. --tie-margin 0 disables
it and `tie_break` records every firing, so it is auditable and reversible.

The label column is reported but never used to choose. Stage 1 measured
label-gated selection as 0.6 points worse than distance selection, because 16.1%
of nearest references lie outside the sample's own label. This prints whether the
same holds on real isolates, which is the first check of that number outside
simulation.
"""
import argparse, collections, csv, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="refbias/cohort.pilot.tsv")
    ap.add_argument("--dir", default="refbias/p1")
    ap.add_argument("--lineages", default="data/tbprof_lineages.csv",
                    help="panel genome lineage labels, to place the chosen reference")
    ap.add_argument("--crossmap-dir", default="is6110/assets/isclean_matched",
                    help="per-reference crossmaps, to count each candidate's "
                         "element intervals. A reference with no crossmap "
                         "leaves the count blank, which is not zero.")
    ap.add_argument("--tie-margin", type=int, default=5,
                    help="two candidates within this many SNPs are a tie, and "
                         "the one carrying more element intervals is chosen. "
                         "0 disables the tie-break and restores pure "
                         "distance selection.")
    ap.add_argument("--out", default="refbias/p1/refmap.tsv")
    ap.add_argument("--allow-missing", action="store_true",
                    help="exit 0 even when some samples have no output. Off by "
                         "default: a short table silently drops those samples "
                         "from every later pass, because P5 sizes itself from "
                         "the refmap")
    a = ap.parse_args()

    panel_lin = {}
    if os.path.exists(a.lineages):
        for r in csv.DictReader(open(a.lineages)):
            panel_lin[r["sample"]] = ((r.get("main_lineage") or "").strip(),
                                      (r.get("sub_lineage") or "").strip())

    # element intervals per candidate, cached: the crossmap has one row per
    # interval the build removed. Defined before the selection loop because
    # the tie-break needs it while choosing, not only when reporting.
    iv_cache = {}

    def n_intervals(ref):
        if not ref:
            return ""
        if ref not in iv_cache:
            q = os.path.join(a.crossmap_dir, f"{ref}.crossmap.tsv")
            if os.path.exists(q):
                with open(q) as fh:
                    iv_cache[ref] = max(sum(1 for _ in fh) - 1, 0)
            else:
                iv_cache[ref] = ""
        return iv_cache[ref]

    rows, missing = [], []
    for c in csv.DictReader(open(a.cohort), delimiter="\t"):
        s = c["sample"]
        p = os.path.join(a.dir, f"{s}.candidates.tsv")
        if not os.path.exists(p):
            missing.append(s); continue
        cand = list(csv.DictReader(open(p), delimiter="\t"))
        if not cand:
            missing.append(s); continue
        best = cand[0]
        nxt = cand[1] if len(cand) > 1 else {}

        # Tie-break on element content. Only the top two are considered: a
        # third candidate further away is not a tie in any useful sense, and
        # widening it would trade distance for content without a measurement
        # to justify the exchange.
        tie = ""
        if a.tie_margin > 0 and nxt.get("snp_distance"):
            gap = int(nxt["snp_distance"]) - int(best["snp_distance"])
            iv_b, iv_n = n_intervals(best["reference"]), n_intervals(nxt["reference"])
            if (gap <= a.tie_margin and iv_b != "" and iv_n != ""
                    and iv_n > iv_b):
                tie = (f"intervals {iv_b}->{iv_n} at margin {gap}")
                best, nxt = nxt, best
        ref = best["reference"]
        rmain, rsub = panel_lin.get(ref, ("", ""))
        slab = (c.get("lineage") or "").strip()
        # "inside the label" means the chosen reference carries the same main
        # lineage the sample was assigned. Unknown either side is not a match.
        if not slab or not rmain:
            inside = "unknown"
        else:
            inside = "yes" if rmain == slab else "no"
        rows.append(dict(
            sample=s, sample_lineage=slab,
            sample_sublineage=(c.get("deepest") or "").strip(),
            meandepth=c.get("meandepth", ""),
            reference=ref, ref_lineage=rmain, ref_sublineage=rsub,
            snp_distance=int(best["snp_distance"]),
            inside_label=inside,
            rank2_reference=nxt.get("reference", ""),
            rank2_distance=nxt.get("snp_distance", ""),
            # ALWAYS the non-negative SNP gap between the two candidates.
            # An earlier version computed rank2 minus best after the tie-break
            # had already swapped them, so margin came out NEGATIVE for every
            # isolate the rule fired on -- and since consumers test
            # `int(margin) <= 5`, a negative value passes every threshold
            # silently. When `tie_break` is set the chosen reference is the
            # FARTHER of the two by this margin; that direction lives in
            # `tie_break`, not in the sign of this column.
            margin=(abs(int(nxt["snp_distance"]) - int(best["snp_distance"])))
                   if nxt.get("snp_distance") else "",
            tie_break=tie))

    for r in rows:
        a_iv, b_iv = n_intervals(r["reference"]), n_intervals(r["rank2_reference"])
        r["ref_intervals"] = a_iv
        r["rank2_intervals"] = b_iv
        r["interval_delta"] = (abs(a_iv - b_iv) if a_iv != "" and b_iv != ""
                               else "")

    if not rows:
        print("no candidates found", file=sys.stderr); return 1
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    print(f"  {len(rows)} isolates with a reference"
          f"{'; MISSING: ' + ', '.join(missing) if missing else ''}\n")

    # the at-risk set: a near-arbitrary choice between candidates whose element
    # content genuinely differs. Neither column alone identifies it.
    def small(r, lim=5):
        return str(r["margin"]).isdigit() and int(r["margin"]) <= lim

    scored = [r for r in rows if r["interval_delta"] != ""]
    differ = [r for r in scored if r["interval_delta"] > 0]
    # Residual risk, not raw risk. A case where the chosen reference is
    # already the richer of the two is not at risk -- there is nothing the
    # tie-break could improve. Counting those was over-reporting: the earlier
    # screen still flagged 4 scale100 isolates after the rule had fired,
    # all of them cases where the selector had picked the richer reference
    # unaided.
    def richer_chosen(r):
        return (r["ref_intervals"] != "" and r["rank2_intervals"] != ""
                and r["ref_intervals"] >= r["rank2_intervals"])
    at_risk = [r for r in differ if small(r) and not richer_chosen(r)]
    nsmall = sum(1 for r in rows if small(r))
    print(f"  reference choice, element content:")
    print(f"    both crossmaps available          {len(scored):4d} of {len(rows)}"
          f"   (blank means no crossmap, not zero intervals)")
    print(f"    candidates carry the same count   {len(scored) - len(differ):4d}")
    print(f"    candidates DIFFER                 {len(differ):4d}")
    print(f"    margin <= 5 SNPs                  {nsmall:4d}")
    fired = [r for r in rows if r["tie_break"]]
    print(f"    tie-break FIRED (chose more intervals) {len(fired):4d}"
          + (f"  at margin <= {a.tie_margin}" if fired else ""))
    for r in fired:
        print(f"      {r['sample']:<17s} -> {r['reference']} "
              f"({r['ref_intervals']} iv) over {r['rank2_reference']} "
              f"({r['rank2_intervals']} iv); {r['tie_break']}")
    print(f"    RESIDUAL RISK after the tie-break {len(at_risk):4d}"
          + ("   <- differ, small margin, and the poorer reference still "
             "chosen" if at_risk
             else "   (every small-margin case now takes the richer reference)"))
    for r in at_risk[:8]:
        print(f"      {r['sample']:<17s} margin {str(r['margin']):>4s}  "
              f"{r['reference']} ({r['ref_intervals']} intervals) vs "
              f"{r['rank2_reference']} ({r['rank2_intervals']})")
    print(f"  {'isolate':<17s}{'sample lin':<11s}{'reference':<16s}"
          f"{'ref lin':<11s}{'dist':>6s}{'margin':>8s}  label")
    for r in sorted(rows, key=lambda x: x["snp_distance"]):
        print(f"  {r['sample']:<17s}{(r['sample_lineage'] or '-'):<11s}"
              f"{r['reference']:<16s}{(r['ref_lineage'] or '-'):<11s}"
              f"{r['snp_distance']:>6d}{str(r['margin']):>8s}  {r['inside_label']}")

    d = sorted(r["snp_distance"] for r in rows)
    print(f"\n  selected-reference distance: min {d[0]}, median {d[len(d)//2]}, "
          f"max {d[-1]}")
    refs = collections.Counter(r["reference"] for r in rows)
    print(f"  distinct references: {len(refs)} for {len(rows)} isolates")
    for ref, n in refs.most_common():
        if n > 1:
            who = [r["sample"] for r in rows if r["reference"] == ref]
            print(f"    {ref} chosen by {n}: {', '.join(who)}")
    lab = collections.Counter(r["inside_label"] for r in rows)
    known = lab["yes"] + lab["no"]
    print(f"\n  chosen reference vs the sample's tb-profiler label:")
    print(f"    inside  {lab['yes']:>3d}")
    print(f"    outside {lab['no']:>3d}"
          + (f"  ({100*lab['no']/known:.1f}% of {known} with both labels known)"
             if known else ""))
    if lab["unknown"]:
        print(f"    unknown {lab['unknown']:>3d}  (sample or reference unlabelled)")
    print(f"    stage 1 on simulated data: 16.1% of nearest references lay "
          f"outside the label")
    print(f"\n  written: {a.out}")
    if missing and not a.allow_missing:
        print(f"FATAL: {len(missing)} samples have no usable output "
              f"({', '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}); "
              f"rerun them, or pass --allow-missing to accept a partial table",
              file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
