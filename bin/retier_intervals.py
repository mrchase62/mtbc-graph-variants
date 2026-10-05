#!/usr/bin/env python3
"""Re-tier the deletion catalogue by measured evidence instead of provenance.

THE PROBLEM WITH THE EXISTING TIER. `support_tier` records WHERE an interval
came from -- A_multi_assembly and B_single_assembly from the 333-genome graph,
C_caller_only from clustering the per-isolate callers. That is provenance, and
the RD cross-reference showed it predicting reliability only crudely: among
known regions of difference with a checkable lineage label, 28 of 45 tier-A
intervals were lineage-concordant against 2 of 8 caller-only ones. Meanwhile
70% of caller-only intervals have ZERO isolates genotyped as deleted in 997 --
the callers proposed them and independent depth never confirmed one.

Provenance is also fixed for the life of the catalogue, while the evidence grows
with every cohort. So this adds an `evidence_tier` beside it, computed from what
was actually measured, and leaves `support_tier` untouched so nothing that reads
it breaks.

WHAT IS MEASURED, and what was rejected.

  depth support    how many isolates are genotyped deleted. The single
                   strongest signal: 5% of tier-A intervals have none against
                   70% of caller-only ones.
  lineage          of the isolates called deleted, the largest share belonging
  precision        to one lineage. A real deletion is inherited and therefore
                   lineage-restricted; a depth artefact scatters. Measured
                   across all 3,658 intervals, known RDs sit at 0.84 to 0.95
                   median precision and non-RDs at 0.45 to 0.50 -- and that
                   separation holds WITHIN every provenance tier, which is why
                   it is worth having.
  context          a measured repeat tract or an element landmark within 500 bp,
                   where depth counts paralogous copies rather than the interval.
                   Repeat content comes from the rebuilt mask, which is measured
                   rather than name-based; the name-based version missed
                   Rv1759c/wag22 entirely.

REJECTED: phylogenetic coherence, the number of independent losses. The
intuition is that a real deletion is inherited once and an artefact scatters,
but measured on this catalogue it runs the other way -- tier A has a median of
8 independent losses and caller-only has 0, because well-supported intervals are
deleted in many lineages while unsupported ones are deleted nowhere. Coherence
is therefore not usable as a discriminator here and is not used.

RECALL IS NOT USED EITHER. The share of a lineage that carries the deletion
looks like the natural partner to precision, but it is low everywhere -- 0.04 to
0.23 median -- because most catalogued intervals are sub-lineage RDs and the
lineage labels available are major lineages. Thresholding on it would reject
almost every real sub-lineage deletion.
"""
import argparse, bisect, collections, csv, os, sys


def load_bed_spans(path):
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[1].isdigit():
                out.append((int(f[1]), int(f[2]),
                            f[3] if len(f) > 3 else ""))
    return sorted(out)


def overlapper(spans):
    st = [s[0] for s in spans]

    def hit(a, b):
        i = bisect.bisect_right(st, b)
        for j in range(max(0, i - 8), min(len(spans), i + 1)):
            if min(spans[j][1], b) > max(spans[j][0], a):
                return spans[j][2] or "masked"
        return ""
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--intervals", default="refbias/assets/sv_intervals.tsv")
    ap.add_argument("--iv-states", required=True,
                    help="svgt_iv_states.tsv for the cohort supplying evidence")
    ap.add_argument("--refmap", required=True)
    ap.add_argument("--mask", default="data/annotation/H37Rv_repeat_mask.bed")
    ap.add_argument("--rds", default="data/annotation/known_RDs.bed")
    ap.add_argument("--min-deleted", type=int, default=2)
    ap.add_argument("--min-precision", type=float, default=0.80)
    ap.add_argument("--min-for-precision", type=int, default=4,
                    help="deleted isolates needed before lineage precision is "
                         "treated as evidence; on two isolates it is 1.0 "
                         "whenever they happen to share a lineage")
    ap.add_argument("--cohort-tag", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ivs = list(csv.DictReader(open(a.intervals, newline=""), delimiter="\t"))
    lin = {r["sample"]: (r.get("sample_lineage") or "?") for r in
           csv.DictReader(open(a.refmap, newline=""), delimiter="\t")}
    dele = collections.defaultdict(list)
    called = collections.Counter()
    with open(a.iv_states, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["state"] in ("ALT", "ABSENT"):
                dele[r["interval"]].append(r["sample"])
            if r["state"] != "NOCALL":
                called[r["interval"]] += 1
    mask_hit = overlapper(load_bed_spans(a.mask))
    rd_hit = overlapper(load_bed_spans(a.rds))
    print(f"  {len(ivs):,} intervals; evidence from {a.iv_states} "
          f"({len(lin):,} isolates)")

    n = collections.Counter()
    for r in ivs:
        s, e = int(r["start"]), int(r["end"])
        d = dele.get(r["interval"], [])
        r["n_deleted"] = len(d)
        r["n_called"] = called.get(r["interval"], 0)
        prec, top = "", ""
        if len(d) >= 2:
            c = collections.Counter(lin.get(x, "?") for x in d)
            top, k = c.most_common(1)[0]
            prec = round(k / len(d), 3)
        r["lineage_precision"] = prec
        r["top_lineage"] = top
        ctx = mask_hit(s, e)
        r["repeat_context"] = ctx.split("|")[0] if ctx else ""
        r["is_known_rd"] = 1 if rd_hit(s, e) else 0

        # PRECISION DECIDES; CONTEXT ONLY EXPLAINS. The first version of this
        # tested context first, which demoted 228 RD-overlapping intervals to
        # E3 for sitting in a repeat region -- and 29 of the 135 known RDs are
        # at least half masked, because real deletions often remove PE/PPE
        # genes. An interval whose deletions are lineage-coherent at 0.95 is
        # more likely real than one scattered at 0.45, whatever the context. So
        # context is recorded as a flag on every row and only separates the
        # scattered intervals into those with an explanation and those without.
        #
        # E1 also requires more than a handful of deleted isolates: precision on
        # two isolates is 1.0 whenever they share a lineage, which is a coin
        # toss dressed as evidence.
        r["context_flag"] = ctx.split("|")[0] if ctx else (
            "element_prox" if str(r.get("is6110_prox", "0")) == "1" else "")
        if len(d) < a.min_deleted:
            t = "E4_unsupported"
        elif prec != "" and float(prec) >= a.min_precision:
            t = ("E1_coherent" if len(d) >= a.min_for_precision
                 else "E2_too_few")
        elif r["context_flag"]:
            t = "E3_scattered_in_repeat"
        else:
            t = "E3_scattered"
        r["evidence_tier"] = t
        r["evidence_cohort"] = a.cohort_tag
        n[t] += 1

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ivs[0]), delimiter="\t")
        w.writeheader(); w.writerows(ivs)

    print(f"\n  EVIDENCE TIERS")
    order = ["E1_coherent", "E2_too_few", "E3_scattered",
             "E3_scattered_in_repeat", "E4_unsupported"]
    for t in order:
        print(f"    {t:<22}{n[t]:>6,}  ({n[t]/len(ivs):.0%})")

    print(f"\n  cross-tabulated against provenance:")
    x = collections.Counter((r["support_tier"], r["evidence_tier"]) for r in ivs)
    prov = sorted({r["support_tier"] for r in ivs})
    print(f"    {'provenance':<20}" + "".join(f"{t.split('_')[0]:>8}" for t in order))
    for p in prov:
        print(f"    {p:<20}" + "".join(f"{x[(p,t)]:>16,}" for t in order))

    print(f"\n  VALIDATION -- known RDs recovered, by evidence tier")
    print(f"  a tier that concentrates the known RDs is doing its job")
    for t in order:
        g = [r for r in ivs if r["evidence_tier"] == t]
        rd = sum(1 for r in g if r["is_known_rd"])
        print(f"    {t:<22}{len(g):>6,} intervals, {rd:>4} overlap a known RD "
              f"({rd/max(1,len(g)):.0%})")
    e1 = [r for r in ivs if r["evidence_tier"] == "E1_coherent"]
    allrd = sum(1 for r in ivs if r["is_known_rd"])
    print(f"    of {allrd:,} RD-overlapping intervals, "
          f"{sum(1 for r in e1 if r['is_known_rd']):,} "
          f"({sum(1 for r in e1 if r['is_known_rd'])/max(1,allrd):.0%}) are E1")
    print(f"\n  -> {a.out}")


if __name__ == "__main__":
    main()
