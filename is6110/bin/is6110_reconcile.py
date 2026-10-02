#!/usr/bin/env python3
"""Stage 2 of P1h: reconcile the two directions a junction can be read from.

WHY THIS IS A JOIN AND NOT A COMPARISON
There are three counts over the same 23 alignments -- 128 chromosome-side on the
normal FASTA, 197 chromosome-side on the IS-clean FASTA, 237 element-side -- and
one depth estimate at 259.7. Each pair of those differs by a single ratio, and a
ratio does not say which rows moved or why. This script replaces the ratios with
a per-row accounting: every element-side stack is traced back through the
chromosome-side scoring to the step that dropped it, and every chromosome-side
site with no element-side stack is reported too.

THE ORDER MATTERS. No threshold is retuned until this table exists, because a
threshold chosen to close a gap is fitted to the answer. That is why
is6110_element_side.py records SA mapping quality as a column instead of
applying it.

THE CHROMOSOME-SIDE STEPS, in the order a position passes through them

  1. min_clips      p1g_isclean.sh only proposes a position carrying >= 5
                    clipped ends, so a junction below that never enters the
                    junction table at all
  2. element SA     is6110_junctions.py sets is6110=1 when the far side lands on
                    the element contig
  3. read-through   is6110_isclean_summary.py drops a position whose reads mostly
                    span it, at a fixed fraction of 0.10
  4. clustering     single-linkage at radius 1000 groups positions into sites
  5. min_peak       a cluster counts only if its peak carries >= 10 element SA
                    records

A stack that the element side proposes and the chromosome side does not was lost
at exactly one of those, and the label says which. Step 4 has two distinct
outcomes and they are reported separately: a stack can be MERGED into a cluster
that is kept under a different peak, which loses a site while keeping a site,
or it can sit in a cluster that is then dropped at step 5.

SELF-CHECK. The chromosome-side scoring is reimplemented here rather than
imported, because is6110_isclean_summary.py holds it inside main(). That is a
duplication risk, so the reproduced per-isolate site count is asserted against
is6110/results/isclean_summary.tsv and the run fails if any isolate disagrees.

THE SECOND JOB, added after P1H_STAGE1.md section 7
Of the 237 element-side stacks, 132 carry target-site duplication geometry and
are independently proposed by ISMapper 94.7% of the time; 69 are one-sided and
are proposed 24.6% of the time. Those 69 are not simply noise -- section 5 shows
one copy in a repeat interval producing one well-placed stack and four one-sided
fragments -- so this script also records, for every stack, the distance to its
nearest neighbour and whether that neighbour carries the complementary element
terminus. Without that split neither 237 nor 132 is a reportable total.

LABELLING NOTE, 2026-09-21.  This script scores on the A/B/C tiers, which have
been replaced by two independent columns -- `evidence` (two_sided/one_sided, a
quality ordering) and `site_class` (ref_shared/ref_lacking, not an ordering).
The rule now lives in is6110_promote_sites.py and is described in
is6110/docs/PROMOTED_TIER.md.  `tierAB` below means `evidence = two_sided`, and
it excludes promoted one-sided calls, so any ratio it forms against a
depth-based estimate reads low.  The numbers this script produced are still the
numbers that were measured; do not read tier A as more confident than tier B.
"""
import argparse, bisect, collections, csv, glob, os, sys


def load_crossmap(path):
    if not path or not os.path.exists(path):
        return []
    return [{k: int(v) for k, v in r.items()}
            for r in csv.DictReader(open(path), delimiter="\t")]


def clean_to_orig(p, rows):
    shift = 0
    for r in rows:
        if p > r["clean_junction"]:
            shift = r["cum_deleted"]
        else:
            break
    return p + shift


def cluster(items, radius):
    """Single-linkage on sorted (pos, ...) tuples. Identical to the grouping in
    is6110_isclean_summary.py, which is the point: this has to reproduce it."""
    out, cur = [], []
    for p in sorted(items, key=lambda x: x[0]):
        if cur and p[0] - cur[-1][0] > radius:
            out.append(cur); cur = []
        cur.append(p)
    if cur:
        out.append(cur)
    return out


def score_chromosome_side(path, max_rt, radius, min_peak):
    """Replay the chromosome-side scoring, keeping the intermediate state.

    Returns (rows_by_pos, element_positions, spanned_positions, kept_peaks,
             pos_to_cluster_index, kept_cluster_indices).
    """
    rows, elem, spanned = {}, [], set()
    for r in csv.DictReader(open(path), delimiter="\t"):
        pos = int(r["pos"])
        rows[pos] = r
        if not int(r["is6110"]):
            continue
        clips = int(r["clips_start"]) + int(r["clips_end"])
        rt = int(r["readthrough"])
        if clips + rt and rt / (clips + rt) > max_rt:
            spanned.add(pos)
            continue
        elem.append((pos, int(r["sa_at_is6110"]),
                     int(r["clips_start"]), int(r["clips_end"])))
    clusters = cluster(elem, radius)
    pos_cl, kept_idx, peaks = {}, set(), {}
    for i, cl in enumerate(clusters):
        pk = max(cl, key=lambda x: x[1])
        for m in cl:
            pos_cl[m[0]] = i
        if pk[1] >= min_peak:
            kept_idx.add(i)
            peaks[i] = pk[0]
    return rows, {p for p, *_ in elem}, spanned, peaks, pos_cl, kept_idx


def nearest(sorted_keys, p):
    if not sorted_keys:
        return None
    i = bisect.bisect_left(sorted_keys, p)
    cands = [sorted_keys[j] for j in (i - 1, i) if 0 <= j < len(sorted_keys)]
    return min(cands, key=lambda k: abs(k - p)) if cands else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--elside-dir", default="refbias/p1h")
    ap.add_argument("--junc-dir", default="refbias/p1g")
    ap.add_argument("--ismapper-dir", default="refbias/p1f")
    ap.add_argument("--crossmap", default="is6110/assets/H37Rv.isclean.crossmap.tsv",
                    help="one crossmap, when every isolate shares a reference")
    ap.add_argument("--crossmap-dir", default=None,
                    help="per-reference crossmaps, for P1i where the 23 isolates "
                         "use 18 different references and therefore 18 coordinate "
                         "systems. Requires --refmap. Setting this DISABLES the "
                         "ISMapper join, see below")
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv",
                    help="sample -> reference, for --crossmap-dir")
    ap.add_argument("--check-against", default="is6110/results/isclean_summary.tsv",
                    help="the chromosome-side counts this run must reproduce. "
                         "Pass an empty string where no such record exists yet, "
                         "which forfeits the self-check and is reported as such")
    ap.add_argument("--radius", type=int, default=1000)
    ap.add_argument("--min-peak", type=int, default=10)
    ap.add_argument("--max-readthrough-frac", type=float, default=0.10)
    ap.add_argument("--min-reads-q", type=int, default=10,
                    help="element-side stacks below this many uniquely-placed "
                         "reads are not carried into the accounting")
    ap.add_argument("--window", type=int, default=20,
                    help="an element-side stack and a chromosome-side position "
                         "are the same junction within this distance. The two "
                         "are computed from the same reads, so they should agree "
                         "exactly; the window absorbs which side of a target-site "
                         "duplication each direction reports")
    ap.add_argument("--ism-window", type=int, default=50)
    ap.add_argument("--out", default="is6110/results/p1h_reconcile.tsv")
    ap.add_argument("--summary-out", default="is6110/results/p1h_reconcile_summary.tsv")
    a = ap.parse_args()

    # COORDINATE FRAMES. Everything joined in this script must be in ONE frame.
    # P1g/P1h are all on IS-clean H37Rv, so one crossmap serves. P1i puts each
    # isolate on its own matched reference -- 18 frames for 23 isolates -- and
    # ISMapper was run against H37Rv, in NC_000962.3 coordinates. Joining those
    # would compare positions from different genomes, which is exactly the slip
    # recorded in P1I_STAGE3.md section 5 and the three before it. So
    # --crossmap-dir turns the ISMapper join OFF rather than computing it in the
    # wrong frame.
    #
    # Restoring it does NOT need new machinery. Each matched reference is a path
    # in the pangenome graph, so `odgi position -r GCF_000195955#1#NC_000962.3`
    # converts an R coordinate to an H37Rv one by path lookup, which is what
    # bin/p4_place.sh already does. Feeding this script projected coordinates is
    # the way to get the column back; inferring them here is not.
    per_sample_xmap = {}
    xmap = []
    if a.crossmap_dir:
        ref_of = {r["sample"]: r["reference"] for r in
                  csv.DictReader(open(a.refmap), delimiter="\t")}
        for smp, ref in ref_of.items():
            cp = os.path.join(a.crossmap_dir, f"{ref}.crossmap.tsv")
            if os.path.exists(cp):
                per_sample_xmap[smp] = load_crossmap(cp)
        if a.ismapper_dir:
            print("  note: --crossmap-dir given, so the ISMapper join is OFF. "
                  "ISMapper positions are in H37Rv coordinates and these stacks "
                  "are not.")
            a.ismapper_dir = ""
    else:
        xmap = load_crossmap(a.crossmap)
    expect = {}
    if a.check_against and os.path.exists(a.check_against):
        for r in csv.DictReader(open(a.check_against), delimiter="\t"):
            expect[r["sample"]] = int(r["n_sites"])
    else:
        print("  note: no --check-against record, so the replayed "
              "chromosome-side scoring is NOT verified against a recorded "
              "count. The accounting below is self-consistent but unaudited.")

    out_rows, mismatches = [], []
    for ep in sorted(glob.glob(os.path.join(a.elside_dir, "*.elstacks.tsv"))):
        s = os.path.basename(ep).split(".")[0]
        jp = os.path.join(a.junc_dir, f"{s}.junctions.tsv")
        if not os.path.exists(jp):
            sys.exit(f"{s}: no chromosome-side table at {jp}")
        rows, elem_pos, spanned, peaks, pos_cl, kept_idx = score_chromosome_side(
            jp, a.max_readthrough_frac, a.radius, a.min_peak)
        # In --crossmap-dir mode a missing crossmap used to fall back to the
        # empty (identity) map without a word, which leaves every stack
        # downstream of an excised copy shifted. Every reference in every cohort
        # so far has one (18, 90 and 150 checked), so a missing one is an error.
        if a.crossmap_dir and s not in per_sample_xmap:
            sys.exit(f"FATAL: {s}: no crossmap for its reference in "
                     f"{a.crossmap_dir}; refusing to use an identity map")
        xm = per_sample_xmap.get(s, xmap)

        if s in expect and len(kept_idx) != expect[s]:
            mismatches.append((s, len(kept_idx), expect[s]))

        stacks = [r for r in csv.DictReader(open(ep), delimiter="\t")
                  if int(r["reads_q"]) >= a.min_reads_q]
        spos = sorted(int(r["clean_pos"]) for r in stacks)
        by_pos = {int(r["clean_pos"]): r for r in stacks}

        ism = []
        tp = os.path.join(a.ismapper_dir, s, s, "IS6110",
                          f"{s}__NC_000962.3_table.txt")
        if os.path.exists(tp):
            for r in csv.DictReader(open(tp), delimiter="\t"):
                try:
                    ism.append((int(r["x"]), int(r["y"]), r.get("call", "")))
                except (ValueError, KeyError):
                    pass

        all_rowpos = sorted(rows)
        elem_sorted = sorted(elem_pos)
        peak_positions = sorted(peaks.values())
        matched_peaks = set()

        for st in stacks:
            p = int(st["clean_pos"])
            orig = int(st["orig_pos"])

            # --- where did the chromosome side lose it, if it did?
            near_peak = nearest(peak_positions, p)
            if near_peak is not None and abs(near_peak - p) <= a.window:
                verdict = "agreed"
                matched_peaks.add(near_peak)
            else:
                ne = nearest(elem_sorted, p)
                nr = nearest(all_rowpos, p)
                if ne is not None and abs(ne - p) <= a.window:
                    ci = pos_cl.get(ne)
                    if ci in kept_idx:
                        verdict = "merged_by_clustering"
                        matched_peaks.add(peaks[ci])
                    else:
                        verdict = "min_peak"
                elif nr is not None and abs(nr - p) <= a.window:
                    verdict = ("readthrough" if nr in spanned
                               else "no_element_sa")
                else:
                    verdict = "not_scanned"

            # --- geometry class
            two = int(st["both_el_termini"])
            span = int(st["span"])
            cls = ("tsd" if two and 2 <= span <= 6 else
                   "two_sided_wide" if two else "one_sided")

            # --- nearest element-side neighbour, and is it complementary?
            i = spos.index(p)
            nbrs = [spos[j] for j in (i - 1, i + 1) if 0 <= j < len(spos)]
            nb = min(nbrs, key=lambda q: abs(q - p)) if nbrs else None
            nb_dist, nb_comp = "", ""
            if nb is not None:
                nb_dist = abs(nb - p)
                o = by_pos[nb]
                mine = ("START" if int(st["el_start"]) and not int(st["el_end"])
                        else "END" if int(st["el_end"]) and not int(st["el_start"])
                        else "BOTH")
                theirs = ("START" if int(o["el_start"]) and not int(o["el_end"])
                          else "END" if int(o["el_end"]) and not int(o["el_start"])
                          else "BOTH")
                nb_comp = int({"START": "END", "END": "START"}.get(mine, "") == theirs)

            ism_hit = int(any(abs(orig - x) <= a.ism_window or
                              abs(orig - y) <= a.ism_window for x, y, _ in ism))

            out_rows.append(dict(
                sample=s, clean_pos=p, orig_pos=orig,
                reads=int(st["reads"]), reads_q=int(st["reads_q"]),
                geometry=cls, span=span,
                el_start=int(st["el_start"]), el_end=int(st["el_end"]),
                sa_mapq_mean=st["sa_mapq_mean"],
                chrom_side=verdict,
                nbr_dist=nb_dist, nbr_complementary=nb_comp,
                ismapper=ism_hit))

        for pk in peak_positions:
            if pk in matched_peaks:
                continue
            out_rows.append(dict(
                sample=s, clean_pos=pk, orig_pos=clean_to_orig(pk, xm),
                reads="", reads_q="", geometry="", span="",
                el_start="", el_end="", sa_mapq_mean="",
                chrom_side="chrom_side_only", nbr_dist="", nbr_complementary="",
                ismapper=int(any(abs(clean_to_orig(pk, xm) - x) <= a.ism_window or
                                 abs(clean_to_orig(pk, xm) - y) <= a.ism_window
                                 for x, y, _ in ism))))

    if mismatches:
        for s, got, want in mismatches:
            print(f"FATAL {s}: reproduced {got} chromosome-side sites, "
                  f"{a.check_against} says {want}", file=sys.stderr)
        sys.exit("the replayed chromosome-side scoring does not match the "
                 "recorded one, so the accounting below would be against the "
                 "wrong baseline")

    out_rows.sort(key=lambda r: (r["sample"], r["clean_pos"]))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(out_rows)

    # ---- the accounting
    piv = collections.Counter()
    ismn = collections.Counter(); ismd = collections.Counter()
    for r in out_rows:
        k = (r["chrom_side"], r["geometry"])
        piv[k] += 1
        ismd[k] += 1; ismn[k] += r["ismapper"]
    order = ["agreed", "merged_by_clustering", "readthrough", "min_peak",
             "no_element_sa", "not_scanned", "chrom_side_only"]
    geos = ["tsd", "two_sided_wide", "one_sided", ""]
    print("\nWhere every element-side stack stands against the chromosome-side "
          "scoring\n")
    hdr = (f"{'chromosome-side outcome':22s} {'TSD':>5s} {'2sided':>7s} "
           f"{'1sided':>7s} {'total':>6s} {'ISMapper':>9s}")
    print(hdr); print("-" * len(hdr))
    with open(a.summary_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom_side", "tsd", "two_sided_wide", "one_sided",
                    "total", "ismapper_confirmed"])
        for k in order:
            cells = [piv[(k, g)] for g in geos]
            tot = sum(cells)
            if not tot:
                continue
            conf = sum(ismn[(k, g)] for g in geos)
            print(f"{k:22s} {cells[0]:5d} {cells[1]:7d} {cells[2]:7d} "
                  f"{tot:6d} {conf:6d} {conf/tot:6.0%}")
            w.writerow([k, cells[0], cells[1], cells[2] + cells[3], tot, conf])
    print("-" * len(hdr))
    tot = len(out_rows); conf = sum(r["ismapper"] for r in out_rows)
    print(f"{'all rows':22s} {sum(piv[(k,'tsd')] for k in order):5d} "
          f"{sum(piv[(k,'two_sided_wide')] for k in order):7d} "
          f"{sum(piv[(k,'one_sided')] for k in order):7d} {tot:6d} "
          f"{conf:6d} {conf/tot:6.0%}")
    print(f"\n  written: {a.out}\n           {a.summary_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
