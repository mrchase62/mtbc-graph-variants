#!/usr/bin/env python3
"""Collate P1g: junction sites on the IS-clean reference, against the other estimates.

WHY POSITIONS MUST BE CLUSTERED
An insertion does not produce one clipped position, it produces a spread. The two
sides of a junction sit 3-4 bp apart across the target-site duplication, and
around a repeat-rich locus the scatter is far wider. SAMEA104061150 carries a
single IS6110 copy by every other measure, and it yields eighteen positions with
element-contig SA support spanning 4.4 kb of original coordinates -- one dominant
peak of 74 records at the true site and a tail of smaller ones around it.

Counting positions would therefore report that isolate as carrying twelve to
eighteen copies. Counting CLUSTERS reports one, which is right. Each cluster is
summarised by its peak position, because the peak is the junction and the tail is
scatter.

THE COMPARISON THIS EXISTS TO MAKE
On an IS-free backbone an isolate's copy number is just its junction-site count,
so this produces a copy-number estimate that shares no machinery with P1d's family
depth. The prediction was stated before the run: site count should track P1d,
about one for lineage 7 and about twenty-four for lineage 2. The element-contig
depth ratio is a third reading from the same alignment -- reads landing on the
element contig against reads landing on the chromosome -- and it should agree too.

Three estimates that agree would settle the copy-number question this arm has
been circling. Three that disagree localise which one is wrong.
"""
import argparse, bisect, csv, glob, os, statistics, sys


def load_crossmap(path):
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


def cluster(positions, radius):
    """Single-linkage on sorted positions. Returns list of lists."""
    out, cur = [], []
    for p in sorted(positions, key=lambda x: x[0]):
        if cur and p[0] - cur[-1][0] > radius:
            out.append(cur); cur = []
        cur.append(p)
    if cur:
        out.append(cur)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default="refbias/p1g")
    ap.add_argument("--crossmap", default="is6110/assets/H37Rv.isclean.crossmap.tsv")
    ap.add_argument("--radius", type=int, default=1000,
                    help="positions within this distance are one insertion site")
    ap.add_argument("--min-peak", type=int, default=10,
                    help="element-contig SA records at a cluster's peak before the "
                         "cluster counts as a site")
    ap.add_argument("--max-readthrough-frac", type=float, default=0.10,
                    help="a position whose reads mostly span it is not a breakpoint")
    ap.add_argument("--family-depth", default="is6110/results/family_depth_h37rv.tsv")
    ap.add_argument("--cohort-lineage", default="refbias/cohort.pilot.tsv",
                    help="supplies the lineage column; the family-depth table "
                         "does not carry one, so without this every lineage "
                         "printed blank and the prediction being tested here is "
                         "stated per lineage")
    ap.add_argument("--normal-dir", default="refbias/p1e",
                    help="junction tables from the NORMAL reference, scored with "
                         "the same filter and clustering, for the sensitivity "
                         "comparison")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sites-out", default=None)
    a = ap.parse_args()

    xmap = load_crossmap(a.crossmap)
    lineage = {}
    if a.cohort_lineage and os.path.exists(a.cohort_lineage):
        for r in csv.DictReader(open(a.cohort_lineage), delimiter="\t"):
            if r.get("sample"):
                lineage[r["sample"]] = r.get("lineage", "") or ""
    p1d = {}
    if os.path.exists(a.family_depth):
        for r in csv.DictReader(open(a.family_depth), delimiter="\t"):
            p1d[r["sample"]] = (float(r["total_copies"]), r.get("lineage", ""))

    def read_hits(path, max_rt):
        """Element-attributed positions surviving the read-through filter.

        Used for BOTH references so the two site counts cannot drift apart on a
        threshold. The normal-reference tables (P1e) carry the same columns in
        original coordinates; the IS-clean tables (P1g) are in clean coordinates
        and are converted only when a site is reported.
        """
        hits = []
        for r in csv.DictReader(open(path), delimiter="\t"):
            if not int(r["is6110"]):
                continue
            # READ-THROUGH FILTER. A real insertion cannot be spanned: reads stop
            # at it. Measured independently on the pilot before this experiment,
            # zero read-through holds for 90% of true IS6110 sites against 31% of
            # non-sites (METHOD_CHOICE.md section 1), so this is a criterion
            # established elsewhere rather than one tuned to this result.
            #
            # It matters here. SAMEA104061150 carries one copy by every other
            # measure and yields two candidate clusters: the true site at the
            # IS6110-11 locus with 161 clips, 72 element records and read-through
            # 0, and one 1.3 kb upstream with 43 clips, 17 records and read-through
            # 84. The second is crossed by most reads and is not a breakpoint.
            clips = int(r["clips_start"]) + int(r["clips_end"])
            rt = int(r["readthrough"])
            if clips + rt and rt / (clips + rt) > max_rt:
                continue
            hits.append((int(r["pos"]), int(r["sa_at_is6110"]),
                         int(r["clips_start"]), int(r["clips_end"])))
        return hits

    def sites_for(hits, radius, min_peak):
        clusters = cluster(hits, radius)
        kept = [(max(cl, key=lambda x: x[1]), cl) for cl in clusters]
        kept = [(pk, cl) for pk, cl in kept if pk[1] >= min_peak]
        return clusters, kept

    rows, sites_rows = [], []
    for path in sorted(glob.glob(os.path.join(a.dir, "*.junctions.tsv"))):
        s = os.path.basename(path).split(".")[0]
        hits = read_hits(path, a.max_readthrough_frac)
        clusters, kept = sites_for(hits, a.radius, a.min_peak)

        # THE SENSITIVITY COMPARISON. The same filter and the same clustering on
        # the normal-reference table for this isolate, so the difference between
        # the two counts is the reference and nothing else.
        n_normal = ""
        if a.normal_dir:
            npath = os.path.join(a.normal_dir, f"{s}.junctions.tsv")
            if os.path.exists(npath):
                _, nkept = sites_for(read_hits(npath, a.max_readthrough_frac),
                                     a.radius, a.min_peak)
                n_normal = len(nkept)
        ed = os.path.join(a.dir, f"{s}.elementdepth.tsv")
        eratio = ""
        if os.path.exists(ed):
            for r in csv.DictReader(open(ed), delimiter="\t"):
                eratio = float(r["element_ratio"])
        fam, lin = p1d.get(s, ("", ""))
        lin = lineage.get(s, "") or lin
        rows.append(dict(sample=s, lineage=lin, n_positions=len(hits),
                         n_clusters=len(clusters), n_sites=len(kept),
                         n_sites_normal_ref=n_normal,
                         element_ratio=eratio, p1d_family_copies=fam))
        for peak, cl in kept:
            sites_rows.append(dict(sample=s, lineage=lin,
                                   clean_pos=peak[0],
                                   orig_pos=clean_to_orig(peak[0], xmap),
                                   element_sa=peak[1],
                                   clips_start=peak[2], clips_end=peak[3],
                                   cluster_positions=len(cl),
                                   cluster_span=cl[-1][0] - cl[0][0]))
    if not rows:
        sys.exit(f"no *.junctions.tsv under {a.dir}")

    rows.sort(key=lambda r: (r["lineage"], r["sample"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    if a.sites_out and sites_rows:
        with open(a.sites_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sites_rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(sites_rows)

    print(f"  {len(rows)} isolates\n")
    print(f"  {'sample':16s} {'lineage':10s} {'sites':>5s} {'normal':>6s} "
          f"{'el.ratio':>8s} {'P1d':>6s}")
    for r in rows:
        fam = f"{r['p1d_family_copies']:.1f}" if r["p1d_family_copies"] != "" else "--"
        er = f"{r['element_ratio']:.2f}" if r["element_ratio"] != "" else "--"
        nr = r["n_sites_normal_ref"]
        print(f"  {r['sample']:16s} {r['lineage']:10s} {r['n_sites']:5d} "
              f"{(str(nr) if nr != '' else '--'):>6s} {er:>8s} {fam:>6s}")

    pairs = [(r["n_sites"], r["p1d_family_copies"]) for r in rows
             if r["p1d_family_copies"] != ""]
    pairs2 = [(r["element_ratio"], r["p1d_family_copies"]) for r in rows
              if r["p1d_family_copies"] != "" and r["element_ratio"] != ""]

    def pearson(v):
        if len(v) < 3:
            return float("nan")
        xa = statistics.mean([x for x, _ in v]); ya = statistics.mean([y for _, y in v])
        num = sum((x - xa) * (y - ya) for x, y in v)
        den = (sum((x - xa) ** 2 for x, _ in v) * sum((y - ya) ** 2 for _, y in v)) ** 0.5
        return num / den if den else float("nan")

    both = [(r["n_sites"], r["n_sites_normal_ref"]) for r in rows
            if r["n_sites_normal_ref"] != ""]
    if both:
        c = sum(x for x, _ in both); n = sum(y for _, y in both)
        gained = sum(1 for x, y in both if x > y)
        lost = sum(1 for x, y in both if x < y)
        print(f"\n  sites on the IS-clean reference : {c}")
        print(f"  sites on the normal reference  : {n}")
        print(f"  isolates where the IS-clean reference reports more : "
              f"{gained} of {len(both)}")
        print(f"  isolates where it reports fewer                   : {lost}")
    print(f"\n  junction sites vs P1d family depth : r = {pearson(pairs):.3f}  (n={len(pairs)})")
    print(f"  element ratio vs P1d family depth  : r = {pearson(pairs2):.3f}  (n={len(pairs2)})")
    print(f"\n  written: {a.out}" + (f"\n           {a.sites_out}" if a.sites_out else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
