#!/usr/bin/env python3
"""Flag genomes carrying an implausible number of SNPs for their sublineage.

WHY THIS IS NOT A GLOBAL THRESHOLD
SNP count against H37Rv is dominated by phylogenetic depth, not by quality: an
M. canettii genome legitimately carries tens of thousands more differences than a
lineage 4 isolate. Comparing every genome to one panel-wide cutoff would flag the
outgroup and miss a contaminated lineage-4 assembly entirely.

So the comparison is WITHIN sublineage, using a robust z-score (median and median
absolute deviation) so that a couple of bad genomes cannot inflate the very
spread they are being judged against. Sublineages too small to estimate a spread
fall back to their major lineage, and genomes with no lineage call are reported
separately rather than silently dropped or silently kept.

An excess of SNPs is the expected signature of contamination, of a mixed
infection that the barcode screen missed, and of systematic assembly error. A
DEFICIT is also reported: unusually few SNPs can mean a collapsed assembly or a
near-duplicate of the reference.

INPUT (audit PGB-8). The counts come from bin/variant_counts.py, which aligns
every candidate directly to H37Rv; the CX333-era counts came from the 2025
412-genome graph, had no producer in bin/, and covered only 254 of the 333
panel genomes. The same screen runs on variant_counts.py's 1 bp indel and
private homopolymer-indel counts (the long-read error checks); --column names
the measurement in the output. With --accessions every candidate must have a
count: a genome with no count was not screened, and the script stops rather
than pass it.

NO SPREAD. A group whose MAD is 0 cannot give a z-score. Such a sublineage
falls back to its major lineage; if that also has MAD 0 the z-score is left
blank and a genome off the group median is flagged NO_SPREAD for review,
rather than given z = 0 and passed.
"""
import argparse, collections, csv, re, statistics


def robust_z(x, med, mad):
    if mad <= 0:
        return None
    return 0.6745 * (x - med) / mad


def med_mad(vals):
    med = statistics.median(vals)
    return med, statistics.median([abs(v - med) for v in vals])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True,
                    help="accession<TAB>count (variant_counts.py counts.*.tsv)")
    ap.add_argument("--accessions",
                    help="the candidates; every one must have a count")
    ap.add_argument("--column", default="snps",
                    help="name of the measurement column in the output")
    ap.add_argument("--lineages", required=True)
    ap.add_argument("--min-group", type=int, default=5,
                    help="sublineages smaller than this fall back to the major lineage")
    ap.add_argument("--z", type=float, default=5.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cnt = {}
    for line in open(a.counts):
        f = line.rstrip("\n").split("\t")
        if len(f) == 2:
            cnt[f[0]] = int(f[1])
    lin = {r["accession"]: (r["strain"] or "") for r in
           csv.DictReader(open(a.lineages), delimiter="\t")}
    if a.accessions:
        want = [l.strip() for l in open(a.accessions) if l.strip()]
        miss = [x for x in want if x not in cnt]
        if miss:
            raise SystemExit(f"{len(miss)} candidate(s) have no count in "
                             f"{a.counts} (not screened is not passed): "
                             f"{' '.join(miss[:10])}")
        cnt = {x: cnt[x] for x in want}

    def sub(x):
        return (lin.get(x, "") or "").split(":")[0] or "no_call"

    def major(x):
        m = re.match(r"(lineage\d+|La\d+|M\.\w+)", lin.get(x, "") or "")
        return m.group(1) if m else "no_call"

    by_sub = collections.defaultdict(list)
    by_maj = collections.defaultdict(list)
    for g in cnt:
        by_sub[sub(g)].append(g)
        by_maj[major(g)].append(g)

    rows = []
    for g, c in cnt.items():
        s, m = sub(g), major(g)
        grp, level = (by_sub[s], f"sublineage {s}") if len(by_sub[s]) >= a.min_group \
            else (by_maj[m], f"major {m}")
        med, mad = med_mad([cnt[x] for x in grp])
        if mad <= 0 and level.startswith("sublineage"):
            grp, level = by_maj[m], f"major {m} (sublineage MAD 0)"
            med, mad = med_mad([cnt[x] for x in grp])
        z = robust_z(c, med, mad)
        if z is None:
            flag = "NO_SPREAD" if c != med else ""
        else:
            flag = "EXCESS" if z >= a.z else ("DEFICIT" if z <= -a.z else "")
        rows.append({"accession": g, "sublineage": s, "major": m,
                     a.column: c, "group": level, "n_group": len(grp),
                     "group_median": int(med), "group_mad": int(mad),
                     "robust_z": "" if z is None else round(z, 2),
                     "flag": flag})

    rows.sort(key=lambda r: -(r["robust_z"] if r["robust_z"] != "" else 0))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"[snp_outlier] {len(rows)} genomes; robust-z threshold {a.z}")
    col = a.column
    print(f"  panel {col}: median {statistics.median([r[col] for r in rows]):,}"
          f"  min {min(r[col] for r in rows):,}  max {max(r[col] for r in rows):,}")
    for f in ("EXCESS", "DEFICIT", "NO_SPREAD"):
        sel = [r for r in rows if r["flag"] == f]
        print(f"\n  {f}: {len(sel)}")
        if sel:
            print(f"    {'accession':<16s} {'sublineage':<18s} {'snps':>8s} "
                  f"{'median':>8s} {'mad':>6s} {'z':>7s} {'n':>4s}")
            for r in (sel if f == "EXCESS" else sel)[:14]:
                print(f"    {r['accession']:<16s} {r['sublineage']:<18s} {r[col]:>8,} "
                      f"{r['group_median']:>8,} {r['group_mad']:>6,} {str(r['robust_z']):>7s} "
                      f"{r['n_group']:>4d}")
    print(f"\n  -> {a.out}")


if __name__ == "__main__":
    main()
