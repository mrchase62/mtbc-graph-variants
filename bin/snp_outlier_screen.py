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
"""
import argparse, collections, csv, re, statistics


def robust_z(x, med, mad):
    if mad <= 0:
        return 0.0
    return 0.6745 * (x - med) / mad


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="accession<TAB>snp_count")
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
        vals = [cnt[x] for x in grp]
        med = statistics.median(vals)
        mad = statistics.median([abs(v - med) for v in vals])
        z = robust_z(c, med, mad)
        rows.append({"accession": g, "sublineage": s, "major": m,
                     "snps": c, "group": level, "n_group": len(grp),
                     "group_median": int(med), "group_mad": int(mad),
                     "robust_z": round(z, 2),
                     "flag": "EXCESS" if z >= a.z else ("DEFICIT" if z <= -a.z else "")})

    rows.sort(key=lambda r: -r["robust_z"])
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"[snp_outlier] {len(rows)} genomes; robust-z threshold {a.z}")
    print(f"  panel SNP counts: median {statistics.median([r['snps'] for r in rows]):,}"
          f"  min {min(r['snps'] for r in rows):,}  max {max(r['snps'] for r in rows):,}")
    for f in ("EXCESS", "DEFICIT"):
        sel = [r for r in rows if r["flag"] == f]
        print(f"\n  {f}: {len(sel)}")
        if sel:
            print(f"    {'accession':<16s} {'sublineage':<18s} {'snps':>8s} "
                  f"{'median':>8s} {'mad':>6s} {'z':>7s} {'n':>4s}")
            for r in (sel if f == "EXCESS" else sel)[:14]:
                print(f"    {r['accession']:<16s} {r['sublineage']:<18s} {r['snps']:>8,} "
                      f"{r['group_median']:>8,} {r['group_mad']:>6,} {r['robust_z']:>7.1f} "
                      f"{r['n_group']:>4d}")
    print(f"\n  -> {a.out}")


if __name__ == "__main__":
    main()
