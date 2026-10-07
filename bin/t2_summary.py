#!/usr/bin/env python3
"""Does the graph's carrier set correspond to genomes that actually contain the sequence?

Not "is the sequence present" -- at a multiallelic locus every genome carries
some version, so non-carriers hit a similar sequence too. The question is whether
identity and coverage in graph-carriers are materially higher than in
graph-non-carriers. An allele whose carriers and non-carriers look the same is
not distinguishing anything and must not enter the accessory panel.
"""
import argparse, csv, glob, gzip, os, statistics, sys


def main():
    ap = argparse.ArgumentParser()
    # no defaults: refbias/t2/ and refbias/T2.accessory_validation.tsv are
    # the pilot's, CX333's
    ap.add_argument("--hits", required=True,
                    help="t2_validate_accessory.sh's OUT")
    ap.add_argument("--meta", required=True,
                    help="t2_extract_candidates.py's --out-meta")
    ap.add_argument("--min-ident", type=float, default=99.0)
    ap.add_argument("--min-cov", type=float, default=0.95)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    meta, carriers = {}, {}
    for r in csv.DictReader(open(a.meta), delimiter="\t"):
        meta[r["allele_id"]] = (int(r["pos"]), int(r["len"]), int(r["AC"]))
        carriers[r["allele_id"]] = set(r["carriers"].split(","))

    # allele -> sample -> (pident, cov)
    hit = {k: {} for k in meta}
    nsamp = 0
    for p in sorted(glob.glob(os.path.join(a.hits, "*.hits.tsv.gz"))):
        nsamp += 1
        for r in csv.DictReader(gzip.open(p, "rt"), delimiter="\t"):
            if r["allele_id"] in hit:
                hit[r["allele_id"]][r["sample"]] = (float(r["pident"]),
                                                    float(r["cov"]))
    print(f"  {nsamp} genomes, {len(meta)} candidate alleles")

    allsamp = set()
    for v in hit.values():
        allsamp |= set(v)

    rows = []
    for aid, (pos, ln, ac) in meta.items():
        car, non = carriers[aid], allsamp - carriers[aid]
        def stats(group):
            vals = [hit[aid].get(s) for s in group]
            full = [1 if (v and v[0] >= a.min_ident and v[1] >= a.min_cov) else 0
                    for v in vals]
            idents = [v[0] for v in vals if v]
            return (statistics.median(idents) if idents else 0.0,
                    sum(full) / len(full) if full else 0.0)
        ci, cf = stats(car)
        ni, nf = stats(non) if non else (0.0, 0.0)
        rows.append(dict(allele_id=aid, pos=pos, len=ln, AC=ac,
                         carrier_median_ident=round(ci, 3),
                         carrier_frac_fulllength=round(cf, 4),
                         noncarrier_median_ident=round(ni, 3),
                         noncarrier_frac_fulllength=round(nf, 4),
                         separation=round(cf - nf, 4)))

    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: -r["separation"]))

    confirmed = [r for r in rows if r["carrier_frac_fulllength"] >= 0.90
                 and r["separation"] >= 0.50]
    present_all = [r for r in rows if r["carrier_frac_fulllength"] >= 0.90
                   and r["separation"] < 0.50]
    absent = [r for r in rows if r["carrier_frac_fulllength"] < 0.90]
    mass = sum(r["AC"] for r in rows)
    print(f"\n  {'class':<44}{'alleles':>8}{'copy mass':>12}")
    for lab, sel in (("CONFIRMED: in carriers, separates from non-carriers", confirmed),
                     ("present in carriers but NOT distinguishing", present_all),
                     ("NOT recovered at full length in carriers", absent)):
        m = sum(r["AC"] for r in sel)
        print(f"  {lab:<44}{len(sel):>8}{m:>12,}")
    print(f"  {'total':<44}{len(rows):>8}{mass:>12,}")
    print(f"\n  thresholds: identity >= {a.min_ident}%, coverage >= {a.min_cov}")
    print(f"  written: {a.out}")
    if absent:
        print(f"\n  worst — not recovered in their own carriers:")
        for r in sorted(absent, key=lambda x: x["carrier_frac_fulllength"])[:10]:
            print(f"    {r['allele_id']:<22} AC {r['AC']:>4}  "
                  f"carrier full-length {r['carrier_frac_fulllength']:.2f}  "
                  f"median ident {r['carrier_median_ident']:.1f}")


if __name__ == "__main__":
    sys.exit(main())
