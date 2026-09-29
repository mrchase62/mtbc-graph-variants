#!/usr/bin/env python3
"""Collate P4: what the routing rule actually routed, across isolates.

The point of reporting this is to check the rule is doing what the measurements
said it would, not to score accuracy (real isolates have no truth). Three things
worth watching:

  - the direct arm should carry the large majority, because core sequence is most
    of the genome and core is routed to direct;
  - the inherited half of the composed arm should be substantial in PE/PPE,
    because that recall is the reason composition is used there at all -- and it
    is also where T16 found 89.2% of composition's false positives, so its size
    is the size of the risk being taken;
  - accessory-scale off-path records should be rare, per section 11's 0.9%. A
    large count would mean the near/accessory split is mis-set.
"""
import argparse, collections, csv, os, sys


def rd(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


def _snp_dist(row):
    """snp_distance, or a value that sorts last when it is blank.

    An arm that PINS every isolate to one reference -- the two-reference
    comparison does exactly that -- has no selection distance to report. Three
    separate summaries assumed the field was always an integer, and each one
    died with the same ValueError and took the rest of the chain with it
    through afterok. Fixed in p2_summary.py first, then here, which is the
    argument for fixing a class of defect rather than its instances.
    """
    v = (row.get("snp_distance") or "").strip()
    try:
        return int(v)
    except ValueError:
        return 1 << 30


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--dir", default="refbias/p4")
    ap.add_argument("--out", default="refbias/p4/p4_summary.tsv")
    a = ap.parse_args()

    rows, missing = [], []
    nodes = collections.Counter()
    builds = set()
    for r in rd(a.refmap):
        s = r["sample"]
        p = os.path.join(a.dir, f"{s}.placed.tsv")
        if not os.path.exists(p):
            missing.append(s); continue
        recs = rd(p)
        c = collections.Counter((x["arm"], x["component"], x["region"]) for x in recs)
        for x in recs:
            builds.add(x["build_id"])
            if x["region"] == "off_path_accessory" and x["node"]:
                nodes[x["node"]] += 1
        rows.append(dict(
            sample=s, reference=r["reference"], snp_distance=r["snp_distance"],
            total=len(recs),
            direct_core=c[("direct", "called", "core")],
            comp_called_pe=c[("composed", "called", "pe_ppe")],
            comp_called_masked=c[("composed", "called", "masked")],
            comp_inherited_pe=c[("composed", "inherited", "pe_ppe")],
            comp_inherited_masked=c[("composed", "inherited", "masked")],
            off_near=c[("composed", "called", "off_path_near")],
            off_accessory=c[("composed", "called", "off_path_accessory")],
            named=sum(1 for x in recs
                      if x["region"] == "off_path_accessory" and x["acc_locus"])))

    if not rows:
        print("no P4 output found", file=sys.stderr); return 1
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"  {len(rows)} isolates"
          f"{'; MISSING: ' + ', '.join(missing) if missing else ''}")
    if len(builds) > 1:
        print(f"  WARNING: records from {len(builds)} graph builds: "
              f"{', '.join(sorted(builds))} -- do not merge these", file=sys.stderr)
    else:
        print(f"  all records from graph build {builds.pop() if builds else '?'}")
    print(f"\n  {'isolate':<17s}{'d':>6s}{'total':>7s}{'core':>7s}{'cPE':>5s}"
          f"{'iPE':>6s}{'iMsk':>6s}{'near':>6s}{'acc':>5s}{'named':>6s}")
    for r in sorted(rows, key=_snp_dist):
        print(f"  {r['sample']:<17s}{r['snp_distance']:>6s}{r['total']:>7d}"
              f"{r['direct_core']:>7d}{r['comp_called_pe']:>5d}"
              f"{r['comp_inherited_pe']:>6d}{r['comp_inherited_masked']:>6d}"
              f"{r['off_near']:>6d}{r['off_accessory']:>5d}{r['named']:>6d}")

    tot = sum(r["total"] for r in rows)
    def share(k):
        n = sum(r[k] for r in rows)
        return n, (100 * n / tot if tot else 0)
    print(f"\n  pooled {tot} records across {len(rows)} isolates")
    for k, lab in (("direct_core", "direct, core"),
                   ("comp_called_pe", "composed called, PE/PPE"),
                   ("comp_called_masked", "composed called, masked"),
                   ("comp_inherited_pe", "composed inherited, PE/PPE"),
                   ("comp_inherited_masked", "composed inherited, masked"),
                   ("off_near", "off-path, small-indel interior"),
                   ("off_accessory", "off-path, accessory scale")):
        n, pc = share(k)
        print(f"    {lab:<34s}{n:>8d}  {pc:>5.1f}%")
    ni, _ = share("comp_inherited_pe")
    nm, _ = share("comp_inherited_masked")
    nc, _ = share("comp_called_pe")
    ncm, _ = share("comp_called_masked")
    comp = ni + nm + nc + ncm
    if comp:
        inh_share = 100 * (ni + nm) / comp
        print(f"\n  within the composed arm, the inherited half is "
              f"{inh_share:.1f}% of records")
        print(f"    T16 found it carried 89.2% of composition's core false "
              f"positives. Those two numbers being equal means NO enrichment:")
        print(f"    the inherited half dominates the errors because it dominates")
        print(f"    the records, not because it is worse per record. Composition")
        print(f"    is uniformly wrong, which is why no routing rule cut it.")
    # Denominators matter here. Section 11's 0.9% was a fraction of
    # SAMPLE-vs-R calls; the pooled total below is dominated by direct-arm core
    # records, which the matched arm never produced. Comparing against the pooled
    # total understates it by the size of the direct arm, so report both.
    na, pa = share("off_accessory")
    nn, _ = share("off_near")
    nc2, _ = share("comp_called_pe")
    nm2, _ = share("comp_called_masked")
    matched_calls = na + nn + nc2 + nm2
    print(f"\n  accessory-scale off-path records: {na}")
    print(f"    {pa:.2f}% of all {tot} placed records (includes the direct arm)")
    if matched_calls:
        print(f"    {100*na/matched_calls:.2f}% of the {matched_calls} matched-arm "
              f"calls, which is the denominator section 11's 0.9% used")
    if na + nn:
        print(f"    {100*na/(na+nn):.1f}% of off-path records; section 11 measured "
              f"3.7% beyond 50 bp")
    shared = sum(1 for n, k in nodes.items() if k > 1)
    print(f"  accessory nodes seen in more than one isolate: {shared} of "
          f"{len(nodes)} distinct")
    print(f"    (this is the cross-sample identity the node key exists to give)")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
