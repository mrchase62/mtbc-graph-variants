#!/usr/bin/env python3
"""How often is one event represented twice in a cohort's merged VCF? Report
only (user, 2026-10-10: measure on scale200, no production code changes).

Input: P5's merged VCF (bin/merge_cohort_vcf.py), whose INFO/CLASS says
which arm a record came from: small (HaplotypeCaller via P4), sv (delly +
dysgu via P4b, CLASS sv and svi), is6110 (the IS6110 arm). One frame
(H37Rv, or node contigs for off-path records).

Per sample, records it carries (GT 1 for ALT 1). Overlaps counted per sample
and per record:
  is6110_ins   SV INS within 200 bp of an IS6110 insertion the same sample
               carries
  is6110_del   SV DEL of 1,200-1,500 bp within 200 bp of an IS6110 site
               (any state): the loss of an element copy the frame carries
  small_sv     SV call within 60 bp of a small indel >= 10 bp the same
               sample carries, sizes within 70% (truvari pctsize)
  small_ge50   small indels of 50 bp or more (in the SV size range by
               definition), and how many have an SV call within 60 bp
Writes <out>.tsv (one row per overlapping pair) and prints a summary.

  overlap.py <merged.vcf.gz> <out prefix>
"""
import bisect
import collections
import csv
import gzip
import re
import sys


def main():
    path, out = sys.argv[1], sys.argv[2]
    recs = collections.defaultdict(lambda: collections.defaultdict(list))  # sample -> class -> [rec]
    is_sites = collections.defaultdict(list)  # chrom -> positions of IS6110 records
    n_rec = collections.Counter()
    samples = []
    for line in gzip.open(path, "rt"):
        if line.startswith("##"):
            continue
        c = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            samples = c[9:]
            continue
        info = dict(x.split("=", 1) if "=" in x else (x, "") for x in c[7].split(";"))
        cls = info.get("CLASS", "")
        alt1 = c[4].split(",")[0]
        pos = int(c[1])
        if cls in ("sv", "svi"):
            kind = "sv_" + info.get("SVTYPE", "?")
            size = abs(int(info.get("SVLEN", "0") or 0))
        elif cls == "is6110":
            kind = "is6110"
            size = 1355
            is_sites[c[0]].append(pos)
        elif cls == "small":
            d = len(alt1) - len(c[3])
            if abs(d) < 10 or alt1.startswith("<"):
                continue
            kind, size = "small", abs(d)
        else:
            continue
        n_rec[kind] += 1
        r = dict(chrom=c[0], pos=pos, end=int(info.get("END", pos)), size=size, id=c[2], kind=kind,
                 sub=(info.get("COMPONENT", "") + ("/" + info["SVTIER"] if "SVTIER" in info else "")) if cls in ("sv", "svi")
                 else info.get("SITECLASS", "") if cls == "is6110" else info.get("REGION", ""))
        for s, f in zip(samples, c[9:]):
            if f.split(":")[0] == "1":
                recs[s][kind].append(r)
    for v in is_sites.values():
        v.sort()
    rows = []
    cnt = collections.Counter()
    carried = collections.Counter()
    for s in samples:
        d = recs[s]
        for k, v in d.items():
            carried[k] += len(v)
        isins = d["is6110"]
        small = d["small"]
        for r in d["sv_INS"]:
            hit = [x for x in isins if x["chrom"] == r["chrom"] and abs(x["pos"] - r["pos"]) <= 200]
            if hit:
                cnt["is6110_ins"] += 1
                rows.append(dict(sample=s, kind="is6110_ins", sv=r["id"], sv_size=r["size"], sv_origin=r["sub"],
                                 other=hit[0]["id"], other_origin=hit[0]["sub"]))
        for r in d["sv_DEL"]:
            if 1200 <= r["size"] <= 1500:
                v = is_sites.get(r["chrom"], [])
                i = bisect.bisect_left(v, r["pos"] - 200)
                if i < len(v) and v[i] <= r["end"] + 200:
                    cnt["is6110_del"] += 1
                    rows.append(dict(sample=s, kind="is6110_del", sv=r["id"], sv_size=r["size"], sv_origin=r["sub"],
                                     other=f"{r['chrom']}:{v[i]}", other_origin=""))
        for kind in ("sv_INS", "sv_DEL"):
            for r in d[kind]:
                hit = [x for x in small if x["chrom"] == r["chrom"] and abs(x["pos"] - r["pos"]) <= 60
                       and r["size"] and min(x["size"], r["size"]) >= 0.7 * max(x["size"], r["size"])]
                if hit:
                    cnt["small_sv"] += 1
                    rows.append(dict(sample=s, kind="small_sv", sv=r["id"], sv_size=r["size"], sv_origin=r["sub"],
                                     other=f"{hit[0]['id'][:60]}({hit[0]['size']})", other_origin=hit[0]["sub"]))
        for x in small:
            if x["size"] >= 50:
                cnt["small_ge50"] += 1
                sv = [r for k in ("sv_INS", "sv_DEL") for r in d[k]
                      if r["chrom"] == x["chrom"] and abs(r["pos"] - x["pos"]) <= 60]
                cnt["small_ge50_with_sv"] += bool(sv)
                rows.append(dict(sample=s, kind="small_ge50", sv=sv[0]["id"] if sv else "", sv_size=sv[0]["size"] if sv else "",
                                 sv_origin=sv[0]["sub"] if sv else "", other=f"{x['id'][:60]}({x['size']})",
                                 other_origin=x["sub"]))
    with open(f"{out}.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=["sample", "kind", "sv", "sv_size", "sv_origin", "other", "other_origin"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"samples {len(samples)}; records (small = indels >= 10 bp): {dict(n_rec)}")
    print(f"sample-level calls carried: {dict(carried)}")
    for k, v in cnt.items():
        print(f"  {k:22s} {v}")
    for kind in ("is6110_ins", "is6110_del", "small_sv"):
        sv = {r["sv"] for r in rows if r["kind"] == kind}
        sizes = sorted(r["sv_size"] for r in rows if r["kind"] == kind)
        if sizes:
            print(f"  {kind}: distinct SV records {len(sv)}; SV size median {sizes[len(sizes) // 2]}, "
                  f"range {sizes[0]}-{sizes[-1]}")

    print("\nby origin (sample-level pairs): kind, SV origin (COMPONENT/SVTIER), other origin (IS6110 SITECLASS or small REGION)")
    for (k, a, b), v in sorted(collections.Counter((r["kind"], r["sv_origin"], r["other_origin"]) for r in rows).items(),
                               key=lambda x: (x[0][0], -x[1])):
        print(f"  {k:12s} {a:28s} {b:14s} {v}")


if __name__ == "__main__":
    main()
