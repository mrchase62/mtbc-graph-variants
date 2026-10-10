#!/usr/bin/env python3
"""Does the tatC-helY short-read excision artifact track the sample's
fold-back chimera rate? Report only; existing tables (run from runroot).

Per isolate carrying tatC-helY (repeat_family.isolates.tsv, found == 1):
  foldback      fold-back chimeras per 1,000 reads (foldback_qc/phaseB_rates.tsv)
  junction      tatC-helY deletion-junction reads per 100x flank depth
  dysgu_af      dysgu DEL AF at tatC-helY
  depth_loss    1 - depth inside the element / flanks
  family_junc   junction reads per 100x at all other members, summed,
                excluding 755253 (a real deletion)
Spearman rank correlation of each against foldback, all isolates and within
each sequencing batch (group: setE NextSeq, novel and in_panel HiSeq), since
batches differ in both. Writes out/qc51/foldback_vs_element.tsv.
"""
import collections
import csv

from scipy.stats import spearmanr

N = "../analysis/novel_events"
Q = f"{N}/phaseB/out/qc51"
FB = "../analysis/foldback_qc/phaseB_rates.tsv"


def main():
    fb = {r["sample"]: float(r["inv_chimera_per_1k"]) for r in csv.DictReader(open(FB), delimiter="\t")}
    rows = collections.defaultdict(dict)
    fam = collections.Counter()
    for r in csv.DictReader(open(f"{Q}/repeat_family.isolates.tsv"), delimiter="\t"):
        if r["found"] != "1" or not float(r["flank_depth"] or 0):
            continue
        per100 = 100 * int(r["junction_reads"]) / float(r["flank_depth"])
        s = r["sample"]
        if r["member"] == "tatC-helY":
            rows[s].update(sample=s, group=r["group"], junction=per100, dysgu_af=float(r["dysgu_af"] or 0),
                           depth_loss=1 - float(r["depth_ratio"]))
        elif not r["member"].startswith("755253"):
            fam[s] += per100
    out = []
    for s, d in rows.items():
        if s in fb:
            out.append(dict(d, foldback=fb[s], family_junc=round(fam[s], 2)))
    missing = [s for s in rows if s not in fb]
    with open(f"{Q}/foldback_vs_element.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=["sample", "group", "foldback", "junction", "dysgu_af", "depth_loss",
                                           "family_junc"], delimiter="\t")
        w.writeheader()
        for r in sorted(out, key=lambda x: (x["group"], x["foldback"])):
            w.writerow({k: round(v, 3) if isinstance(v, float) else v for k, v in r.items()})
    print(f"carriers with fold-back rate: {len(out)}; missing: {missing}")
    sets = [("all", out)] + [(g, [r for r in out if r["group"] == g]) for g in sorted({r["group"] for r in out})]
    for name, rs in sets:
        f = [r["foldback"] for r in rs]
        print(f"\n{name}: n={len(rs)}, fold-back per 1k median {sorted(f)[len(f) // 2]:.2f} "
              f"(range {min(f):.2f}-{max(f):.2f})")
        for k in ("junction", "dysgu_af", "depth_loss", "family_junc"):
            v = [r[k] for r in rs]
            rho, p = spearmanr(f, v)
            print(f"  {k:12s} median {sorted(v)[len(v) // 2]:7.3f}   rho {rho:+.2f}  p {p:.3g}")


if __name__ == "__main__":
    main()
