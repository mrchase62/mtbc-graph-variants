#!/usr/bin/env python3
"""Part D (../FOLDBACK_SV_EFFECT_PLAN.md): do fold-back chimeras hide real
SVs? Report only. Run from runroot with $MTB_PY after 02_post.sbatch.

Scores (novel_events/score.py, callable sequence, strict precision) the
fbrecall calls against out/truth.tsv, pooled by fold-back group:
  kept      selfcal54 isolates with v2 rate <= 4.9 (the cut-off, old <= 2.5)
  excluded  selfcal54 isolates with v2 rate > 4.9
  chimeric  the 9 chimeric Phase B controls (all reads)
  chimeric_fb100  the same 9 with every fold-back pair removed (cohort
            fbrecall100; same reference and truth as the all-reads sample),
            also compared isolate by isolate with chimeric (the causal arm)
and per isolate. Rates are the own-genome v2 rates (selfcal/out/calib_points.tsv).
Callers: proto_asm (our caller v6 + assembly), dysgu_af20, and the current
combination proto_asm+dysgu_af20. Writes out/recall_by_group.tsv and
out/recall_by_isolate.tsv and prints a summary.
"""
import collections
import csv
import os
import subprocess
import sys

from scipy.stats import spearmanr

D = "../analysis/foldback_qc/recall"
N = "../analysis/novel_events"
RM = "refbias/fbrecall/run/p1/refmap.tsv"
RM100 = "refbias/fbrecall100/run/p1/refmap.tsv"
CUT = 4.9
CALLERS = ("proto_asm", "dysgu_af20", "proto_asm+dysgu_af20")


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def score(samples, truth, tag, p2="refbias/fbrecall/run/p2", refmap=RM):
    t = f"{D}/out/score_work/{tag}.truth.tsv"
    os.makedirs(os.path.dirname(t), exist_ok=True)
    with open(t, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(truth[0]), delimiter="\t")
        w.writeheader()
        w.writerows(r for r in truth if r["sample"] in samples)
    out = f"{D}/out/score_work/{tag}.score.tsv"
    subprocess.run([sys.executable, f"{N}/score.py", "--truth", t, "--proto-dir", f"{D}/out/proto",
                    "--asm-dir", f"{D}/out/asm", "--p2-dir", p2,
                    "--refmap", refmap, "--mask-dir", f"{N}/out/masks", "--out", out],
                   check=True, stdout=subprocess.DEVNULL)
    return {r["caller"]: r for r in tsv(out) if r["stratum"] == "all" and r["caller"] in CALLERS}


def main():
    truth = tsv(f"{D}/out/truth.tsv")
    rate = {p["sample"]: float(p["rate"]) for p in tsv("../analysis/foldback_qc/selfcal/out/calib_points.tsv")}
    samples = [r["sample"] for r in tsv(RM)]
    grp = {s: ("chimeric" if s.startswith("mar_") else "kept" if rate[s] <= CUT else "excluded") for s in samples}
    ntrue = collections.Counter(r["sample"] for r in truth)

    rows = []
    for g in ("kept", "excluded", "chimeric"):
        ss = {s for s in samples if grp[s] == g}
        sc = score(ss, truth, g)
        for c in CALLERS:
            r = sc.get(c, {})
            rows.append(dict(group=g, isolates=len(ss), with_truth=sum(1 for s in ss if ntrue[s]),
                             caller=c, n_true=r.get("n_true", 0), typed_recall=r.get("typed_recall", ""),
                             breakpoint_recall=r.get("breakpoint_recall", ""),
                             typed_calls=r.get("typed_calls", ""), precision=r.get("precision", ""),
                             precision_lenient=r.get("precision_lenient", "")))
    chim = sorted(s for s in samples if grp[s] == "chimeric")
    t100 = [dict(r, sample=r["sample"] + "_fb100") for r in truth if r["sample"] in chim]
    f100 = dict(p2="refbias/fbrecall100/run/p2", refmap=RM100)
    sc = score({s + "_fb100" for s in chim}, t100, "chimeric_fb100", **f100)
    for c in CALLERS:
        r = sc.get(c, {})
        rows.append(dict(group="chimeric_fb100", isolates=len(chim),
                         with_truth=sum(1 for s in chim if ntrue[s]),
                         caller=c, n_true=r.get("n_true", 0), typed_recall=r.get("typed_recall", ""),
                         breakpoint_recall=r.get("breakpoint_recall", ""),
                         typed_calls=r.get("typed_calls", ""), precision=r.get("precision", ""),
                         precision_lenient=r.get("precision_lenient", "")))
    with open(f"{D}/out/recall_by_group.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("group\tisolates\twith_truth\tcaller\tn_true\ttyped_recall\tbreakpoint_recall\ttyped_calls\tprecision\tprecision_lenient")
    for r in rows:
        print("\t".join(str(v) for v in r.values()))

    per = []
    for s in samples:
        if not ntrue[s]:
            continue
        sc = score({s}, truth, s)
        r = sc.get("proto_asm+dysgu_af20", {})
        per.append(dict(sample=s, group=grp[s], rate=rate[s], n_true=ntrue[s],
                        typed_recall=float(r.get("typed_recall") or 0),
                        precision=r.get("precision", ""), typed_calls=r.get("typed_calls", "")))
    with open(f"{D}/out/recall_by_isolate.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(per[0]), delimiter="\t")
        w.writeheader()
        w.writerows(per)
    print("\nCausal arm, isolate by isolate (all reads -> fold-back pairs removed):")
    print("isolate\tcallable_true\tcaller\trecall_all\trecall_fb100\tcalls_all\tcalls_fb100\tprec_all\tprec_fb100")
    for s in chim:
        if not ntrue[s]:
            continue
        a = score({s}, truth, s)
        b = score({s + "_fb100"}, t100, s + "_fb100", **f100)
        for c in CALLERS:
            x, y = a.get(c, {}), b.get(c, {})
            print(f"{s}\t{x.get('n_true')}\t{c}\t{x.get('typed_recall')}\t{y.get('typed_recall')}\t"
                  f"{x.get('typed_calls')}\t{y.get('typed_calls')}\t{x.get('precision')}\t{y.get('precision')}")
    rho, p = spearmanr([x["rate"] for x in per], [x["typed_recall"] for x in per])
    print(f"\nper isolate (n {len(per)} with truth): Spearman recall vs rate {rho:+.2f} (p {p:.2g})")
    print(f"isolates without truth events (no recall possible): {sum(1 for s in samples if not ntrue[s])}")


if __name__ == "__main__":
    main()
