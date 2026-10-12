#!/usr/bin/env python3
"""Real-data check of our caller wired into production (user, 2026-10-12).
Run from runroot with $MTB_PY after 01_fbrecall_check.sbatch.

1. Reproduction. P2's new step against the evaluation outputs it should
   reproduce, sample by sample:
     caller   work/fbrecall_p2/<s>.ours.proto.tsv  vs recall/out/proto
     events   run/p2/<s>.ours.events.tsv           vs recall/out/asm_af/keep_na_0.1
   (same BAMs, settings and masks), and the two masks P2 built afresh
   against novel_events/out/masks.
2. P4b scored. The P4b rows a production user would get from P2's calls
   (component called, FILTER PASS: P5's --keep-filter) are written as an
   events table in R coordinates and scored by score.py (callable sequence,
   strict precision) per group, beside the evaluation's combined caller with
   the same rules (proto_asm+dysgu_af20, --hc-dir, --element-gt from P2's
   genotypes). P4b merges the two callers' calls of one event into one row,
   so it has fewer calls than the evaluation's sum.
Writes out/compare_report.txt (printed) and out/score_*.tsv.
"""
import csv
import filecmp
import os
import subprocess
import sys

N = "../analysis/novel_events"
R = "../analysis/foldback_qc/recall"
O = "../analysis/sv_wiring/out"
RM = "refbias/fbrecall/run/p1/refmap.tsv"
P2 = "refbias/fbrecall/run/p2"
P4B = "refbias/fbrecall/run/p4b_test"
CUT = 4.9


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def same_rows(a, b):
    """Equal as tables (row order and columns as written)."""
    return open(a).read() == open(b).read()


def score(samples, tag, asm_dir, extra):
    t = f"{O}/truth.{tag}.tsv"
    truth = tsv(f"{R}/out/truth.tsv")
    with open(t, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(truth[0]), delimiter="\t")
        w.writeheader()
        w.writerows(r for r in truth if r["sample"] in samples)
    out = f"{O}/score_{tag}.tsv"
    subprocess.run([sys.executable, f"{N}/score.py", "--truth", t, "--proto-dir", asm_dir, "--asm-dir", asm_dir,
                    "--p2-dir", P2, "--refmap", RM, "--mask-dir", f"{N}/out/masks", "--out", out] + extra,
                   check=True, stdout=subprocess.DEVNULL)
    return {r["caller"]: r for r in tsv(out) if r["stratum"] == "all"}


def main():
    os.makedirs(f"{O}/p4b_events", exist_ok=True)
    samples = [r["sample"] for r in tsv(RM)]
    lines = []

    def say(x=""):
        print(x)
        lines.append(x)

    say("1. Reproduction of the evaluation outputs")
    n_p = n_e = 0
    diff = []
    for s in samples:
        a = same_rows(f"refbias/work/fbrecall_p2/{s}.ours.proto.tsv", f"{R}/out/proto/{s}.events.tsv")
        b = same_rows(f"{P2}/{s}.ours.events.tsv", f"{R}/out/asm_af/keep_na_0.1/{s}.events.tsv")
        n_p += a
        n_e += b
        if not (a and b):
            diff.append(s)
    say(f"  caller tables identical: {n_p} of {len(samples)}; events after assembly and AF floor: {n_e} of {len(samples)}")
    if diff:
        say(f"  differing: {', '.join(diff)}")
    for g in ("GCF_002357935", "GCF_002357955"):
        new = f"refbias/callable_masks/7713a8d71d8e-fix1/{g}.bed"
        say(f"  mask {g} built in P2: {'identical' if filecmp.cmp(new, f'{N}/out/masks/{g}.bed', shallow=False) else 'DIFFERENT'}")
    gts = [tsv(f"{P2}/{s}.element_gt.tsv")[0] for s in samples]
    with open(f"{O}/element_gt.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(gts[0]), delimiter="\t")
        w.writeheader()
        w.writerows(gts)
    n = {}
    for g in gts:
        n[g["genotype"]] = n.get(g["genotype"], 0) + 1
    say(f"  tatC-helY genotypes: {n}")

    # P4b rows as events tables
    rule_log = {}
    for s in samples:
        rows = [r for r in tsv(f"{P4B}/{s}.sv_placed.tsv") if r["component"] == "called" and r["filter"] == "PASS"]
        with open(f"{O}/p4b_events/{s}.events.tsv", "w") as fo:
            fo.write("sample\ttype\tstart\tend\tsize\tsupport\tleft\tright\tdepth_ratio\tsource\n")
            for r in rows:
                fo.write(f"{s}\t{r['svtype']}\t{r['r_pos']}\t{r['r_end']}\t{r['svlen']}\t\t\t\t\t{r['src']}\n")
        rule_log[s] = rows

    rate = {p["sample"]: float(p["rate"]) for p in tsv("../analysis/foldback_qc/selfcal/out/calib_points.tsv")}
    groups = {
        "kept, without GCF_035581225 (26)": [s for s in samples if s.startswith("ps_") and rate[s] <= CUT
                                             and s != "ps_GCF_035581225"],
        "excluded (12)": [s for s in samples if s.startswith("ps_") and rate[s] > CUT],
        "chimeric controls (9)": [s for s in samples if s.startswith("mar_")],
    }
    say("\n2. P4b rows (called, PASS) scored, against the evaluation's combined caller with the same rules")
    say("group\tset\ttyped_recall\ttyped_calls\tprecision_strict\tprecision_lenient")
    for g, ss in groups.items():
        tag = g.split()[0].rstrip(",")
        ev = score(set(ss), f"eval_{tag}", f"{R}/out/asm_af/keep_na_0.1",
                   ["--hc-dir", P2, "--element-gt", f"{O}/element_gt.tsv"])["proto_asm+dysgu_af20"]
        pb = score(set(ss), f"p4b_{tag}", f"{O}/p4b_events", [])["proto_asm"]
        for name, r in (("evaluation", ev), ("P4b", pb)):
            say(f"{g}\t{name}\t{r['typed_recall']}\t{r['typed_calls']}\t{r['precision']}\t{r['precision_lenient']}")
    with open(f"{O}/compare_report.txt", "w") as fo:
        fo.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
