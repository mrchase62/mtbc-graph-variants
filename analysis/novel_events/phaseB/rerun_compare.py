#!/usr/bin/env python3
"""Before/after for the Phase B rerun on clean data (run from runroot after
11_call_score.sbatch). Scores both runs on the SAME isolates (those in the
marinQC51 refmap), callable sequence only:

  before  raw reads (marinB63 P2; our caller + assembly from out/v3b)
  after   reads QC'd, fold-back <= 2 per 1,000, fixed P2 (marinQC51; out/qc51)

Isolates whose P1 reference changed between the runs are left out of both:
the truth set is in the coordinates of the reference chosen before.

For novel and setE_polished: score.py's typed recall and precision per
caller. For the in-panel controls (no callable true events): typed calls in
callable sequence per isolate, all of them false by construction, counted
directly because score.py skips isolates without truth rows.

Writes out/qc51/compare.tsv and prints a summary.
"""
import collections
import csv
import os
import subprocess
import sys

N = "../analysis/novel_events"
sys.path.insert(0, N)
from score import load_masks, masked, read_proto, read_vcf  # noqa: E402

RM_AFTER = "refbias/marinQC51/run/p1/refmap.tsv"
RM_BEFORE = "refbias/marinB63/run/p1/refmap.tsv"
CALLERS = ("proto_asm", "dysgu_pass", "proto_asm+dysgu_pass")
RUNS = {"before": dict(proto=f"{N}/phaseB/out/v3b/proto", asm=f"{N}/phaseB/out/v3b/asm",
                       p2="refbias/marinB63/run/p2", rm=RM_BEFORE),
        "after": dict(proto=f"{N}/phaseB/out/qc51/proto", asm=f"{N}/phaseB/out/qc51/asm",
                      p2="refbias/marinQC51/run/p2", rm=RM_AFTER)}


def main():
    keep = [r["sample"] for r in csv.DictReader(open(RM_AFTER), delimiter="\t")]
    groups = {r["sample"]: r["group"] for r in csv.DictReader(open(f"{N}/phaseB/out/samples.tsv"), delimiter="\t")}
    ref_b = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM_BEFORE), delimiter="\t")}
    ref_a = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM_AFTER), delimiter="\t")}
    changed = [s for s in keep if ref_a[s] != ref_b.get(s)]
    keep = [s for s in keep if s not in changed]
    out_dir = f"{N}/phaseB/out/qc51"
    rows = []
    for grp in ("novel", "setE_polished"):
        iso = {s for s in keep if groups[s] == grp}
        t = f"{out_dir}/truth_same_ref.{grp}.tsv"
        with open(f"{out_dir}/truth.{grp}.tsv") as fi, open(t, "w") as fo:  # 11_call_score.sbatch's truth
            for k, line in enumerate(fi):
                if k == 0 or line.split("\t", 1)[0] in iso:
                    fo.write(line)
        for run, p in RUNS.items():
            o = f"{out_dir}/score_{run}.{grp}.tsv"
            subprocess.run([sys.executable, f"{N}/score.py", "--truth", t, "--proto-dir", p["proto"],
                            "--asm-dir", p["asm"], "--p2-dir", p["p2"], "--refmap", p["rm"],
                            "--mask-dir", f"{N}/out/masks", "--out", o],
                           check=True, stdout=subprocess.DEVNULL)
            for r in csv.DictReader(open(o), delimiter="\t"):
                if r["stratum"] == "all" and r["caller"] in CALLERS:
                    rows.append(dict(group=grp, run=run, isolates=len(iso), caller=r["caller"],
                                     n_true=r["n_true"], typed_recall=r["typed_recall"],
                                     breakpoint_recall=r["breakpoint_recall"],
                                     typed_calls=r["typed_calls"], precision=r["precision"]))
    ctl = [s for s in keep if groups[s] == "in_panel"]
    for run, p in RUNS.items():
        masks = load_masks(p["rm"], f"{N}/out/masks")
        per = collections.defaultdict(list)
        for s in ctl:
            calls = dict(proto_asm=read_proto(f"{p['asm']}/{s}.events.tsv"),
                         dysgu_pass=read_vcf(f"{p['p2']}/{s}.dysgu.vcf", True))
            for c, cs in calls.items():
                per[c].append(sum(1 for x in cs if x["type"] != "BND" and not masked(masks[s], x["start"], x["end"])))
        for c, v in per.items():
            rows.append(dict(group="in_panel", run=run, isolates=len(ctl), caller=c, n_true=0,
                             typed_calls=sum(v), precision="", typed_recall="", breakpoint_recall="",
                             per_isolate=",".join(map(str, v))))
    with open(f"{out_dir}/compare.tsv", "w") as fo:
        cols = ["group", "run", "isolates", "caller", "n_true", "typed_recall", "breakpoint_recall",
                "typed_calls", "precision", "per_isolate"]
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"isolates compared {len(keep)}; left out, P1 reference changed: {', '.join(changed) or 'none'}")
    for r in rows:
        print("\t".join(str(r.get(k, "")) for k in ("group", "run", "caller", "n_true", "typed_recall",
                                                       "typed_calls", "precision", "per_isolate")))


if __name__ == "__main__":
    main()
