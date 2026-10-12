#!/usr/bin/env python3
"""Excluded-isolates check (user, 2026-10-12): does removing every fold-back
read pair bring the excluded selfcal54 isolates (v2 > 4.9) to the kept
group's precision, without losing recall? If so, filtering reads could
replace most of the sample exclusion. Report only. Run from runroot with
$MTB_PY after 05_post_excl.sbatch.

Groups, scored as analyse.py (score.py, callable sequence, strict precision;
truth = own genome vs the part D reference, out/truth.tsv):
  kept           v2 <= 4.9, all reads (fbrecall)
  excluded       v2 > 4.9, all reads (fbrecall)
  excluded_fb100 the same 12 with every fold-back pair removed (fbrecallx100)
Each at two settings of our caller: no AF floor (as part D) and the
production floor 0.1 (assemble_candidates.af_floor, events without a read
count kept); dysgu at AF >= 0.2 throughout. Then isolate by isolate,
excluded vs excluded_fb100 (combined caller, floor 0.1).
Writes out/excl_by_group.tsv and prints a summary.
"""
import csv
import os
import subprocess
import sys

import pysam

sys.path.insert(0, "../analysis/novel_events")
from assemble_candidates import af_floor  # noqa: E402

D = "../analysis/foldback_qc/recall"
N = "../analysis/novel_events"
RM = "refbias/fbrecall/run/p1/refmap.tsv"
RMX = "refbias/fbrecallx100/run/p1/refmap.tsv"
CUT = 4.9
CALLERS = ("proto_asm", "dysgu_af20", "proto_asm+dysgu_af20")


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def floored(samples, coh, floor):
    """asm events of samples with the AF floor applied, in a folder of their own."""
    if not floor:
        return f"{D}/out/asm"
    od = f"{D}/out/asm_af/excl_{floor}"
    os.makedirs(od, exist_ok=True)
    for s in samples:
        ev = tsv(f"{D}/out/asm/{s}.events.tsv")
        bam = pysam.AlignmentFile(f"refbias/work/{coh}_p2/{s}.bam")
        keep = af_floor(ev, bam, bam.references[0], floor)
        cols = list(ev[0]) if ev else "sample type start end size support left right depth_ratio source".split()
        with open(f"{od}/{s}.events.tsv", "w") as fo:
            w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(keep)
    return od


def score(samples, truth, tag, asm, p2, refmap):
    t = f"{D}/out/score_work/x_{tag}.truth.tsv"
    os.makedirs(os.path.dirname(t), exist_ok=True)
    with open(t, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(truth[0]), delimiter="\t")
        w.writeheader()
        w.writerows(r for r in truth if r["sample"] in samples)
    out = f"{D}/out/score_work/x_{tag}.score.tsv"
    subprocess.run([sys.executable, f"{N}/score.py", "--truth", t, "--proto-dir", f"{D}/out/proto",
                    "--asm-dir", asm, "--p2-dir", p2, "--refmap", refmap,
                    "--mask-dir", f"{N}/out/masks", "--out", out], check=True, stdout=subprocess.DEVNULL)
    return {r["caller"]: r for r in tsv(out) if r["stratum"] == "all" and r["caller"] in CALLERS}


def main():
    truth = tsv(f"{D}/out/truth.tsv")
    rate = {p["sample"]: float(p["rate"]) for p in tsv("../analysis/foldback_qc/selfcal/out/calib_points.tsv")}
    ps = [r["sample"] for r in tsv(RM) if r["sample"].startswith("ps_")]
    kept = {s for s in ps if rate[s] <= CUT}
    excl = sorted(s for s in ps if rate[s] > CUT)
    tx = [dict(r, sample=r["sample"] + "_fb100") for r in truth if r["sample"] in excl]
    arms = {"kept": (kept, truth, "fbrecall", RM),
            "excluded": (set(excl), truth, "fbrecall", RM),
            "excluded_fb100": ({s + "_fb100" for s in excl}, tx, "fbrecallx100", RMX)}
    rows = []
    for floor in (0, 0.1):
        for g, (ss, tr, coh, rm) in arms.items():
            asm = floored(ss, coh, floor)
            sc = score(ss, tr, f"{g}_{floor}", asm, f"refbias/{coh}/run/p2", rm)
            for c in CALLERS:
                r = sc.get(c, {})
                rows.append(dict(floor=floor, group=g, isolates=len(ss), caller=c, n_true=r.get("n_true", 0),
                                 typed_recall=r.get("typed_recall", ""), typed_calls=r.get("typed_calls", ""),
                                 precision=r.get("precision", ""), precision_lenient=r.get("precision_lenient", "")))
    with open(f"{D}/out/excl_by_group.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("floor\tgroup\tisolates\tcaller\tn_true\ttyped_recall\ttyped_calls\tprecision\tprecision_lenient")
    for r in rows:
        print("\t".join(str(v) for v in r.values()))

    print("\nIsolate by isolate, combined caller, floor 0.1 (all reads -> fold-back pairs removed)")
    print("isolate\tv2_rate\tcallable_true\trecall_all\trecall_fb100\tcalls_all\tcalls_fb100\tprec_all\tprec_fb100")
    a_asm, b_asm = f"{D}/out/asm_af/excl_0.1", f"{D}/out/asm_af/excl_0.1"
    for s in excl:
        a = score({s}, truth, s + "_x", a_asm, "refbias/fbrecall/run/p2", RM).get("proto_asm+dysgu_af20", {})
        b = score({s + "_fb100"}, tx, s + "_fb100_x", b_asm, "refbias/fbrecallx100/run/p2",
                  RMX).get("proto_asm+dysgu_af20", {})
        print(f"{s}\t{rate[s]:.1f}\t{a.get('n_true')}\t{a.get('typed_recall')}\t{b.get('typed_recall')}\t"
              f"{a.get('typed_calls')}\t{b.get('typed_calls')}\t{a.get('precision')}\t{b.get('precision')}")


if __name__ == "__main__":
    main()
