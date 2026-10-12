#!/usr/bin/env python3
"""AF floor on our caller (parts A-B suggested 0.2): what it costs in recall
and gains in precision, on part D's truth (report only). Run from runroot
with $MTB_PY after analyse.py.

Our caller's AF, as ../fb_effect.py: supporting reads / reads overlapping the
start position in the P2 BAM. Events with support 0 (typed by assembly, no
read count) have no AF; the floor is applied two ways:
  keep_na  events without an AF are kept (the floor acts only where measured)
  drop_na  events without an AF are dropped
dysgu stays at AF >= 0.2. Groups as analyse.py (kept, excluded, chimeric,
chimeric_fb100). Scoring: score.py as in analyse.py (callable, strict).
Writes out/affloor.tsv.
"""
import csv
import os
import subprocess
import sys

import pysam

D = "../analysis/foldback_qc/recall"
N = "../analysis/novel_events"
CUT = 4.9
FLOORS = (0.0, 0.1, 0.15, 0.2, 0.25)
CALLERS = ("proto_asm", "proto_asm+dysgu_af20")


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def main():
    rate = {p["sample"]: float(p["rate"]) for p in tsv("../analysis/foldback_qc/selfcal/out/calib_points.tsv")}
    cohorts = {"fbrecall": tsv("refbias/fbrecall/run/p1/refmap.tsv"),
               "fbrecall100": tsv("refbias/fbrecall100/run/p1/refmap.tsv")}
    af = {}
    for coh, rows in cohorts.items():
        for r in rows:
            s = r["sample"]
            bam = pysam.AlignmentFile(f"refbias/work/{coh}_p2/{s}.bam")
            ctg = bam.references[0]
            evs = tsv(f"{D}/out/asm/{s}.events.tsv") if os.path.exists(f"{D}/out/asm/{s}.events.tsv") else []
            out = []
            for e in evs:
                st = int(e["start"])
                sup = int(e["support"]) if e.get("support", "").isdigit() else 0
                dp = bam.count(ctg, max(0, st - 1), st)
                out.append((e, min(1.0, sup / dp) if sup and dp else None))
            af[s] = out
    truth = tsv(f"{D}/out/truth.tsv")
    groups = {}
    for r in cohorts["fbrecall"]:
        s = r["sample"]
        groups.setdefault("chimeric" if s.startswith("mar_") else "kept" if rate[s] <= CUT else "excluded", set()).add(s)
    groups["chimeric_fb100"] = {r["sample"] for r in cohorts["fbrecall100"]}
    res = []
    for mode in ("keep_na", "drop_na"):
        for f in FLOORS:
            if f == 0 and mode == "drop_na":
                continue
            od = f"{D}/out/asm_af/{mode}_{f}"
            os.makedirs(od, exist_ok=True)
            for s, evs in af.items():
                keep = [e for e, a in evs if (a is None and mode == "keep_na") or (a is not None and a >= f) or f == 0]
                with open(f"{od}/{s}.events.tsv", "w") as fo:
                    if evs:
                        w = csv.DictWriter(fo, fieldnames=list(evs[0][0]), delimiter="\t")
                        w.writeheader()
                        w.writerows(keep)
                    else:
                        fo.write("sample\ttype\tstart\tend\tsize\tsupport\tleft\tright\tdepth_ratio\tsource\n")
            for g, ss in groups.items():
                fb = g == "chimeric_fb100"
                t = f"{D}/out/score_work/af_{g}.truth.tsv"
                with open(t, "w") as fo:
                    w = csv.DictWriter(fo, fieldnames=list(truth[0]), delimiter="\t")
                    w.writeheader()
                    for r in truth:
                        name = r["sample"] + "_fb100" if fb else r["sample"]
                        if name in ss:
                            w.writerow(dict(r, sample=name))
                o = f"{D}/out/score_work/af_{g}_{mode}_{f}.score.tsv"
                subprocess.run([sys.executable, f"{N}/score.py", "--truth", t, "--proto-dir", f"{D}/out/proto",
                                "--asm-dir", od,
                                "--p2-dir", "refbias/fbrecall100/run/p2" if fb else "refbias/fbrecall/run/p2",
                                "--refmap", "refbias/fbrecall100/run/p1/refmap.tsv" if fb else "refbias/fbrecall/run/p1/refmap.tsv",
                                "--mask-dir", f"{N}/out/masks", "--out", o], check=True, stdout=subprocess.DEVNULL)
                for r in tsv(o):
                    if r["stratum"] == "all" and r["caller"] in CALLERS:
                        res.append(dict(mode=mode, floor=f, group=g, caller=r["caller"], n_true=r["n_true"],
                                        typed_recall=r["typed_recall"], typed_calls=r["typed_calls"],
                                        precision=r["precision"]))
    with open(f"{D}/out/affloor.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(res[0]), delimiter="\t")
        w.writeheader()
        w.writerows(res)
    nna = sum(a is None for evs in af.values() for _, a in evs)
    print(f"our events: {sum(len(v) for v in af.values())}, without AF (support 0): {nna}")
    for c in CALLERS:
        print(f"\n{c}: typed recall / calls / precision")
        print("mode\tfloor\t" + "\t".join(("kept", "excluded", "chimeric", "chimeric_fb100")))
        for mode in ("keep_na", "drop_na"):
            for f in FLOORS:
                cells = []
                for g in ("kept", "excluded", "chimeric", "chimeric_fb100"):
                    r = next((x for x in res if x["mode"] == mode and x["floor"] == f and x["group"] == g
                              and x["caller"] == c), None)
                    cells.append(f"{r['typed_recall']}/{r['typed_calls']}/{r['precision']}" if r else "")
                if any(cells):
                    print(f"{mode}\t{f}\t" + "\t".join(cells))


if __name__ == "__main__":
    main()
