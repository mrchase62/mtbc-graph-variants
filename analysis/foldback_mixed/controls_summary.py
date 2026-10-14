#!/usr/bin/env python3
"""Step 6 summary: Mutect2 PASS calls on each control's real and synthetic
reads, unmasked only, by type and AF class, with the control's fold-back
rate. Real minus synthetic is what the library adds (chimeras, error profile
and any true within-host diversity); comparing that excess between chimeric
(> 4.8) and clean (<= 0.31) controls isolates the chimera effect.

    controls_summary.py --dir controls --out out/controls.tsv
"""
import argparse, csv, gzip, statistics as S
import numpy as np

REPO = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
MASK = f"{REPO}/runroot/refbias/build/7713a8d71d8e-fix1/assets/repeat_mask.bed"
RATES = f"{REPO}/analysis/foldback_qc/phaseB_rates.tsv"
CTRL = f"{REPO}/analysis/novel_events/controls/controls.txt"
mask = np.zeros(4411532 + 200, bool)
for l in open(MASK):
    c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True


def counts(path):
    n = {f"{t}_{c}": 0 for t in ("snv", "indel") for c in ("fixed", "mixed", "trace")}
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l[0] == "#": continue
            f = l.rstrip("\n").split("\t")
            if f[6] != "PASS": continue
            pos, ref, alts = int(f[1]), f[3], f[4].split(",")
            if mask[pos - 1:pos - 1 + len(ref)].any(): continue
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            for alt, af in zip(alts, fmt["AF"].split(",")):
                af = float(af)
                if max(len(ref), len(alt)) > 50: continue
                t = "snv" if len(ref) == 1 and len(alt) == 1 else "indel"
                c = "fixed" if af >= 0.9 else ("mixed" if af >= 0.1 else "trace")
                n[f"{t}_{c}"] += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    fb = {r[0]: float(r[3]) for r in csv.reader(open(RATES), delimiter="\t") if r[0] != "sample"}
    rows = []
    for s in (l.strip() for l in open(CTRL) if l.strip()):
        r = {"sample": s, "foldback": fb.get(s, float("nan"))}
        for arm in ("real", "syn"):
            for k, v in counts(f"{a.dir}/{arm}_{s}.vcf.gz").items():
                r[f"{arm}_{k}"] = v
        rows.append(r)
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print("sample\tfoldback\tmixed SNV real/syn\tmixed indel real/syn\ttrace SNV real/syn")
    for r in sorted(rows, key=lambda r: r["foldback"]):
        print(f"{r['sample']}\t{r['foldback']:.2f}\t{r['real_snv_mixed']}/{r['syn_snv_mixed']}\t"
              f"{r['real_indel_mixed']}/{r['syn_indel_mixed']}\t{r['real_snv_trace']}/{r['syn_snv_trace']}")
    for lab, grp in (("clean (<= 0.31)", [r for r in rows if r["foldback"] <= 0.31]),
                     ("chimeric (> 4.8)", [r for r in rows if r["foldback"] > 4.8])):
        if not grp: continue
        ex = lambda k: S.median(r[f"real_{k}"] - r[f"syn_{k}"] for r in grp)
        print(f"{lab}: n={len(grp)}; median real-minus-synthetic: mixed SNV {ex('snv_mixed'):g}, "
              f"mixed indel {ex('indel_mixed'):g}, trace SNV {ex('snv_trace'):g}")


if __name__ == "__main__":
    main()
