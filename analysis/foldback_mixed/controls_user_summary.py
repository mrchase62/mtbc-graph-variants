#!/usr/bin/env python3
"""Step 6 rerun summary: the user's Mutect2 settings in two frames (H37Rv and
each control's matched panel genome), beside the default-settings run.

Matched-frame calls are lifted to H37Rv (minimap2 -cx asm5 of the matched
genome against H37Rv, then paftools liftover) so the same repeat mask
applies; calls that do not lift are counted as "unplaced", not dropped
silently. PASS only; classes as in the default run (fixed >= 0.9, mixed
0.1-0.9, trace < 0.1).

    controls_user_summary.py --dir controls_user --default out/controls.tsv --out out/controls_user.tsv
"""
import argparse, csv, glob, gzip, os, statistics as S, subprocess
import numpy as np

REPO = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
BUILD = f"{REPO}/runroot/refbias/build/7713a8d71d8e-fix1"
MASK = f"{BUILD}/assets/repeat_mask.bed"
H37RV = f"{BUILD}/refs/GCF_000195955.fasta"
RATES = f"{REPO}/analysis/foldback_qc/phaseB_rates.tsv"
mask = np.zeros(4411532 + 200, bool)
for l in open(MASK):
    c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True


def calls(path):
    out = []
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l[0] == "#": continue
            f = l.rstrip("\n").split("\t")
            if f[6] != "PASS": continue
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            for alt, af in zip(f[4].split(","), fmt["AF"].split(",")):
                if max(len(f[3]), len(alt)) > 50: continue
                out.append((f[0], int(f[1]), len(f[3]), "snv" if len(f[3]) == len(alt) == 1 else "indel", float(af)))
    return out


def lift(ref, rows, work):
    """{(contig, pos): h37rv_pos} for matched-frame positions."""
    paf = f"{work}/{ref}.paf"
    if not os.path.exists(paf):
        subprocess.run(f"{os.environ['MTB_SV_MINIMAP2']} -cx asm5 --cs -t 4 {H37RV} {BUILD}/refs/{ref}.fasta > {paf}.tmp "
                       f"&& mv {paf}.tmp {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
    bed = "".join(f"{c}\t{p - 1}\t{p}\t{c}:{p}\n" for c, p, *_ in rows)
    r = subprocess.run([os.environ["MTB_K8"], os.environ["MTB_PAFTOOLS"], "liftover", "-l", "1", paf, "-"],
                       input=bed, capture_output=True, text=True, check=True)
    m = {}
    for l in r.stdout.splitlines():
        f = l.split("\t")
        c, _, p = f[3].rsplit("_", 2)  # paftools renames to contig_start_end
        m[(c, int(p))] = int(f[2])
    return m


def count(rows, frame, ref, work):
    n = {f"{t}_{c}": 0 for t in ("snv", "indel") for c in ("fixed", "mixed", "trace")}
    n["unplaced"] = 0
    lm = lift(ref, rows, work) if frame == "matched" else None
    for c, p, rl, t, af in rows:
        hp, span = (p, rl) if lm is None else (lm.get((c, p)), 1)
        if hp is None:
            n["unplaced"] += 1; continue
        if mask[hp - 1:hp - 1 + span].any(): continue
        k = "fixed" if af >= 0.9 else ("mixed" if af >= 0.1 else "trace")
        n[f"{t}_{k}"] += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--default", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    work = f"{a.dir}/lift"; os.makedirs(work, exist_ok=True)
    fb = {r[0]: float(r[3]) for r in csv.reader(open(RATES), delimiter="\t") if r[0] != "sample"}
    dflt = {r["sample"]: r for r in csv.DictReader(open(a.default), delimiter="\t")}
    res = {}
    for info in sorted(glob.glob(f"{a.dir}/*.info.tsv")):
        s, arm, frame, ref, depth, frac = open(info).read().split()
        v = info.replace(".info.tsv", ".vcf.gz")
        if not os.path.exists(v): continue
        for k, x in count(calls(v), frame, ref, work).items():
            res.setdefault(s, {"sample": s, "foldback": fb.get(s)})[f"{frame}_{arm}_{k}"] = x
        res[s][f"{frame}_ref"] = ref
    rows = list(res.values())
    keys = sorted({k for r in rows for k in r} - {"sample", "foldback"})
    with open(a.out, "w") as fo:
        fo.write("\t".join(["sample", "foldback"] + keys) + "\n")
        for r in rows:
            fo.write("\t".join(str(r.get(k, "")) for k in ["sample", "foldback"] + keys) + "\n")
    print("mixed SNV / mixed indel, real-syn, per control (default H37Rv | user H37Rv | user matched)")
    for r in sorted(rows, key=lambda r: r["foldback"]):
        d = dflt[r["sample"]]
        g = lambda f, t, arm: r.get(f"{f}_{arm}_{t}_mixed", "NA")
        print(f"{r['sample']:16s} {r['foldback']:6.2f}  default {d['real_snv_mixed']}/{d['syn_snv_mixed']}, {d['real_indel_mixed']}/{d['syn_indel_mixed']}"
              f"  | user-H37Rv {g('h37rv', 'snv', 'real')}/{g('h37rv', 'snv', 'syn')}, {g('h37rv', 'indel', 'real')}/{g('h37rv', 'indel', 'syn')}"
              f"  | user-matched {g('matched', 'snv', 'real')}/{g('matched', 'snv', 'syn')}, {g('matched', 'indel', 'real')}/{g('matched', 'indel', 'syn')}")
    for lab, grp in (("clean (<= 0.31)", [r for r in rows if r["foldback"] <= 0.31]),
                     ("chimeric (> 4.8)", [r for r in rows if r["foldback"] > 4.8])):
        out = [f"{lab}: n={len(grp)}; median real-minus-synthetic"]
        for f in ("h37rv", "matched"):
            for t in ("snv_mixed", "indel_mixed", "snv_trace"):
                v = [r[f"{f}_real_{t}"] - r[f"{f}_syn_{t}"] for r in grp if f"{f}_real_{t}" in r and f"{f}_syn_{t}" in r]
                if v: out.append(f"user-{f} {t} {S.median(v):g}")
        print("; ".join(out))


if __name__ == "__main__":
    main()
