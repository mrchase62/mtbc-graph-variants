#!/usr/bin/env python3
"""Score the matched-reference settings grid against assembly truth.

Truth for each control: its own complete assembly aligned to its matched
reference (minimap2 -cx asm5 --cs, paftools call -l 1000 -L 10000). Only
reference positions callable by that alignment (paftools R lines) and, once
lifted to H37Rv, outside the build's repeat mask are scored.

Callers (all on the matched reference, PASS only):
  A, B, C, D  Mutect2 settings (A = the user's; controls_user/, controls_grid/)
  HC          production's haploid HaplotypeCaller (P2: refbias/marinB63/run/p2,
              refbias/marinSyn25/run/p2), full depth; every call counts as fixed

Per caller x reads (real / synthetic) x type (SNV / indel):
  fixed_tp     fixed call (AF >= 0.9) matching a truth variant
  fixed_fp     fixed call with no truth variant  (artifactual fixed call)
  demoted      truth variant called only at AF 0.1-0.9
  missed       truth variant not called at AF >= 0.1
  mixed_false  call at AF 0.1-0.9 with no truth variant (artifactual mixed call)
A SNV matches on position and base; an indel on length change within 10 bp
(left-alignment differs between paftools and GATK).

    grid_score.py --out out/grid_score.tsv
"""
import argparse, bisect, collections, csv, gzip, os, statistics as S, subprocess
import numpy as np

REPO = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
RR = f"{REPO}/runroot/refbias"
BUILD = f"{RR}/build/7713a8d71d8e-fix1"
FM = f"{REPO}/analysis/foldback_mixed"
CTRL = f"{REPO}/analysis/novel_events/controls/controls.txt"
TRUTH = f"{REPO}/analysis/novel_events/phaseB/truth_genomes.tsv"
ASM = f"{REPO}/analysis/external_assemblies/genomes"
RATES = f"{REPO}/analysis/foldback_qc/phaseB_rates.tsv"
MM2, K8, PT = os.environ["MTB_SV_MINIMAP2"], os.environ["MTB_K8"], os.environ["MTB_PAFTOOLS"]
mask = np.zeros(4411532 + 200, bool)
for l in open(f"{BUILD}/assets/repeat_mask.bed"):
    c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, stderr=subprocess.DEVNULL)


def truth(s, ref, work):
    fa = next(r["file"] for r in csv.DictReader(open(TRUTH), delimiter="\t") if r["sample"] == s)
    var = f"{work}/{s}.{ref}.var"
    if not os.path.exists(var):
        sh(f"{MM2} -cx asm5 --cs -t 4 {BUILD}/refs/{ref}.fasta {ASM}/{fa} | sort -k6,6 -k8,8n "
           f"| {K8} {PT} call -l 1000 -L 10000 - > {var}.tmp && mv {var}.tmp {var}")
    callable_, vs = [], []
    for l in open(var):
        f = l.rstrip("\n").split("\t")
        if f[0] == "R":
            callable_.append((int(f[2]), int(f[3])))
        elif f[0] == "V" and f[4] == "1":
            r, a = f[6].upper().replace("-", ""), f[7].upper().replace("-", "")
            vs.append((int(f[2]) + 1, r, a))  # 1-based
    return callable_, vs


def h37_lifter(ref, positions, work):
    paf = f"{FM}/controls_user/lift/{ref}.paf"
    bed = "".join(f"{c}\t{p - 1}\t{p}\t{c}_{p - 1}_{p}\n" for c, p in positions)
    r = subprocess.run([K8, PT, "liftover", "-l", "1", paf, "-"], input=bed, capture_output=True, text=True, check=True)
    m = {}
    for l in r.stdout.splitlines():
        f = l.split("\t"); c, _, p = f[3].rsplit("_", 2); m[(c, int(p))] = int(f[2])
    return m


def vcf_calls(path, hc=False):
    out = []
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l[0] == "#": continue
            f = l.rstrip("\n").split("\t")
            if not hc and f[6] != "PASS": continue
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = fmt.get("GT", "1")
            alts = f[4].split(",")
            if hc:
                if gt in (".", "0"): continue
                k = int(gt.split("/")[0].split("|")[0])
                if k == 0 or alts[k - 1] == "*": continue
                out.append((f[0], int(f[1]), f[3], alts[k - 1], 1.0))
            else:
                for alt, af in zip(alts, fmt["AF"].split(",")):
                    out.append((f[0], int(f[1]), f[3], alt, float(af)))
    return out


def score(calls, tv, ok_pos):
    """tv: list of (pos, ref, alt) truth; ok_pos(pos) -> scorable."""
    snv_t = {(p, a): k for k, (p, r, a) in enumerate(tv) if len(r) == len(a) == 1}
    ind_t = sorted((p, len(a) - len(r), k) for k, (p, r, a) in enumerate(tv) if not (len(r) == len(a) == 1))
    ind_p = [x[0] for x in ind_t]
    best = {}  # truth index -> max AF seen
    n = collections.Counter()
    for c, p, r, a, af in calls:
        if not ok_pos(p) or af < 0.1 or max(len(r), len(a)) > 50: continue
        t = "snv" if len(r) == len(a) == 1 else "indel"
        hit = None
        if t == "snv":
            hit = snv_t.get((p, a))
        else:
            dl = len(a) - len(r)
            i = bisect.bisect_left(ind_p, p - 10)
            while i < len(ind_t) and ind_t[i][0] <= p + 10:
                if ind_t[i][1] == dl: hit = ind_t[i][2]; break
                i += 1
        fixed = af >= 0.9
        if hit is None:
            n[f"{t}_{'fixed_fp' if fixed else 'mixed_false'}"] += 1
        else:
            best[hit] = max(best.get(hit, 0), af)
    for k, (p, r, a) in enumerate(tv):
        if not ok_pos(p): continue
        t = "snv" if len(r) == len(a) == 1 else "indel"
        n[f"{t}_truth"] += 1
        af = best.get(k)
        n[f"{t}_{'missed' if af is None else ('fixed_tp' if af >= 0.9 else 'demoted')}"] += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    work = f"{FM}/controls_grid/truth"; os.makedirs(work, exist_ok=True)
    fb = {r[0]: float(r[3]) for r in csv.reader(open(RATES), delimiter="\t") if r[0] != "sample"}
    rows = []
    for s in (l.strip() for l in open(CTRL) if l.strip()):
        ref = [l.split("\t")[1].strip() for l in open(f"{RR}/marinB63/run/p1/{s}.p1.done") if l.startswith("reference")][0]
        cal, tv = truth(s, ref, work)
        contig = open(f"{BUILD}/refs/{ref}.fasta").readline()[1:].split()[0]
        cs = sorted(cal); starts = [x[0] for x in cs]
        def in_cal(p):
            i = bisect.bisect_right(starts, p - 1) - 1
            return i >= 0 and cs[i][0] <= p - 1 < cs[i][1]
        files = {("HC", "real"): (f"{RR}/marinB63/run/p2/{s}.vcf.gz", True),
                 ("HC", "syn"): (f"{RR}/marinSyn25/run/p2/syn_{s}.vcf.gz", True)}
        for arm in ("real", "syn"):
            files[("A", arm)] = (f"{FM}/controls_user/matched_{arm}_{s}.vcf.gz", False)
            for v in "BCD":
                files[(v, arm)] = (f"{FM}/controls_grid/{v}_{arm}_{s}.vcf.gz", False)
        allcalls = {k: vcf_calls(p, hc) for k, (p, hc) in files.items() if os.path.exists(p)}
        pos = {(contig, p) for p, _, _ in tv} | {(contig, c[1]) for v in allcalls.values() for c in v}
        lm = h37_lifter(ref, sorted(pos), work)
        ok = lambda p: in_cal(p) and (contig, p) in lm and not mask[lm[(contig, p)] - 1]
        for (caller, arm), calls in allcalls.items():
            n = score(calls, tv, ok)
            rows.append({"sample": s, "foldback": fb.get(s), "group": "chimeric" if fb.get(s, 0) > 4.8 else "clean",
                         "caller": caller, "reads": arm, **n})
    keys = ["sample", "foldback", "group", "caller", "reads"] + sorted({k for r in rows for k in r} - {"sample", "foldback", "group", "caller", "reads"})
    with open(a.out, "w") as fo:
        fo.write("\t".join(keys) + "\n")
        for r in rows: fo.write("\t".join(str(r.get(k, 0)) for k in keys) + "\n")
    print("medians per control: real reads (synthetic in brackets)")
    print("group     caller  type   truth  fixed_fp  mixed_false  demoted  missed")
    for g in ("clean", "chimeric"):
        for caller in ("HC", "A", "B", "C", "D"):
            for t in ("snv", "indel"):
                def med(arm, k):
                    v = [r.get(f"{t}_{k}", 0) for r in rows if r["group"] == g and r["caller"] == caller and r["reads"] == arm]
                    return f"{S.median(v):g}" if v else "NA"
                print(f"{g:9s} {caller:6s}  {t:5s}  {med('real', 'truth'):>5s}  {med('real', 'fixed_fp'):>4s} ({med('syn', 'fixed_fp')})"
                      f"  {med('real', 'mixed_false'):>4s} ({med('syn', 'mixed_false')})  {med('real', 'demoted'):>3s} ({med('syn', 'demoted')})"
                      f"  {med('real', 'missed'):>3s} ({med('syn', 'missed')})")


if __name__ == "__main__":
    main()
