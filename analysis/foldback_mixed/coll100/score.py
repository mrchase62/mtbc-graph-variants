#!/usr/bin/env python3
"""Score the 50 clean + 50 chimeric collection isolates (run.sbatch outputs).

There is no assembly truth for collection isolates, so fixed SNVs on the
matched reference are checked against an independent frame:

  expected  positions where the isolate's own H37Rv-frame call (P1's haploid
            HaplotypeCaller, work/<cohort>_p1/<s>.h37rv.vcf.gz; REF where it
            has no record) differs from the matched genome's direct-alignment
            state (analysis/panel_checks/states/<ref>.npz), both called,
            outside the repeat mask. These are the SNVs the matched-reference
            caller should find.
  unsupported  a fixed SNV on the matched reference, lifted to H37Rv, that is
            not an expected difference (candidate artifact; can also be a
            frame difference)
  missed    an expected difference with no fixed call (AF >= 0.9); "demoted"
            when the caller has it at AF 0.1-0.9 instead

Per caller (HaplotypeCaller and the user's Mutect2 on the matched reference;
the user's Mutect2 on H37Rv for counts only) and isolate: fixed / mixed
(0.1-0.9) / trace (0.05-0.1) SNV and indel counts, unmasked, PASS (HC: all
non-REF calls are fixed). Group medians for clean vs chimeric are printed.

    score.py --out out/coll100_score.tsv
"""
import argparse, bisect, collections, csv, gzip, os, statistics as S, subprocess
import numpy as np

R = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
RR = f"{R}/runroot/refbias"
BUILD = f"{RR}/build/7713a8d71d8e-fix1"
C = f"{R}/analysis/foldback_mixed/coll100"
STATES = f"{R}/analysis/panel_checks/states"
K8, PT, MM2 = os.environ["MTB_K8"], os.environ["MTB_PAFTOOLS"], os.environ["MTB_SV_MINIMAP2"]
L = 4411532
mask = np.zeros(L + 200, bool)
for l in open(f"{BUILD}/assets/repeat_mask.bed"):
    c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True
CODE = {b: i for i, b in enumerate("ACGT", 1)}


def vcf(path, hc):
    out = []
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l[0] == "#": continue
            f = l.rstrip("\n").split("\t")
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            alts = f[4].split(",")
            if hc:
                gt = fmt.get("GT", ".").replace("|", "/").split("/")[0]
                if not gt.isdigit() or gt == "0" or alts[int(gt) - 1] == "*": continue
                out.append((f[0], int(f[1]), f[3], alts[int(gt) - 1], 1.0))
            elif f[6] == "PASS":
                for alt, af in zip(alts, fmt["AF"].split(",")):
                    out.append((f[0], int(f[1]), f[3], alt, float(af)))
    return [x for x in out if max(len(x[2]), len(x[3])) <= 50]


def lifter(ref, positions, work):
    paf = f"{work}/{ref}.paf"
    if not os.path.exists(paf):
        subprocess.run(f"{MM2} -cx asm5 --cs -t 4 {BUILD}/refs/GCF_000195955.fasta {BUILD}/refs/{ref}.fasta "
                       f"> {paf}.tmp && mv {paf}.tmp {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
    bed = "".join(f"{c}\t{p - 1}\t{p}\t{c}_{p - 1}_{p}\n" for c, p in positions)
    r = subprocess.run([K8, PT, "liftover", "-l", "1", paf, "-"], input=bed, capture_output=True, text=True, check=True)
    m = {}
    for l in r.stdout.splitlines():
        f = l.split("\t"); c, _, p = f[3].rsplit("_", 2); m[(c, int(p))] = (int(f[2]), f[5])
    return m


RC = str.maketrans("ACGT", "TGCA")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    work = f"{C}/lift"; os.makedirs(work, exist_ok=True)
    seq = "".join(l.strip() for l in open(f"{BUILD}/refs/GCF_000195955.fasta") if not l.startswith(">"))
    rc = np.frombuffer(seq.encode(), np.uint8)
    lut = np.zeros(256, np.uint8)
    for b, i in CODE.items(): lut[ord(b)] = i
    rc = lut[rc]
    rows = []
    for smp in csv.DictReader(open(f"{C}/samples.tsv"), delimiter="\t"):
        s, grp = smp["sample"], smp["group"]
        info = f"{C}/calls/m2_matched_{s}.info.tsv"
        if not os.path.exists(f"{C}/calls/m2_h37rv_{s}.vcf.gz"):
            print("missing", s); continue
        ref = open(info).read().split("\t")[4]
        coh = "new1k" if grp == "clean" else "fbmix50"
        dist = next(r["snp_distance"] for r in csv.DictReader(open(f"{RR}/{coh}/p1/refmap.tsv"), delimiter="\t") if r["sample"] == s)
        # isolate in H37Rv frame (P1 HC): base per position, unknown at non-SNV records
        iso = rc.copy()
        for c, p, r_, al, _ in vcf(f"{RR}/work/{coh}_p1/{s}.h37rv.vcf.gz", True):
            if len(r_) == len(al) == 1 and al in CODE: iso[p - 1] = CODE[al]
            else: iso[p - 1:p - 1 + max(len(r_), 1)] = 0
        g = np.load(f"{STATES}/{ref}.npz")["state"]
        exp = np.where((iso > 0) & (g > 0) & (iso != g) & ~mask[:L])[0] + 1  # 1-based H37Rv positions
        exp_set = set(exp.tolist())
        callers = {"HC_matched": (vcf(f"{C}/calls/hc_matched_{s}.vcf.gz", True), True),
                   "M2_matched": (vcf(f"{C}/calls/m2_matched_{s}.vcf.gz", False), True),
                   "M2_h37rv": (vcf(f"{C}/calls/m2_h37rv_{s}.vcf.gz", False), False)}
        contig = open(f"{BUILD}/refs/{ref}.fasta").readline()[1:].split()[0]
        pos = {(contig, p) for k, (cl, m) in callers.items() if m for _, p, *_ in cl}
        lm = lifter(ref, sorted(pos), work)
        for name, (calls, matched) in callers.items():
            n = collections.Counter(); fixed_exp, low_exp = set(), set()
            for c, p, r_, al, af in calls:
                if matched:
                    if (c, p) not in lm: n["unplaced"] += 1; continue
                    hp, strand = lm[(c, p)]
                else:
                    hp, strand = p, "+"
                if mask[hp - 1: hp - 1 + (len(r_) if not matched else 1)].any(): continue
                t = "snv" if len(r_) == len(al) == 1 else "indel"
                k = "fixed" if af >= 0.9 else ("mixed" if af >= 0.1 else "trace")
                n[f"{t}_{k}"] += 1
                if matched and t == "snv" and af >= 0.1:
                    if hp in exp_set:
                        (fixed_exp if k == "fixed" else low_exp).add(hp)
                    elif k == "fixed":
                        n["snv_fixed_unsupported"] += 1
            if matched:
                n["snv_expected"] = len(exp_set)
                n["snv_missed"] = len(exp_set - fixed_exp - low_exp)
                n["snv_demoted"] = len(low_exp - fixed_exp)
            rows.append({"sample": s, "group": grp, "lineage": smp["lineage"], "foldback": smp["foldback"],
                         "matched_ref": ref, "p1_snp_distance": dist, "caller": name, **n})
    keys = ["sample", "group", "lineage", "foldback", "matched_ref", "p1_snp_distance", "caller"]
    keys += sorted({k for r in rows for k in r} - set(keys))
    with open(a.out, "w") as fo:
        fo.write("\t".join(keys) + "\n")
        for r in rows: fo.write("\t".join(str(r.get(k, 0)) for k in keys) + "\n")
    cols = ["snv_expected", "snv_fixed", "snv_fixed_unsupported", "snv_missed", "snv_demoted", "snv_mixed",
            "snv_trace", "indel_fixed", "indel_mixed", "indel_trace"]
    print("median per isolate (unmasked)\ncaller      group     n   " + "  ".join(c.replace("snv_", "S:").replace("indel_", "I:") for c in cols))
    for name in ("HC_matched", "M2_matched", "M2_h37rv"):
        for g in ("clean", "chimeric"):
            x = [r for r in rows if r["caller"] == name and r["group"] == g]
            print(f"{name:11s} {g:9s} {len(x):3d}  " + "  ".join(
                (f"{S.median(r.get(c, 0) for r in x):g}" if x else "NA") if (name != "M2_h37rv" or not c.endswith(("expected", "unsupported", "missed", "demoted"))) else "-"
                for c in cols))


if __name__ == "__main__":
    main()
