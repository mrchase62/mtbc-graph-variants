#!/usr/bin/env python3
"""How much does mapping an isolate's Illumina reads to its own finished
genome call? (report only; run from runroot with $MTB_PY after
11_call_score.sbatch)

Controls: Phase B in_panel isolates (out/samples.tsv) whose P1 matched
reference R is the panel genome of the same isolate (genomes.tsv
same_isolate_in_cx333 == refmap reference), and not in exclude.tsv. Every
difference called against R is then noise: read errors, mapping in repeats,
or a difference between the assembly and the DNA that was sequenced.

Counted per isolate, split by R's callable mask (callable_mask.py; a call
whose span +-50 bp touches it is masked):
  HaplotypeCaller (P2 <sample>.vcf.gz), PASS: SNPs; indels 1-9 bp; 10-49 bp
  our SV caller after assembly (out/qc51/asm): typed calls (DEL, INS, REPL,
      INV, REARR) and unpaired breakpoints (BND)
  dysgu (P2 <sample>.dysgu.vcf): PASS; PASS with AF >= 0.2 (dysgu_af20)
The tatC-helY and HaplotypeCaller hand-off filters are not applied: this is
the raw output of each step.
Writes out/qc51/self_noise.tsv and prints totals.
"""
import csv
import gzip
import os
import sys

N = "../analysis/novel_events"
sys.path.insert(0, N)
from score import load_masks, masked, read_proto, read_vcf  # noqa: E402

RM = "refbias/marinQC51/run/p1/refmap.tsv"
P2 = "refbias/marinQC51/run/p2"
O = f"{N}/phaseB/out/qc51"


def hc_counts(path, m):
    n = dict.fromkeys(("snp", "ind1_9", "ind10_49"), 0)
    nm = dict(n)
    for line in gzip.open(path, "rt"):
        if line.startswith("#"):
            continue
        c = line.split("\t", 8)
        if c[6] not in ("PASS", "."):
            continue
        pos, ref, alts = int(c[1]), c[3], c[4].split(",")
        d = max(abs(len(a) - len(ref)) for a in alts)
        k = "snp" if d == 0 else "ind1_9" if d < 10 else "ind10_49" if d < 50 else None
        if k is None:
            continue
        (nm if masked(m, pos, pos + len(ref) - 1) else n)[k] += 1
    return n, nm


def main():
    rm = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM), delimiter="\t")}
    own = {r["cohort_sample"]: r["same_isolate_in_cx333"]
           for r in csv.DictReader(open(f"{N}/phaseB/genomes.tsv"), delimiter="\t")}
    grp = {r["sample"]: r["group"] for r in csv.DictReader(open(f"{N}/phaseB/out/samples.tsv"), delimiter="\t")}
    excl = {r["sample"] for r in csv.DictReader(open(f"{N}/phaseB/exclude.tsv"), delimiter="\t")}
    masks = load_masks(RM, f"{N}/out/masks")
    rows = []
    for s in sorted(rm):
        if grp.get(s) != "in_panel" or s in excl or own.get(s) != rm[s]:
            continue
        m = masks[s]
        hc, hcm = hc_counts(f"{P2}/{s}.vcf.gz", m)
        r = dict(sample=s, reference=rm[s])
        for k in hc:
            r[f"hc_{k}"], r[f"hc_{k}_masked"] = hc[k], hcm[k]
        ours = read_proto(f"{O}/asm/{s}.events.tsv")
        for lab, keep in (("sv_typed", lambda c: c["type"] != "BND"), ("sv_bnd", lambda c: c["type"] == "BND")):
            cs = [c for c in ours if keep(c)]
            r[lab] = sum(1 for c in cs if not masked(m, c["start"], c["end"]))
            r[f"{lab}_masked"] = len(cs) - r[lab]
        for lab, af in (("dysgu_pass", None), ("dysgu_af20", 0.2)):
            cs = read_vcf(f"{P2}/{s}.dysgu.vcf", True, af)
            r[lab] = sum(1 for c in cs if not masked(m, c["start"], c["end"]))
            r[f"{lab}_masked"] = len(cs) - r[lab]
        rows.append(r)
    with open(f"{O}/self_noise.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} isolates mapped to their own genome")
    for k in rows[0]:
        if k in ("sample", "reference"):
            continue
        v = [r[k] for r in rows]
        print(f"  {k:22s} total {sum(v):5d}  per isolate {','.join(map(str, v))}")


if __name__ == "__main__":
    main()
