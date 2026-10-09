#!/usr/bin/env python3
"""Fold-back chimera rates for the cohorts in use (run from runroot after
01_measure.sbatch): per cohort, how many isolates a cutoff would remove, at
0.5, 1 and 4.8 per 1,000 reads (1 is the proposed cutoff: the 15 clean
in-panel controls are at or below 0.31 and the 10 chimeric ones at or above
4.8). For gwas1000, also by lineage and by RRDR status, since that cohort was
stratified on both and losses that differ by stratum would unbalance it.
Read length (avglen_1) and study come from Peter's results table and the ENA
table where available. Writes cohorts_in_use.summary.txt."""
import collections
import csv
import statistics as st

Q = "../analysis/foldback_qc"
COHORTS = [("pilot", "refbias/cohort.crams.tsv"), ("scale100", "refbias/scale100.crams.tsv"),
           ("l49", "refbias/l49.crams.tsv"), ("l7", "refbias/l7.crams.tsv"),
           ("scale200", "refbias/scale200.crams.tsv"), ("gwas1000", "refbias/gwas1000.crams.tsv")]
CUTS = (0.5, 1.0, 4.8)


def line(name, v):
    if not v:
        return f"{name:28s} n=0"
    s = f"{name:28s} n={len(v):4d} median={st.median(v):6.2f}"
    for c in CUTS:
        k = sum(x > c for x in v)
        s += f"  >{c:g}: {k:4d} ({k / len(v):4.0%})"
    return s


def main():
    rate = {r["sample"]: float(r["inv_chimera_per_1k"])
            for r in csv.DictReader(open(f"{Q}/cohorts_in_use_rates.tsv"), delimiter="\t")}
    out = []
    out.append(f"measured: {len(rate)} isolates")
    for name, path in COHORTS:
        ids = [r["sample"] for r in csv.DictReader(open(path), delimiter="\t")]
        out.append(line(name, [rate[s] for s in ids if s in rate]))
    ph = list(csv.DictReader(open("refbias/gwas1000.phenotype.tsv"), delimiter="\t"))
    out.append("\ngwas1000 by lineage")
    by = collections.defaultdict(list)
    for r in ph:
        if r["sample"] in rate:
            by[r["lineage"] or "(none)"].append(rate[r["sample"]])
    for k in sorted(by):
        out.append(line("  " + k, by[k]))
    out.append("\ngwas1000 by RRDR status (1 = RRDR variant)")
    by = collections.defaultdict(list)
    for r in ph:
        if r["sample"] in rate:
            by[r["rrdr"]].append(rate[r["sample"]])
    for k in sorted(by):
        out.append(line("  rrdr=" + k, by[k]))
    T = "/n/netscratch/sfortune_lab/Lab/pculviner/notebooks/260728_sv_exploration/completed_results_all.csv"
    rd = csv.DictReader(open(T))
    idc = rd.fieldnames[0]
    rl = {r[idc]: float(r["avglen_1"] or 0) for r in rd if r[idc] in rate}
    out.append("\nall measured, by read length (avglen_1)")
    for lo, hi in ((0, 110), (110, 160), (160, 400)):
        out.append(line(f"  {lo}-{hi} bp", [rate[s] for s in rate if lo <= rl.get(s, 0) < hi]))
    text = "\n".join(out)
    open(f"{Q}/cohorts_in_use.summary.txt", "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
