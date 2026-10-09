#!/usr/bin/env python3
"""Runner tables for a QC'd Phase B cohort (after 08_qc.sbatch), as built by
hand for marinQC5: qc_summary.py per sample -> qc_table.tsv; PASS samples ->
cohort.tsv (meandepth = depth after QC and downsampling; coverage,
mapping_rate, error_rate and duplicate_rate are the fixed placeholders
marinQC5 used, since the runner does not filter on them here) and crams.tsv
(unaligned BAM paths relative to MTB_CRAM_ROOT). Run from runroot:
  python3 make_qc_cohort.py --list <samples> --cohort <C> --cram-root <MTB_CRAM_ROOT>
"""
import argparse
import csv
import os
import subprocess
import sys

H = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--cram-root", required=True)
    a = ap.parse_args()
    base = os.path.join("refbias", a.cohort)
    rows = []
    for s in (x.strip() for x in open(a.list) if x.strip()):
        q = os.path.join(base, "qc", s)
        kr = os.path.join(q, "kraken2.tsv")
        subprocess.run([sys.executable, os.path.join(H, "qc_summary.py"), "--sample", s,
                        "--fastp", f"{q}/fastp.json", "--tbprofiler", f"{q}/tbprofiler.json",
                        "--final-bases", open(f"{q}/final_bases.txt").read().strip(),
                        "--out", f"{q}/qc.tsv"] + (["--kraken", kr] if os.path.exists(kr) else []),
                       check=True, stdout=subprocess.DEVNULL)
        rows.append(next(csv.DictReader(open(f"{q}/qc.tsv"), delimiter="\t")))
    with open(os.path.join(base, "qc_table.tsv"), "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    ok = [r for r in rows if r["status"] == "PASS"]
    with open(os.path.join(base, "cohort.tsv"), "w") as fo:
        fo.write("sample\tlineage\tdeepest\tmeandepth\tcoverage\tmapping_rate\terror_rate\tduplicate_rate\n")
        for r in ok:
            fo.write(f"{r['sample']}\t{r['lineage']}\t{r['sub_lineage']}\t{r['final_depth']}\t100\t1\t0.002\t0\n")
    with open(os.path.join(base, "crams.tsv"), "w") as fo:
        fo.write("sample\tcram_relpath\n")
        for r in ok:
            p = os.path.realpath(os.path.join(base, "reads", f"{r['sample']}.bam"))
            fo.write(f"{r['sample']}\t{os.path.relpath(p, a.cram_root)}\n")
    for r in rows:
        print(r["sample"], r["status"], r["final_depth"], r["reason"])


if __name__ == "__main__":
    main()
