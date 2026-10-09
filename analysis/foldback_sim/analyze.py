#!/usr/bin/env python3
"""Summary of the fold-back chimera study, one row per condition (run from
runroot after 02_call_score.sbatch):

  injected        fraction of read pairs replaced, model
  measured        fold-back chimeras per 1,000 reads on the P2 BAM (the QC metric)
  SV callers      typed recall, typed calls and typed false calls in callable
                  sequence (score.py; false = typed calls x (1 - precision))
  small variants  P2 HaplotypeCaller calls (non-reference genotype) compared
                  with the same genome's 0% condition, which has the same base
                  reads: calls added and calls lost, SNPs and indels apart.
                  Any difference is caused by the chimeric pairs.

Writes out/summary.tsv and prints it.
"""
import csv
import gzip
import os

F = "../analysis/foldback_sim"
P2 = "refbias/foldsim/run/p2"
CALLERS = ("proto_asm", "dysgu_pass", "delly_pass", "proto_asm+dysgu_pass")


def small_variants(path):
    out = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            gt = c[9].split(":")[0] if len(c) > 9 else ""
            if gt in ("0", ".", "0/0", "./."):
                continue
            for alt in c[4].split(","):
                if alt != "*":
                    out.add((c[1], c[3], alt))
    return out


def kind(v):
    return "snp" if len(v[1]) == len(v[2]) == 1 else "indel"


def main():
    conds = list(csv.DictReader(open(f"{F}/conditions.tsv"), delimiter="\t"))
    base = {}
    for r in conds:
        if r["model"] == "none":
            base[r["source"]] = small_variants(f"{P2}/{r['sample']}.vcf.gz")
    rows = []
    for r in conds:
        s = r["sample"]
        sp = next(csv.DictReader(open(f"{F}/out/split/{s}.tsv"), delimiter="\t"))
        row = dict(sample=s, genome=r["source"][4:], model=r["model"], frac=r["frac"],
                   foldback_per_1k=sp["inv_chimera_per_1k"])
        sc = {x["caller"]: x for x in csv.DictReader(open(f"{F}/out/score/{s}.tsv"), delimiter="\t")
              if x["stratum"] == "all"}
        for c in CALLERS:
            x = sc.get(c)
            if not x:
                continue
            n, p = int(x["typed_calls"]), float(x["precision"]) if x["precision"] != "nan" else 0.0
            row[f"{c}_recall"] = x["typed_recall"]
            row[f"{c}_calls"] = n
            row[f"{c}_false"] = round(n * (1 - p))
        sv = small_variants(f"{P2}/{s}.vcf.gz")
        b = base[r["source"]]
        for k in ("snp", "indel"):
            row[f"{k}_added"] = sum(1 for v in sv - b if kind(v) == k)
            row[f"{k}_lost"] = sum(1 for v in b - sv if kind(v) == k)
        rows.append(row)
    os.makedirs(f"{F}/out", exist_ok=True)
    with open(f"{F}/out/summary.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    main()
