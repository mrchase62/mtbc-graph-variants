#!/usr/bin/env python3
"""Draw the collection comparison for the fold-back vs variant-call question
(user, 2026-10-14: "Run the 50 clean and 50 chimeric too").

clean:    50 new1k isolates with fold-back <= 0.31 (P1 already run, build fix1)
chimeric: 50 standard-pass lineage 1-4 isolates with fold-back > 4.8, not in
          new1k (P1 runs for them as registry row fbmix50, build fix1)
Both 12-13 per lineage (L1, L2 get 13), random with a fixed seed, meandepth
>= 60 (the read standard). Writes samples.tsv and, for the chimeric 50, the
cohort and CRAM tables refbias_run.sh needs.
"""
import csv, gzip, json, os, random
R = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
FB = f"{R}/analysis/foldback_qc"
OUT = f"{R}/analysis/foldback_mixed/coll100"
LIN = ["lineage1", "lineage2", "lineage3", "lineage4"]
rnd = random.Random(20261014)
fb = {}
for f in (f"{FB}/screen50k/rates_all.tsv", f"{FB}/cohorts_in_use_rates.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        fb[r["sample"]] = float(r["inv_chimera_per_1k"])
qc = {r["sample"]: r for r in csv.DictReader(open(f"{FB}/screen50k/pass_standard.tsv"), delimiter="\t")}
new1k = {r["sample"]: r for r in csv.DictReader(open(f"{R}/runroot/refbias/cohort.new1k.tsv"), delimiter="\t")}
meta = json.load(gzip.open(f"{os.environ['MTB_CRAM_ROOT']}/metadata/samples.json.gz"))
pick = []
for k, l in enumerate(LIN):
    n = 13 if k < 2 else 12
    clean = sorted(s for s, r in new1k.items() if r["lineage"] == l and fb.get(s, 9) <= 0.31 and float(r["meandepth"]) >= 60)
    chim = sorted(s for s, r in qc.items() if r["lineage"] == l and s not in new1k and fb.get(s, 0) > 4.8
                  and float(r["meandepth"]) >= 60 and s in meta)
    pick += [(s, "clean", l) for s in rnd.sample(clean, n)] + [(s, "chimeric", l) for s in rnd.sample(chim, n)]
with open(f"{OUT}/samples.tsv", "w") as fo:
    fo.write("sample\tgroup\tlineage\tfoldback\tcram_relpath\n")
    for s, g, l in pick:
        fo.write(f"{s}\t{g}\t{l}\t{fb[s]}\t{meta[s]['path']}\n")
cols = ["sample", "lineage", "deepest", "meandepth", "coverage", "mapping_rate", "error_rate", "duplicate_rate"]
with open(f"{R}/runroot/refbias/cohort.fbmix50.tsv", "w") as fc, open(f"{R}/runroot/refbias/fbmix50.crams.tsv", "w") as fr:
    fc.write("\t".join(cols) + "\n"); fr.write("sample\tcram_relpath\n")
    for s, g, l in pick:
        if g != "chimeric": continue
        fc.write("\t".join(qc[s][c] for c in cols) + "\n"); fr.write(f"{s}\t{meta[s]['path']}\n")
print("drawn:", len(pick))
