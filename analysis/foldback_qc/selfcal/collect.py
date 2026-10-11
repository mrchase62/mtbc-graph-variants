#!/usr/bin/env python3
"""One row per isolate for the fold-back cut-off calibration and the
read-quality table (user, 2026-10-11: set sample filters from the
distribution in real data). Run from runroot after 03_post.sbatch.

Rows: the 54 selfcal54 isolates (those failing QC keep their QC fields and
nothing after) and the 11 Phase B own-genome controls (marinQC51), all with
reads QC'd the same way and mapped to their own genome.
Columns:
  metadata   study, instrument, library fields (ENA), first_public
  QC         status, reason, raw/final depth, MTBC fraction, top other taxon
  reads      samtools stats on the P2 BAM: read length mean/max, insert mean
             and sd, error rate (mismatches per aligned base), properly
             paired fraction; duplicate fraction (foldback_v2)
  library    adapter read-through on the raw FASTQ (adapters/adapter_count.py:
             Nextera/Tn5 vs TruSeq-style, reads per million)
  fold-back  v2 per 1,000 (split, clip) and the old count
  noise      self_noise.py counts in callable sequence (and masked)
Writes selfcal/out/selfcal_table.tsv.
"""
import csv
import os

S = "../analysis/foldback_qc/selfcal"
O = f"{S}/out"
PB = "../analysis/novel_events/phaseB"
STATS = {"average length": "read_len_mean", "maximum length": "read_len_max",
         "insert size average": "insert_mean", "insert size standard deviation": "insert_sd",
         "error rate": "error_rate", "reads properly paired": "_proper", "sequences": "_seqs"}


def tsv(path, key):
    if not os.path.exists(path):
        return {}
    return {r[key]: r for r in csv.DictReader(open(path), delimiter="\t")}


def stats(s):
    p = f"{O}/stats/{s}.tsv"
    if not os.path.exists(p):
        return {}
    v = {}
    for line in open(p):
        k, _, x = line.rstrip("\n").partition(":\t")
        if k in STATS:
            v[STATS[k]] = x.split("\t")[0]
    if v.get("_seqs"):
        v["proper_frac"] = round(int(v["_proper"]) / int(v["_seqs"]), 4)
    return {k: x for k, x in v.items() if not k.startswith("_")}


def main():
    meta = tsv(f"{S}/genomes.tsv", "cohort_sample")
    pbm = {"mar_" + r["sample"].replace("-", "_"): r for r in
           csv.DictReader(open(f"{PB}/runs.tsv"), delimiter="\t")}
    qc = {**tsv("refbias/marinQC51/qc_table.tsv", "sample"), **tsv("refbias/selfcal54/qc_table.tsv", "sample")}
    noise = {**tsv(f"{PB}/out/qc51/self_noise.tsv", "sample"), **tsv(f"{O}/self_noise.tsv", "sample")}
    ad = {**tsv("../analysis/foldback_qc/adapters/phaseB.tsv", "sample"),
          **tsv("../analysis/foldback_qc/adapters/selfcal.tsv", "sample")}
    controls = list(tsv(f"{PB}/out/qc51/self_noise.tsv", "sample"))
    rows = []
    for s in list(meta) + controls:
        m = meta.get(s)
        r = dict(sample=s, set="selfcal54" if m else "phaseB_control")
        if m:
            r.update(genome=m["genome"], strain=m["strain"], run=m["run"], study=m["study_accession"],
                     instrument=m["instrument_model"], library_selection=m["library_selection"],
                     library_name=m["library_name"], library_protocol=m["library_construction_protocol"],
                     nominal_length=m["nominal_length"], first_public=m["first_public"])
        else:
            p = pbm.get(s, {})
            r.update(genome=noise[s]["reference"], strain=p.get("sample", ""), run=p.get("run", ""),
                     study="PRJNA719670", instrument=p.get("instrument", ""))
        q = qc.get(s, {})
        for k in ("status", "reason", "raw_depth", "final_depth", "mtbc_frac", "top_other", "top_other_frac", "lineage"):
            r[f"qc_{k}" if k in ("status", "reason") else k] = q.get(k, "")
        x = ad.get(s, {})
        for k in ("adapter_call", "nextera_per_M", "truseq_per_M"):
            r[k] = x.get(k, "")
        r.update(stats(s))
        fb = tsv(f"{O}/foldback/{s}.tsv", "sample").get(s, {})
        for k in ("foldback_per_1k", "split_per_1k", "clip_per_1k", "old_per_1k", "dup_frac", "reads"):
            r[f"fb_{k}" if k == "reads" else k] = fb.get(k, "")
        n = noise.get(s, {})
        for k, x in n.items():
            if k not in ("sample", "reference"):
                r[f"noise_{k}"] = x
        rows.append(r)
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(f"{O}/selfcal_table.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows ({sum(r['set'] == 'selfcal54' for r in rows)} selfcal54); "
          f"QC PASS {sum(r['qc_status'] == 'PASS' for r in rows)}; "
          f"fold-back measured {sum(bool(r['foldback_per_1k']) for r in rows)}; "
          f"noise counted {sum(bool(r.get('noise_hc_snp', '')) for r in rows)}")


if __name__ == "__main__":
    main()
