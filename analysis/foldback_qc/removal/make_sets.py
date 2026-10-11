#!/usr/bin/env python3
"""Runner tables for the read-removal experiment (cohort fbremove;
../FOLDBACK_SV_EFFECT_PLAN.md, part C). Run from runroot with $MTB_PY.

  --stage 1  after QC (phaseB/make_qc_cohort.py has written qc_table.tsv):
             each QC PASS isolate as <isolate>_all (its QC'd unaligned BAM),
             reference pinned to its own panel genome (phaseB/genomes.tsv
             same_isolate_in_cx333); writes cohort.tsv, crams.tsv,
             run/p1/refmap.tsv and the output root's build record.
  --stage 2  after 01_sets.sbatch: adds <isolate>_fb50, _fb100 and _rand
             (reads/<isolate>_<set>.bam) to the three tables. Samples already
             there are kept, so P2 skips the finished _all samples.
"""
import argparse
import csv
import datetime
import os

C = "refbias/fbremove"
BUILD = "7713a8d71d8e-fix1"
COLS = ("sample sample_lineage sample_sublineage meandepth reference ref_lineage ref_sublineage "
        "snp_distance inside_label rank2_reference rank2_distance margin tie_break ref_intervals "
        "rank2_intervals interval_delta").split()
SETS = ("all", "fb50", "fb100", "rand")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, choices=(1, 2), required=True)
    ap.add_argument("--cram-root", default=os.environ.get("MTB_CRAM_ROOT"))
    a = ap.parse_args()
    own = {r["cohort_sample"]: r["same_isolate_in_cx333"] for r in
           csv.DictReader(open("../analysis/novel_events/phaseB/genomes.tsv"), delimiter="\t")}
    qc = [r for r in csv.DictReader(open(f"{C}/qc_table.tsv"), delimiter="\t") if r["status"] == "PASS"]
    sets = SETS[:1] if a.stage == 1 else SETS
    rows = []
    for r in qc:
        s = r["sample"]
        for k in sets:
            bam = f"{C}/reads/{s}.bam" if k == "all" else f"{C}/reads/{s}_{k}.bam"
            assert os.path.exists(bam), bam
            rows.append((f"{s}_{k}", r, own[s], os.path.relpath(os.path.realpath(bam), a.cram_root)))
    with open(f"{C}/cohort.tsv", "w") as fo:
        fo.write("sample\tlineage\tdeepest\tmeandepth\tcoverage\tmapping_rate\terror_rate\tduplicate_rate\n")
        for s, r, _, _ in rows:
            fo.write(f"{s}\t{r['lineage']}\t{r['sub_lineage']}\t{r['final_depth']}\t100\t1\t0.002\t0\n")
    with open(f"{C}/crams.tsv", "w") as fo:
        fo.write("sample\tcram_relpath\n")
        for s, _, _, p in rows:
            fo.write(f"{s}\t{p}\n")
    os.makedirs(f"{C}/run/p1", exist_ok=True)
    if not os.path.exists(f"{C}/run/.mtb_build"):
        with open(f"{C}/run/.mtb_build", "w") as fo:
            fo.write(f"build_id\t{BUILD}\nrecorded\t{datetime.datetime.now().astimezone().isoformat(timespec='seconds')}\n")
    with open(f"{C}/run/p1/refmap.tsv", "w") as fo:
        fo.write("\t".join(COLS) + "\n")
        for s, r, ref, _ in rows:
            v = dict.fromkeys(COLS, "")
            v.update(sample=s, sample_lineage=r["lineage"], sample_sublineage=r["sub_lineage"],
                     meandepth=r["final_depth"], reference=ref, ref_lineage=r["lineage"],
                     ref_sublineage=r["sub_lineage"], snp_distance="0", inside_label="pinned_own_genome")
            fo.write("\t".join(v[c] for c in COLS) + "\n")
    print(f"stage {a.stage}: {len(rows)} samples ({len(qc)} isolates passed QC)")


if __name__ == "__main__":
    main()
