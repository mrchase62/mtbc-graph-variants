#!/usr/bin/env python3
"""Pinned refmap for selfcal54: every sample that passed QC (cohort.tsv, from
phaseB/make_qc_cohort.py) mapped to its own panel genome (genomes.tsv
genome), in P1's refmap layout, as the hand-made refmaps of foldsim and
marinSyn25. Also writes the output root's build record (.mtb_build), which
the runner requires before it accepts a folder that already holds a file.
Run from runroot:  $MTB_PY ../analysis/foldback_qc/selfcal/make_refmap.py
"""
import csv
import datetime
import os

C = "refbias/selfcal54"
BUILD = "7713a8d71d8e-fix1"
COLS = ("sample sample_lineage sample_sublineage meandepth reference ref_lineage ref_sublineage "
        "snp_distance inside_label rank2_reference rank2_distance margin tie_break ref_intervals "
        "rank2_intervals interval_delta").split()


def main():
    g = {r["cohort_sample"]: r["genome"] for r in
         csv.DictReader(open("../analysis/foldback_qc/selfcal/genomes.tsv"), delimiter="\t")}
    rows = list(csv.DictReader(open(f"{C}/cohort.tsv"), delimiter="\t"))
    os.makedirs(f"{C}/run/p1", exist_ok=True)
    with open(f"{C}/run/.mtb_build", "w") as fo:
        fo.write(f"build_id\t{BUILD}\nrecorded\t{datetime.datetime.now().astimezone().isoformat(timespec='seconds')}\n")
    with open(f"{C}/run/p1/refmap.tsv", "w") as fo:
        fo.write("\t".join(COLS) + "\n")
        for r in rows:
            ref = g[r["sample"]]
            assert os.path.exists(f"refbias/build/{BUILD}/refs/{ref}.fasta"), ref
            v = dict.fromkeys(COLS, "")
            v.update(sample=r["sample"], sample_lineage=r["lineage"], sample_sublineage=r["deepest"],
                     meandepth=r["meandepth"], reference=ref, ref_lineage=r["lineage"],
                     ref_sublineage=r["deepest"], snp_distance="0", inside_label="pinned_own_genome")
            fo.write("\t".join(v[c] for c in COLS) + "\n")
    print(f"{len(rows)} samples pinned to their own genome")


if __name__ == "__main__":
    main()
