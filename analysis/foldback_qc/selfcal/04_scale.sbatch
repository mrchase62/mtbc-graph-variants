#!/bin/bash
# Scale check (CALIBRATION_PLAN.md limit 6): the cut-off was set on QC'd
# reads mapped to the isolate's own genome; production rates come from
# H37Rv-aligned production CRAMs. For the 2 selfcal isolates that also have
# a production CRAM, measure on both:
#   old count  ../novel_events/controls/split_reads.awk (as 01_measure.sbatch)
#   v2         ../foldback_v2.py (all reads)
# and record which SRA runs the CRAM's reads come from (read-name prefixes).
# CRAMs are read in place, not copied. Submit from runroot.
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -J fb_scale
#SBATCH -o ../analysis/foldback_qc/selfcal/logs/scale_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
mtb_require_cram_root
O=../analysis/foldback_qc/selfcal/out/scale
mkdir -p $O
B=refbias/build/7713a8d71d8e-fix1/refs
while read -r bs g; do
    C=$(ls "$MTB_CRAM_ROOT"/data/*/$bs.cram)
    "$MTB_SAMTOOLS" view -T "$MTB_H37RV" -F 0x904 "$C" \
        | awk -v s=prod_$g -f ../analysis/novel_events/controls/split_reads.awk > $O/prod_$g.old.tsv
    "$MTB_SAMTOOLS" view -T "$MTB_H37RV" -F 0x904 "$C" | cut -f1 | sed 's/\..*//' | sort | uniq -c > $O/prod_$g.runs.txt
    "$MTB_PY" ../analysis/foldback_qc/foldback_v2.py --bam "$C" --ref "$MTB_H37RV" --sample prod_$g --out $O/prod_$g.v2.tsv
    "$MTB_SAMTOOLS" view -F 0x904 refbias/work/selfcal54_p2/ps_$g.bam \
        | awk -v s=self_$g -f ../analysis/novel_events/controls/split_reads.awk > $O/self_$g.old.tsv
done <<'L'
SAMD00089643 GCF_037076265
SAMEA104027311 GCF_900520315
L
head -n 50 $O/*
