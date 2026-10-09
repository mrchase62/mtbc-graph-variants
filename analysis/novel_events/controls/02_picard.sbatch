#!/bin/bash
# Picard CollectAlignmentSummaryMetrics (GATK 4.1.3 container) on the 25
# in-panel control P2 BAMs, to compare its PCT_CHIMERAS with the fold-back
# count in 01_artifacts.sbatch. Picard's definition (tool help): a read pair
# is chimeric when the mates map to different chromosomes, the insert is over
# MAX_INSERT_SIZE (100 kb), or the pair is not in the expected orientation
# (FR). Defaults kept. Reference links come from 01_artifacts.sbatch.
# Submit from runroot:
#   sbatch --array=1-25 ../analysis/novel_events/controls/02_picard.sbatch
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 00:30:00
#SBATCH -J ctl_picard
#SBATCH -o ../analysis/novel_events/logs/ctl_picard_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
N=../analysis/novel_events
O=$N/controls/out
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $N/controls/controls.txt)
R=$(awk -F'\t' -v s="$S" '$1==s{print $5}' refbias/marinB63/run/p1/refmap.tsv)
singularity exec -B /n/netscratch,/n/boslfs02 /n/boslfs02/LABS/sfortune_lab/Lab/containers/gatk_4.1.3.0.sif \
    gatk --java-options "-Xmx4g" CollectAlignmentSummaryMetrics -I refbias/work/marinB63_p2/$S.bam \
    -R $O/ref/$R/$R.fasta -O $O/$S/alignment_summary.txt 2> $O/$S/picard.log
grep -A4 "^CATEGORY" $O/$S/alignment_summary.txt
