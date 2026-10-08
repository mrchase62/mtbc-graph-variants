#!/usr/bin/env bash
#SBATCH -J marinB63_prep
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 0-01:00
#SBATCH -o ../analysis/novel_events/logs/prepB_%A_%a.out
# Phase B: one Marin 2024 Illumina pair -> an unaligned BAM the standard
# runner reads through the CRAM table (as novelA40). Runs deeper than about
# 120x are subsampled to about 120x (samtools view -s, seed 7), so that 28
# runs of up to 1,150x do not dominate the cost and every sample is at a
# similar depth. Row = array index + 2 of analysis/novel_events/phaseB/genomes.tsv.
# Submit from runroot:  sbatch --array=0-62 ../analysis/novel_events/phaseB/02_prepare.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ../config/project_env.sh
G=../analysis/novel_events/phaseB/genomes.tsv
n=$((SLURM_ARRAY_TASK_ID + 2))
col() { awk -F'\t' -v n="$n" -v c="$1" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} NR==n{print $h[c]}' "$G"; }
S=$(col cohort_sample); R1=$(col r1); R2=$(col r2); F=$(col downsample)
IN=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_illumina
O=refbias/marinB63/reads
"$MTB_SAMTOOLS" import -@ 2 -1 "$IN/$R1" -2 "$IN/$R2" -o "$O/${S}.full.tmp.bam"
if awk -v f="$F" 'BEGIN{exit !(f < 1)}'; then
    "$MTB_SAMTOOLS" view -@ 2 -b -s "7${F#0}" -o "$O/${S}.tmp.bam" "$O/${S}.full.tmp.bam"
    rm -f "$O/${S}.full.tmp.bam"
else
    mv -f "$O/${S}.full.tmp.bam" "$O/${S}.tmp.bam"
fi
mv -f "$O/${S}.tmp.bam" "$O/${S}.bam"
echo "[marinB63] $S: $R1 $R2 fraction $F"
