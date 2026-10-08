#!/usr/bin/env bash
#SBATCH -J novelA40_sim
#SBATCH -p serial_requeue
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 0-01:00
#SBATCH -o ../analysis/novel_events/logs/sim_%A_%a.out
# Phase A, step 1: simulated reads from one external assembly
# (analysis/novel_events/genomes.tsv, row = array index + 1), as an unaligned
# BAM that the standard runner reads through the CRAM table. Settings are the
# IS6110 benchmark's (analysis/is6110_simbench/prepare.sbatch): wgsim, 150 bp
# pairs, fragment 450 +/- 50, 0.2% base error, no mutations (the assembly is
# the truth), 1,170,000 pairs (about 80x), fixed seed per genome.
# Submit from runroot:  sbatch --array=0-39 ../analysis/novel_events/01_simulate.sbatch
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source ../config/project_env.sh
G=../analysis/novel_events/genomes.tsv
read -r S ID FILE _ < <(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 2))" 'NR==n' "$G")
O=refbias/novelA40/reads
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
zcat "../analysis/external_assemblies/genomes/${FILE}" > "$T/asm.fa"
"$MTB_WGSIM" -N 1170000 -1 150 -2 150 -d 450 -s 50 -e 0.002 -r 0 -R 0 -X 0 \
    -S "$((SLURM_ARRAY_TASK_ID + 501))" "$T/asm.fa" "$T/r1.fq" "$T/r2.fq" > /dev/null
"$MTB_SAMTOOLS" import -@ 2 -1 "$T/r1.fq" -2 "$T/r2.fq" -o "$O/${S}.tmp.bam"
mv -f "$O/${S}.tmp.bam" "$O/${S}.bam"
echo "[novelA40] ${S} ${ID}: reads written"
