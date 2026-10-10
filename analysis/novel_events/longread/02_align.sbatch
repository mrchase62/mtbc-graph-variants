#!/bin/bash
# Align one setE Nanopore run to the isolate's matched reference (marinQC51
# P1 refmap, column 6 of runs.tsv) with minimap2 -ax map-ont; sorted, indexed
# BAM. Row = array index + 2 of runs.tsv. User approved 2026-10-10.
# Submit from the repository root, after 01_download.sbatch:
#   sbatch --array=0-16%6 --dependency=aftercorr:<download job> analysis/novel_events/longread/02_align.sbatch
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-02:00
#SBATCH -J ont_aln
#SBATCH -o analysis/novel_events/logs/ont_aln_%A_%a.out
set -euo pipefail
source config/project_env.sh
D=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_ont
REFS=runroot/refbias/build/7713a8d71d8e-fix1/refs
T=analysis/novel_events/longread/runs.tsv
n=$((SLURM_ARRAY_TASK_ID + 2))
S=$(awk -F'\t' -v n="$n" 'NR==n{print $1}' "$T")
R=$(awk -F'\t' -v n="$n" 'NR==n{print $6}' "$T")
FQ="$D/$(basename "$(awk -F'\t' -v n="$n" 'NR==n{print $7}' "$T")")"
[[ -s "$FQ" && -s "$REFS/$R.fasta" ]] || { echo "FATAL: missing $FQ or $REFS/$R.fasta" >&2; exit 1; }
mkdir -p "$D/bam"
O="$D/bam/$S.bam"
"$MTB_MINIMAP2" -ax map-ont -t 4 "$REFS/$R.fasta" "$FQ" \
    | "$MTB_SAMTOOLS" sort -@ 1 -m 1G -o "$O.tmp.bam" -
mv -f "$O.tmp.bam" "$O"
"$MTB_SAMTOOLS" index "$O"
echo "$S -> $R: $("$MTB_SAMTOOLS" view -c -F 0x904 "$O") primary alignments"
