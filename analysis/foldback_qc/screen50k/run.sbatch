#!/bin/bash
# Fold-back screen of the whole collection (user, 2026-10-13: "run the
# screen; wait on copying"). One array task per chunk of ~1,300 isolates from
# chunks/cNNN.tsv (sample<TAB>path relative to MTB_CRAM_ROOT): the isolates
# that pass select_isolates.py and are lineage 1-4 or 7, minus the 1,200
# already measured in ../cohorts_in_use_rates.tsv. Same count as
# ../01_measure.sbatch (split_reads.awk; CRAMs read in place, nothing copied).
# Submit from runroot:
#   sbatch --array=0-N%10 ../analysis/foldback_qc/screen50k/run.sbatch
# Output: analysis/foldback_qc/screen50k/rates/cNNN.tsv
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=6G
#SBATCH -t 01:30:00
#SBATCH -J fb_screen50k
#SBATCH -o ../analysis/foldback_qc/screen50k/logs/run_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
mtb_require_cram_root
S=../analysis/foldback_qc/screen50k
C=$(printf "c%03d" "$SLURM_ARRAY_TASK_ID")
LIST=$S/chunks/$C.tsv
mkdir -p $S/rates
O=$(mktemp -d /tmp/fbscreen.$C.XXXXXX)
trap '[[ -n "$O" && "$O" == /tmp/fbscreen.* ]] && rm -rf -- "$O"' EXIT
export O MTB_SAMTOOLS MTB_H37RV MTB_CRAM_ROOT
cut -f1,2 "$LIST" | xargs -P 8 -L 1 bash -c '
    "$MTB_SAMTOOLS" view -T "$MTB_H37RV" -F 0x904 "$MTB_CRAM_ROOT/$1" \
        | awk -v s=$0 -f ../analysis/novel_events/controls/split_reads.awk > $O/$0.tsv'
{ printf "sample\treads\tdup_like\tinv_chimera_per_1k\tother_split_per_1k\n"; cat $O/*.tsv; } > $S/rates/$C.tsv
echo "$C: $(($(wc -l < $S/rates/$C.tsv) - 1)) of $(wc -l < "$LIST") measured"
