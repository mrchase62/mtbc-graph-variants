#!/bin/bash
# Fold-back chimera rate per isolate for a list of production CRAMs
# (sample<TAB>path relative to MTB_CRAM_ROOT), with the count defined in
# analysis/novel_events/controls/split_reads.awk: primary alignments whose
# first supplementary alignment (SA tag) is on the same contig, the opposite
# strand and within 1 kb, per 1,000 primary alignments. CRAMs are aligned to
# H37Rv (bwa mem 0.7.19) and keep the SA tag. Submit from runroot:
#   sbatch ../analysis/foldback_qc/01_measure.sbatch <list> <name>
# Output: analysis/foldback_qc/<name>.tsv
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=6G
#SBATCH -t 01:30:00
#SBATCH -J foldback_qc
#SBATCH -o ../analysis/foldback_qc/logs/measure_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
mtb_require_cram_root
LIST="${1:?list}"; NAME="${2:?name}"
Q=../analysis/foldback_qc
O=$(mktemp -d /tmp/foldback.$NAME.XXXXXX)
trap '[[ -n "$O" && "$O" == /tmp/foldback.* ]] && rm -rf -- "$O"' EXIT
export O MTB_SAMTOOLS MTB_H37RV MTB_CRAM_ROOT
cut -f1,2 "$LIST" | xargs -P 8 -L 1 bash -c '
    "$MTB_SAMTOOLS" view -T "$MTB_H37RV" -F 0x904 "$MTB_CRAM_ROOT/$1" \
        | awk -v s=$0 -f ../analysis/novel_events/controls/split_reads.awk > $O/$0.tsv'
{ printf "sample\treads\tdup_like\tinv_chimera_per_1k\tother_split_per_1k\n"; cat $O/*.tsv; } > $Q/$NAME.tsv
echo "$(($(wc -l < $Q/$NAME.tsv) - 1)) of $(wc -l < "$LIST") measured"
