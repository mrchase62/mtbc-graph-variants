#!/bin/bash
# Fold-back chimera rate (split_reads.awk, as in 01_artifacts.sbatch) on 500
# production CRAMs drawn at random (seed 7) from the 54,461 in
# $MTB_CRAM_ROOT/metadata/samples.json.gz -> production_sample.tsv. The
# CRAMs are aligned to H37Rv (bwa mem 0.7.19) and keep the SA tag.
# Submit from runroot:
#   sbatch ../analysis/novel_events/controls/03_production.sbatch
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -t 00:45:00
#SBATCH -J prod_chimera
#SBATCH -o ../analysis/novel_events/logs/prod_chimera_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
mtb_require_cram_root
N=../analysis/novel_events
O=$N/controls/out/production
mkdir -p $O
export N O MTB_SAMTOOLS MTB_H37RV MTB_CRAM_ROOT
awk -F'\t' 'NR>1{print $1"\t"$2}' $N/controls/production_sample.tsv | xargs -P 4 -L 1 bash -c '
    "$MTB_SAMTOOLS" view -T "$MTB_H37RV" -F 0x904 "$MTB_CRAM_ROOT/$1" \
        | awk -v s=$0 -f $N/controls/split_reads.awk > $O/$0.tsv'
{ printf "sample\treads\tdup_like\tinv_chimera_per_1k\tother_split_per_1k\n"; cat $O/*.tsv; } > $N/controls/production_chimera.tsv
wc -l $N/controls/production_chimera.tsv
