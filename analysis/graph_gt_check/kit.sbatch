#!/bin/bash
# Graph genotype check, step 1 (user, 2026-10-14: "yes, run the graph
# genotype check"): outgroup kit route A (bin/outgroup_kit.py: direct
# minimap2 asm5 alignment to H37Rv, paftools call, SVs, IS6110) on every
# non-H37Rv P200 panel genome, independent of the graph. One task per 20
# genomes. Submit from the repo root:
#   sbatch --array=0-9 analysis/graph_gt_check/kit.sbatch
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 01:30:00
#SBATCH -J ggt_kit
#SBATCH -o analysis/graph_gt_check/logs/kit_%A_%a.out
set -euo pipefail
source config/project_env.sh
O=analysis/graph_gt_check/kit
A=runroot/panel/P200/data/assemblies
tail -n +2 analysis/panel_select/panel_P200.tsv | cut -f1 | grep -v GCF_000195955 \
  | awk -v t=$SLURM_ARRAY_TASK_ID 'NR % 10 == t' | while read -r g; do
    [[ -s $O/$g.sv.tsv ]] && continue
    "$MTB_PY" bin/outgroup_kit.py --genome $A/$g.fna.gz --name $g --h37rv "$MTB_H37RV" \
        --is6110 "$MTB_DATA/annotation/IS6110.fasta" --outdir $O > $O/$g.log 2>&1 || echo "FAILED $g"
done
echo "task $SLURM_ARRAY_TASK_ID done"
