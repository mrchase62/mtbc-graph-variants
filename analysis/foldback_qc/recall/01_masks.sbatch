#!/bin/bash
# Uncallable masks (novel_events/callable_mask.py, as 07_masks.sbatch) for the
# fbrecall matched references that have none yet (part D). Submit from runroot:
#   sbatch ../analysis/foldback_qc/recall/01_masks.sbatch
#SBATCH -p shared
#SBATCH -c 5
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -J recall_masks
#SBATCH -o ../analysis/foldback_qc/recall/logs/masks_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=$N/out/masks
export N O MTB_PY MTB_H37RV MTB_K8 MTB_PAFTOOLS
export M2=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin/minimap2
export MASK=refbias/build/7713a8d71d8e-fix1/assets/repeat_mask.bed
awk -F'\t' 'NR>1{print $5}' refbias/fbrecall/run/p1/refmap.tsv | sort -u |
  xargs -P 5 -I{} bash -c '[[ -s $O/{}.bed ]] || { "$MTB_PY" $N/callable_mask.py \
      --ref refbias/build/7713a8d71d8e-fix1/refs/{}.fasta --h37rv "$MTB_H37RV" --h37rv-mask $MASK \
      --minimap2 $M2 --k8 "$MTB_K8" --paftools "$MTB_PAFTOOLS" --out $O/{}.bed.tmp && mv $O/{}.bed.tmp $O/{}.bed; }'
echo "$(awk -F'\t' 'NR>1{print $5}' refbias/fbrecall/run/p1/refmap.tsv | sort -u | while read g; do [[ -s $O/$g.bed ]] && echo $g; done | wc -l) of 29 masks present"
