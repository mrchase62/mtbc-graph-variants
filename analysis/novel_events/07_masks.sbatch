#!/bin/bash
# Uncallable masks (callable_mask.py) for every matched reference used in
# Phase A (novelA40) and Phase B (marinB63), 5 at a time. Approved
# 2026-10-09 at about 5 billing-hours. Submit from runroot.
#SBATCH -p shared
#SBATCH -c 5
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -J novel_masks
#SBATCH -o ../analysis/novel_events/logs/novel_masks_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=$N/out/masks
mkdir -p $O
export N O MTB_PY MTB_H37RV MTB_K8 MTB_PAFTOOLS
export M2=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin/minimap2
export MASK=refbias/build/7713a8d71d8e-fix1/assets/repeat_mask.bed
cat refbias/{marinB63,novelA40}/run/p1/refmap.tsv | awk -F'\t' '$5!="reference"{print $5}' | sort -u |
  xargs -P 5 -I{} bash -c '[[ -s $O/{}.bed ]] || "$MTB_PY" $N/callable_mask.py \
      --ref refbias/build/7713a8d71d8e-fix1/refs/{}.fasta --h37rv "$MTB_H37RV" --h37rv-mask $MASK \
      --minimap2 $M2 --k8 "$MTB_K8" --paftools "$MTB_PAFTOOLS" --out $O/{}.bed.tmp && mv $O/{}.bed.tmp $O/{}.bed'
ls $O/*.bed | wc -l
