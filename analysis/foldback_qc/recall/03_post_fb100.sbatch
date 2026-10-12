#!/bin/bash
# Part D causal arm, after P2 of cohort fbrecall100: our caller v6 + local
# assembly on each <isolate>_fb100 BAM (as 02_post.sbatch), into the same
# out/proto and out/asm. The truth is the all-reads sample's (analyse.py).
# Submit from runroot:
#   sbatch --dependency=afterok:<p2 summary>:<01_masks> ../analysis/foldback_qc/recall/03_post_fb100.sbatch
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -J recall_post100
#SBATCH -o ../analysis/foldback_qc/recall/logs/post100_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=../analysis/foldback_qc/recall/out
B=refbias/build/7713a8d71d8e-fix1/refs
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
mkdir -p $O/proto $O/asm $O/asm_work
export N O B E MTB_PY
awk -F'\t' 'NR>1{print $1"\t"$5}' refbias/fbrecall100/run/p1/refmap.tsv | xargs -P 8 -L 1 bash -c '
    BAM=refbias/work/fbrecall100_p2/$0.bam
    "$MTB_PY" $N/breakpoint_caller.py --bam $BAM --sample $0 --ref $B/$1.fasta \
        --mask-bed $N/out/masks/$1.bed --out $O/proto/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $BAM \
        --ref $B/$1.fasta --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv'
