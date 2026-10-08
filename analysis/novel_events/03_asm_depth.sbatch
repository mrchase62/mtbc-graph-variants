#!/bin/bash
# Phase A follow-up: step 1 (local assembly of prototype v2 candidates) and
# step 2 (depth scan) over the novelA40 P2 BAMs, 4 samples at a time, then
# scoring. Submit from runroot:  sbatch ../analysis/novel_events/03_asm_depth.sbatch
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:30:00
#SBATCH -J novelA40_asm
#SBATCH -o ../analysis/novel_events/logs/asm_depth_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
mkdir -p $N/out/proto_v2_asm $N/out/depth $N/out/asm_work
export N E MTB_PY MTB_SAMTOOLS
awk -F'\t' 'NR>1{print $1"\t"$5}' refbias/novelA40/run/p1/refmap.tsv |
  xargs -P 4 -L 1 bash -c '
    "$MTB_PY" $N/assemble_candidates.py --events $N/out/proto_v2/$0.events.tsv \
        --bam refbias/work/novelA40_p2/$0.bam --ref refbias/build/7713a8d71d8e-fix1/refs/$1.fasta \
        --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --workdir $N/out/asm_work/$0 --out $N/out/proto_v2_asm/$0.events.tsv
    "$MTB_PY" $N/depth_scan.py --bam refbias/work/novelA40_p2/$0.bam --sample $0 \
        --samtools "$MTB_SAMTOOLS" --out $N/out/depth/$0.events.tsv'
"$MTB_PY" $N/score.py --truth $N/out/truth.tsv --proto-dir $N/out/proto_v2 \
    --asm-dir $N/out/proto_v2_asm --depth-dir $N/out/depth \
    --p2-dir refbias/novelA40/run/p2 --out $N/out/score_v2_asm_depth.tsv
