#!/bin/bash
# After the marinSyn25 P2 summary: our caller v3b + local assembly on the
# synthetic controls (Phase B settings, as 06_full_v3b), then compare.py.
# Submit from runroot:
#   sbatch --dependency=afterok:<p2sum> ../analysis/novel_events/realsyn/02_call_compare.sbatch
#SBATCH -p shared
#SBATCH -c 6
#SBATCH --mem=20G
#SBATCH -t 01:00:00
#SBATCH -J realsyn_call
#SBATCH -o ../analysis/novel_events/realsyn/logs/call_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=$N/realsyn/out
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
mkdir -p $O/proto $O/asm $O/asm_work
export N O E MTB_PY
awk -F'\t' 'NR>1{print $1"\t"$5}' refbias/marinSyn25/run/p1/refmap.tsv | xargs -P 6 -L 1 bash -c '
    B=refbias/work/marinSyn25_p2/$0.bam
    R=refbias/build/7713a8d71d8e-fix1/refs/$1.fasta
    "$MTB_PY" $N/breakpoint_caller.py --bam $B --sample $0 --out $O/proto/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $B \
        --ref $R --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv'
"$MTB_PY" $N/realsyn/compare.py
column -t $O/compare.tsv | cut -c1-250
