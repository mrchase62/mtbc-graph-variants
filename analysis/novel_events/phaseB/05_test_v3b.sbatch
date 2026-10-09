#!/bin/bash
# Phase B step 4, second pass: the caller with sliding-window clip trimming
# (v3b), on the 5 test isolates and the 40 Phase A simulations; caller and
# local assembly only (the depth scan is unchanged from 04_test_v3).
# Submit from runroot.
#SBATCH -p shared
#SBATCH -c 5
#SBATCH --mem=16G
#SBATCH -t 00:40:00
#SBATCH -J marinB63_test_v3b
#SBATCH -o ../analysis/novel_events/logs/marinB63_test_v3b_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=$N/phaseB/out/test_v3b
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
mkdir -p $O/proto $O/asm $O/asm_work $N/out/proto_v3b
export N O E MTB_PY MTB_SAMTOOLS
awk -F'\t' 'NR>1{print $1}' $N/genomes.tsv |
  xargs -P 5 -I{} "$MTB_PY" $N/breakpoint_caller.py \
      --bam refbias/work/novelA40_p2/{}.bam --sample {} --out $N/out/proto_v3b/{}.events.tsv
"$MTB_PY" $N/score.py --truth $N/out/truth.tsv --proto-dir $N/out/proto_v3b \
    --p2-dir refbias/novelA40/run/p2 --out $N/out/proto_v3b/score.tsv
for s in mar_QC_6 mar_TB1236 mar_TB1612 mar_01_R1134 mar_02_R1708; do
    awk -F'\t' -v s=$s '$1==s{print $1"\t"$5}' refbias/marinB63/run/p1/refmap.tsv
done | xargs -P 5 -L 1 bash -c '
    B=refbias/work/marinB63_p2/$0.bam
    "$MTB_PY" $N/breakpoint_caller.py --bam $B --sample $0 --out $O/proto/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $B \
        --ref refbias/build/7713a8d71d8e-fix1/refs/$1.fasta --sample $0 \
        --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv'
for s in mar_QC_6 mar_TB1236 mar_TB1612; do
    "$MTB_PY" $N/score.py --truth $N/phaseB/out/test_v3/truth.$s.tsv --proto-dir $O/proto \
        --asm-dir $O/asm --depth-dir $N/phaseB/out/test_v3/depth \
        --p2-dir refbias/marinB63/run/p2 --out $O/score.$s.tsv
done
