#!/bin/bash
# Phase B step 4: test the real-read changes on 5 isolates before the full run.
#   - prototype v3 (quality-trimmed clips, main-position and clip-sequence
#     agreement filters), on the 5 test isolates and, as a regression check,
#     on the 40 Phase A simulated isolates (caller only);
#   - local assembly of typed candidates and strong one-sided clusters
#     (support >= 20% of median depth), at most 40 windows per isolate;
#   - the depth scan with GC correction.
# Test isolates: mar_QC_6 (setE_polished, has v2 results), mar_TB1236 and
# mar_TB1612 (novel), mar_01_R1134 and mar_02_R1708 (in_panel controls).
# Submit from runroot.
#SBATCH -p shared
#SBATCH -c 5
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J marinB63_test_v3
#SBATCH -o ../analysis/novel_events/logs/marinB63_test_v3_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=$N/phaseB/out/test_v3
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
mkdir -p $O/proto $O/asm $O/depth $O/asm_work $N/out/proto_v3
export N O E MTB_PY MTB_SAMTOOLS
# regression check on the simulated isolates
awk -F'\t' 'NR>1{print $1}' $N/genomes.tsv |
  xargs -P 5 -I{} "$MTB_PY" $N/breakpoint_caller.py \
      --bam refbias/work/novelA40_p2/{}.bam --sample {} --out $N/out/proto_v3/{}.events.tsv
"$MTB_PY" $N/score.py --truth $N/out/truth.tsv --proto-dir $N/out/proto_v3 \
    --p2-dir refbias/novelA40/run/p2 --out $N/out/proto_v3/score.tsv
date
for s in mar_QC_6 mar_TB1236 mar_TB1612 mar_01_R1134 mar_02_R1708; do
    awk -F'\t' -v s=$s '$1==s{print $1"\t"$5}' refbias/marinB63/run/p1/refmap.tsv
done | xargs -P 5 -L 1 bash -c '
    B=refbias/work/marinB63_p2/$0.bam
    R=refbias/build/7713a8d71d8e-fix1/refs/$1.fasta
    "$MTB_PY" $N/breakpoint_caller.py --bam $B --sample $0 --out $O/proto/$0.events.tsv
    "$MTB_PY" $N/depth_scan.py --bam $B --sample $0 --samtools "$MTB_SAMTOOLS" --ref $R \
        --out $O/depth/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $B \
        --ref $R --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv
    date'
for s in mar_QC_6 mar_TB1236 mar_TB1612; do
    awk -F'\t' -v s=$s 'NR==1 || $1==s' $N/phaseB/out/truth.tsv > $O/truth.$s.tsv
    "$MTB_PY" $N/score.py --truth $O/truth.$s.tsv --proto-dir $O/proto --asm-dir $O/asm \
        --depth-dir $O/depth --p2-dir refbias/marinB63/run/p2 --out $O/score.$s.tsv
done
