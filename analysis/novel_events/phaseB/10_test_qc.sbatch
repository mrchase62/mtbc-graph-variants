#!/bin/bash
# The 5 test isolates after read QC (cohort marinQC5): prototype v3b, local
# assembly (same settings as 06_full_v3b) and scoring on callable sequence,
# for comparison with the same isolates without QC (out/v3b). The depth scan
# was dropped 2026-10-09 (too many false calls on real reads).
# Arguments (default marinQC5 qc5): cohort, output folder under phaseB/out.
# Submit from runroot after the cohort's P2 summary:
#   sbatch ../analysis/novel_events/phaseB/10_test_qc.sbatch marinQC5b qc5b
#SBATCH -p shared
#SBATCH -c 5
#SBATCH --mem=16G
#SBATCH -t 00:40:00
#SBATCH -J qc5_call
#SBATCH -o ../analysis/novel_events/logs/qc5_call_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
C=${1:-marinQC5}
O=$N/phaseB/out/${2:-qc5}
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
RM=refbias/$C/run/p1/refmap.tsv
mkdir -p $O/proto $O/asm $O/asm_work
export C N O E MTB_PY MTB_SAMTOOLS
awk -F'\t' 'NR>1{print $1"\t"$5}' $RM | xargs -P 5 -L 1 bash -c '
    B=refbias/work/${C}_p2/$0.bam
    R=refbias/build/7713a8d71d8e-fix1/refs/$1.fasta
    "$MTB_PY" $N/breakpoint_caller.py --bam $B --sample $0 --out $O/proto/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $B \
        --ref $R --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv'
awk -F'\t' 'NR==FNR{k[$1]=1;next} FNR==1||($1 in k)' $N/phaseB/qc_test5.txt $N/phaseB/out/truth.tsv > $O/truth.tsv
"$MTB_PY" $N/score.py --truth $O/truth.tsv --proto-dir $O/proto --asm-dir $O/asm \
    --p2-dir refbias/$C/run/p2 --refmap $RM --mask-dir $N/out/masks --out $O/score_masked.tsv
