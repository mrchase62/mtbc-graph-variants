#!/bin/bash
# After the foldsim P2 summary: per condition, (1) the fold-back count on the
# P2 BAM (novel_events/controls/split_reads.awk), (2) our breakpoint caller
# v3b and local assembly with the Phase B settings, (3) score.py on callable
# sequence against the source genome's Phase A truth (renamed to the
# condition), all callers. analyze.py then builds the summary.
# Submit from runroot:
#   sbatch --dependency=afterok:<p2sum> ../analysis/foldback_sim/02_call_score.sbatch
#SBATCH -p shared
#SBATCH -c 7
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -J foldsim_call
#SBATCH -o ../analysis/foldback_sim/logs/call_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
F=../analysis/foldback_sim
N=../analysis/novel_events
O=$F/out
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
RM=refbias/foldsim/run/p1/refmap.tsv
mkdir -p $O/proto $O/asm $O/asm_work $O/split $O/truth $O/score
export F N O E MTB_PY MTB_SAMTOOLS RM
awk -F'\t' 'NR>1{print $1"\t"$2"\t"$4}' $F/conditions.tsv | xargs -P 7 -L 1 bash -c '
    S=$0; SRC=$1; R=refbias/build/7713a8d71d8e-fix1/refs/$2.fasta
    B=refbias/work/foldsim_p2/$S.bam
    { printf "sample\treads\tdup_like\tinv_chimera_per_1k\tother_split_per_1k\n"
      "$MTB_SAMTOOLS" view -F 0x904 $B | awk -v s=$S -f $N/controls/split_reads.awk; } > $O/split/$S.tsv
    "$MTB_PY" $N/breakpoint_caller.py --bam $B --sample $S --out $O/proto/$S.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$S.events.tsv --bam $B \
        --ref $R --sample $S --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$S --out $O/asm/$S.events.tsv
    awk -F"\t" -v s=$SRC -v n=$S "NR==1{print;next} \$1==s{\$1=n; print}" OFS="\t" $N/out/truth.tsv > $O/truth/$S.tsv
    "$MTB_PY" $N/score.py --truth $O/truth/$S.tsv --proto-dir $O/proto --asm-dir $O/asm \
        --p2-dir refbias/foldsim/run/p2 --refmap $RM --mask-dir $N/out/masks \
        --out $O/score/$S.tsv > /dev/null'
"$MTB_PY" $F/analyze.py
