#!/bin/bash
# Part D (../FOLDBACK_SV_EFFECT_PLAN.md), after P2 of cohort fbrecall: each
# isolate's QC'd reads on the nearest other panel genome R (make_cohort.py).
#   1. our caller v6 + local assembly on the P2 BAM (settings as
#      phaseB/11_call_score.sbatch), mask = R's callable mask;
#   2. truth: novel_events/truth.py, the isolate's own genome against R
#      (every assembly difference of >= 50 bp is a real event).
# Scoring by fold-back band: analyse.py. Submit from runroot:
#   sbatch --dependency=afterok:<p2 summary> ../analysis/foldback_qc/recall/02_post.sbatch
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -J recall_post
#SBATCH -o ../analysis/foldback_qc/recall/logs/post_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=../analysis/foldback_qc/recall/out
B=refbias/build/7713a8d71d8e-fix1/refs
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
RM=refbias/fbrecall/run/p1/refmap.tsv
G=refbias/work/fbrecall_genomes
mkdir -p $O/proto $O/asm $O/asm_work $O/truth_work $G
export N O B E MTB_PY
awk -F'\t' 'NR>1{print $1"\t"$5}' $RM | xargs -P 8 -L 1 bash -c '
    BAM=refbias/work/fbrecall_p2/$0.bam
    "$MTB_PY" $N/breakpoint_caller.py --bam $BAM --sample $0 --ref $B/$1.fasta \
        --mask-bed $N/out/masks/$1.bed --out $O/proto/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $BAM \
        --ref $B/$1.fasta --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv'
awk -F'\t' 'NR>1{print $2}' ../analysis/foldback_qc/recall/genomes.tsv | sort -u | while read g; do
    [[ -s $G/$g.fasta.gz ]] || gzip -c $B/$g.fasta > $G/$g.fasta.gz
done
"$MTB_PY" $N/truth.py --genomes ../analysis/foldback_qc/recall/genomes.tsv --refmap $RM --refs $B \
    --asm-dir $G --minimap2 $E/minimap2 --k8 "$MTB_K8" --paftools "$MTB_PAFTOOLS" \
    --workdir $O/truth_work --out $O/truth.tsv
echo "truth events: $(($(wc -l < $O/truth.tsv) - 1))"
