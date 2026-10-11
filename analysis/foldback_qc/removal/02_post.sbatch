#!/bin/bash
# Read-removal experiment, step 3: after P2 of all four read sets. Per
# sample (as selfcal/03_post.sbatch): fold-back scan v2 (all reads; checks
# the removal and gives junction sites), samtools stats, our caller v6 +
# local assembly; then phaseB/self_noise.py --all-own. Submit from runroot:
#   sbatch --dependency=afterok:<p2 summary> ../analysis/foldback_qc/removal/02_post.sbatch
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH -J fbremove_post
#SBATCH -o ../analysis/foldback_qc/removal/logs/post_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=../analysis/foldback_qc/removal/out
B=refbias/build/7713a8d71d8e-fix1/refs
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
RM=refbias/fbremove/run/p1/refmap.tsv
mkdir -p $O/foldback $O/sites $O/stats $O/proto $O/asm $O/asm_work
export N O B E MTB_PY MTB_SAMTOOLS
awk -F'\t' 'NR>1{print $1"\t"$5}' $RM | xargs -P 8 -L 1 bash -c '
    BAM=refbias/work/fbremove_p2/$0.bam
    "$MTB_PY" ../analysis/foldback_qc/foldback_v2.py --bam $BAM --ref $B/$1.fasta --sample $0 \
        --out $O/foldback/$0.tsv --sites $O/sites/$0.tsv.gz
    "$MTB_SAMTOOLS" stats $BAM | grep "^SN" | cut -f2- > $O/stats/$0.tsv
    "$MTB_PY" $N/breakpoint_caller.py --bam $BAM --sample $0 --ref $B/$1.fasta \
        --mask-bed $N/out/masks/$1.bed --out $O/proto/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $BAM \
        --ref $B/$1.fasta --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv'
"$MTB_PY" $N/phaseB/self_noise.py --all-own --refmap $RM --p2 refbias/fbremove/run/p2 \
    --asm-dir $O/asm --out $O/self_noise.tsv
