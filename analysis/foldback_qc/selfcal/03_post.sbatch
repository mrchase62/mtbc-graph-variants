#!/bin/bash
# Fold-back calibration, after P2 of selfcal54 (each isolate on its own
# genome). Per sample, on the P2 BAM (QC'd reads, downsampled to 100x):
#   1. fold-back scan v2 (../foldback_v2.py, 25% of reads by name hash, as
#      the v2 validation);
#   2. read metrics: samtools stats SN lines (read length, insert size
#      mean/sd, error rate, duplicates, pairing);
#   3. our SV caller v6 + local assembly (as phaseB/11_call_score.sbatch).
# Steps 1-2 also for the 11 Phase B own-genome controls on their QC'd
# marinQC51 BAMs, so every point uses the same read preparation.
# Then phaseB/self_noise.py --all-own and collect.py. Submit from runroot:
#   sbatch --dependency=afterok:<p2 summary> ../analysis/foldback_qc/selfcal/03_post.sbatch
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -J selfcal_post
#SBATCH -o ../analysis/foldback_qc/selfcal/logs/post_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
S=../analysis/foldback_qc/selfcal
O=$S/out
B=refbias/build/7713a8d71d8e-fix1/refs
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
RM=refbias/selfcal54/run/p1/refmap.tsv
mkdir -p $O/foldback $O/sites $O/stats $O/proto $O/asm $O/asm_work
export N S O B E MTB_PY MTB_SAMTOOLS
# sample, reference, BAM folder, whether to run the SV caller
{ awk -F'\t' 'NR>1{print $1"\t"$5"\trefbias/work/selfcal54_p2\t1"}' $RM
  awk -F'\t' 'NR>1{print $1"\t"$2"\trefbias/work/marinQC51_p2\t0"}' $N/phaseB/out/qc51/self_noise.tsv
} | xargs -P 8 -L 1 bash -c '
    BAM=$2/$0.bam
    "$MTB_PY" ../analysis/foldback_qc/foldback_v2.py --bam $BAM --ref $B/$1.fasta --sample $0 \
        --frac 0.25 --out $O/foldback/$0.tsv --sites $O/sites/$0.tsv.gz
    "$MTB_SAMTOOLS" stats $BAM | grep "^SN" | cut -f2- > $O/stats/$0.tsv
    if [[ $3 == 1 ]]; then
        "$MTB_PY" $N/breakpoint_caller.py --bam $BAM --sample $0 --ref $B/$1.fasta \
            --mask-bed $N/out/masks/$1.bed --out $O/proto/$0.events.tsv
        "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $BAM \
            --ref $B/$1.fasta --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
            --bnd-min-frac 0.2 --max-windows 40 \
            --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv
    fi'
"$MTB_PY" $N/phaseB/self_noise.py --all-own --refmap $RM --p2 refbias/selfcal54/run/p2 \
    --asm-dir $O/asm --out $O/self_noise.tsv
"$MTB_PY" $S/collect.py
