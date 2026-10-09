#!/bin/bash
# Phase B rerun on clean data (cohort marinQC51: read QC, kraken2, fold-back
# cutoff 2 per 1,000, fixed P2): our caller v3b + local assembly on every
# isolate in the cohort's refmap (settings as 06_full_v3b), then score.py on
# callable sequence for each group (novel, setE_polished, in_panel), with the
# truth restricted to isolates in the cohort. Submit from runroot after the
# P2 summary:
#   sbatch --dependency=afterok:<p2sum> ../analysis/novel_events/phaseB/11_call_score.sbatch
#SBATCH -p shared
#SBATCH -c 6
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -J qc51_call
#SBATCH -o ../analysis/novel_events/logs/qc51_call_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
C=marinQC51
O=$N/phaseB/out/qc51
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
RM=refbias/$C/run/p1/refmap.tsv
mkdir -p $O/proto $O/asm $O/asm_work
export C N O E MTB_PY
awk -F'\t' 'NR>1{print $1"\t"$5}' $RM | xargs -P 6 -L 1 bash -c '
    B=refbias/work/${C}_p2/$0.bam
    R=refbias/build/7713a8d71d8e-fix1/refs/$1.fasta
    "$MTB_PY" $N/breakpoint_caller.py --bam $B --sample $0 --out $O/proto/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $B \
        --ref $R --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv'
for grp in novel setE_polished in_panel; do
    awk -F'\t' -v g=$grp 'FILENAME==ARGV[1]{if(FNR>1) inc[$1]=1; next}
        FILENAME==ARGV[2]{if($5==g && ($1 in inc)) keep[$1]=1; next}
        FNR==1 || ($1 in keep)' $RM $N/phaseB/truth_genomes.tsv $N/phaseB/out/truth.tsv > $O/truth.$grp.tsv
    "$MTB_PY" $N/score.py --truth $O/truth.$grp.tsv --proto-dir $O/proto --asm-dir $O/asm \
        --p2-dir refbias/$C/run/p2 --refmap $RM --mask-dir $N/out/masks --out $O/score_masked.$grp.tsv
done
