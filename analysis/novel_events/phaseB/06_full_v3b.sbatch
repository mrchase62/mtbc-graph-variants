#!/bin/bash
# Phase B full run on the v3b settings (approved 2026-10-09, about 9
# billing-hours): the 57 isolates not in the 5-isolate test (their outputs were
# copied into out/v3b from test_v3b, and test_v3 for the depth scan).
#   - prototype v3b (sliding-window clip trim, main-position and clip-sequence
#     agreement filters);
#   - local assembly of typed candidates and one-sided clusters with support
#     >= 20% of median depth, at most 40 windows per isolate;
#   - the GC-corrected depth scan, scored separately (not trusted on real
#     reads: user's option 1).
# Then scoring for each group; in_panel gives the false-call check.
# Submit from runroot.
#SBATCH -p shared
#SBATCH -c 5
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J marinB63_full_v3b
#SBATCH -o ../analysis/novel_events/logs/marinB63_full_v3b_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=$N/phaseB/out/v3b
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
mkdir -p $O/proto $O/asm $O/depth $O/asm_work
export N O E MTB_PY MTB_SAMTOOLS
awk -F'\t' 'NR>1{print $1"\t"$5}' refbias/marinB63/run/p1/refmap.tsv |
  while IFS=$'\t' read -r s r; do [[ -s $O/asm/$s.events.tsv ]] || printf '%s\t%s\n' "$s" "$r"; done |
  xargs -P 5 -L 1 bash -c '
    B=refbias/work/marinB63_p2/$0.bam
    R=refbias/build/7713a8d71d8e-fix1/refs/$1.fasta
    "$MTB_PY" $N/breakpoint_caller.py --bam $B --sample $0 --out $O/proto/$0.events.tsv
    "$MTB_PY" $N/depth_scan.py --bam $B --sample $0 --samtools "$MTB_SAMTOOLS" --ref $R \
        --out $O/depth/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/$0.events.tsv --bam $B \
        --ref $R --sample $0 --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/$0 --out $O/asm/$0.events.tsv'
for grp in novel setE_polished in_panel; do
    awk -F'\t' -v g=$grp 'NR==FNR{if($5==g) keep[$1]=1; next} FNR==1 || ($1 in keep)' \
        $N/phaseB/truth_genomes.tsv $N/phaseB/out/truth.tsv > $O/truth.$grp.tsv
    "$MTB_PY" $N/score.py --truth $O/truth.$grp.tsv --proto-dir $O/proto \
        --asm-dir $O/asm --depth-dir $O/depth \
        --p2-dir refbias/marinB63/run/p2 --out $O/score.$grp.tsv
done
