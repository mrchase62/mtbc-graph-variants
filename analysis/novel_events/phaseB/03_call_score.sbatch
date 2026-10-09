#!/bin/bash
# Phase B: the prototype (v2), local assembly of its candidates, and the
# depth scan on the 62 marinB63 P2 BAMs, 4 samples at a time; then scoring
# by group (novel, setE_polished; in_panel as the false-call check).
# Submit from runroot.
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -J marinB63_call
#SBATCH -o ../analysis/novel_events/logs/marinB63_call_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=$N/phaseB/out
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
mkdir -p $O/proto_v2 $O/proto_v2_asm $O/depth $O/asm_work
export N O E MTB_PY MTB_SAMTOOLS
awk -F'\t' 'NR>1{print $1"\t"$5}' refbias/marinB63/run/p1/refmap.tsv |
  xargs -P 4 -L 1 bash -c '
    B=refbias/work/marinB63_p2/$0.bam
    "$MTB_PY" $N/breakpoint_caller.py --bam $B --sample $0 --out $O/proto_v2/$0.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto_v2/$0.events.tsv --bam $B \
        --ref refbias/build/7713a8d71d8e-fix1/refs/$1.fasta --sample $0 \
        --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --workdir $O/asm_work/$0 --out $O/proto_v2_asm/$0.events.tsv
    "$MTB_PY" $N/depth_scan.py --bam $B --sample $0 --samtools "$MTB_SAMTOOLS" \
        --out $O/depth/$0.events.tsv'
for grp in novel setE_polished; do
    awk -F'\t' -v g=$grp 'NR==FNR{if($5==g) keep[$1]=1; next} FNR==1 || ($1 in keep)' \
        $N/phaseB/truth_genomes.tsv $O/truth.tsv > $O/truth.$grp.tsv
    "$MTB_PY" $N/score.py --truth $O/truth.$grp.tsv --proto-dir $O/proto_v2 \
        --asm-dir $O/proto_v2_asm --depth-dir $O/depth \
        --p2-dir refbias/marinB63/run/p2 --out $O/score.$grp.tsv
done
