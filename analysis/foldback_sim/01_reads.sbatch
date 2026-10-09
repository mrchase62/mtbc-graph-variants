#!/bin/bash
# Synthetic reads for the fold-back chimera study, one array task per row of
# conditions.tsv (1-21). Submit from runroot:
#   sbatch --array=1-21 ../analysis/foldback_sim/01_reads.sbatch
#   1. wgsim from the source assembly with the Phase A settings (150 bp pairs,
#      fragment 450 +/- 50, 0.2% base error, no mutations), 1,470,000 pairs
#      (about 100x). The seed is fixed per source genome, so every condition
#      of a genome starts from the same base reads.
#   2. foldback_sim.py replaces the condition's fraction of pairs with
#      fold-back chimeric pairs (model ir or random; none at 0).
#   3. Unaligned BAM for the runner (refbias/foldsim/reads), as Phase A.
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -J foldsim_reads
#SBATCH -o ../analysis/foldback_sim/logs/reads_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
F=../analysis/foldback_sim
read -r S SRC FILE REF MODEL FRAC SEED < <(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==n' $F/conditions.tsv)
O=refbias/foldsim
mkdir -p $O/reads $O/junctions
W=$(mktemp -d /tmp/foldsim.$S.XXXXXX)
trap '[[ -n "$W" && "$W" == /tmp/foldsim.* ]] && rm -rf -- "$W"' EXIT
zcat ../analysis/external_assemblies/genomes/$FILE > $W/g.fa
"$MTB_WGSIM" -N 1470000 -1 150 -2 150 -d 450 -s 50 -e 0.002 -r 0 -R 0 -X 0 -S $SEED \
    $W/g.fa $W/b1.fq $W/b2.fq > /dev/null
if [[ "$MODEL" == none ]]; then
    mv $W/b1.fq $W/r1.fq; mv $W/b2.fq $W/r2.fq
else
    "$MTB_PY" $F/foldback_sim.py --genome $W/g.fa --r1 $W/b1.fq --r2 $W/b2.fq --frac $FRAC \
        --model $MODEL --seed $SEED --out1 $W/r1.fq --out2 $W/r2.fq --junctions $O/junctions/$S.tsv
fi
"$MTB_SAMTOOLS" import -1 $W/r1.fq -2 $W/r2.fq -o $O/reads/$S.tmp.bam
mv -f $O/reads/$S.tmp.bam $O/reads/$S.bam
echo "$S done"
