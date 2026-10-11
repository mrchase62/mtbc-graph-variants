#!/bin/bash
# Fold-back scan v2 (foldback_v2.py) validation: the 21 synthetic foldsim
# samples (chimeras injected at known rates, refbias/work/foldsim_p2) and the
# 62 Phase B isolates on raw reads (refbias/work/marinB63_p2), each against
# its matched reference, 25% of reads by name hash. Submit from runroot:
#   sbatch ../analysis/foldback_qc/03_v2_validate.sbatch
# Output: analysis/foldback_qc/v2/{rates.tsv, sites/<sample>.tsv.gz}
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -t 01:00:00
#SBATCH -J foldback_v2
#SBATCH -o ../analysis/foldback_qc/logs/v2_validate_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
Q=../analysis/foldback_qc
V=$Q/v2
mkdir -p $V/sites $V/rows
export Q V MTB_PY
# the reference column is found by name (refmap layouts differ)
REFCOL='NR==1{for(i=1;i<=NF;i++) if($i=="reference") c=i; next} {print $1"\t"$c"\t"d}'
{ awk -F'\t' -v d=refbias/work/foldsim_p2 "$REFCOL" refbias/foldsim/run/p1/refmap.tsv
  awk -F'\t' -v d=refbias/work/marinB63_p2 "$REFCOL" refbias/marinB63/run/p1/refmap.tsv
} | xargs -P 8 -L 1 bash -c '
    "$MTB_PY" $Q/foldback_v2.py --bam $2/$0.bam --ref refbias/build/7713a8d71d8e-fix1/refs/$1.fasta \
        --sample $0 --frac 0.25 --out $V/rows/$0.tsv --sites $V/sites/$0.tsv.gz'
{ head -1 $(ls $V/rows/*.tsv | head -1); for f in $V/rows/*.tsv; do tail -n +2 $f; done; } > $V/rates.tsv
echo "$(($(wc -l < $V/rates.tsv) - 1)) samples"
