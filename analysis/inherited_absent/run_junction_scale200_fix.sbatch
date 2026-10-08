#!/bin/bash
# junction_check.py over every scale200_fix sample, 4 at a time.
# Submit from runroot, after run_scale200_fix.sbatch. About 2 s per sample.
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 00:45:00
#SBATCH -J junction_s200fix
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
O=../analysis/inherited_absent/out/scale200_fix
mkdir -p "$O/junction"
export O MTB_PY
awk -F'\t' 'NR>1{print $1}' "$O/summary.tsv" |
  xargs -P 4 -I{} bash -c '"$MTB_PY" ../analysis/inherited_absent/junction_check.py \
      --records "$O/records.tsv" --summary "$O/summary.tsv" \
      --refmap refbias/scale200_fix/p1/refmap.tsv \
      --graph-vcf refbias/build/7713a8d71d8e-fix1/assets/graph_collapsed.vcf.gz \
      --bam refbias/work/scale200fix_p1/{}.h37rv.bam --sample {} \
      --out "$O/junction/{}.tsv"'
f=( "$O"/junction/*.tsv )
{ head -1 "${f[0]}"; for x in "${f[@]}"; do tail -n +2 "$x"; done; } > "$O/junction.tsv"
wc -l "$O/junction.tsv"
