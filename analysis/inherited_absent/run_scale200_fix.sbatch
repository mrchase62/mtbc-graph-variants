#!/bin/bash
# count_false_absent.py over every scale200_fix sample, 4 at a time.
# Submit from runroot. About 16 s and 1 GB per sample.
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 01:30:00
#SBATCH -J false_absent_s200fix
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
O=../analysis/inherited_absent/out/scale200_fix
mkdir -p "$O/per_sample"
export O MTB_PY MTB_SAMTOOLS
awk -F'\t' 'NR>1{print $1}' refbias/cohort.scale200.tsv |
  xargs -P 4 -I{} bash -c '"$MTB_PY" ../analysis/inherited_absent/count_false_absent.py \
      --vcf refbias/scale200_fix/p5/merged.vcf.gz \
      --bam refbias/work/scale200fix_p1/{}.h37rv.bam --sample {} \
      --samtools "$MTB_SAMTOOLS" \
      --out "$O/per_sample/{}.tsv" --summary "$O/per_sample/{}.summary.tsv"'
f=( "$O"/per_sample/*.summary.tsv )
{ head -1 "${f[0]}"; for x in "${f[@]}"; do tail -n +2 "$x"; done; } > "$O/summary.tsv"
{ head -1 "$O/per_sample/$(basename "${f[0]}" .summary.tsv).tsv"
  for x in "$O"/per_sample/*.tsv; do [[ "$x" == *.summary.tsv ]] || tail -n +2 "$x"; done; } > "$O/records.tsv"
wc -l "$O/summary.tsv" "$O/records.tsv"
