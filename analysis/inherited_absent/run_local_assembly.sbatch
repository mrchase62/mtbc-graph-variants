#!/bin/bash
# local_assembly.py over one event set, one sample per process, 4 at a time.
# Submit from runroot:  sbatch run_local_assembly.sbatch calibration|unresolved
# About 7 s per event.
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J local_asm_s200fix
set -euo pipefail
SET="${1:?calibration or unresolved}"
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
O=../analysis/inherited_absent/out/scale200_fix/asm
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
mkdir -p "$O/$SET/per_sample"
export O E SET MTB_PY
awk -F'\t' 'NR>1{print $2}' "$O/events.$SET.tsv" | sort -u |
  xargs -P 4 -I{} bash -c '"$MTB_PY" ../analysis/inherited_absent/local_assembly.py \
      --events "$O/events.$SET.tsv" --sample {} \
      --bam refbias/work/scale200fix_p1/{}.h37rv.bam \
      --h37rv-mmi "$O/H37Rv.asm5.mmi" --spades "$E/python $E/spades.py" \
      --minimap2 "$E/minimap2" --workdir "$O/$SET/work" \
      --out "$O/$SET/per_sample/{}.tsv"'
f=( "$O/$SET"/per_sample/*.tsv )
{ head -1 "${f[0]}"; for x in "${f[@]}"; do tail -n +2 "$x"; done; } > "$O/$SET/results.tsv"
wc -l "$O/$SET/results.tsv"
