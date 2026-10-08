#!/bin/bash
# Prototype breakpoint caller over the novelA40 P2 BAMs, 4 at a time.
# Submit from runroot:  sbatch ../analysis/novel_events/02_call.sbatch <version>
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -J novelA40_call
#SBATCH -o ../analysis/novel_events/logs/call_%j.out
set -euo pipefail
V="${1:?version label}"
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
O=../analysis/novel_events/out/proto_$V
mkdir -p "$O"
export O MTB_PY
awk -F'\t' 'NR>1{print $1}' ../analysis/novel_events/genomes.tsv |
  xargs -P 4 -I{} "$MTB_PY" ../analysis/novel_events/breakpoint_caller.py \
      --bam refbias/work/novelA40_p2/{}.bam --sample {} --out "$O/{}.events.tsv"
"$MTB_PY" ../analysis/novel_events/score.py --truth ../analysis/novel_events/out/truth.tsv \
    --proto-dir "$O" --p2-dir refbias/novelA40/run/p2 --out "$O/score.tsv"
