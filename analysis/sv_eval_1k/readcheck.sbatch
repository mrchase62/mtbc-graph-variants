#!/bin/bash
#SBATCH -p shared -c 4 --mem=8G -t 2:00:00 -J sveval_reads
#SBATCH -o /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/sv_eval_1k/out/readcheck.%j.log
set -euo pipefail
cd /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
source config/project_env.sh
export TMPDIR=${TMPDIR:-/tmp}
$MTB_PY analysis/sv_eval_1k/readcheck.py --acc-threads 4
