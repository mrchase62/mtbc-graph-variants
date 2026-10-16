#!/bin/bash
#SBATCH -p shared -c 4 --mem=16G -t 1:00:00 -J rdanc_pairs
#SBATCH -o /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/rd_ancestry_1k/out/pairs.%j.log
set -euo pipefail
cd /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
source config/project_env.sh
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd analysis/rd_ancestry_1k
$MTB_PY pairs.py
