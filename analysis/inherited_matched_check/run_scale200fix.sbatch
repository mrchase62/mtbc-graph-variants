#!/bin/bash
# Matched-genome check of inherited SVs on scale200_fix (user, 2026-10-13:
# "yes, run the matched-genome check"). Analysis only: reads the cohort's P4b
# rows and P2 BAMs, writes under analysis/inherited_matched_check/out/.
#   sbatch analysis/inherited_matched_check/run_scale200fix.sbatch   (from the repo root)
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 04:00:00
#SBATCH -J inh_matched
#SBATCH -o analysis/inherited_matched_check/logs/%x_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source config/project_env.sh
R=runroot/refbias
"$MTB_PY" analysis/inherited_matched_check/check.py \
    --rows "$R/scale200_fix/p4b/SAM*.sv_placed.tsv" \
    --bam-dir "$R/work/scale200fix_p2" --build "$R/build/7713a8d71d8e-fix1" \
    --minimap2 "$MTB_SV_MINIMAP2" --out analysis/inherited_matched_check/out/scale200fix.tsv
