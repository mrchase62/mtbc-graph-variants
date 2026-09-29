#!/usr/bin/env bash
# LEVEL 1 across a cohort. One task per sample, 9.5 CPU-s each measured, so a
# 200-sample cohort is ~32 CPU-min and a 1000-sample cohort ~2.6 CPU-hours.
set -euo pipefail
source config/project_env.sh 2>/dev/null || true
COHORT_TAG="${COHORT_TAG:?set COHORT_TAG}"
LIST="${LIST:?set LIST}"
S="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST")"
[[ -n "$S" ]] || { echo "no sample at index ${SLURM_ARRAY_TASK_ID}"; exit 0; }
OUT="accessory/${COHORT_TAG}/${S}.presence.tsv"
[[ -s "$OUT" ]] && { echo "already done: $OUT"; exit 0; }
exec bash accessory/bin/locus_presence_one.sh "$S" "$COHORT_TAG"
