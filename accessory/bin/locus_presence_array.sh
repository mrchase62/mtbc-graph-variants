#!/usr/bin/env bash
#SBATCH --job-name=acc_presence
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-00:30
#SBATCH -p shared
#SBATCH --mem=4000
#SBATCH --output=slurm/accpres_%A_%a.out
#SBATCH --error=slurm/accpres_%A_%a.err
#
# LEVEL 1 across a cohort. One task per sample, 9.5 CPU-s each measured, so a
# 200-sample cohort is ~32 CPU-min and a 1000-sample cohort ~2.6 CPU-hours.
#
# Submitted by bin/refbias_run.sh as part of pass p3, so every cohort has the
# accessory presence tables p5_finish.sh merges into the VCF. They used to be
# run by hand, and only scale200 and gwas1000 had them.
#
#   sbatch --array=1-N accessory/bin/locus_presence_array.sh
#
# The sample is the task's refmap row (REFMAP, from the runner), or the task's
# line of LIST when LIST is set. An existing table is kept: the presence call
# reads only the sample's reads and the accessory catalogue, and its output is
# byte-identical between the code before and after the 2026-09-29 review.
set -euo pipefail
_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _env="$_c"; break; }
done
[[ -n "$_env" ]] || { echo "FATAL: no project_env.sh" >&2; exit 1; }
source "$_env"
cd "${SLURM_SUBMIT_DIR:-.}"

COHORT_TAG="${COHORT_TAG:?set COHORT_TAG}"
I="${SLURM_ARRAY_TASK_ID:?an array step}"
if [[ -n "${LIST:-}" ]]; then
    S="$(sed -n "$((I + 1))p" "$LIST")"
else
    S="$(awk -F'\t' -v n="$((I + 1))" 'NR==n{print $1}' "${REFMAP:?set REFMAP or LIST}")"
fi
[[ -n "$S" ]] || { echo "FATAL: no sample at index ${I}" >&2; exit 1; }
OUTDIR="${OUTDIR:-accessory/${COHORT_TAG}}"
OUT="${OUTDIR}/${S}.presence.tsv"
[[ -s "$OUT" ]] && { echo "already done: $OUT"; exit 0; }
OUTDIR="$OUTDIR" exec bash accessory/bin/locus_presence_one.sh "$S" "$COHORT_TAG"
