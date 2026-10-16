#!/usr/bin/env bash
#SBATCH -J p5gate
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 0-02:00
#SBATCH -o slurm/p5gate_%j.out
# The cohort acceptance gate, the last pass of every chain (2026-10-16):
# bin/cohort_gate.py on the merged VCF; a FAIL exits 1, which the chain
# watchdog reports. GATE_BASELINE (a previous merged VCF) and GATE_TBPROFILER
# (a TB-Profiler JSON path pattern with {s}) are optional; GATE_ARGS passes
# further options (thresholds) through.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
for _c in "${MTB_ENV_FILE:-}" config/project_env.sh "$(dirname "${BASH_SOURCE[0]}")/../config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { source "$_c"; break; }
done
OUTDIR="${OUTDIR:?OUTDIR unset; run through bin/refbias_run.sh}"
BUILD="${MTB_BUILD_DIR:?MTB_BUILD_DIR unset}"
args=(--vcf "${OUTDIR}/merged.vcf.gz" --build "$BUILD" --out "${OUTDIR}/gate.tsv")
[[ -n "${GATE_BASELINE:-}" ]] && args+=(--baseline "$GATE_BASELINE")
[[ -n "${GATE_TBPROFILER:-}" ]] && args+=(--tbprofiler "$GATE_TBPROFILER")
read -r -a _extra <<< "${GATE_ARGS:-}"
exec "$MTB_PY" bin/cohort_gate.py "${args[@]}" "${_extra[@]}"
