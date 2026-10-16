#!/usr/bin/env bash
#SBATCH -J assoc_tail
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-12:00
#SBATCH -o slurm/assoc_tail_%j.out
# Run cohort_assoc_tail.sh unattended (2026-10-16): the tail submits its tree
# as a job and stops, to be re-run when the tree exists. This wrapper runs the
# tail and, if it submitted a tree job, resubmits itself to start after that
# job succeeds; otherwise the tail has finished (or failed, with its status).
#   cd runroot && sbatch [--dependency=afterok:<gate job>] \
#       ../assoc/bin/cohort_assoc_tail_job.sh <cohort> <phenotype>
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
C="${1:?usage: $0 <cohort> <phenotype>}"; PH="${2:?phenotype}"
REPO=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
log="$(mktemp)"
bash "${REPO}/assoc/bin/cohort_assoc_tail.sh" "$C" "$PH" 2>&1 | tee "$log"
rc=${PIPESTATUS[0]}
J="$(grep -o 'submitted [0-9]*; re-run this script' "$log" | awk '{print $2}' | tail -1)"
rm -f "$log"
if [[ "$rc" -eq 0 && -n "$J" ]]; then
    N="$(sbatch --parsable --dependency=afterok:"$J" -J "assoc_tail_${C}" \
         -o "slurm/assoc_tail_${C}_%j.out" --export=ALL "${REPO}/assoc/bin/cohort_assoc_tail_job.sh" "$C" "$PH")"
    echo "[assoc_tail_job] tree job ${J}; this wrapper resubmitted as ${N} to finish after it"
fi
exit "$rc"
