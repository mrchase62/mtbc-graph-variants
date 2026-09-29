#!/usr/bin/env bash
# Print resolved paths and check that everything the pipeline needs is present.
set -uo pipefail
# --- locate config/project_env.sh -----------------------------------------
# Under sbatch, BASH_SOURCE[0] is Slurm's spool copy of this script rather than
# the file in bin/, so the relative lookup alone is not enough. Fall back to the
# submit directory and then to an already-exported MTB_WORK.
_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/../config/project_env.sh" \
          "${MTB_WORK:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || {
    echo "FATAL: cannot locate config/project_env.sh." >&2
    echo "       Submit from the project root, or export MTB_ENV_FILE." >&2
    exit 1
}
source "$_mtb_env"

mtb_show_config
echo
echo "--- presence check"
rc=0
for f in "$MTB_PGGB_SIF" "$MTB_VG_SIF" "$MTB_REF_FASTA" "$MTB_H37RV"; do
    if [[ -e "$f" ]]; then printf '  OK      %s\n' "$f"
    else printf '  MISSING %s\n' "$f"; rc=1; fi
done
for f in "$MTB_SNPEFF_JAR" "$MTB_JAVA"; do
    # Optional: only needed for bin/run_snpEff.sh, so warn rather than fail.
    if [[ -e "$f" ]]; then printf '  OK      %s\n' "$f"
    else printf '  WARN    %s (snpEff annotation unavailable)\n' "$f"; fi
done
for d in "$MTB_DATA" "$MTB_GRAPHS" "$MTB_FASTAS"; do
    if [[ -d "$d" ]]; then printf '  OK      %s/\n' "$d"
    else printf '  MISSING %s/\n' "$d"; rc=1; fi
done
echo
if [[ $rc -ne 0 ]]; then
    echo "Some inputs are absent — run: ${MTB_BIN}/stage_in.sh"
else
    echo "All required inputs present."
fi
exit $rc
