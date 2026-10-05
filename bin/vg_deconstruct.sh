#!/usr/bin/env bash
#SBATCH --job-name=vg_deconstruct
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 0-04:00
#SBATCH -p shared
#SBATCH --mem=64G
#SBATCH --output=slurm/vg_deconstruct_%j.out
#SBATCH --error=slurm/vg_deconstruct_%j.err
#
# Deconstruct a graph to VCF relative to the H37Rv reference path.
#
#   sbatch bin/vg_deconstruct.sh <graph-dir-name> [graph-file] [out.vcf]
#
# If <graph-file> is omitted the single *.smooth.final.gfa in the directory is used.
# Note: -e and -d are deprecated in vg 1.69 and are deliberately not passed.
set -euo pipefail
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

if [[ $# -lt 1 ]]; then mtb_usage "${BASH_SOURCE[0]}"; exit 2; fi

GRAPH_DIR="${MTB_GRAPHS}/$1"
[[ -d "$GRAPH_DIR" ]] || { echo "no such graph dir: $GRAPH_DIR" >&2; exit 1; }

if [[ -n "${2:-}" ]]; then
    GRAPH_HOST="${GRAPH_DIR}/$2"
else
    mapfile -t cands < <(find "$GRAPH_DIR" -maxdepth 1 -name '*.smooth.final.gfa' | sort)
    (( ${#cands[@]} == 1 )) || {
        echo "expected exactly one *.smooth.final.gfa in $GRAPH_DIR, found ${#cands[@]}" >&2
        printf '  %s\n' "${cands[@]}" >&2; exit 1; }
    GRAPH_HOST="${cands[0]}"
fi

OUT_HOST="${GRAPH_DIR}/${3:-variants.vcf}"

mtb_require_work
mtb_require_file "$MTB_VG_SIF" "$GRAPH_HOST"

GRAPH_C="$(mtb_in_graphs "$GRAPH_HOST")"
echo "[vg_deconstruct] graph: ${GRAPH_HOST}"
echo "[vg_deconstruct] ref path: ${MTB_REF_PATH}"

singularity exec $(mtb_bind) "$MTB_VG_SIF" \
    vg deconstruct -a -t "${SLURM_CPUS_PER_TASK:-8}" \
        -P "$MTB_REF_PATH" \
        "$GRAPH_C" > "$OUT_HOST"

echo "[vg_deconstruct] wrote ${OUT_HOST} ($(grep -vc '^#' "$OUT_HOST") records)"
