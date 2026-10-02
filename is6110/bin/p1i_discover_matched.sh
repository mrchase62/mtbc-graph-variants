#!/usr/bin/env bash
#SBATCH --job-name=P1i_discover
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 0-01:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH --output=slurm/P1idisc_%A_%a.out
#SBATCH --error=slurm/P1idisc_%A_%a.err
#
# P1i stage 1, as a driver. Discover the IS6110 intervals in each distinct
# SNP-matched reference and write them as a GFF, which stage 2 excises.
#
#   REFMAP=refbias/scale/p1/refmap.tsv bash is6110/bin/p1i_discover_matched.sh
#
# WHY THIS EXISTS NOW AND NOT BEFORE. The pilot's eighteen references were
# discovered one at a time by hand, which was tractable at eighteen and is not
# at sixty-one. The step itself is unchanged -- one minimap2 alignment of the
# canonical element into the assembly per reference -- so this is packaging,
# not a new method. Idempotent: a reference whose GFF exists is skipped.
#
# NOTE ON THE GRADE. Stage 1 returns 14 of 16 boundaries base-exact on H37Rv,
# the two misses being copies longer than the 1355 bp query, which leaves 22 bp
# of 21,700. That gate was written as no-go and the pilot proceeded past it by
# decision. The same decision is inherited here and is not re-taken silently;
# see is6110/docs/P1I_STAGE1.md section 4.
set -euo pipefail

for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
source "$_mtb_env"
cd "${SLURM_SUBMIT_DIR:-.}"

REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
REFS="${REFS:-refbias/build/7713a8d71d8e/refs}"
OUTDIR="${OUTDIR:-is6110/assets/matched_gff}"
ELEMENT="${ELEMENT:-is6110/assets/IS6110.query.fasta}"
mkdir -p "$OUTDIR"

n_new=0; n_skip=0; n_fail=0
for R in $(awk -F'\t' 'NR>1{print $5}' "$REFMAP" | sort -u); do
    GFF="${OUTDIR}/${R}.is6110.gff"
    if [[ -s "$GFF" ]]; then n_skip=$((n_skip+1)); continue; fi
    REF="${REFS}/${R}.fasta"
    [[ -s "$REF" ]] || { echo "MISSING reference ${REF}" >&2; n_fail=$((n_fail+1)); continue; }
    if "$MTB_PY" is6110/bin/is6110_discover_elements.py \
            --assembly "$REF" --name "$R" --query "$ELEMENT" \
            --minimap2 "$MTB_MINIMAP2" --out-gff "${GFF}.tmp" > /dev/null 2>&1 \
       && [[ -s "${GFF}.tmp" ]]; then
        mv -f "${GFF}.tmp" "$GFF"
        N=$(grep -vc '^#' "$GFF" || true)
        echo "[P1i-disc] ${R}: ${N} intervals$([[ "$N" -eq 0 ]] && echo '  (none: a genome carrying no copies is legitimate, not a failure)')"
        n_new=$((n_new+1))
    else
        rm -f "${GFF}.tmp"
        echo "[P1i-disc] ${R}: FAILED" >&2
        n_fail=$((n_fail+1))
    fi
done
echo "[P1i-disc] ${n_new} discovered, ${n_skip} already present, ${n_fail} failed"
[[ "$n_fail" -eq 0 ]]
