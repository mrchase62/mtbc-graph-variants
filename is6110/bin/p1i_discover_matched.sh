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
#   MTB_BUILD_DIR=<build> bash is6110/bin/p1i_discover_matched.sh
#                       every panel genome of that build, for P0's
#                       is6110_intervals step (then p1i_build_matched.sh)
#   REFMAP=<refmap>     only that cohort's matched references
#
# WHY THIS EXISTS NOW AND NOT BEFORE. The pilot's eighteen references were
# discovered one at a time by hand, which was tractable at eighteen and is not
# at sixty-one. The step itself is unchanged -- one minimap2 alignment of the
# canonical element into the assembly per reference -- so this is packaging,
# not a new method. Idempotent: a reference whose GFF exists and was found in
# the same reference bytes is skipped.
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

# THE BUILD'S REFERENCES, AND EVERY PANEL GENOME BY DEFAULT -- the same
# change as p1i_build_matched.sh. REFS was build 7713a8d71d8e's refs/ (CX333),
# hard-coded, and the list was the pilot's refmap, so P0's is6110_intervals
# step could not be fed on a new build. The build is MTB_BUILD_DIR (or the
# only completed build); the references are every accession in its
# assets/accessions.txt, or the refmap's when REFMAP is set.
BUILD="$(mtb_resolve_build)" || exit 1
REFS="${REFS:-${BUILD}/refs}"
OUTDIR="${OUTDIR:-is6110/assets/matched_gff}"
ELEMENT="${ELEMENT:-is6110/assets/IS6110.query.fasta}"
if [[ -n "${REFMAP:-}" ]]; then
    mapfile -t REFLIST < <(awk -F'\t' 'NR>1{print $5}' "$REFMAP" | sort -u)
    echo "[P1i-disc] ${#REFLIST[@]} references from ${REFMAP}"
else
    mtb_require_file "${BUILD}/assets/accessions.txt" || exit 1
    mapfile -t REFLIST < <(sort -u "${BUILD}/assets/accessions.txt")
    echo "[P1i-disc] all ${#REFLIST[@]} panel genomes of ${BUILD}"
fi
mkdir -p "$OUTDIR"

n_new=0; n_skip=0; n_fail=0
for R in "${REFLIST[@]}"; do
    GFF="${OUTDIR}/${R}.is6110.gff"
    REF="${REFS}/${R}.fasta"
    [[ -s "$REF" ]] || { echo "MISSING reference ${REF}" >&2; n_fail=$((n_fail+1)); continue; }
    # which reference bytes a GFF was found in: OUTDIR is shared by builds,
    # and an accession's refs FASTA can differ between them (a new rotation,
    # a new assembly version) under the same name. A GFF with no record, or
    # from other bytes, is redone; p1i_build_matched.sh then rebuilds too.
    REFSHA="${OUTDIR}/${R}.ref.sha256"
    _sha="$(sha256sum "$REF" | cut -d' ' -f1)"
    if [[ -s "$GFF" && "$(cat "$REFSHA" 2>/dev/null)" == "$_sha" ]]; then
        n_skip=$((n_skip+1)); continue
    fi
    if "$MTB_PY" is6110/bin/is6110_discover_elements.py \
            --assembly "$REF" --name "$R" --query "$ELEMENT" \
            --minimap2 "$MTB_MINIMAP2" --out-gff "${GFF}.tmp" > /dev/null 2>&1 \
       && [[ -s "${GFF}.tmp" ]]; then
        mv -f "${GFF}.tmp" "$GFF"
        printf '%s\n' "$_sha" > "$REFSHA"
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
