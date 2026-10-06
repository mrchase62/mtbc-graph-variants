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

# AN EXISTING TABLE IS KEPT ONLY IF IT WAS MADE AGAINST THIS BUILD'S
# CATALOGUE (audit P3IS-8). Locus ids are ACC_<H37Rv anchor>, and on a new
# graph the catalogue is rebuilt: an old table merged by id would describe a
# different sequence, with no warning. The catalogue locus_presence.py reads
# (accessory/assets/, its default) must be byte-identical to the build's copy
# (P0 step assets), and each table's sidecar <table>.build records the build
# id and catalogue checksum it was made with. A table from another build or
# catalogue is refused, not overwritten: use a new OUTDIR or move it aside.
BUILD="$(mtb_resolve_build)" || exit 1
BUILD_ID="$(awk -F'\t' '$1=="build_id"{print $2}' "${BUILD}/build_info.tsv")"
CAT="${ACC_CATALOGUE:-accessory/assets/accessory_catalogue}"
for ext in tsv fasta; do
    [[ -s "${BUILD}/assets/accessory_catalogue.${ext}" ]] || {
        echo "FATAL: build ${BUILD_ID} has no accessory_catalogue.${ext}; run bin/p0_prepare.sh --step assets" >&2; exit 1; }
    cmp -s "${CAT}.${ext}" "${BUILD}/assets/accessory_catalogue.${ext}" || {
        echo "FATAL: ${CAT}.${ext}, which locus_presence.py reads, differs from build ${BUILD_ID}'s copy" >&2; exit 1; }
done
CAT_SHA="$(cat "${CAT}.tsv" "${CAT}.fasta" | sha256sum | cut -c1-16)"
SIDE="${OUT}.build"
if [[ -s "$OUT" ]]; then
    _b="$(mtb_kv "$SIDE" build_id)"; _c="$(mtb_kv "$SIDE" catalogue_sha)"
    if [[ "$_b" == "$BUILD_ID" && "$_c" == "$CAT_SHA" ]]; then
        echo "already done: $OUT (build ${BUILD_ID})"; exit 0
    fi
    # From before the sidecar: kept only if written after the current
    # catalogue file, so it cannot predate a catalogue rebuild.
    if [[ ! -e "$SIDE" && "$OUT" -nt "${CAT}.tsv" && "$OUT" -nt "${CAT}.fasta" ]]; then
        printf 'build_id\t%s\ncatalogue_sha\t%s\nsource\tinferred_newer_than_catalogue\n' \
            "$BUILD_ID" "$CAT_SHA" > "$SIDE"
        echo "already done: $OUT (newer than build ${BUILD_ID}'s catalogue; recorded)"; exit 0
    fi
    echo "FATAL: ${OUT} was made against build '${_b:-unrecorded}', catalogue" \
         "'${_c:-unrecorded}', not ${BUILD_ID}/${CAT_SHA}. Use a new OUTDIR or" \
         "move it aside." >&2
    exit 1
fi
OUTDIR="$OUTDIR" bash accessory/bin/locus_presence_one.sh "$S" "$COHORT_TAG"
[[ -s "$OUT" ]] || { echo "FATAL: no ${OUT} written" >&2; exit 1; }
printf 'build_id\t%s\ncatalogue_sha\t%s\ncompleted\t%s\n' "$BUILD_ID" "$CAT_SHA" "$(date -Is)" > "$SIDE"
