#!/usr/bin/env bash
#SBATCH --job-name=P1_select
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=16000
#SBATCH --output=slurm/P1_%A_%a.out
#SBATCH --error=slurm/P1_%A_%a.err
#
# P1: choose the matched reference for each cohort isolate.
#
#   sbatch --array=1-23 bin/p1_select_reference.sh
#   bash bin/p1_select_reference.sh <sample>
#   bash bin/p1_select_reference.sh --summary
#
# This is pass one of the two-pass design. It aligns to H37Rv to obtain a SNP
# profile, then picks the nearest of the 333 panel genomes by profile distance.
# The H37Rv alignment is not waste: the profile it produces is the selection
# input, and its VCF is the direct-calling arm the pilot compares against.
#
# Selection is by DISTANCE, not by lineage label. Stage 1 measured label-gating
# as 0.6 points worse because 16.1% of nearest references lie outside the
# sample's label, so the tb-profiler call is carried as an annotation and is not
# allowed to constrain the choice.
#
# No near-identity exclusion here. Stages 1 and 2 excluded panel genomes within
# ~50 SNPs of the sample, because there the "sample" WAS a panel genome and could
# select itself, making the arm meaninglessly easy. A real isolate is not in the
# panel, so there is nothing to exclude.
#
# Everything reads from a P0 build: the H37Rv reference and its indices, and the
# panel SNP matrix. Nothing here derives a per-graph asset, and outputs inherit
# the build stamp through simulate_and_call.sh.
#
# The CRAM tree is READ-ONLY. Nothing is written outside refbias/.
set -euo pipefail

_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: no project_env.sh" >&2; exit 1; }
source "$_mtb_env"

BUILD_ROOT="${BUILD_ROOT:-refbias/build}"
BUILD="${MTB_BUILD_DIR:-}"
if [[ -z "$BUILD" ]]; then
    mapfile -t _c < <(find "$BUILD_ROOT" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
    [[ "${#_c[@]}" -eq 1 ]] || { echo "FATAL: set MTB_BUILD_DIR (${#_c[@]} builds found)" >&2; exit 1; }
    BUILD="${_c[0]}"
fi
export MTB_BUILD_DIR="$BUILD"
BUILD_ID="$(awk -F'\t' '$1=="build_id"{print $2}' "${BUILD}/build_info.tsv")"
[[ -s "${BUILD}/logs/refs.done" && -s "${BUILD}/logs/assets.done" ]] \
    || { echo "FATAL: P0 build ${BUILD_ID} is incomplete (refs/assets)" >&2; exit 1; }

COHORT="${COHORT:-refbias/cohort.pilot.tsv}"
CRAMMAP="${CRAMMAP:-refbias/cohort.crams.tsv}"
CRAMROOT="${CRAMROOT:-${MTB_CRAM_ROOT}}"
CRAMREF="${CRAMREF:-${CRAMROOT}/metadata/reference.fasta}"
H37RV="${BUILD}/refs/GCF_000195955.fasta"
PANEL_SNPS="${BUILD}/assets/panel_snps.vcf.gz"
OUTDIR="${OUTDIR:-refbias/p1}"
WORK="${WORK:-refbias/work/p1}"

for f in "$H37RV" "$PANEL_SNPS" "$COHORT" "$CRAMMAP" "$CRAMREF"; do
    [[ -r "$f" ]] || { echo "FATAL: unreadable: $f" >&2; exit 1; }
done
mkdir -p "$OUTDIR" "$WORK" slurm

# --- summary ------------------------------------------------------------------
# The step is normally $1. bin/refbias_run.sh submits through sbatch,
# which passes environment but not positional arguments, so P1STEP is
# accepted as an equivalent. $1 wins if both are given.
if [[ "${1:-${P1STEP:-}}" == "--summary" ]]; then
    # Python, not shell: IFS=$'\t' read collapses consecutive tabs because tab is
    # IFS whitespace, which silently shifted every column after an empty lineage
    # field. See bin/p1_summary.py.
    exec "$MTB_PY" bin/p1_summary.py --cohort "$COHORT" --dir "$OUTDIR" \
        --out "${OUTDIR}/refmap.tsv"
fi

# --- per-sample ---------------------------------------------------------------
if [[ $# -ge 1 ]]; then
    SAMPLE="$1"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "FATAL: no SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    SAMPLE="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$COHORT")"
fi
[[ -n "${SAMPLE:-}" ]] || { echo "FATAL: no cohort row" >&2; exit 1; }
grep -q "^${SAMPLE}	" "$COHORT" \
    || { echo "FATAL: ${SAMPLE} is not in the screened cohort ${COHORT}" >&2; exit 1; }

RELPATH="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $2; exit}' "$CRAMMAP")"
[[ -n "$RELPATH" ]] || { echo "FATAL: ${SAMPLE}: no CRAM path in ${CRAMMAP}" >&2; exit 1; }
CRAM="${CRAMROOT}/${RELPATH}"
[[ -r "$CRAM" ]] || { echo "FATAL: ${SAMPLE}: CRAM unreadable: ${CRAM}" >&2; exit 1; }

if [[ -s "${OUTDIR}/${SAMPLE}.candidates.tsv" && -s "${WORK}/${SAMPLE}.h37rv.vcf.gz" ]]; then
    echo "[P1] ${SAMPLE}: already done"; exit 0
fi

echo "[P1] ${SAMPLE}: build ${BUILD_ID}"
FQ1="${WORK}/${SAMPLE}_1.fq"; FQ2="${WORK}/${SAMPLE}_2.fq"
if [[ ! -s "$FQ1" ]]; then
    "$MTB_SAMTOOLS" collate -@ 2 -u -O -T "${WORK}/${SAMPLE}.collate" \
        --reference "$CRAMREF" "$CRAM" \
      | "$MTB_SAMTOOLS" fastq -@ 2 -n -1 "$FQ1" -2 "$FQ2" -0 /dev/null -s /dev/null -
fi
NR1=$(( $(wc -l < "$FQ1") / 4 ))
[[ "$NR1" -gt 0 ]] || { echo "FATAL: ${SAMPLE}: no reads extracted" >&2; exit 1; }
echo "[P1] ${SAMPLE}: ${NR1} read pairs"

# pass one: align to H37Rv and call, to get the SNP profile
SIM_FQ_PREFIX="${WORK}/${SAMPLE}" SIM_KEEP_FQ=1 SIM_BAM_SUFFIX=".h37rv" \
    bash bin/simulate_and_call.sh "$H37RV" "$H37RV" "$WORK" "$SAMPLE" 1 150
mv -f "${WORK}/${SAMPLE}.vcf.gz" "${WORK}/${SAMPLE}.h37rv.vcf.gz"
mv -f "${WORK}/${SAMPLE}.vcf.gz.tbi" "${WORK}/${SAMPLE}.h37rv.vcf.gz.tbi" 2>/dev/null || true

# selection: nearest panel genome by SNP-profile distance
REFID="$("$MTB_PY" bin/t8_select_reference.py \
    --vcf "${WORK}/${SAMPLE}.h37rv.vcf.gz" \
    --panel-snps "$PANEL_SNPS" \
    --out "${OUTDIR}/${SAMPLE}.candidates.tsv" | tail -1)"
[[ -n "$REFID" ]] || { echo "FATAL: ${SAMPLE}: selection produced no reference" >&2; exit 1; }
[[ -s "${BUILD}/refs/${REFID}.fasta.bwt" ]] \
    || { echo "FATAL: ${SAMPLE}: chose ${REFID} but P0 has no index for it" >&2; exit 1; }
echo "[P1] ${SAMPLE}: reference ${REFID} (d=$(awk -F'\t' 'NR==2{print $3}' "${OUTDIR}/${SAMPLE}.candidates.tsv"))"
rm -f "$FQ1" "$FQ2"
echo "[P1] ${SAMPLE}: done"
