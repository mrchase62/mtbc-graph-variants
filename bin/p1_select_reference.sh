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
# One resolver for every pass: MTB_BUILD_DIR, else the single COMPLETED build.
# "The only directory under refbias/build" accepted a half-built one (P0P2-13).
BUILD="$(mtb_resolve_build)" || exit 1
export MTB_BUILD_DIR="$BUILD"
BUILD_ID="$(awk -F'\t' '$1=="build_id"{print $2}' "${BUILD}/build_info.tsv")"
[[ -n "$BUILD_ID" ]] || { echo "FATAL: no build_id in ${BUILD}/build_info.tsv" >&2; exit 1; }
[[ -s "${BUILD}/logs/refs.done" && -s "${BUILD}/logs/assets.done" ]] \
    || { echo "FATAL: P0 build ${BUILD_ID} is incomplete (refs/assets)" >&2; exit 1; }

COHORT="${COHORT:-refbias/cohort.pilot.tsv}"
CRAMMAP="${CRAMMAP:-refbias/cohort.crams.tsv}"
# The read collection is site configuration (config/site.local.sh), and the
# reference it was encoded against is MTB_CRAM_REF. This used to default to
# ${CRAMROOT}/metadata/reference.fasta -- one collection's layout -- while P1g
# and P1i already read MTB_CRAM_REF, so a new collection either failed here or
# was decoded against whatever file happened to sit at that path.
CRAMROOT="${CRAMROOT:-${MTB_CRAM_ROOT:-}}"
CRAMREF="${CRAMREF:-${MTB_CRAM_REF:-}}"
MTB_CRAM_ROOT="$CRAMROOT" mtb_require_cram_root || exit 1
[[ -n "$CRAMREF" ]] || { echo "FATAL: MTB_CRAM_REF is unset; set it to the reference the CRAMs were encoded against" >&2; exit 1; }
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
    # The tie-break's interval counts are a build asset (P0 step
    # is6110_intervals), passed explicitly; the summary refuses a candidates
    # file whose .p1.done marker names another build.
    exec "$MTB_PY" bin/p1_summary.py --cohort "$COHORT" --dir "$OUTDIR" \
        --intervals "${BUILD}/assets/is6110_intervals.tsv" \
        --build-id "$BUILD_ID" \
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

# DONE MEANS DONE AGAINST THIS BUILD (audit P0P2-1). Existence alone kept the
# candidates of whichever build's panel_snps.vcf.gz made them, so a rerun on a
# new graph into the same folders silently selected references from the old
# panel. The H37Rv VCF carries the build stamp, and <s>.p1.done records the
# build the candidates were selected against; both must name this build.
# Products from another build are refused, never overwritten or reused: rerun
# into new output folders, or move the old ones aside.
CAND="${OUTDIR}/${SAMPLE}.candidates.tsv"
H37VCF="${WORK}/${SAMPLE}.h37rv.vcf.gz"
P1DONE="${OUTDIR}/${SAMPLE}.p1.done"
_m_build="$(mtb_kv "$P1DONE" build_id)"
_v_build="$(mtb_vcf_build_id "$H37VCF")"
_foreign=""
[[ -s "$P1DONE" && "$_m_build" != "$BUILD_ID" ]] && _foreign="${P1DONE} names build '${_m_build:-none}'"
[[ -s "$H37VCF" && "$_v_build" != "$BUILD_ID" ]] && _foreign="${H37VCF} is stamped '${_v_build:-unstamped}'"
[[ -s "$CAND" && ! -s "$P1DONE" && "$_v_build" != "$BUILD_ID" ]] \
    && _foreign="${CAND} has no build marker and no H37Rv VCF of build ${BUILD_ID}"
if [[ -n "$_foreign" ]]; then
    echo "FATAL: ${SAMPLE}: P1 output from another build: ${_foreign}; this" \
         "is build ${BUILD_ID}. Rerun into new output folders (OUTDIR, WORK)" \
         "or move the old outputs aside." >&2
    exit 1
fi
if [[ -s "$CAND" && -s "$H37VCF" ]]; then
    if [[ "$_m_build" == "$BUILD_ID" ]]; then
        echo "[P1] ${SAMPLE}: already done (build ${BUILD_ID})"; exit 0
    fi
    # From before the marker: the candidates were selected in the same task
    # that wrote the H37Rv VCF, so that VCF's stamp is their build. Recorded,
    # not assumed again.
    printf 'build_id\t%s\nreference\t%s\nsource\tinferred_from_h37rv_vcf_stamp\n' \
        "$BUILD_ID" "$(awk -F'\t' 'NR==2{print $2}' "$CAND")" > "$P1DONE"
    echo "[P1] ${SAMPLE}: already done (build ${BUILD_ID}, from the H37Rv VCF stamp)"; exit 0
fi

echo "[P1] ${SAMPLE}: build ${BUILD_ID}"
FQ1="${WORK}/${SAMPLE}_1.fq"; FQ2="${WORK}/${SAMPLE}_2.fq"
# Extract to temporary names and rename both only on success. Writing straight
# into the final names meant a task killed mid-extraction (these run under 2-3
# hour limits) left a truncated FQ1 that the rerun's `-s` test accepted, and the
# isolate was then aligned from part of its reads.
if [[ ! -s "$FQ1" || ! -s "$FQ2" ]]; then
    rm -f "${FQ1}.tmp" "${FQ2}.tmp"
    "$MTB_SAMTOOLS" collate -@ 2 -u -O -T "${WORK}/${SAMPLE}.collate" \
        --reference "$CRAMREF" "$CRAM" \
      | "$MTB_SAMTOOLS" fastq -@ 2 -n -1 "${FQ1}.tmp" -2 "${FQ2}.tmp" -0 /dev/null -s /dev/null -
    mv -f "${FQ1}.tmp" "$FQ1"; mv -f "${FQ2}.tmp" "$FQ2"
fi
[[ $(( $(wc -l < "$FQ1") )) -eq $(( $(wc -l < "$FQ2") )) ]] \
    || { echo "FATAL: ${SAMPLE}: ${FQ1} and ${FQ2} have different read counts" >&2; exit 1; }
NR1=$(( $(wc -l < "$FQ1") / 4 ))
[[ "$NR1" -gt 0 ]] || { echo "FATAL: ${SAMPLE}: no reads extracted" >&2; exit 1; }
echo "[P1] ${SAMPLE}: ${NR1} read pairs"

# pass one: align to H37Rv and call, to get the SNP profile
SIM_REQUIRE_FQ=1 SIM_FQ_PREFIX="${WORK}/${SAMPLE}" SIM_KEEP_FQ=1 SIM_BAM_SUFFIX=".h37rv" \
    bash bin/simulate_and_call.sh "$H37RV" "$H37RV" "$WORK" "$SAMPLE" 1 150
mv -f "${WORK}/${SAMPLE}.vcf.gz" "${WORK}/${SAMPLE}.h37rv.vcf.gz"
mv -f "${WORK}/${SAMPLE}.vcf.gz.tbi" "${WORK}/${SAMPLE}.h37rv.vcf.gz.tbi" 2>/dev/null || true

# selection: nearest panel genome by SNP-profile distance, over the panel
# sites this isolate covers (its H37Rv BAM; D20, D24)
H37BAM="${WORK}/${SAMPLE}.h37rv.bam"
[[ -s "$H37BAM" ]] || { echo "FATAL: ${SAMPLE}: no H37Rv BAM ${H37BAM} for the coverage the selection needs" >&2; exit 1; }
REFID="$("$MTB_PY" bin/t8_select_reference.py \
    --vcf "${WORK}/${SAMPLE}.h37rv.vcf.gz" \
    --bam "$H37BAM" --samtools "$MTB_SAMTOOLS" \
    --panel-snps "$PANEL_SNPS" \
    --out "${OUTDIR}/${SAMPLE}.candidates.tsv" | tail -1)"
[[ -n "$REFID" ]] || { echo "FATAL: ${SAMPLE}: selection produced no reference" >&2; exit 1; }
[[ -s "${BUILD}/refs/${REFID}.fasta.bwt" ]] \
    || { echo "FATAL: ${SAMPLE}: chose ${REFID} but P0 has no index for it" >&2; exit 1; }
echo "[P1] ${SAMPLE}: reference ${REFID} (d=$(awk -F'\t' 'NR==2{print $3}' "${OUTDIR}/${SAMPLE}.candidates.tsv"))"
rm -f "$FQ1" "$FQ2"
printf 'build_id\t%s\nreference\t%s\ncompleted\t%s\n' "$BUILD_ID" "$REFID" "$(date -Is)" \
    > "${P1DONE}.tmp"
mv -f "${P1DONE}.tmp" "$P1DONE"
echo "[P1] ${SAMPLE}: done"
