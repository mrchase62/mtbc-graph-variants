#!/usr/bin/env bash
#SBATCH --job-name=P2_call
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -t 0-03:00
#SBATCH -p shared
#SBATCH --mem=16000
#SBATCH --output=slurm/P2_%A_%a.out
#SBATCH --error=slurm/P2_%A_%a.err
#
# P2: align each isolate to ITS matched reference and call, in R coordinates.
#
#   sbatch --array=1-23 bin/p2_call.sh
#   bash bin/p2_call.sh <sample>
#   bash bin/p2_call.sh --summary
#
# Pass two of the two-pass design. P1 chose the reference; this aligns to it and
# produces the callsets the later stages place into the common frame:
#
#   small variants   GATK HaplotypeCaller, haploid, via simulate_and_call.sh
#   structural       delly and dysgu on the same BAM
#
# Both callsets stay in R coordinates. Nothing is lifted or composed here --
# T16 showed composition is the defect and that no region rule separates its
# errors, so the frame work belongs in P4 under T15's repeat-mask routing, not
# smuggled into the caller stage.
#
# A GVCF is produced as well. Three-state genotypes need to distinguish "reads
# covered this and matched the reference" from "nothing was observed here", and a
# plain VCF cannot: stage 4 found 464 of 1131 inherited differences rejected
# spuriously when MIN_DP was read instead of DP. The GVCF is the only record that
# keeps the distinction, and it cannot be reconstructed later.
#
# The BAM is retained. SV calling needs it, depth genotyping of the 674
# copy-number loci is still untested and will need it, and re-aligning 23
# isolates to recover it costs more than storing it.
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
# One resolver for every pass: MTB_BUILD_DIR, else the single COMPLETED build
# (audit P0P2-13).
BUILD="$(mtb_resolve_build)" || exit 1
export MTB_BUILD_DIR="$BUILD"
BUILD_ID="$(awk -F'\t' '$1=="build_id"{print $2}' "${BUILD}/build_info.tsv")"
[[ -n "$BUILD_ID" ]] || { echo "FATAL: no build_id in ${BUILD}/build_info.tsv" >&2; exit 1; }
[[ -s "${BUILD}/logs/refs.done" ]] || { echo "FATAL: P0 ${BUILD_ID} refs incomplete" >&2; exit 1; }

REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
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
OUTDIR="${OUTDIR:-refbias/p2}"
WORK="${WORK:-refbias/work/p2}"
# P1's own H37Rv VCFs, for the summary's burden comparison. refbias_run.sh
# exports P1WORK per cohort; without it the summary read the pilot's P1 folder
# and left h37rv_small blank for every isolate not in the pilot (audit P0P2-8).
P1WORK="${P1WORK:-refbias/work/p1}"
DELLY_ENV="${DELLY_ENV:-${MTB_DELLY_ENV:?MTB_DELLY_ENV is unset; see config/project_env.sh}}"
DYSGU="${DYSGU:-refbias/work/svvenv/bin/dysgu}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

[[ -s "$REFMAP" ]] || { echo "FATAL: no P1 refmap at ${REFMAP}; run P1 first" >&2; exit 1; }
# A refmap selected against another build's panel is refused. p1_summary.py
# writes <refmap>.build; a hand-made refmap (pinned-reference arms) has none.
_rm_build="$(mtb_kv "${REFMAP}.build" build_id)"
if [[ -n "$_rm_build" && "$_rm_build" != "$BUILD_ID" ]]; then
    echo "FATAL: ${REFMAP} was selected against build ${_rm_build}, not ${BUILD_ID}" >&2
    exit 1
fi
mkdir -p "$OUTDIR" "$WORK" slurm

# The step is normally $1. bin/refbias_run.sh submits through sbatch,
# which passes environment but not positional arguments, so P2STEP is
# accepted as an equivalent. $1 wins if both are given.
if [[ "${1:-${P2STEP:-}}" == "--summary" ]]; then
    exec "$MTB_PY" bin/p2_summary.py --refmap "$REFMAP" --dir "$OUTDIR" \
        --work "$WORK" --p1-work "$P1WORK" --out "${OUTDIR}/p2_summary.tsv"
fi

if [[ $# -ge 1 ]]; then
    SAMPLE="$1"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "FATAL: no SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    SAMPLE="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$REFMAP")"
fi
[[ -n "${SAMPLE:-}" ]] || { echo "FATAL: no refmap row" >&2; exit 1; }

# awk, not `IFS=$'\t' read`: tab is IFS whitespace so consecutive tabs collapse
# and an empty field shifts every later column. That silently put a mean depth in
# a sub-lineage column in P1's first summary.
REFID="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $5; exit}' "$REFMAP")"
[[ -n "$REFID" ]] || { echo "FATAL: ${SAMPLE}: no reference in ${REFMAP}" >&2; exit 1; }
# REFSDIR lets an arm point at a different set of prepared references without
# forking this script. The accessory arm builds R plus the novel accessory
# contigs into ${BUILD}/refs_acc and sets REFSDIR to it; everything downstream
# -- the asset check, the caller, the output naming -- is unchanged, which is
# the point: a variant inside an accessory contig should come out of the same
# caller by the same path as one on the chromosome.
REFSDIR="${REFSDIR:-${BUILD}/refs}"
REF="${REFSDIR}/${REFID}.fasta"
for f in "$REF" "${REF}.bwt" "${REF}.fai" "${REF%.*}.dict"; do
    [[ -s "$f" ]] || { echo "FATAL: ${SAMPLE}: P0 asset missing: $f" >&2; exit 1; }
done

RELPATH="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $2; exit}' "$CRAMMAP")"
# An empty RELPATH used to slip through: CRAM became the collection ROOT, the
# `-r` test passed because a directory is readable, and the run died inside
# samtools with "Is a directory" instead of naming the real problem, which is
# that the sample is not in $CRAMMAP. p1_select_reference.sh checks this; this
# script did not.
[[ -n "$RELPATH" ]] || { echo "FATAL: ${SAMPLE}: no CRAM path in ${CRAMMAP}" >&2; exit 1; }
CRAM="${CRAMROOT}/${RELPATH}"
[[ -f "$CRAM" && -r "$CRAM" ]] \
    || { echo "FATAL: ${SAMPLE}: CRAM is not a readable file: ${CRAM}" >&2; exit 1; }

# DONE MEANS THE MARKER, written last. Testing only for vcf.gz and delly.vcf
# accepted a sample whose job died during dysgu (no dysgu VCF, delly never
# stamped). Samples finished before the marker existed are recognised by their
# complete outputs -- including a dysgu VCF with a #CHROM line -- and marked.
#
# AND DONE MEANS DONE AGAINST THIS BUILD AND THIS REFERENCE (audit P0P2-1).
# The marker held only a date, so after a rebuild, or a refmap that changed
# its choice, P2 kept calls made against the previous reference. The marker
# now records both, and outputs that name another build or reference are
# refused: rerun into new output folders, or move them aside.
DONE="${OUTDIR}/${SAMPLE}.p2.done"
_write_done() {
    printf 'build_id\t%s\nreference\t%s\ncompleted\t%s\n' \
        "$BUILD_ID" "$REFID" "$(date -Is)" > "${DONE}.tmp"
    mv -f "${DONE}.tmp" "$DONE"
}
# What the small-variant VCF says it was made from: its build stamp, and the
# reference in HaplotypeCaller's command line.
_vcf_ref() {
    gzip -cd "$1" 2>/dev/null | awk '/^#CHROM/{exit}
        /^##GATKCommandLine=<ID=HaplotypeCaller/{
            n=split($0,w," ")
            for(i=1;i<n;i++) if(w[i]=="--reference"||w[i]=="-R"){print w[i+1]; exit}}' || true
}
_m_build="$(mtb_kv "$DONE" build_id)"; _m_ref="$(mtb_kv "$DONE" reference)"
_v="${OUTDIR}/${SAMPLE}.vcf.gz"
_v_build="$(mtb_vcf_build_id "$_v")"
_v_ref=""; [[ -s "$_v" ]] && _v_ref="$(_vcf_ref "$_v")"
_foreign=""
if [[ -n "$_m_build" ]]; then
    [[ "$_m_build" == "$BUILD_ID" && "$_m_ref" == "$REFID" ]] \
        || _foreign="${DONE} records build ${_m_build}, reference ${_m_ref:-?}"
elif [[ -s "$_v" ]]; then
    [[ "$_v_build" == "$BUILD_ID" ]] \
        || _foreign="${_v} is stamped '${_v_build:-unstamped}'"
    [[ -z "$_v_ref" || "$(basename "$_v_ref")" == "${REFID}.fasta" ]] \
        || _foreign="${_foreign:+${_foreign}; }${_v} was called against ${_v_ref}"
fi
if [[ -n "$_foreign" ]]; then
    echo "FATAL: ${SAMPLE}: P2 output from another build or reference: ${_foreign};" \
         "this is build ${BUILD_ID}, reference ${REFID}. Rerun into new output" \
         "folders (OUTDIR, WORK) or move the old outputs aside." >&2
    exit 1
fi
if [[ -n "$_m_build" ]]; then
    echo "[P2] ${SAMPLE}: already done (build ${BUILD_ID}, reference ${REFID})"; exit 0
fi
# Outputs from before the marker recorded build and reference (a date-only
# marker, or none): accepted only when the VCF's stamp and caller reference,
# checked above, are this build's and this reference.
_legacy_complete() {
    [[ -s "${OUTDIR}/${SAMPLE}.vcf.gz" && -s "${OUTDIR}/${SAMPLE}.delly.vcf" ]] || return 1
    [[ -n "$_v_ref" ]] || return 1
    [[ -x "$DYSGU" ]] || return 0                     # dysgu not part of this site
    [[ -s "${OUTDIR}/${SAMPLE}.dysgu.vcf" ]] || return 1
    grep -q '^#CHROM' "${OUTDIR}/${SAMPLE}.dysgu.vcf"
}
if _legacy_complete; then
    _write_done
    echo "[P2] ${SAMPLE}: already done (outputs of build ${BUILD_ID}, reference" \
         "${REFID}, from before the marker recorded them)"; exit 0
fi
echo "[P2] ${SAMPLE}: reference ${REFID}, build ${BUILD_ID}"

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
echo "[P2] ${SAMPLE}: ${NR1} read pairs"

# --- small variants, in R coordinates ----------------------------------------
SIM_REQUIRE_FQ=1 SIM_FQ_PREFIX="${WORK}/${SAMPLE}" SIM_KEEP_FQ=1 SIM_GVCF=1 MTB_THREADS="$THREADS" \
    bash bin/simulate_and_call.sh "$REF" "$REF" "$WORK" "$SAMPLE" 1 150
mv -f "${WORK}/${SAMPLE}.vcf.gz"     "${OUTDIR}/${SAMPLE}.vcf.gz"
mv -f "${WORK}/${SAMPLE}.vcf.gz.tbi" "${OUTDIR}/${SAMPLE}.vcf.gz.tbi" 2>/dev/null || true
[[ -s "${WORK}/${SAMPLE}.g.vcf.gz" ]] && {
    mv -f "${WORK}/${SAMPLE}.g.vcf.gz" "${OUTDIR}/${SAMPLE}.g.vcf.gz"
    mv -f "${WORK}/${SAMPLE}.g.vcf.gz.tbi" "${OUTDIR}/${SAMPLE}.g.vcf.gz.tbi" 2>/dev/null || true
}

BAM="${WORK}/${SAMPLE}.bam"
[[ -s "$BAM" ]] || { echo "FATAL: ${SAMPLE}: no BAM after alignment" >&2; exit 1; }

# --- structural variants, same BAM, same coordinates -------------------------
export LD_LIBRARY_PATH="${DELLY_ENV}/lib:${LD_LIBRARY_PATH:-}"
"${DELLY_ENV}/bin/delly" call -g "$REF" -o "${WORK}/${SAMPLE}.delly.bcf" "$BAM" \
    > "${WORK}/${SAMPLE}.delly.log" 2>&1 || echo "[P2] ${SAMPLE}: delly nonzero exit" >&2
if [[ -s "${WORK}/${SAMPLE}.delly.bcf" ]]; then
    "$MTB_BCFTOOLS" view "${WORK}/${SAMPLE}.delly.bcf" > "${OUTDIR}/${SAMPLE}.delly.vcf.tmp"
    mv -f "${OUTDIR}/${SAMPLE}.delly.vcf.tmp" "${OUTDIR}/${SAMPLE}.delly.vcf"
else
    echo "FATAL: ${SAMPLE}: delly produced no BCF" >&2; exit 1
fi

if [[ -x "$DYSGU" ]]; then
    # A nonzero exit used to be logged and the partial stdout kept and stamped
    # as a complete callset. Keep the output only on success; otherwise fail
    # the sample so the gap is visible (once in the project's history so far).
    _dy_rc=0
    "$DYSGU" run --clean -x -p "$THREADS" "$REF" \
        "${WORK}/${SAMPLE}.dysgu_tmp" "$BAM" > "${OUTDIR}/${SAMPLE}.dysgu.vcf.tmp" \
        2> "${WORK}/${SAMPLE}.dysgu.log" || _dy_rc=$?
    rm -rf "${WORK}/${SAMPLE}.dysgu_tmp"
    if [[ "$_dy_rc" -ne 0 ]]; then
        rm -f "${OUTDIR}/${SAMPLE}.dysgu.vcf.tmp"
        echo "FATAL: ${SAMPLE}: dysgu exited ${_dy_rc}; see ${WORK}/${SAMPLE}.dysgu.log" >&2
        exit 1
    fi
    mv -f "${OUTDIR}/${SAMPLE}.dysgu.vcf.tmp" "${OUTDIR}/${SAMPLE}.dysgu.vcf"
else
    echo "[P2] ${SAMPLE}: dysgu not installed at ${DYSGU}" >&2
fi

# Stamp the SV VCFs too. simulate_and_call.sh already stamped the small-variant
# VCF and GVCF, but delly and dysgu write their own headers.
_sv_vcfs=("${OUTDIR}/${SAMPLE}.delly.vcf")
[[ -s "${OUTDIR}/${SAMPLE}.dysgu.vcf" ]] && _sv_vcfs+=("${OUTDIR}/${SAMPLE}.dysgu.vcf")
bash bin/stamp_build_id.sh "${_sv_vcfs[@]}"

rm -f "$FQ1" "$FQ2"
echo "[P2] ${SAMPLE}: small $(zcat "${OUTDIR}/${SAMPLE}.vcf.gz" | grep -vc '^#'), "\
"delly $(grep -vc '^#' "${OUTDIR}/${SAMPLE}.delly.vcf"), "\
"dysgu $([[ -s "${OUTDIR}/${SAMPLE}.dysgu.vcf" ]] && grep -vc '^#' "${OUTDIR}/${SAMPLE}.dysgu.vcf" || echo 0)"
_write_done
echo "[P2] ${SAMPLE}: done"
