#!/usr/bin/env bash
#SBATCH --job-name=P4b_sv
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH -t 0-00:40
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH --output=slurm/P4b_%A_%a.out
#SBATCH --error=slurm/P4b_%A_%a.err
#
# P4b: place the SV calls P2 has been producing and nothing has been reading.
#
#   sbatch --array=1-23 bin/p4b_place_sv.sh
#   bash bin/p4b_place_sv.sh <sample>
#   bash bin/p4b_place_sv.sh --summary
#
# P2 emits delly and dysgu VCFs per sample in the matched reference's
# coordinates. P4 places only the small-variant VCF, so the matrix has contained
# zero SV records. This projects both breakpoints of every merged event through
# the graph and emits them in the same frame as everything else.
set -euo pipefail
_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
source "$_mtb_env"
cd "${SLURM_SUBMIT_DIR:-.}"

# Exactly one build, or MTB_BUILD_DIR. `find | head -1` picked whichever build
# the filesystem listed first when there were several, while p4_place.sh
# refused -- so two passes of one chain could run against different builds.
BUILD="${MTB_BUILD_DIR:-}"
if [[ -z "$BUILD" ]]; then
    mapfile -t _c < <(find "${BUILD_ROOT:-refbias/build}" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
    [[ "${#_c[@]}" -eq 1 ]] || { echo "FATAL: set MTB_BUILD_DIR (${#_c[@]} builds)" >&2; exit 1; }
    BUILD="${_c[0]}"
fi
export MTB_BUILD_DIR="$BUILD"
BUILD_ID="$(awk -F'\t' '$1=="build_id"{print $2}' "${BUILD}/build_info.tsv")"
REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
P2DIR="${P2DIR:-refbias/p2}"
OUTDIR="${OUTDIR:-refbias/p4b}"
WORK="${WORK:-refbias/work/p4b}"
OG="${OG:-$(ls graphs/CX333.s10k.k23.K15/*.smooth.final.og 2>/dev/null | head -1)}"
ODGI="${MTB_ODGI:?MTB_ODGI is unset; see config/project_env.sh}"
H37RV_PATH="${H37RV_PATH:-GCF_000195955#1#NC_000962.3}"
PATHS="${BUILD}/assets/paths.txt"
MASK="${BUILD}/assets/repeat_mask.bed"
NT="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "$OUTDIR" "$WORK" slurm

# The step is normally $1. bin/refbias_run.sh submits through sbatch,
# which passes environment but not positional arguments, so P4BSTEP is
# accepted as an equivalent. $1 wins if both are given.
if [[ "${1:-${P4BSTEP:-}}" == "--summary" ]]; then
    exec "$MTB_PY" bin/p4b_summary.py --refmap "$REFMAP" --dir "$OUTDIR" \
        --out "${OUTDIR}/p4b_summary.tsv"
fi

if [[ $# -ge 1 ]]; then
    SAMPLE="$1"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "FATAL: no task id" >&2; exit 1; }
    SAMPLE="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$REFMAP")"
fi
[[ -n "${SAMPLE:-}" ]] || { echo "FATAL: no refmap row" >&2; exit 1; }
REFID="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $5; exit}' "$REFMAP")"
RPATH="$(grep -m1 "^${REFID}#" "$PATHS" || true)"
[[ -n "$RPATH" ]] || { echo "FATAL: ${SAMPLE}: ${REFID} not a graph path" >&2; exit 1; }
DELLY="${P2DIR}/${SAMPLE}.delly.vcf"
DYSGU="${P2DIR}/${SAMPLE}.dysgu.vcf"
OUT="${OUTDIR}/${SAMPLE}.sv_placed.tsv"
[[ -s "$OUT" ]] && { echo "[P4b] ${SAMPLE}: already done"; exit 0; }

# WHICH FRAME ARE THE CALLER VCFs IN? Measured, not assumed. delly and dysgu
# run on the P2 alignment to ${BUILD}/refs, so their breakpoints are refs-frame,
# while `odgi position` reads the panel frame the graph was built in. For 110 of
# the 333 accessions those differ by a rotation and for 22 also by strand.
# graphframe/docs/GRAPH_FRAME_RESOLUTION.md
#
# A header-only caller VCF passes [[ -s ]] and has nothing to test; frame_detect
# then exits 2 and a bare $(...) under `set -e` killed the task. n=0 is "nothing
# to test", not a mismatch, so it is let through.
for _v in "$DELLY" "$DYSGU"; do
    [[ -s "$_v" ]] || continue
    _fd_rc=0
    _fd="$("$MTB_PY" graphframe/bin/frame_detect.py --vcf "$_v" \
        --accession "$REFID" --refs "${BUILD}/refs")" || _fd_rc=$?
    VFRAME="$(cut -f2 <<< "$_fd")"
    _fd_n="$(grep -o 'n=[0-9]*' <<< "$_fd" | cut -d= -f2)"
    case "$VFRAME" in
        refs|either) ;;
        *) [[ "${_fd_n:-}" == "0" ]] && continue
           echo "FATAL: ${SAMPLE}: ${_v} is in the '${VFRAME}' frame, not refs" \
                "(frame_detect exit ${_fd_rc}: ${_fd})" >&2
           exit 1 ;;
    esac
done

# every breakpoint of every caller record, projected in one odgi call
{ for f in "$DELLY" "$DYSGU"; do
    [[ -s "$f" ]] || continue
    awk -v p="$RPATH" '!/^#/ {
        pos=$2; end=pos;
        if (match($8, /(^|;)END=[0-9]+/)) { e=substr($8, RSTART, RLENGTH); sub(/.*END=/, "", e); end=e }
        print p","(pos-1)",+"; print p","(end-1)",+"
    }' "$f"
  done; } | sort -u \
  | "$MTB_PY" graphframe/bin/frame_convert.py to-panel \
  > "${WORK}/${SAMPLE}.bp.txt"
NBP=$(wc -l < "${WORK}/${SAMPLE}.bp.txt")
if [[ "$NBP" -eq 0 ]]; then
    echo "[P4b] ${SAMPLE}: no SV breakpoints to project"
    : > "${WORK}/${SAMPLE}.bp.pos"
else
    "$ODGI" position -i "$OG" -F "${WORK}/${SAMPLE}.bp.txt" -r "$H37RV_PATH" \
        -t "$NT" 2> "${WORK}/${SAMPLE}.odgi.log" \
      | "$MTB_PY" graphframe/bin/frame_convert.py from-panel \
      > "${WORK}/${SAMPLE}.bp.pos"
    [[ -s "${WORK}/${SAMPLE}.bp.pos" ]] \
        || { echo "FATAL: ${SAMPLE}: odgi produced nothing for ${NBP} breakpoints" >&2; exit 1; }
fi
echo "[P4b] ${SAMPLE}: ${NBP} distinct breakpoints projected"

GRAPH_VCF="${GRAPH_VCF:-$(dirname "$OG")/all_variants.nolab.vcf.gz}"
[[ -s "$GRAPH_VCF" ]] || { echo "FATAL: no graph VCF at ${GRAPH_VCF}" >&2; exit 1; }
"$MTB_PY" bin/p4b_place_sv.py --sample "$SAMPLE" --reference "$REFID" \
    --build-id "$BUILD_ID" --delly "$DELLY" --dysgu "$DYSGU" \
    --graph-vcf "$GRAPH_VCF" \
    --positions "${WORK}/${SAMPLE}.bp.pos" --mask "$MASK" --out "$OUT"
rm -f "${WORK}/${SAMPLE}.bp.txt"
echo "[P4b] ${SAMPLE}: done"
