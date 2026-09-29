#!/usr/bin/env bash
#SBATCH --job-name=P4_place
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 0-01:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH --output=slurm/P4_%A_%a.out
#SBATCH --error=slurm/P4_%A_%a.err
#
# P4: place each isolate's calls into the common frame, routing by region.
#
#   sbatch --array=1-23 bin/p4_place.sh
#   bash bin/p4_place.sh <sample>
#   bash bin/p4_place.sh --summary
#
# The routing rule is not a choice made here; it is what the tests measured.
#
#   region              arm chosen        core PPV / PE-PPE PPV
#   core sequence       DIRECT (H37Rv)    0.9546   vs composed 0.9135
#   PE/PPE and masked   MATCHED (R)       0.8664   vs direct    0.7217
#
# T15's hybrid on those two rules scored PPV 0.9366 with 34% fewer false
# positives than either arm alone, and nothing found since displaces it. T16
# tried to remove the need for it by adopting the graph frame and failed:
# composition is uniformly ~13% wrong on inherited differences and no region rule
# separates those errors, so routing by repeat mask -- which does separate them --
# remains the answer.
#
# Both arms already exist, which is the third use pass one earns. P1's H37Rv VCF
# is the direct arm; P2's R VCF is the matched arm. Nothing is re-aligned.
#
# Off-path records are keyed on the graph NODE, not on an H37Rv anchor and not on
# an accessory contig. Two samples matched to different references that carry the
# same accessory locus land on the same node, which is exact; the H37Rv anchor is
# shared only to within tens of bp, and the accessory catalogue has no locus
# within 50 bp for 29% of accessory-scale calls.
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
    [[ "${#_c[@]}" -eq 1 ]] || { echo "FATAL: set MTB_BUILD_DIR (${#_c[@]} builds)" >&2; exit 1; }
    BUILD="${_c[0]}"
fi
export MTB_BUILD_DIR="$BUILD"
BUILD_ID="$(awk -F'\t' '$1=="build_id"{print $2}' "${BUILD}/build_info.tsv")"

REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
P1WORK="${P1WORK:-refbias/work/p1}"
P2DIR="${P2DIR:-refbias/p2}"
OUTDIR="${OUTDIR:-refbias/p4}"
WORK="${WORK:-refbias/work/p4}"
OG="${OG:-$(ls graphs/CX333.s10k.k23.K15/*.smooth.final.og 2>/dev/null | head -1)}"
ODGI="${MTB_ODGI:-/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/odgi/bin/odgi}"
H37RV_PATH="${H37RV_PATH:-GCF_000195955#1#NC_000962.3}"
PATHS="${BUILD}/assets/paths.txt"
MASK="${BUILD}/assets/repeat_mask.bed"
LOCI="${BUILD}/assets/accessory_loci.tsv"
NT="${SLURM_CPUS_PER_TASK:-4}"

for f in "$REFMAP" "$OG" "$ODGI" "$PATHS" "$MASK" "$LOCI"; do
    [[ -e "$f" ]] || { echo "FATAL: missing: $f" >&2; exit 1; }
done
mkdir -p "$OUTDIR" "$WORK" slurm

# The step is normally $1. bin/refbias_run.sh submits through sbatch,
# which passes environment but not positional arguments, so P4STEP is
# accepted as an equivalent. $1 wins if both are given.
if [[ "${1:-${P4STEP:-}}" == "--summary" ]]; then
    exec "$MTB_PY" bin/p4_summary.py --refmap "$REFMAP" --dir "$OUTDIR" \
        --out "${OUTDIR}/p4_summary.tsv"
fi

if [[ $# -ge 1 ]]; then
    SAMPLE="$1"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "FATAL: no SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    SAMPLE="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$REFMAP")"
fi
[[ -n "${SAMPLE:-}" ]] || { echo "FATAL: no refmap row" >&2; exit 1; }
REFID="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $5; exit}' "$REFMAP")"
[[ -n "$REFID" ]] || { echo "FATAL: ${SAMPLE}: no reference" >&2; exit 1; }
RPATH="$(grep -m1 "^${REFID}#" "$PATHS" || true)"
[[ -n "$RPATH" ]] || { echo "FATAL: ${SAMPLE}: ${REFID} is not a graph path" >&2; exit 1; }
DIRECT="${P1WORK}/${SAMPLE}.h37rv.vcf.gz"
MATCHED="${P2DIR}/${SAMPLE}.vcf.gz"
for f in "$DIRECT" "$MATCHED"; do
    [[ -s "$f" ]] || { echo "FATAL: ${SAMPLE}: missing input ${f}" >&2; exit 1; }
done

if [[ -s "${OUTDIR}/${SAMPLE}.placed.tsv" ]]; then
    echo "[P4] ${SAMPLE}: already done"; exit 0
fi
echo "[P4] ${SAMPLE}: reference ${REFID}, build ${BUILD_ID}"


# WHICH FRAME IS THIS SAMPLE'S INPUT IN? Measured, not assumed. P1 and P2 align
# to ${BUILD}/refs, so their calls are refs-frame and the conversions below are
# required; the stage 2 simulation aligned to the rotated panel sequence and its
# calls are panel-frame, where the same conversions would be the defect rather
# than the fix. The two are indistinguishable from a VCF header -- same contig,
# same length -- so frame_detect.py decides it from the REF alleles and this
# stops rather than guesses. graphframe/docs/GRAPH_FRAME_RESOLUTION.md section 6.
VFRAME="$("$MTB_PY" graphframe/bin/frame_detect.py --vcf "$MATCHED" \
    --accession "$REFID" --refs "${BUILD}/refs" | cut -f2)"
case "$VFRAME" in
    refs|either) ;;
    *) echo "FATAL: ${SAMPLE}: $MATCHED is in the '${VFRAME}' frame, not refs; the
   conversion in this script would be wrong for it" >&2; exit 1 ;;
esac
# R-coordinate positions of the matched arm's calls, once. These come out of a
# VCF called against refbias/build/<id>/refs/<R>.fasta, so they are in the REFS
# frame, while `odgi position` reads the PANEL frame the graph was built in. For
# 110 of the 333 accessions those differ by a rotation, and for 22 also by
# strand, so the two conversions below are not optional bookkeeping: without
# them the projection lands the offset distance away and scores 29.7% against
# H37Rv's own base instead of 100%. See graphframe/docs/GRAPH_FRAME_RESOLUTION.md.
zcat "$MATCHED" | awk -v p="$RPATH" '!/^#/ {print p","($2-1)",+"}' \
    > "${WORK}/${SAMPLE}.rpos.txt"
NPOS=$(wc -l < "${WORK}/${SAMPLE}.rpos.txt")
# AN EMPTY MATCHED VCF IS A LEGITIMATE STATE, NOT A FAILURE. It means the
# matched reference is so close to the isolate that there is nothing to call
# against it -- the limiting case, seen for real when a truth genome that is
# itself a panel path selects its own assembly and P2 emits 0 records with a
# 12,954-block gVCF proving the alignment worked. Exiting here lost the whole
# sample, direct arm included, over the absence of something optional. The
# composed arm is simply empty and the odgi passes are skipped.
if [[ "$NPOS" -eq 0 ]]; then
    echo "[P4] ${SAMPLE}: matched VCF has no records -- the reference is that "
    echo "     close; the composed CALLED arm will be empty and only the "
    echo "     direct arm and the inherited half contribute"
    : > "${WORK}/${SAMPLE}.h37rv.pos"
    : > "${WORK}/${SAMPLE}.node.pos"
fi
if [[ "$NPOS" -gt 0 ]]; then
"$MTB_PY" graphframe/bin/frame_convert.py to-panel \
    < "${WORK}/${SAMPLE}.rpos.txt" > "${WORK}/${SAMPLE}.rpos.panel.txt"

# two odgi passes: H37Rv projection (with dist.to.ref), and graph node positions
"$ODGI" position -i "$OG" -F "${WORK}/${SAMPLE}.rpos.panel.txt" -r "$H37RV_PATH" -t "$NT" \
    2> "${WORK}/${SAMPLE}.odgi.ref.log" \
  | "$MTB_PY" graphframe/bin/frame_convert.py from-panel \
    > "${WORK}/${SAMPLE}.h37rv.pos"
"$ODGI" position -i "$OG" -F "${WORK}/${SAMPLE}.rpos.panel.txt" -v -t "$NT" \
    2> "${WORK}/${SAMPLE}.odgi.node.log" \
  | "$MTB_PY" graphframe/bin/frame_convert.py from-panel \
    > "${WORK}/${SAMPLE}.node.pos"
for f in "${WORK}/${SAMPLE}.h37rv.pos" "${WORK}/${SAMPLE}.node.pos"; do
    [[ -s "$f" ]] || { echo "FATAL: ${SAMPLE}: odgi produced nothing in ${f}" >&2; exit 1; }
done
fi

"$MTB_PY" bin/p4_place.py \
    --sample "$SAMPLE" --reference "$REFID" --build-id "$BUILD_ID" \
    --h37rv "${BUILD}/refs/GCF_000195955.fasta" \
    --ref-fasta "${BUILD}/refs/${REFID}.fasta" \
    --direct "$DIRECT" --matched "$MATCHED" \
    --h37rv-pos "${WORK}/${SAMPLE}.h37rv.pos" \
    --node-pos "${WORK}/${SAMPLE}.node.pos" \
    --mask "$MASK" --loci "$LOCI" \
    --out "${OUTDIR}/${SAMPLE}.placed.tsv"
rm -f "${WORK}/${SAMPLE}.rpos.txt" "${WORK}/${SAMPLE}.rpos.panel.txt"
echo "[P4] ${SAMPLE}: done"
