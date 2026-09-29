#!/usr/bin/env bash
#SBATCH --job-name=P5_merge
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 0-06:00
# Six hours, not one. The --states step's cost tracks the KEY count, not the
# isolate's own data, because each task walks the whole key set to emit a dense
# row per key. At 100 isolates the worst task took 41.7 min of a one-hour
# request, which is not a margin to carry into a larger cohort; the failure mode
# is a wall-clock kill part way through a per-sample state file, which leaves a
# truncated file that looks like output. Memory is not the constraint: peak RSS
# was 0.65 GB against the 8 GB requested, and keys grow sublinearly with cohort
# size (40,445 at 23 isolates, 61,050 at 100).
#SBATCH -p shared
#SBATCH --mem=32000
# 32 GB, not 8. gwas1000's matrix step died OUT_OF_MEMORY at the 8 GB
# ceiling 98 seconds in: p5_matrix.py held every sample's whole state
# table as a dict of strings, which at 997 isolates and 173,685 keys is
# 173 million Python dict entries. That is fixed -- the matrix is now one
# uint8 array of 173 MB -- but the keys table and the per-row write still
# want room, and a job that dies here costs the whole tail of the chain.
#SBATCH --output=slurm/P5_%A_%a.out
#SBATCH --error=slurm/P5_%A_%a.err
#
# P5: merge the cohort into one matrix in the common frame.
#
#   bash bin/p5_merge.sh --keys                        # step 1, once
#   sbatch --array=1-23 bin/p5_merge.sh --states       # step 2, per sample
#   bash bin/p5_merge.sh --matrix                      # step 3, once
#
# This is where section 2.2.2's two open questions get answered against real
# output rather than in the abstract:
#
#   canonical REF   H37Rv's own base at the position, so a row is identical
#                   whichever reference a sample was called against. Measured:
#                   only 174 of 38,764 keys (0.45%) had samples disagreeing on
#                   REF, so the allele-normalisation problem is real but small.
#
#   ABSENT vs       four states, not three: ALT, REF, ABSENT, NOCALL. ABSENT
#   no-call         means the key does not project into that sample's own
#                   reference, so there is no sequence there to call; NOCALL
#                   means it projects but nothing was observed.
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
    [[ "${#_c[@]}" -eq 1 ]] || { echo "FATAL: set MTB_BUILD_DIR" >&2; exit 1; }
    BUILD="${_c[0]}"
fi
export MTB_BUILD_DIR="$BUILD"

REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
P4DIR="${P4DIR:-refbias/p4}"
P2DIR="${P2DIR:-refbias/p2}"
OUTDIR="${OUTDIR:-refbias/p5}"
WORK="${WORK:-refbias/work/p5}"
OG="${OG:-$(ls graphs/CX333.s10k.k23.K15/*.smooth.final.og 2>/dev/null | head -1)}"
ODGI="${MTB_ODGI:?MTB_ODGI is unset; see config/project_env.sh}"
H37RV_PATH="${H37RV_PATH:-GCF_000195955#1#NC_000962.3}"
H37RV_FA="${BUILD}/refs/GCF_000195955.fasta"
PATHS="${BUILD}/assets/paths.txt"
KEYS="${OUTDIR}/keys.tsv"
NT="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "$OUTDIR" "$WORK" slurm

# The step is normally $1. bin/refbias_run.sh submits through sbatch, which
# can pass environment but not positional arguments, so P5STEP is accepted as
# an equivalent. $1 wins if both are given.
STEP="${1:-${P5STEP:-}}"
case "$STEP" in
  --keys)
    # Write to a temporary file and move it into place ONLY if the key set
    # actually changed. p5_finish.sh's freshness guard compares the state
    # files' mtimes against keys.tsv, on the reasoning that a state file older
    # than the keys was computed against a different key set. That reasoning is
    # right, but it is defeated by an unconditional rewrite: re-running a
    # complete cohort regenerates a byte-identical keys.tsv with a new mtime
    # while every --states task skips with "already done" and keeps its old
    # one, so all 23 states read as stale and the matrix step refuses a run
    # whose inputs are in fact correct. That happened on the pilot re-run of
    # 2026-09-23. Preserving the mtime when the content is unchanged makes the
    # guard mean what it says: keys.tsv is newer only when the keys differ.
    "$MTB_PY" bin/p5_keys.py --refmap "$REFMAP" --dir "$P4DIR" \
        --h37rv "$H37RV_FA" --out "${KEYS}.tmp"
    if [[ -s "$KEYS" ]] && cmp -s "${KEYS}.tmp" "$KEYS"; then
        echo "[P5] keys unchanged (${KEYS}); keeping the existing file and its mtime"
        rm -f "${KEYS}.tmp"
    else
        mv -f "${KEYS}.tmp" "$KEYS"
        echo "[P5] keys written: $(wc -l < "$KEYS") lines"
    fi
    exit 0
    ;;
  --matrix)
    exec "$MTB_PY" bin/p5_matrix.py --refmap "$REFMAP" --keys "$KEYS" \
        --dir "$OUTDIR" --out "${OUTDIR}/matrix.tsv"
    ;;
  --states) ;;
  *) echo "usage: $0 --keys | --states | --matrix   (or P5STEP=...)" >&2; exit 2 ;;
esac

[[ -s "$KEYS" ]] || { echo "FATAL: run --keys first" >&2; exit 1; }
if [[ $# -ge 2 ]]; then
    SAMPLE="$2"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "FATAL: no SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    SAMPLE="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$REFMAP")"
fi
[[ -n "${SAMPLE:-}" ]] || { echo "FATAL: no refmap row" >&2; exit 1; }
REFID="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $5; exit}' "$REFMAP")"
RPATH="$(grep -m1 "^${REFID}#" "$PATHS" || true)"
[[ -n "$RPATH" ]] || { echo "FATAL: ${SAMPLE}: ${REFID} not a graph path" >&2; exit 1; }
GVCF="${P2DIR}/${SAMPLE}.g.vcf.gz"
PLACED="${P4DIR}/${SAMPLE}.placed.tsv"
for f in "$GVCF" "$PLACED"; do
    [[ -s "$f" ]] || { echo "FATAL: ${SAMPLE}: missing ${f}" >&2; exit 1; }
done
# "Already done" means done AGAINST THIS KEY SET. The skip used to test only
# that a states file existed, so after the keys changed every task reported
# done and p5_finish.sh then refused the stale files -- to be deleted by hand.
# The sparse header carries the checksum of the key list it was written
# against; compare it with keys.tsv and redo the sample on a mismatch.
_st="${OUTDIR}/${SAMPLE}.states.tsv"
if [[ -s "$_st" ]]; then
    _have="$(awk -F'\t' '$1=="#keys_sha1"{print $2; exit} !/^#/{exit}' "$_st")"
    _want="$("$MTB_PY" -c 'import csv,sys; sys.path.insert(0,"bin"); import p5_states_io as io; print(io.keys_sha1(list(csv.DictReader(open(sys.argv[1]),delimiter="\t"))))' "$KEYS")"
    if [[ -n "$_have" && "$_have" == "$_want" ]]; then
        echo "[P5] ${SAMPLE}: already done against this key set"; exit 0
    fi
    echo "[P5] ${SAMPLE}: states file is from a different key set" \
         "(${_have:-no checksum} vs ${_want}); recomputing"
fi


# WHICH FRAME IS THIS SAMPLE'S INPUT IN? Measured, not assumed. P1 and P2 align
# to ${BUILD}/refs, so their calls are refs-frame and the conversions below are
# required; the stage 2 simulation aligned to the rotated panel sequence and its
# calls are panel-frame, where the same conversions would be the defect rather
# than the fix. The two are indistinguishable from a VCF header -- same contig,
# same length -- so frame_detect.py decides it from the REF alleles and this
# stops rather than guesses. graphframe/docs/GRAPH_FRAME_RESOLUTION.md section 6.
VFRAME="$("$MTB_PY" graphframe/bin/frame_detect.py --vcf "$GVCF" \
    --accession "$REFID" --refs "${BUILD}/refs" | cut -f2)"
case "$VFRAME" in
    refs|either) ;;
    *) echo "FATAL: ${SAMPLE}: $GVCF is in the '${VFRAME}' frame, not refs; the
   conversion in this script would be wrong for it" >&2; exit 1 ;;
esac
# PROJECT ONCE PER REFERENCE, NOT ONCE PER SAMPLE. Every H37Rv-framed key is
# projected into the sample's reference so REF can be told from ABSENT -- a key
# that does not project has no sequence in that reference to call. But the
# projection is a function of (key set, reference) and of nothing else, and
# isolates share references heavily: gwas1000 is 997 isolates over 150
# references, so 847 of its 997 odgi runs recomputed an answer another task had
# already produced, and odgi is the dominant cost of this pass. scale200 is
# 2.2x redundant, gwas1000 6.6x, and a 10,000-isolate cohort over the
# 333-genome panel would be about 30x.
#
# The cache is keyed on the reference AND on a checksum of the key column, so a
# changed key set cannot be answered from a stale file -- the same discipline
# the sparse state files use. Written to a temporary name and renamed into
# place, because fifty tasks run at once and a half-written projection read as
# complete is the failure that cost this project a whole array.
# THE STORE PERSISTS BETWEEN RUNS AND IS KEYED BY POSITION. The earlier cache
# lived in the run's work directory and was keyed on a checksum of the key
# column, so gwas1000 got no benefit at all from the 150 references scale200
# had already projected -- a different cohort has a different key set and so a
# different checksum, although the positions overlap almost entirely. Keyed on
# (reference, position) and kept under the build, a new cohort pays only for
# positions nobody has projected yet. bin/proj_store.py says how it stays
# safe under fifty concurrent writers.
PROJSTORE="${PROJSTORE:-${BUILD}/proj}"
PROJDIR="${WORK}/proj"; mkdir -p "$PROJDIR" "$PROJSTORE"
KSHA="$(awk -F'\t' 'NR>1{print $1}' "$KEYS" | sha1sum | cut -c1-16)"
PROJ="${PROJDIR}/${REFID}.${KSHA}.pos"
NPOS="$(awk -F'\t' 'NR>1 && $2=="h37rv"' "$KEYS" | wc -l)"
if [[ -s "$PROJ" ]] && [[ "$(grep -vc '^#' "$PROJ")" -eq "$NPOS" ]]; then
    echo "[P5] ${SAMPLE}: reusing the ${REFID} projection ($NPOS positions)"
else
    PANEL="${PROJDIR}/anchors.${KSHA}.panel.txt"
    if [[ ! -s "$PANEL" ]]; then
        tr -d '\r' < "$KEYS" \
          | awk -F'\t' -v p="$H37RV_PATH" 'NR>1 && $2=="h37rv" {print p","($3-1)",+"}' \
          > "${PANEL}.raw.$$"
        [[ -s "${PANEL}.raw.$$" ]] \
            || { echo "FATAL: ${SAMPLE}: no H37Rv keys extracted" >&2; exit 1; }
        # The H37Rv anchors are in the refs frame and `odgi position` reads the
        # panel frame the graph was built in; the target coordinates come back
        # in the panel frame while the depth and gVCF files they index are in
        # the refs frame. For 110 of 333 accessions those frames differ by a
        # rotation and for 22 also by strand, so both conversions are
        # load-bearing. graphframe/docs/GRAPH_FRAME_RESOLUTION.md
        "$MTB_PY" graphframe/bin/frame_convert.py to-panel \
            < "${PANEL}.raw.$$" > "${PANEL}.$$"
        mv -f "${PANEL}.$$" "$PANEL"; rm -f "${PANEL}.raw.$$"
    fi
    # Ask the store what it already has, project only the remainder, fold the
    # answer back in, then emit the whole set from the store.
    MISS="${PROJDIR}/${REFID}.${KSHA}.miss.$$"
    "$MTB_PY" bin/proj_store.py missing --store "$PROJSTORE" --ref "$REFID" \
        --need "$PANEL" --out "$MISS"
    if [[ -s "$MISS" ]]; then
        "$ODGI" position -i "$OG" -F "$MISS" -r "$RPATH" -t "$NT" \
            2> "${WORK}/${SAMPLE}.odgi.log" \
          | "$MTB_PY" graphframe/bin/frame_convert.py from-panel \
            > "${MISS}.pos"
        # ODGI RETURNING NOTHING IS NOT AN ERROR HERE. Some H37Rv positions
        # have no equivalent in a given panel path, so `missing` keeps asking
        # for them, odgi keeps returning nothing, and there is nothing to add --
        # the deadlock bin/proj_store.py documents for GCF_965121955. It cost
        # 9 of gwas1000's 997 tasks on the 2026-09-28 run. The store already
        # held 151,603 of 151,941 positions for that reference, 99.78%, so the
        # right response is to proceed on what projects and let the floor below
        # decide, not to fail the sample.
        if [[ ! -s "${MISS}.pos" ]]; then
            echo "[P5] ${SAMPLE}: none of the $(grep -vc '^#' "$MISS" || true)" \
                 "outstanding positions project into ${REFID}; proceeding" >&2
            rm -f "${MISS}.pos"
        fi
        if [[ -s "${MISS}.pos" ]]; then
            "$MTB_PY" bin/proj_store.py add --store "$PROJSTORE" --ref "$REFID" \
                --result "${MISS}.pos"
            rm -f "${MISS}.pos"
        fi
    fi
    rm -f "$MISS"
    "$MTB_PY" bin/proj_store.py emit --store "$PROJSTORE" --ref "$REFID" \
        --need "$PANEL" --out "${PROJ}.$$" --allow-missing
    got="$(grep -vc '^#' "${PROJ}.$$" || true)"
    # A FLOOR, NOT EQUALITY. Demanding all NPOS positions makes an unprojectable
    # remainder fatal forever, which is the same mistake p5_svgt.sh already
    # corrected by requiring 95% of its probe set. p5_states.py reports how many
    # keys had no projection and states them NOCALL, so a shortfall is recorded
    # in the output rather than hidden -- but a LARGE shortfall means the store
    # is wrong, not the graph, and must still stop the run.
    FLOOR_PCT="${PROJ_FLOOR_PCT:-99}"
    need_min=$(( NPOS * FLOOR_PCT / 100 ))
    if [[ "$got" -lt "$need_min" ]]; then
        echo "FATAL: ${SAMPLE}: the store returned ${got} of ${NPOS} positions" \
             "for ${REFID}, below the ${FLOOR_PCT}% floor" >&2
        rm -f "${PROJ}.$$"; exit 1
    fi
    if [[ "$got" -ne "$NPOS" ]]; then
        echo "[P5] ${SAMPLE}: ${REFID} projection is ${got} of ${NPOS}" \
             "($(( got * 100 / NPOS ))%); the rest have no equivalent in that path"
    fi
    mv -f "${PROJ}.$$" "$PROJ"
    echo "[P5] ${SAMPLE}: ${REFID} projection assembled from the store"
fi

"$MTB_PY" bin/p5_states.py --sample "$SAMPLE" --reference "$REFID" \
    --keys "$KEYS" --placed "$PLACED" --gvcf "$GVCF" --h37rv "$H37RV_FA" \
    --projected "$PROJ" \
    --out "${OUTDIR}/${SAMPLE}.states.tsv"
# nothing per-sample to remove: the projection is a shared cache
echo "[P5] ${SAMPLE}: done"
