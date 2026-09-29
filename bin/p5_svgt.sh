#!/usr/bin/env bash
#SBATCH --job-name=P5_svgt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 0-04:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH --output=slurm/P5svgt_%A_%a.out
#SBATCH --error=slurm/P5svgt_%A_%a.err
#
# P5 step 4: genotype deletion ABSENCE, so the SV block is not presence-only.
#
#   bash   bin/p5_svgt.sh --probes                  # step 1, once
#   sbatch --array=1-200 bin/p5_svgt.sh --states    # step 2, per sample
#   bash   bin/p5_svgt.sh --merge                   # step 3, once
#
# WHY THIS PASS EXISTS. bin/p5_sv_matrix.py writes ALT, NOCALL and ABSENT and
# deliberately never REF, because delly's measured sensitivity is 0.29 for
# 50-499 bp deletions, so silence is mostly a missed call rather than an
# absence. The cost of that correctness showed up the first time the merged VCF
# met a phylogenetic method: over scale200 the SV block is 1.80% ALT and 98.20%
# NOCALL with no REF at all, and with no leaf ever in the reference state a
# parsimony reconstruction has nothing to contrast, so every one of its 6,415
# SV records scored zero independent gains. The fix is not to reinterpret
# silence -- it is to go and measure the absence, which for a deletion is
# directly observable as read depth across the interval it would remove.
#
# The projection block below is the same one p5_merge.sh uses for key
# positions, and for the same reason: the gVCF is in the sample's own refs
# frame while the SV matrix is in H37Rv's, and for 110 of 333 accessions those
# frames differ by a rotation and for 22 also by strand.
#
# Insertions are NOT genotyped here and stay NOCALL. Depth across an insertion
# point is the same whether or not the sample carries the insertion; testing
# that needs spanning reads from the BAM, which is a different instrument.
set -euo pipefail

_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: no project_env.sh" >&2; exit 1; }
source "$_mtb_env"
cd "${SLURM_SUBMIT_DIR:-.}"

BUILD_ROOT="${BUILD_ROOT:-refbias/build}"
BUILD="${MTB_BUILD_DIR:-}"
if [[ -z "$BUILD" ]]; then
    mapfile -t _c < <(find "$BUILD_ROOT" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
    [[ "${#_c[@]}" -eq 1 ]] || { echo "FATAL: set MTB_BUILD_DIR" >&2; exit 1; }
    BUILD="${_c[0]}"
fi
export MTB_BUILD_DIR="$BUILD"

REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
P2DIR="${P2DIR:-refbias/p2}"
OUTDIR="${OUTDIR:-refbias/p5}"
WORK="${WORK:-refbias/work/p5svgt}"
SVMATRIX="${SVMATRIX:-${OUTDIR}/sv_matrix.tsv}"
# INTERVAL MODE. With IVTAB set, the pass genotypes every sample at every
# interval of a catalogue built once by bin/sv_intervals.py, and the caller's
# per-sample breakpoints define nothing. That is the rebuilt arm; the
# sv_matrix path below is the original one and is kept so the two can be
# compared. assoc/SV_SCATTER_RESULTS.md is the evidence for the change.
IVTAB="${IVTAB:-}"
if [[ -n "$IVTAB" ]]; then
    SVDIR="${SVDIR:-${OUTDIR}/svgt_iv}"
    PROBES_NAME="probes.iv.txt"
    OUTSUF="svgt_iv"
else
    SVDIR="${SVDIR:-${OUTDIR}/svgt}"
    PROBES_NAME="probes.h37rv.txt"
    OUTSUF="svgt"
fi
OG="${OG:-$(ls graphs/CX333.s10k.k23.K15/*.smooth.final.og 2>/dev/null | head -1)}"
ODGI="${MTB_ODGI:?MTB_ODGI is unset; see config/project_env.sh}"
H37RV_PATH="${H37RV_PATH:-GCF_000195955#1#NC_000962.3}"
PATHS="${BUILD}/assets/paths.txt"
PROBES="${WORK}/${PROBES_NAME}"
NT="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "$OUTDIR" "$WORK" "$SVDIR" slurm

STEP="${1:-${SVGTSTEP:-}}"
case "$STEP" in
  --probes)
    if [[ -n "$IVTAB" ]]; then
        [[ -s "$IVTAB" ]] || { echo "FATAL: no ${IVTAB}" >&2; exit 1; }
        exec "$MTB_PY" bin/p5_sv_genotype.py --intervals "$IVTAB" \
            --h37rv-path "$H37RV_PATH" --probes-out "$PROBES"
    fi
    [[ -s "$SVMATRIX" ]] || { echo "FATAL: no ${SVMATRIX}" >&2; exit 1; }
    exec "$MTB_PY" bin/p5_sv_genotype.py --sv-matrix "$SVMATRIX" \
        --h37rv-path "$H37RV_PATH" --probes-out "$PROBES"
    ;;
  --merge)
    exec "$MTB_PY" - "$SVDIR" "${OUTDIR}/${OUTSUF}_states.tsv" <<'PY'
import csv, glob, os, sys, collections
d, out = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(os.path.join(d, "*.tsv")))
if not files:
    sys.exit(f"FATAL: no per-sample tables in {d}")
n = collections.Counter()
with open(out, "w", newline="") as fh:
    w = None
    for f in files:
        for r in csv.DictReader(open(f, newline=""), delimiter="\t"):
            if w is None:
                w = csv.DictWriter(fh, fieldnames=list(r), delimiter="\t")
                w.writeheader()
            w.writerow(r); n[r["state"]] += 1
t = sum(n.values())
print(f"  {len(files)} samples, {t:,} SV cells  " +
      "  ".join(f"{k} {v:,} ({v/t:.1%})" for k, v in sorted(n.items())))
print(f"  -> {out}")
PY
    ;;
  --states) ;;
  *) echo "usage: $0 --probes | --states | --merge   (or SVGTSTEP=...)" >&2; exit 2 ;;
esac

[[ -s "$PROBES" ]] || { echo "FATAL: run --probes first" >&2; exit 1; }
if [[ $# -ge 2 ]]; then
    SAMPLE="$2"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "FATAL: no SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    SAMPLE="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$REFMAP")"
fi
[[ -n "${SAMPLE:-}" ]] || { echo "FATAL: no refmap row" >&2; exit 1; }
OUT="${SVDIR}/${SAMPLE}.${OUTSUF}.tsv"
# Done only if written AFTER the current probe list. The probes are rebuilt
# from the SV matrix or interval catalogue by --probes, and SV keys shift when
# cluster membership changes, so an older per-sample table either fails to join
# or -- worse -- joins to the wrong cluster. The skip used to test existence.
if [[ -s "$OUT" && "$OUT" -nt "$PROBES" ]]; then
    echo "[P5svgt] ${SAMPLE}: already done against the current probes"; exit 0
fi
[[ -s "$OUT" ]] && echo "[P5svgt] ${SAMPLE}: output predates ${PROBES}; recomputing"

REFID="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $5; exit}' "$REFMAP")"
RPATH="$(grep -m1 "^${REFID}#" "$PATHS" || true)"
[[ -n "$RPATH" ]] || { echo "FATAL: ${SAMPLE}: ${REFID} not a graph path" >&2; exit 1; }
GVCF="${P2DIR}/${SAMPLE}.g.vcf.gz"
[[ -s "$GVCF" ]] || { echo "FATAL: ${SAMPLE}: missing ${GVCF}" >&2; exit 1; }
# The alignment, for mapping-quality-filtered depth in element-proximal
# intervals. P2WORK is where p2_call.sh leaves the BAMs.
BAMDIR="${BAMDIR:-${P2WORK:-refbias/work/p2}}"
BAM="${BAMDIR}/${SAMPLE}.bam"

# Same frame check p5_merge.sh makes, and for the same reason: a panel-frame
# gVCF would make the conversions below the defect rather than the fix.
VFRAME="$("$MTB_PY" graphframe/bin/frame_detect.py --vcf "$GVCF" \
    --accession "$REFID" --refs "${BUILD}/refs" | cut -f2)"
case "$VFRAME" in
    refs|either) ;;
    *) echo "FATAL: ${SAMPLE}: $GVCF is in the '${VFRAME}' frame" >&2; exit 1 ;;
esac

# REUSE A PROJECTION THAT IS STILL VALID, AND ONLY THEN. odgi is the whole cost
# of this pass and a task that dies after it leaves a usable projection behind,
# so reuse is worth having -- but the first version of this guard checked only
# that the file existed and was newer than the probe list, and a run cancelled
# mid-odgi leaves a file that is both. 111 truncated projections were reused,
# the shortest 19% complete, and every probe missing from one reads as "does
# not project", which this script's own logic then calls ABSENT. One sample
# came back 80% ABSENT and the whole array had to be discarded.
#
# So: the projection is written to .tmp and moved into place only once it is
# complete, and reuse additionally requires the row count to match the probe
# count, which catches anything left from before that discipline.
# ONE PROJECTION PER REFERENCE, NOT PER SAMPLE. odgi is the whole cost of this
# pass -- about an hour a task -- and the projection depends only on (probe
# list, reference), not on the isolate. gwas1000 is 997 isolates over 150
# references, so 847 of its tasks were recomputing an answer another task had
# already produced. Cached under proj/ and keyed on the reference and a
# checksum of the probe list.
NPROBE="$(grep -c . "$PROBES")"
# NOT an exact match. odgi legitimately drops a small number of positions that
# do not project into a given path -- 177,209 of 177,553 on gwas1000, 0.19% --
# and requiring equality failed 8 tasks outright. A truncated file from a
# killed task is nothing like that: the one this guard was written for was 19%
# complete. So the floor is a fraction, which separates the two.
MINFRAC="${MINFRAC:-95}"
NMIN=$(( NPROBE * MINFRAC / 100 ))

# THE STORE IS KEYED ON (REFERENCE, POSITION) AND LIVES UNDER THE BUILD, so it
# is cohort-independent by construction. The per-cohort cache this replaced was
# keyed on the reference plus a checksum of the probe LIST, under
# ${WORK}/proj -- which meant scale200 and gwas1000 could not share a single
# projection even where they select the same reference and use the same
# catalogue. Measured: 81 of gwas1000's 150 references were already projected
# for scale200 with the identical interval-mode probe set, and the old cache
# would have rebuilt all 150 at about 134 CPU-minutes each. 181 CPU-hours of
# duplicated odgi.
#
# bin/proj_store.py says how it stays consistent: append-only part files, and
# readers take the first value for each position, so two tasks projecting the
# same reference at once cannot corrupt it.
PROJSTORE="${PROJSTORE:-${BUILD}/proj}"
PROJDIR="${WORK}/proj"; mkdir -p "$PROJDIR" "$PROJSTORE"
POS="${PROJDIR}/${REFID}.$(sha1sum "$PROBES" | cut -c1-16).pos"

# Ask the store what it is missing, project only that, fold it back, then emit
# the whole probe set from the store.
MISS="${PROJDIR}/${REFID}.miss.$$"
"$MTB_PY" bin/proj_store.py missing --store "$PROJSTORE" --ref "$REFID" \
    --need "$PROBES" --out "$MISS"
if [[ -s "$MISS" ]]; then
    echo "[P5svgt] ${SAMPLE}: projecting $(grep -c . "$MISS") new positions for ${REFID}"
    "$MTB_PY" graphframe/bin/frame_convert.py to-panel \
        < "$MISS" > "${MISS}.panel"
    "$ODGI" position -i "$OG" -F "${MISS}.panel" -r "$RPATH" -t "$NT" \
        2> "${WORK}/${SAMPLE}.odgi.log" \
      | "$MTB_PY" graphframe/bin/frame_convert.py from-panel \
        > "${MISS}.pos"
    # AN EMPTY RESULT IS NOT A FAILURE HERE, and treating it as one cost 12
    # tasks. odgi legitimately drops positions that do not project into a given
    # path -- 0.19% of them -- and the store asks for exactly the positions it
    # does not yet have, so once every remaining one is undroppable the request
    # comes back with nothing. GCF_965121955 sat at 39,620 of 39,685 stored and
    # died on the last 65 every time. The real gate is the emit floor below,
    # which asks whether the store can now answer the whole probe set well
    # enough; that is the question that matters, and 39,620 of 39,685 is 99.8%.
    if [[ -s "${MISS}.pos" ]]; then
        "$MTB_PY" bin/proj_store.py add --store "$PROJSTORE" --ref "$REFID" \
            --result "${MISS}.pos"
    else
        echo "[P5svgt] ${SAMPLE}: none of those positions project into ${REFID};"\
             " the emit floor decides whether that matters"
    fi
    rm -f "${MISS}.panel" "${MISS}.pos"
else
    echo "[P5svgt] ${SAMPLE}: ${REFID} fully covered by the store"
fi
rm -f "$MISS"
"$MTB_PY" bin/proj_store.py emit --store "$PROJSTORE" --ref "$REFID" \
    --need "$PROBES" --allow-missing --out "${POS}.$$"
_got="$(grep -vc '^#' "${POS}.$$" || true)"
[[ "$_got" -ge "$NMIN" ]] || {
    echo "FATAL: ${SAMPLE}: the store emitted ${_got} of ${NPROBE} probes for "\
         "${REFID}, below the ${MINFRAC}% floor" >&2
    rm -f "${POS}.$$"; exit 1; }
mv -f "${POS}.$$" "$POS"
echo "[P5svgt] ${SAMPLE}: ${REFID} projection ready (${_got} rows)"

if [[ -n "$IVTAB" ]]; then
    [[ -s "$BAM" ]] || { echo "FATAL: ${SAMPLE}: missing ${BAM}; the interval
   mode needs it for mapping-quality-filtered depth. Pass BAMDIR, or
   MAPQSCOPE=never to fall back to the gVCF everywhere." >&2; exit 1; }
    # The inherited class needs positive evidence from the second frame before
    # it is promoted to ALT; see bin/p5_sv_genotype.py and
    # sv2frame/TWOFRAME_RESULTS.md. TWOFRAME points at that table; without it
    # inherited candidates come back NOCALL, which is the honest state.
    # TRUST_INHERITED=1 restores the pre-2026-09-27 behaviour.
    TWOFRAME="${TWOFRAME:-sv2frame/${COHORT_TAG:-$(basename "$(dirname "$OUTDIR")")}/per_sample/${SAMPLE}.2frame.tsv}"
    TF_ARGS=()
    [[ -s "$TWOFRAME" ]] && TF_ARGS+=(--twoframe "$TWOFRAME")
    [[ -n "${TRUST_INHERITED:-}" ]] && TF_ARGS+=(--trust-inherited)
    "$MTB_PY" bin/p5_sv_genotype.py --intervals "$IVTAB" --sample "$SAMPLE" \
        --reference "$REFID" --gvcf "$GVCF" --bam "$BAM" \
        --mapq-scope "${MAPQSCOPE:-is6110}" --min-mapq "${MINMAPQ:-30}" \
        "${TF_ARGS[@]}" \
        --projected "$POS" --out "${OUT}.tmp"
else
    "$MTB_PY" bin/p5_sv_genotype.py --sv-matrix "$SVMATRIX" --sample "$SAMPLE" \
        --gvcf "$GVCF" --projected "$POS" --out "${OUT}.tmp"
fi
mv -f "${OUT}.tmp" "$OUT"
rm -f "${WORK}/${SAMPLE}.probes.panel.txt"
echo "[P5svgt] ${SAMPLE}: done"
