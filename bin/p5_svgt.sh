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


# Fill the projection store for one reference: ask it what it is missing
# from the probe list, project only that with odgi, fold it back.
#   fill_store <refid> <rpath> <log label>
fill_store() {
    local REFID="$1" RPATH="$2" _lab="$3"
    # Ask the store what it is missing, project only that, fold it back, then emit
    # the whole probe set from the store.
    MISS="${PROJDIR}/${REFID}.miss.$$"
    "$MTB_PY" bin/proj_store.py missing --store "$PROJSTORE" --ref "$REFID" \
        --need "$PROBES" --out "$MISS"
    if [[ -s "$MISS" ]]; then
        echo "[P5svgt] ${_lab}: projecting $(grep -c . "$MISS") new positions for ${REFID}"
        "$MTB_PY" graphframe/bin/frame_convert.py to-panel \
            < "$MISS" > "${MISS}.panel"
        "$ODGI" position -i "$OG" -F "${MISS}.panel" -r "$RPATH" -t "$NT" \
            2> "${WORK}/${_lab}.odgi.log" \
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
            echo "[P5svgt] ${_lab}: none of those positions project into ${REFID};"\
                 " the emit floor decides whether that matters"
        fi
        rm -f "${MISS}.panel" "${MISS}.pos"
    else
        echo "[P5svgt] ${_lab}: ${REFID} fully covered by the store"
    fi
    rm -f "$MISS"
}

# ---- the interval catalogue and the second frame, for interval mode --------
# Both used to be built by hand -- the catalogue once, from scale200's caller
# matrix, into refbias/assets/sv_intervals.tsv, and the two-frame tables by
# sv2frame/bin/sv_twoframe.sh -- so nothing tied them to the cohort they were
# used for. gwas1000 was genotyped against scale200's catalogue, and the merge
# then dropped every gwas1000 caller deletion that catalogue did not hold.
# They are steps of this pass now, written under the cohort's own p5 directory.
#
#   SVCAT=cohort   graph deletions plus this cohort's caller deletions called
#                  in at least SVCAT_MIN_CARRIERS isolates (default 2, decided
#                  2026-10-01); the rest stay in the VCF presence-only, flagged
#                  UNCATALOGUED. A singleton cannot show repeated origins, and
#                  caller false positives concentrate there; genotyping all of
#                  gwas1000's would have cost about 240 CPU-hours of odgi
#                  against about 85 for the 4,408 recurrent ones
#   SVCAT=graph    graph deletions only; caller deletions stay presence-only
SVCAT="${SVCAT:-cohort}"
SVCAT_MIN_CARRIERS="${SVCAT_MIN_CARRIERS:-2}"
TWOFRAMEDIR="${TWOFRAMEDIR:-${OUTDIR}/twoframe}"
COHORT_NAME="${COHORT_NAME:-}"
case "$COHORT_NAME" in
  ""|pilot|pilot_rerun) _ISTAG="" ;;
  scale100)             _ISTAG="scale_" ;;
  *)                    _ISTAG="${COHORT_NAME}_" ;;
esac
IS6110KEYS="${IS6110KEYS:-is6110/results/${_ISTAG}p1i_cohort_keys.tsv}"

STEP="${1:-${SVGTSTEP:-}}"
case "$STEP" in
  --catalogue)
    [[ -n "$IVTAB" ]] || { echo "FATAL: --catalogue needs IVTAB, the path to write" >&2; exit 1; }
    GRAPH_VCF="${GRAPH_VCF:-$(dirname "$OG")/all_variants.decomposed.vcf.gz}"
    [[ -s "$GRAPH_VCF" ]] || { echo "FATAL: no graph VCF at ${GRAPH_VCF}" >&2; exit 1; }
    CAT_ARGS=(--graph-vcf "$GRAPH_VCF")
    case "$SVCAT" in
      cohort) [[ -s "$SVMATRIX" ]] || { echo "FATAL: no ${SVMATRIX}; run p5 --pre" >&2; exit 1; }
              CAT_ARGS+=(--sv-matrix "$SVMATRIX" --min-caller-carriers "$SVCAT_MIN_CARRIERS") ;;
      graph)  ;;
      *) echo "FATAL: SVCAT must be cohort or graph, not '${SVCAT}'" >&2; exit 1 ;;
    esac
    # IS6110 landmarks from this cohort's own key table, when it has one
    if [[ -s "$IS6110KEYS" ]]; then CAT_ARGS+=(--is6110-keys "$IS6110KEYS")
    else echo "  note: no IS6110 key table at ${IS6110KEYS}; H37Rv copies only"; fi
    "$MTB_PY" bin/sv_intervals.py "${CAT_ARGS[@]}" --out "${IVTAB}.tmp"
    mv -f "${IVTAB}.tmp" "$IVTAB"
    exit 0 ;;
  --twoframe)
    # The H37Rv-frame half of the two-frame genotyper, one sample per task.
    # About 8 s a sample: the catalogue is already in H37Rv coordinates, so it
    # needs only the read collection's own CRAM. The matched-frame column is
    # left empty -- it would be this pass's own output, and the genotyper reads
    # only the H37Rv-frame state and clips from this table.
    [[ -s "$IVTAB" ]] || { echo "FATAL: no catalogue at ${IVTAB}; run --catalogue" >&2; exit 1; }
    if [[ $# -ge 2 ]]; then S="$2"
    else S="$(awk -F'\t' -v n="$(( ${SLURM_ARRAY_TASK_ID:?--twoframe is an array step} + 1 ))" 'NR==n{print $1}' "$REFMAP")"; fi
    [[ -n "$S" ]] || { echo "FATAL: no refmap row" >&2; exit 1; }
    mkdir -p "$TWOFRAMEDIR"
    TF="${TWOFRAMEDIR}/${S}.2frame.tsv"
    if [[ -s "$TF" && "$TF" -nt "$IVTAB" ]]; then
        echo "[P5svgt] ${S}: two-frame table current"; exit 0
    fi
    CRAMTAB="${CRAMTAB:-${CRAMS:?CRAMS is unset; run through bin/refbias_run.sh}}"
    REL="$(awk -F'\t' -v s="$S" '$1==s{print $2; exit}' "$CRAMTAB")"
    [[ -n "$REL" ]] || { echo "FATAL: ${S} not in ${CRAMTAB}" >&2; exit 1; }
    CRAM="${MTB_CRAM_ROOT:?MTB_CRAM_ROOT is unset}/${REL}"
    [[ -s "$CRAM" ]] || { echo "FATAL: ${S}: no CRAM at ${CRAM}" >&2; exit 1; }
    # The IS-clean H37Rv alignment from pass p1g, for element-proximal
    # intervals: on plain H37Rv a read carrying element sequence has sixteen
    # places to align. Required, so a missing one cannot silently change the
    # instrument; TWOFRAME_NO_ISCLEAN=1 for a cohort that never ran p1g.
    TF_EXTRA=()
    IC="${P1GDIR:-refbias/p1g}/${S}.isclean.bam"
    if [[ -s "$IC" ]]; then TF_EXTRA+=(--isclean-bam "$IC")
    elif [[ -z "${TWOFRAME_NO_ISCLEAN:-}" ]]; then
        echo "FATAL: ${S}: no p1g alignment at ${IC} (restore it with" \
             "bin/archive_alignments.sh --restore, or set TWOFRAME_NO_ISCLEAN=1)" >&2
        exit 1
    fi
    "${MTB_PY_VT:-$MTB_PY}" sv2frame/bin/sv_twoframe.py --intervals "$IVTAB" \
        --sample "$S" --cram "$CRAM" --h37rv-fasta "${MTB_CRAM_REF:?MTB_CRAM_REF is unset}" \
        "${TF_EXTRA[@]}" --out "${TF}.tmp.$$"
    mv -f "${TF}.tmp.$$" "$TF"
    exit 0 ;;
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
  --project)
    # ONE TASK PER REFERENCE, BEFORE THE PER-SAMPLE ARRAY. Filling the store
    # from inside the per-sample tasks made every sample of an as-yet
    # unprojected reference project it at the same time: on scale200 49 of 139
    # projections repeated one another, about a third of the pass's CPU. Task i
    # fills the store for the i-th distinct reference of the refmap (sorted);
    # tasks past the last reference exit at once. The --states step still
    # fills anything missing, so it is a no-op there when this step ran.
    [[ -s "$PROBES" ]] || { echo "FATAL: run --probes first" >&2; exit 1; }
    I="${SLURM_ARRAY_TASK_ID:?--project is an array step}"
    mapfile -t _refs < <(awk -F'\t' 'NR>1 && $5!=""{print $5}' "$REFMAP" | sort -u)
    if [[ "$I" -gt "${#_refs[@]}" ]]; then
        echo "[P5svgt] task ${I}: only ${#_refs[@]} references; nothing to do"; exit 0
    fi
    REFID="${_refs[$((I - 1))]}"
    RPATH="$(grep -m1 "^${REFID}#" "$PATHS" || true)"
    [[ -n "$RPATH" ]] || { echo "FATAL: ${REFID} not a graph path" >&2; exit 1; }
    fill_store "$REFID" "$RPATH" "ref_${REFID}"
    exit 0 ;;
  --states) ;;
  *) echo "usage: $0 --catalogue | --twoframe | --probes | --project | --states | --merge   (or SVGTSTEP=...)" >&2; exit 2 ;;
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

POS="${PROJDIR}/${REFID}.$(sha1sum "$PROBES" | cut -c1-16).pos"
fill_store "$REFID" "$RPATH" "$SAMPLE"
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
    TWOFRAME="${TWOFRAME:-${TWOFRAMEDIR}/${SAMPLE}.2frame.tsv}"
    TF_ARGS=()
    if [[ -n "${TRUST_INHERITED:-}" ]]; then
        TF_ARGS+=(--trust-inherited)
    elif [[ -s "$TWOFRAME" && "$TWOFRAME" -nt "$IVTAB" ]]; then
        TF_ARGS+=(--twoframe "$TWOFRAME")
    else
        # Missing or older than the catalogue: a stale table keys on interval
        # ids that may no longer exist, which reads as "no second frame" for
        # every interval and silently turns every inherited ALT into NOCALL.
        echo "FATAL: ${SAMPLE}: no two-frame table newer than ${IVTAB} at" \
             "${TWOFRAME}; run --twoframe (or TRUST_INHERITED=1)" >&2
        exit 1
    fi
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
