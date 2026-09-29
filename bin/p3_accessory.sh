#!/usr/bin/env bash
#SBATCH --job-name=P3_acc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 0-01:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH --output=slurm/P3_%A_%a.out
#SBATCH --error=slurm/P3_%A_%a.err
#
# P3: accessory content, from the reads the MATCHED reference could not place.
#
#   sbatch --array=1-23 bin/p3_accessory.sh
#   bash bin/p3_accessory.sh <sample>
#   bash bin/p3_accessory.sh --summary
#
# Two arms, because the accessory catalogue has two kinds of locus and they need
# opposite treatments.
#
#   arm A, presence/absence  49 novel + 79 mosaic contigs. This sequence is
#                            largely absent from H37Rv and from most references,
#                            so a carrier shows up as reads that will not place
#                            on R but do place on the contig. Reuses T10's
#                            mechanics, with one change: the unplaced reads come
#                            from the R-aligned BAM rather than the H37Rv CRAM,
#                            so what is recovered is what the MATCHED reference
#                            could not explain -- a much smaller and more
#                            interesting set.
#
#   arm B, copy number       674 loci, of which IS6110 is 501. Here the sequence
#                            IS present in H37Rv (union coverage > 0.90); what
#                            differs is how many copies. Unplaced reads say
#                            nothing, because the reads place fine -- they just
#                            pile onto R's copies. The signal is DEPTH: a sample
#                            with three copies where R has one shows roughly
#                            three times median depth over that locus.
#
# Arm B needs each locus in R coordinates, and the catalogue anchors them on
# H37Rv. `odgi position` is run in the reverse direction to the rest of this
# pipeline -- H37Rv path positions in, R path positions out -- once per sample for
# all 674 loci. A locus whose anchor does not project is reported as such rather
# than silently dropped, since "no copy-number call" and "zero copies" are
# different statements and conflating them is the encoding hazard T6 measured.
#
# Depth uses `samtools depth -aa`, not -a: -a omits contigs with no coverage
# entirely, which in T10 made 8 of 62 contigs look observed when the real span
# was 128.
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
P2WORK="${P2WORK:-refbias/work/p2}"
OUTDIR="${OUTDIR:-refbias/p3}"
WORK="${WORK:-refbias/work/p3}"
LOCI="${BUILD}/assets/accessory_loci.tsv"
PATHS="${BUILD}/assets/paths.txt"
OG="${OG:-$(ls graphs/CX333.s10k.k23.K15/*.smooth.final.og 2>/dev/null | head -1)}"
H37RV_PATH="${H37RV_PATH:-GCF_000195955#1#NC_000962.3}"
BWA="${MTB_BWA:-${MTB_QC_BIN}/bwa}"
ODGI="${MTB_ODGI:?MTB_ODGI is unset; see config/project_env.sh}"
NT="${SLURM_CPUS_PER_TASK:-4}"

for f in "$REFMAP" "$LOCI" "$PATHS" "$OG" "$BWA" "$ODGI"; do
    [[ -e "$f" ]] || { echo "FATAL: missing: $f" >&2; exit 1; }
done
mkdir -p "$OUTDIR" "$WORK" slurm

# The step is normally $1. bin/refbias_run.sh submits through sbatch,
# which passes environment but not positional arguments, so P3STEP is
# accepted as an equivalent. $1 wins if both are given.
if [[ "${1:-${P3STEP:-}}" == "--summary" ]]; then
    exec "$MTB_PY" bin/p3_summary.py --refmap "$REFMAP" --dir "$OUTDIR" \
        --loci "$LOCI" --out "${OUTDIR}/p3_summary.tsv"
fi

# --- the combined accessory panel, indexed once per build --------------------
# Build-scoped, so it belongs beside P0's other assets rather than in a sample
# work dir. Race-guarded the same way P0 indexes references, since array tasks
# start together.
PANEL="${BUILD}/assets/accessory_panel.fasta"
# COMPLETE MEANS EVERY FILE, and .bwt is moved into place LAST. Other tasks
# used to take .bwt alone as "indexed", and it was moved before .pac and .sa,
# so a task starting in that window ran bwa mem on a half-installed index and
# failed -- and one failed task cancels everything downstream under afterok.
_panel_ready() {
    local e
    [[ -s "$PANEL" ]] || return 1
    for e in amb ann pac sa fai bwt; do [[ -s "${PANEL}.${e}" ]] || return 1; done
}
if ! _panel_ready; then
    lock="${PANEL}.build.$$"
    cat "${BUILD}/assets/accessory_novel.fasta" \
        "${BUILD}/assets/accessory_mosaic.fasta" > "$lock"
    "$BWA" index "$lock" > "${WORK}/panel_index.log" 2>&1
    "$MTB_SAMTOOLS" faidx "$lock"
    [[ -s "$PANEL" ]] || cp -f "$lock" "$PANEL"
    for ext in amb ann pac sa fai bwt; do       # bwt last: it is the signal
        [[ -f "${lock}.${ext}" ]] && mv -f "${lock}.${ext}" "${PANEL}.${ext}"
    done
    rm -f "$lock"
fi
_panel_ready || { echo "FATAL: accessory panel index incomplete at ${PANEL}" >&2; exit 1; }

if [[ $# -ge 1 ]]; then
    SAMPLE="$1"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "FATAL: no SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    SAMPLE="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$REFMAP")"
fi
[[ -n "${SAMPLE:-}" ]] || { echo "FATAL: no refmap row" >&2; exit 1; }
REFID="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $5; exit}' "$REFMAP")"
[[ -n "$REFID" ]] || { echo "FATAL: ${SAMPLE}: no reference" >&2; exit 1; }
BAM="${P2WORK}/${SAMPLE}.bam"
[[ -s "$BAM" ]] || { echo "FATAL: ${SAMPLE}: no P2 BAM at ${BAM}; run P2 first" >&2; exit 1; }
RPATH="$(grep -m1 "^${REFID}#" "$PATHS" || true)"
[[ -n "$RPATH" ]] || { echo "FATAL: ${SAMPLE}: ${REFID} is not a graph path" >&2; exit 1; }

if [[ -s "${OUTDIR}/${SAMPLE}.accessory.tsv" && -s "${OUTDIR}/${SAMPLE}.copynumber.tsv" ]]; then
    echo "[P3] ${SAMPLE}: already done"; exit 0
fi
echo "[P3] ${SAMPLE}: reference ${REFID}, build ${BUILD_ID}"

# --- arm A: reads R could not place --------------------------------------------
TOTAL=$("$MTB_SAMTOOLS" view -c "$BAM")
"$MTB_SAMTOOLS" view -@ 2 -u -f 4 "$BAM" \
  | "$MTB_SAMTOOLS" fastq -@ 1 -n - > "${WORK}/${SAMPLE}.unplaced.fq" 2>/dev/null
NUM=$(( $(wc -l < "${WORK}/${SAMPLE}.unplaced.fq") / 4 ))
echo "[P3] ${SAMPLE}: ${NUM} unplaced of ${TOTAL} reads against ${REFID}"

ACC="${OUTDIR}/${SAMPLE}.accessory.tsv"
if [[ "$NUM" -eq 0 ]]; then
    # Emit every contig at zero rather than an empty file: a carrier-free sample
    # and an unrun sample must not look the same downstream.
    { printf 'contig\tlen\treads\tcovered_bp\tbreadth\tmean_depth\n'
      awk '/^>/{name=substr($1,2)} /^>/{next} {l[name]+=length($0)} END{for(n in l) printf "%s\t%d\t0\t0\t0.0000\t0.0000\n",n,l[n]}' \
          "$PANEL"
    } > "$ACC"
else
    "$BWA" mem -t "$NT" "$PANEL" "${WORK}/${SAMPLE}.unplaced.fq" \
        2> "${WORK}/${SAMPLE}.panel.bwa.log" \
      | "$MTB_SAMTOOLS" sort -@ 2 -o "${WORK}/${SAMPLE}.panel.bam" -
    "$MTB_SAMTOOLS" index "${WORK}/${SAMPLE}.panel.bam"
    "$MTB_SAMTOOLS" depth -aa "${WORK}/${SAMPLE}.panel.bam" > "${WORK}/${SAMPLE}.panel.depth"
    "$MTB_SAMTOOLS" idxstats "${WORK}/${SAMPLE}.panel.bam" > "${WORK}/${SAMPLE}.panel.idxstats"
    "$MTB_PY" bin/p3_panel_support.py \
        --depth "${WORK}/${SAMPLE}.panel.depth" \
        --idxstats "${WORK}/${SAMPLE}.panel.idxstats" --out "$ACC"
fi
rm -f "${WORK}/${SAMPLE}.unplaced.fq"

# --- arm B: copy number by depth ----------------------------------------------
# project each locus anchor from H37Rv into R, one odgi call for all of them
# tr -d '\r' is load-bearing. Every TSV this project writes has CRLF line
# terminators, because Python's csv.writer on a file opened without newline=""
# emits \r\n. Python readers never notice, since text-mode universal newlines
# strips it, so this has been invisible for the whole project -- but awk sees the
# last field as "copy_number\r", matches nothing, and writes an empty anchor
# file. The downstream guard caught it; a silent zero would have been worse.
tr -d '\r' < "$LOCI" \
  | awk -F'\t' -v p="$H37RV_PATH" 'NR>1 && $12=="copy_number" {print p","($2-1)",+"}' \
  > "${WORK}/${SAMPLE}.anchors.txt"
[[ -s "${WORK}/${SAMPLE}.anchors.txt" ]] \
    || { echo "FATAL: ${SAMPLE}: no copy_number anchors extracted from ${LOCI}" >&2; exit 1; }
NANCH=$(wc -l < "${WORK}/${SAMPLE}.anchors.txt")
# WHICH FRAME IS THIS SAMPLE'S BAM IN? Measured, not assumed. The depth this
# script reads comes from the P2 alignment, so the P2 VCF beside it settles the
# question: P1 and P2 align to ${BUILD}/refs and are refs-frame, while the stage
# 2 simulation aligned to the rotated panel sequence and is panel-frame, and
# nothing in a VCF header tells them apart. If this BAM were panel-frame the
# conversion below would be the defect rather than the fix, so stop rather than
# guess. graphframe/docs/GRAPH_FRAME_RESOLUTION.md section 6.
P2VCF="${P2DIR:-refbias/p2}/${SAMPLE}.vcf.gz"
if [[ -s "$P2VCF" ]]; then
    VFRAME="$("$MTB_PY" graphframe/bin/frame_detect.py --vcf "$P2VCF" \
        --accession "$REFID" --refs "${BUILD}/refs" | cut -f2)"
    case "$VFRAME" in
        refs|either) ;;
        *) echo "FATAL: ${SAMPLE}: the P2 calls are in the '${VFRAME}' frame, not
   refs; the conversion in this script would be wrong for it" >&2; exit 1 ;;
    esac
fi

# The H37Rv anchors below are in the refs frame and `odgi position` reads the
# panel frame the graph was built in; the target R coordinates come back in the
# panel frame while the depth and gVCF files they index are in the refs frame.
# For 110 of 333 accessions those frames differ by a rotation and for 22 also by
# strand, so both conversions are load-bearing. graphframe/docs/GRAPH_FRAME_RESOLUTION.md
"$MTB_PY" graphframe/bin/frame_convert.py to-panel \
    < "${WORK}/${SAMPLE}.anchors.txt" > "${WORK}/${SAMPLE}.anchors.panel.txt"
"$ODGI" position -i "$OG" -F "${WORK}/${SAMPLE}.anchors.panel.txt" -r "$RPATH" -t "$NT" \
    2> "${WORK}/${SAMPLE}.odgi.log" \
  | "$MTB_PY" graphframe/bin/frame_convert.py from-panel \
    > "${WORK}/${SAMPLE}.anchors.pos"
"$MTB_SAMTOOLS" depth -aa "$BAM" > "${WORK}/${SAMPLE}.r.depth"
cp -f "${WORK}/${SAMPLE}.r.depth" "${WORK}/${SAMPLE}.r.depth.keep"
"$MTB_PY" bin/p3_copynumber.py \
    --loci "$LOCI" --positions "${WORK}/${SAMPLE}.anchors.pos" \
    --depth "${WORK}/${SAMPLE}.r.depth" --sample "$SAMPLE" --reference "$REFID" \
    --out "${OUTDIR}/${SAMPLE}.copynumber.tsv"

# --- IS6110 family copy number, from the element spans themselves -------------
# Depth at ONE copy of a repeated element reflects the whole family, because
# reads from every copy map onto each. That makes the per-locus table useless for
# total copy number and makes this the right measurement: median depth across
# H37Rv's 16 element spans, projected into the sample's reference, divided by the
# sample's genome median.
ELGFF="${BUILD}/assets/is6110_elements.gff"
if [[ -s "$ELGFF" ]]; then
    awk -F'\t' -v p="$H37RV_PATH" '!/^#/ && NF>4 {print p","($4-1)",+\n"p","($5-1)",+"}' \
        "$ELGFF" > "${WORK}/${SAMPLE}.elements.txt"
    "$MTB_PY" graphframe/bin/frame_convert.py to-panel \
        < "${WORK}/${SAMPLE}.elements.txt" > "${WORK}/${SAMPLE}.elements.panel.txt"
    "$ODGI" position -i "$OG" -F "${WORK}/${SAMPLE}.elements.panel.txt" -r "$RPATH" \
        -t "$NT" 2>> "${WORK}/${SAMPLE}.odgi.log" \
      | "$MTB_PY" graphframe/bin/frame_convert.py from-panel \
        > "${WORK}/${SAMPLE}.elements.pos"
    "$MTB_PY" bin/p3_element_depth.py --gff "$ELGFF" \
        --positions "${WORK}/${SAMPLE}.elements.pos" \
        --depth "${WORK}/${SAMPLE}.r.depth.keep" --sample "$SAMPLE" \
        --reference "$REFID" --out "${OUTDIR}/${SAMPLE}.elementdepth.tsv"
fi
rm -f "${WORK}/${SAMPLE}.r.depth" "${WORK}/${SAMPLE}.r.depth.keep" \
      "${WORK}/${SAMPLE}.panel.depth" "${WORK}/${SAMPLE}.elements.txt" \
      "${WORK}/${SAMPLE}.elements.panel.txt" "${WORK}/${SAMPLE}.anchors.panel.txt"
echo "[P3] ${SAMPLE}: ${NANCH} copy-number loci projected; done"
