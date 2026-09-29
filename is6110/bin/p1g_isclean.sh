#!/usr/bin/env bash
#SBATCH --job-name=P1g_isclean
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH --output=slurm/P1g_%A_%a.out
#SBATCH --error=slurm/P1g_%A_%a.err
#
# P1g: realign to the IS6110-clean reference and read junctions off it.
#
#   sbatch --array=1-23 is6110/bin/p1g_isclean.sh
#   bash is6110/bin/p1g_isclean.sh <sample>
#   bash is6110/bin/p1g_isclean.sh --summary
#
# THE POINT
# On normal H37Rv a read carrying IS6110 sequence has sixteen places it can
# align, so the aligner can tuck a junction read wholly inside a reference copy
# and the junction is lost rather than mis-scored. On is6110/assets/
# H37Rv.isclean.fasta the chromosome carries no element sequence at all, so such
# a read MUST clip at its junction. H37Rv's own copies then show up as ordinary
# insertions exactly like novel ones, which is the thing that would replace the
# broken depth-based reference-copy call (is6110/README.md section 4).
#
# The reference also carries the canonical element as a SEPARATE contig, so the
# clipped part still has somewhere to realign and the SA tag still resolves the
# junction -- and attribution becomes a contig-name test with no tolerance and no
# ambiguity about which copy. See is6110/bin/is6110_build_isclean.py.
#
# THREE THINGS THIS PRODUCES that the arm does not otherwise have:
#   1. junctions visible here but not on normal H37Rv = a sensitivity number
#   2. a working reference-copy presence call
#   3. a third copy-number estimate -- on an IS-free backbone the junction count
#      IS the copy number. Prediction: it tracks P1d family depth, ~1 for
#      lineage 7 and ~24 for lineage 2. If it does not, one of them is wrong.
#
# NOT A REPLACEMENT FOR P1d. Family depth needs the copies present in the
# reference; it cannot be computed here. Run both.
#
# Resources match P1f, which does the same CRAM -> FASTQ -> align work.
set -euo pipefail

for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
source "$_mtb_env"
cd "${SLURM_SUBMIT_DIR:-.}"

COHORT="${COHORT:-refbias/cohort.pilot.tsv}"
OUTDIR="${OUTDIR:-refbias/p1g}"
CRAMS="${CRAMS:-refbias/cohort.crams.tsv}"
COLLECTION="${COLLECTION:-${MTB_CRAM_ROOT}}"
CREF="${CREF:-${MTB_CRAM_REF}}"
CLEAN="${CLEAN:-is6110/assets/H37Rv.isclean.fasta}"
CLEAN_CONTIG="${CLEAN_CONTIG:-NC_000962.3_isclean}"
ELEMENT_CONTIG="${ELEMENT_CONTIG:-IS6110}"
BWA="${MTB_BWA:-${MTB_QC_BIN}/bwa}"
NT="${NT:-4}"
# TEMPORARY FILES. Deliberately NOT named SCRATCH: this cluster exports
# SCRATCH=/n/netscratch site-wide, so "${SCRATCH:-...}" never reached its
# fallback. Every task then wrote its temporary files to the root of the shared
# filesystem -- which failed on permissions -- and the EXIT trap below ran as
# "rm -rf /n/netscratch". Forty-four tasks spent 45 minutes recursing the whole
# 3.8 PB filesystem before being cancelled. TMPBASE is a base directory only;
# the path actually removed comes from mktemp -d, so it is always one this task
# created, as bin/assemble_from_cram.sh does it.
TMPBASE="${MTB_TMPBASE:-${TMPDIR:-/tmp}}"
WORKTMP=""
HELD_LOCK=""
mkdir -p "$OUTDIR" slurm

# The step is normally $1. bin/refbias_run.sh submits through sbatch, which
# passes environment but not positional arguments, so P1GSTEP is accepted as an
# equivalent -- the same accommodation p2_call.sh through p5_merge.sh already
# make. Without it the summary step died with "FATAL: no task id", because a
# non-array job has no SLURM_ARRAY_TASK_ID and the script fell through to the
# per-sample branch. That is what the first back-filled l7 run hit.
if [[ "${1:-${P1GSTEP:-}}" == "--summary" ]]; then
    exec "$MTB_PY" is6110/bin/is6110_isclean_summary.py \
        --dir "$OUTDIR" --crossmap is6110/assets/H37Rv.isclean.crossmap.tsv \
        --out "${OUTDIR}/isclean_summary.tsv"
fi

for f in "$CLEAN" "${CLEAN}.bwt" "$BWA"; do
    [[ -e "$f" ]] || { echo "FATAL: missing ${f}" >&2; exit 1; }
done

if [[ $# -ge 1 ]]; then
    SAMPLE="$1"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "FATAL: no task id" >&2; exit 1; }
    SAMPLE="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$COHORT")"
fi
[[ -n "${SAMPLE:-}" ]] || { echo "FATAL: no cohort row" >&2; exit 1; }

BAM="${OUTDIR}/${SAMPLE}.isclean.bam"
JUNC="${OUTDIR}/${SAMPLE}.junctions.tsv"
[[ -s "$JUNC" ]] && { echo "[P1g] ${SAMPLE}: already done"; exit 0; }

REL="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $2}' "$CRAMS" 2>/dev/null || true)"
[[ -n "$REL" ]] || { echo "FATAL: ${SAMPLE}: no CRAM path in ${CRAMS}" >&2; exit 1; }
CRAM="${COLLECTION}/${REL}"
[[ -s "$CRAM" ]] || { echo "FATAL: ${SAMPLE}: no CRAM at ${CRAM}" >&2; exit 1; }

# CONCURRENCY. Two tasks on the same sample corrupt each other's work: they share
# every output path, so one mv's its sorted alignment over a file the other has
# already indexed. That happened to SAMN13568032 on 2026-09-18 when a duplicate
# probe task was submitted five minutes after the first, and the resulting
# alignment had to be discarded. The completion test above only catches a
# FINISHED sample, never one in progress, which is precisely the gap. mkdir is
# atomic, so it is the lock.
#
# A task killed with SIGKILL skips its EXIT trap and leaves the lock behind. The
# message below names the path so it can be removed by hand.
LOCK="${OUTDIR}/.lock.${SAMPLE}"
if ! mkdir "$LOCK" 2>/dev/null; then
    echo "[P1g] ${SAMPLE}: ${LOCK} is held by another task, leaving this sample to it"
    echo "[P1g] ${SAMPLE}: if no such task is running, remove ${LOCK} and re-run"
    exit 0
fi
HELD_LOCK="$LOCK"
trap 'if [[ -n "${WORKTMP:-}" && -d "$WORKTMP" ]]; then rm -rf "$WORKTMP"; fi
      if [[ -n "${HELD_LOCK:-}" && -d "$HELD_LOCK" ]]; then rmdir "$HELD_LOCK"; fi' EXIT

[[ -d "$TMPBASE" && -w "$TMPBASE" ]] || { echo "FATAL: TMPBASE ${TMPBASE} is not a writable directory" >&2; exit 1; }
WORKTMP="$(mktemp -d "${TMPBASE}/p1g.${SLURM_JOB_ID:-$$}.XXXXXX")"

if [[ ! -s "$BAM" ]]; then
    echo "[P1g] ${SAMPLE}: CRAM -> FASTQ -> IS-clean alignment"
    "$MTB_SAMTOOLS" collate -@ "$NT" -u -O --reference "$CREF" "$CRAM" "${WORKTMP}/c" \
      | "$MTB_SAMTOOLS" fastq -@ 2 -n \
            -1 "${WORKTMP}/r1.fq.gz" -2 "${WORKTMP}/r2.fq.gz" -0 /dev/null -s /dev/null -
    "$BWA" mem -t "$NT" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
        "$CLEAN" "${WORKTMP}/r1.fq.gz" "${WORKTMP}/r2.fq.gz" 2> "${OUTDIR}/${SAMPLE}.bwa.log" \
      | "$MTB_SAMTOOLS" sort -@ 2 -o "${BAM}.${SLURM_JOB_ID:-$$}.tmp" -
    mv -f "${BAM}.${SLURM_JOB_ID:-$$}.tmp" "$BAM"
    "$MTB_SAMTOOLS" index "$BAM"
fi

# Every position on the clean chromosome carrying a clipped end, then the
# junctions at those positions. --element-contig makes attribution exact.
#
# COORDINATES: both numbers below are 1-BASED, because is6110_junctions.py takes
# 1-based sites and converts internally. SAM POS ($4) is already the 1-based
# leftmost aligned base, so a left clip is $4 and a right clip is $4 + reflen - 1.
# Emitting $4-1 for one and $4+reflen-1 for the other mixes conventions within a
# single file; that bug was written here first and caught only because position 0
# made pysam raise. It is the third coordinate slip in this arm -- see
# STAGE_D_RESULTS.md section 7.
echo "[P1g] ${SAMPLE}: scanning junctions"
"$MTB_SAMTOOLS" view -@ 2 "$BAM" "$CLEAN_CONTIG" \
  | awk 'BEGIN{OFS="\t"} $6 ~ /^[0-9]+[SH]/ {print $4} $6 ~ /[0-9]+[SH]$/ {
         n=0; c=$6; while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
           l=substr(c,RSTART,RLENGTH); op=substr(l,length(l));
           if (op ~ /[MDN=X]/) n+=substr(l,1,length(l)-1);
           c=substr(c,RLENGTH+1) } print $4+n-1 }' \
  | sort -n | uniq -c | awk -v m="${MIN_CLIPS:-5}" '$1>=m{print $2}' \
  > "${WORKTMP}/pos.txt"
{ echo pos; cat "${WORKTMP}/pos.txt"; } > "${WORKTMP}/sites.tsv"
echo "[P1g] ${SAMPLE}: $(wc -l < "${WORKTMP}/pos.txt") candidate positions"

"$MTB_PY" is6110/bin/is6110_junctions.py \
    --alignment "$BAM" --contig "$CLEAN_CONTIG" --sample "$SAMPLE" \
    --sites "${WORKTMP}/sites.tsv" --element-contig "$ELEMENT_CONTIG" \
    --min-clips "${MIN_CLIPS:-5}" --out "${JUNC}.tmp"
mv -f "${JUNC}.tmp" "$JUNC"

# family abundance from the same alignment: element contig against chromosome
EL=$("$MTB_SAMTOOLS" depth -a -Q 0 -q 0 -r "$ELEMENT_CONTIG" "$BAM" | awk '{s+=$3;n++} END{print (n?s/n:0)}')
CH=$("$MTB_SAMTOOLS" depth -a -Q 0 -q 0 -r "$CLEAN_CONTIG" "$BAM" | awk '{s+=$3;n++} END{print (n?s/n:0)}')
printf 'sample\telement_depth\tchrom_depth\telement_ratio\n%s\t%.3f\t%.3f\t%.4f\n' \
    "$SAMPLE" "$EL" "$CH" "$(awk -v a="$EL" -v b="$CH" 'BEGIN{print (b?a/b:0)}')" \
    > "${OUTDIR}/${SAMPLE}.elementdepth.tsv"
echo "[P1g] ${SAMPLE}: element/chrom depth ratio $(awk -v a="$EL" -v b="$CH" 'BEGIN{printf "%.3f", (b?a/b:0)}')"
