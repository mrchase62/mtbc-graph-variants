#!/usr/bin/env bash
# LABELLING NOTE, 2026-09-21.  This script scores on the A/B/C tiers, which have
# been replaced by two independent columns -- `evidence` (two_sided/one_sided, a
# quality ordering) and `site_class` (ref_shared/ref_lacking, not an ordering).
# The rule now lives in is6110_promote_sites.py and is described in
# is6110/docs/PROMOTED_TIER.md.  `tierAB` below means `evidence = two_sided`, and
# it excludes promoted one-sided calls, so any ratio it forms against a
# depth-based estimate reads low.  The numbers this script produced are still the
# numbers that were measured; do not read tier A as more confident than tier B.

#SBATCH --job-name=P1i_matched
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH --output=slurm/P1i_%A_%a.out
#SBATCH --error=slurm/P1i_%A_%a.err
#
# P1i stage 3: align each isolate to ITS OWN IS-clean matched reference and read
# junctions off it.
#
#   sbatch --array=1-23 is6110/bin/p1i_matched.sh
#   bash is6110/bin/p1i_matched.sh <sample>
#
# WHAT THIS CHANGES FROM P1g
# P1g aligns every isolate to one IS-clean H37Rv. This aligns each isolate to an
# IS-clean build of the reference the SNP-distance selector chose for it
# (refbias/p1/refmap.tsv, 18 distinct references for 23 isolates, distance 78 to
# 1657). Everything else -- the detector, the read-through filter, the
# clustering, the element-side scan -- is held identical, so a difference
# between P1g and P1i is attributable to the reference and nothing else.
#
# WHY IT SHOULD MATTER, measured rather than assumed
# is6110/docs/P1I_ACCESSORY_SITES.md section 6: across the 23 isolates, 1,829
# reads are unmapped against H37Rv while their mate sits on the element contig,
# and 89.8% of them map to the isolate's own matched reference against 0.2% to
# H37Rv. They cluster into 29 candidate loci. That sequence is ordinary
# chromosome in the matched reference, so those loci need no accessory panel and
# no extra contigs -- only this alignment.
#
# PREDICTIONS, stated before the run so they cannot be chosen afterwards
# (P1I_PLAN.md stage 4):
#   1. the one-sided stack fraction falls from 69/237
#   2. the 29 loci in refbias/p1h/unplaced/loci.tsv appear as junction sites
#   3. tier A+B rises from 161
#   4. isolates whose reference is far away gain least
# NOT predicted: the four low-copy isolates gain. Their flanks sit in a direct
# repeat, not outside the reference, and they contributed 0-4 of those reads.
#
# CAVEAT CARRIED FORWARD. 18 references means 18 coordinate systems. Per-isolate
# counting is unaffected; comparing a site BETWEEN isolates is not, and that
# projection is deliberately not attempted here.
#
# Resources match p1g_isclean.sh: same CRAM -> FASTQ -> align work.
set -euo pipefail

for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
source "$_mtb_env"
cd "${SLURM_SUBMIT_DIR:-.}"

COHORT="${COHORT:-refbias/cohort.pilot.tsv}"
OUTDIR="${OUTDIR:-refbias/p1i}"
CRAMS="${CRAMS:-refbias/cohort.crams.tsv}"
COLLECTION="${COLLECTION:-${MTB_CRAM_ROOT}}"
CREF="${CREF:-${MTB_CRAM_REF}}"
REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
CLEANDIR="${CLEANDIR:-is6110/assets/isclean_matched}"
ELEMENT_CONTIG="${ELEMENT_CONTIG:-IS6110}"
BWA="${MTB_BWA:-${MTB_QC_BIN}/bwa}"
NT="${NT:-4}"
# TEMPORARY FILES: base directory only, path from mktemp -d, so the argument to
# rm -rf is always one this task created. is6110/README.md section 7.
TMPBASE="${MTB_TMPBASE:-${TMPDIR:-/tmp}}"
WORKTMP=""
HELD_LOCK=""
mkdir -p "$OUTDIR" slurm

if [[ $# -ge 1 ]]; then
    SAMPLE="$1"
else
    # This script has NO --summary mode, unlike its fixed-arm sibling: the
    # matched arm's per-cohort aggregation is pass p1iv (bin/p1i_vcf.sh), which
    # labels, reconciles and writes the VCFs from the whole cohort's stacks.
    # bin/refbias_run.sh therefore submits p1i as an array only, and p1iv
    # depends on that array directly.
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "FATAL: no task id" >&2; exit 1; }
    SAMPLE="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$COHORT")"
fi
[[ -n "${SAMPLE:-}" ]] || { echo "FATAL: no cohort row" >&2; exit 1; }

REFID="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $5}' "$REFMAP")"
[[ -n "$REFID" ]] || { echo "FATAL: ${SAMPLE}: no matched reference in ${REFMAP}" >&2; exit 1; }
CLEAN="${CLEANDIR}/${REFID}.isclean.fasta"
[[ -s "${CLEAN}.bwt" ]] || { echo "FATAL: ${CLEAN} not built or not indexed -- run is6110/bin/p1i_build_matched.sh" >&2; exit 1; }
# The contig name comes from the manifest when there is a row, and otherwise
# from the build's own FASTA header, which is where the manifest got it. The
# manifest covers 6 references against 65 builds on disk -- it was written for
# the pilot and never extended -- so requiring a row made every reference built
# later fail with "no manifest row" even though its index was present and
# usable. The header is the more direct source and cannot go stale.
CLEAN_CONTIG="$(awk -F'\t' -v r="$REFID" '$1==r{print $2"_isclean"}' "${CLEANDIR}/manifest.tsv" 2>/dev/null)"
if [[ -z "$CLEAN_CONTIG" ]]; then
    CLEAN_CONTIG="$(head -1 "$CLEAN" | sed 's/^>//; s/[[:space:]].*//')"
    echo "[P1i] ${REFID}: no manifest row; contig ${CLEAN_CONTIG} read from ${CLEAN}"
fi
[[ -n "$CLEAN_CONTIG" ]] || { echo "FATAL: ${REFID}: cannot determine the clean contig name" >&2; exit 1; }

BAM="${OUTDIR}/${SAMPLE}.isclean.bam"
JUNC="${OUTDIR}/${SAMPLE}.junctions.tsv"
# Done only when EVERY per-sample output exists. The junction table alone
# was the marker, but it is written first; a task that died in the later steps
# was then skipped on rerun and left without them.
if [[ -s "${OUTDIR}/${SAMPLE}.junctions.tsv" && -s "${OUTDIR}/${SAMPLE}.elstacks.tsv" && -s "${OUTDIR}/${SAMPLE}.elementdepth.tsv" ]]; then
    echo "[P1i] ${SAMPLE}: already done"; exit 0
fi

REL="$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $2}' "$CRAMS" 2>/dev/null || true)"
[[ -n "$REL" ]] || { echo "FATAL: ${SAMPLE}: no CRAM path" >&2; exit 1; }
CRAM="${COLLECTION}/${REL}"
[[ -s "$CRAM" ]] || { echo "FATAL: ${SAMPLE}: no CRAM at ${CRAM}" >&2; exit 1; }

# CONCURRENCY: two tasks on one sample share every output path. mkdir is atomic.
LOCK="${OUTDIR}/.lock.${SAMPLE}"
if ! mkdir "$LOCK" 2>/dev/null; then
    echo "[P1i] ${SAMPLE}: ${LOCK} held by another task, leaving it"
    echo "[P1i] ${SAMPLE}: if no such task runs, remove ${LOCK} and re-run"
    exit 0
fi
HELD_LOCK="$LOCK"
trap 'if [[ -n "${WORKTMP:-}" && -d "$WORKTMP" ]]; then rm -rf "$WORKTMP"; fi
      if [[ -n "${HELD_LOCK:-}" && -d "$HELD_LOCK" ]]; then rmdir "$HELD_LOCK"; fi' EXIT

[[ -d "$TMPBASE" && -w "$TMPBASE" ]] || { echo "FATAL: TMPBASE ${TMPBASE} not writable" >&2; exit 1; }
WORKTMP="$(mktemp -d "${TMPBASE}/p1i.${SLURM_JOB_ID:-$$}.XXXXXX")"

if [[ ! -s "$BAM" ]]; then
    echo "[P1i] ${SAMPLE} -> ${REFID}: CRAM -> FASTQ -> IS-clean matched alignment"
    "$MTB_SAMTOOLS" collate -@ "$NT" -u -O --reference "$CREF" "$CRAM" "${WORKTMP}/c" \
      | "$MTB_SAMTOOLS" fastq -@ 2 -n \
            -1 "${WORKTMP}/r1.fq.gz" -2 "${WORKTMP}/r2.fq.gz" -0 /dev/null -s /dev/null -
    "$BWA" mem -t "$NT" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
        "$CLEAN" "${WORKTMP}/r1.fq.gz" "${WORKTMP}/r2.fq.gz" 2> "${OUTDIR}/${SAMPLE}.bwa.log" \
      | "$MTB_SAMTOOLS" sort -@ 2 -o "${BAM}.${SLURM_JOB_ID:-$$}.tmp" -
    mv -f "${BAM}.${SLURM_JOB_ID:-$$}.tmp" "$BAM"
    "$MTB_SAMTOOLS" index "$BAM"
fi

# COORDINATES: both numbers are 1-BASED. SAM POS is the 1-based leftmost aligned
# base, so a left clip is $4 and a right clip is $4 + reflen - 1. Mixing the two
# conventions in one file was the third coordinate slip in this arm.
echo "[P1i] ${SAMPLE}: scanning junctions on ${CLEAN_CONTIG}"
"$MTB_SAMTOOLS" view -@ 2 "$BAM" "$CLEAN_CONTIG" \
  | awk 'BEGIN{OFS="\t"} $6 ~ /^[0-9]+[SH]/ {print $4} $6 ~ /[0-9]+[SH]$/ {
         n=0; c=$6; while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
           l=substr(c,RSTART,RLENGTH); op=substr(l,length(l));
           if (op ~ /[MDN=X]/) n+=substr(l,1,length(l)-1);
           c=substr(c,RLENGTH+1) } print $4+n-1 }' \
  | sort -n | uniq -c | awk -v m="${MIN_CLIPS:-5}" '$1>=m{print $2}' \
  > "${WORKTMP}/pos.txt"
{ echo pos; cat "${WORKTMP}/pos.txt"; } > "${WORKTMP}/sites.tsv"
echo "[P1i] ${SAMPLE}: $(wc -l < "${WORKTMP}/pos.txt") candidate positions"

"$MTB_PY" is6110/bin/is6110_junctions.py \
    --alignment "$BAM" --contig "$CLEAN_CONTIG" --sample "$SAMPLE" \
    --sites "${WORKTMP}/sites.tsv" --element-contig "$ELEMENT_CONTIG" \
    --min-clips "${MIN_CLIPS:-5}" --out "${JUNC}.tmp"
mv -f "${JUNC}.tmp" "$JUNC"

# the element-side scan, same script and same settings as P1h stage 1
"$MTB_PY" is6110/bin/is6110_element_side.py \
    --alignment "$BAM" --sample "$SAMPLE" \
    --element-contig "$ELEMENT_CONTIG" --chrom-contig "$CLEAN_CONTIG" \
    --crossmap "${CLEANDIR}/${REFID}.crossmap.tsv" \
    --out "${OUTDIR}/${SAMPLE}.elstacks.tsv.tmp"
mv -f "${OUTDIR}/${SAMPLE}.elstacks.tsv.tmp" "${OUTDIR}/${SAMPLE}.elstacks.tsv"

EL=$("$MTB_SAMTOOLS" depth -a -Q 0 -q 0 -r "$ELEMENT_CONTIG" "$BAM" | awk '{s+=$3;n++} END{print (n?s/n:0)}')
CH=$("$MTB_SAMTOOLS" depth -a -Q 0 -q 0 -r "$CLEAN_CONTIG" "$BAM" | awk '{s+=$3;n++} END{print (n?s/n:0)}')
printf 'sample\treference\telement_depth\tchrom_depth\telement_ratio\n%s\t%s\t%.3f\t%.3f\t%.4f\n' \
    "$SAMPLE" "$REFID" "$EL" "$CH" "$(awk -v a="$EL" -v b="$CH" 'BEGIN{print (b?a/b:0)}')" \
    > "${OUTDIR}/${SAMPLE}.elementdepth.tsv.tmp"
mv -f "${OUTDIR}/${SAMPLE}.elementdepth.tsv.tmp" "${OUTDIR}/${SAMPLE}.elementdepth.tsv"
echo "[P1i] ${SAMPLE}: done, element/chrom ratio $(awk -v a="$EL" -v b="$CH" 'BEGIN{printf "%.3f", (b?a/b:0)}')"
