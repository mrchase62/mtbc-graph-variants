#!/usr/bin/env bash
# LEVEL 1 for one sample: presence/absence of each accessory locus.
#
# Thin wrapper so the two-instrument call is reproducible and so the CRAM and
# matched-reference lookups come from the same tables every other arm uses:
# the cohort CRAM table for the alignment and p1/refmap.tsv for the reference
# whose accessory content the reference route consults.
set -euo pipefail
# Found relative to this script, and a failure is fatal: `|| true` on a
# relative path let the script run with every tool variable empty.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/config/project_env.sh"

S="${1:?usage: locus_presence_one.sh <sample> [cohort]}"
COHORT_TAG="${2:-${COHORT_TAG:-gwas1000}}"
CRAMTAB="${CRAMTAB:-refbias/${COHORT_TAG}.crams.tsv}"
CRAMROOT="${CRAMROOT:-${MTB_CRAM_ROOT}}"
REFMAP="${REFMAP:-refbias/${COHORT_TAG}/p1/refmap.tsv}"
OUTDIR="${OUTDIR:-accessory/${COHORT_TAG}}"

REL="$(awk -F'\t' -v s="$S" '$1==s{print $2; exit}' "$CRAMTAB")"
[[ -n "$REL" ]] || { echo "FATAL: ${S} not in ${CRAMTAB}" >&2; exit 1; }
REF="$(awk -F"\t" -v s="$S" 'NR==1{for(i=1;i<=NF;i++)if($i=="reference")c=i;next} $1==s{print $c; exit}' "$REFMAP")"

# the catalogue: the build's (locus_presence_array.sh passes it), never
# locus_presence.py's default under accessory/assets/
CAT="${ACC_CATALOGUE:?ACC_CATALOGUE is unset: <build>/assets/accessory_catalogue (run through locus_presence_array.sh)}"
# THE P1 ALIGNMENT TO H37Rv (refbias_run.sh exports P1WORK; p3acc waits for
# P1): the junction route reads the loci H37Rv carries at their anchor in the
# preprocessed reads every other arm reads. Without P1WORK that route reads
# the CRAM, and says so; with it a missing alignment is FATAL.
H37BAM=""
if [[ -n "${P1WORK:-}" ]]; then
    H37BAM="${P1WORK}/${S}.h37rv.bam"
    [[ -s "$H37BAM" ]] || { echo "FATAL: ${S}: no P1 alignment at ${H37BAM}" \
        "(bin/archive_alignments.sh --restore brings back an archived one)" >&2; exit 1; }
else
    echo "NOTE: P1WORK unset; the junction route reads ${CRAMROOT}/${REL}"
fi
mkdir -p "$OUTDIR"
# --index-dir: the per-catalogue index and unique positions, built by the
# first task and shared by the rest (keyed on the catalogue, H37Rv and the
# presence code version)
"$MTB_PY_VT" accessory/bin/locus_presence.py \
  --sample "$S" --reference "$REF" \
  --catalogue "${CAT}.tsv" --fasta "${CAT}.fasta" \
  --cram "${CRAMROOT}/${REL}" --h37rv "$MTB_H37RV" \
  ${H37BAM:+--h37rv-bam "$H37BAM"} --index-dir "${OUTDIR}/.presence_index" \
  --out "${OUTDIR}/${S}.presence.tsv"
