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
mkdir -p "$OUTDIR"
"$MTB_PY_VT" accessory/bin/locus_presence.py \
  --sample "$S" --reference "$REF" \
  --catalogue "${CAT}.tsv" --fasta "${CAT}.fasta" \
  --cram "${CRAMROOT}/${REL}" --h37rv "$MTB_H37RV" \
  --out "${OUTDIR}/${S}.presence.tsv"
