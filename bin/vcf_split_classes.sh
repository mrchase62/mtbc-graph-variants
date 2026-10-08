#!/usr/bin/env bash
# Split the collapsed graph VCF into SNP / indel / SV views that share one
# coordinate system and one set of genotypes.
#
#   bin/vcf_split_classes.sh <graph-dir-name> [--sv-min 50] [--in FILE] [--outdir DIR]
#
# --outdir writes the four files somewhere other than the graph's directory,
# so a rebuild on an existing graph does not replace the files beside it.
#
# Because all three come from the same graph VCF, an SV and a SNP in the
# same genome are directly comparable and can go into the same downstream
# analysis (tree, PastML, GWAS). That was not true of the 2025 outputs, where
# SVs lived in a genotype-free TSV.
set -euo pipefail
# --- locate config/project_env.sh -----------------------------------------
# Under sbatch, BASH_SOURCE[0] is Slurm's spool copy of this script rather than
# the file in bin/, so the relative lookup alone is not enough. Fall back to the
# submit directory and then to an already-exported MTB_WORK.
_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/../config/project_env.sh" \
          "${MTB_WORK:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || {
    echo "FATAL: cannot locate config/project_env.sh." >&2
    echo "       Submit from the project root, or export MTB_ENV_FILE." >&2
    exit 1
}
source "$_mtb_env"

[[ $# -ge 1 ]] || { mtb_usage "${BASH_SOURCE[0]}"; exit 2; }
GRAPH_NAME="$1"
GRAPH_DIR="${MTB_GRAPHS}/${GRAPH_NAME}"; shift
SV_MIN=50
# the collapsed product, never the decomposed intermediate, whose duplicate
# keys split each allele's carriers (audit PGB-9)
IN="${GRAPH_DIR}/all_variants.collapsed.vcf.gz"
OUTDIR="$GRAPH_DIR"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sv-min) SV_MIN="$2"; shift 2 ;;
        --in)     IN="$2";     shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
mtb_require_file "$IN"
IN="$(readlink -f "$IN")"     # absolute before the cd below
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# What counts as "large" needs care. Sizing a record by the LENGTH DIFFERENCE
# alone -- abs(strlen(REF)-strlen(ALT)) -- misses every length-preserving change:
# a 1,379 bp block substitution has a difference of zero and was silently
# classified as a small variant, alongside blocks up to 13 kb. Those are divergent
# haplotypes that vcfwave could not decompose, overwhelmingly in PE_PGRS, and
# snpEff labels them "missense_variant" or "frameshift_variant" as if they were
# point changes.
#
# So a record is large when EITHER side is big, not only when the two differ:
#   difference >= SV_MIN                      (indels and unbalanced events)
#   OR both REF and ALT are >= SV_MIN         (block substitutions, MNPs)
# The second clause is min(strlen(REF),strlen(ALT)) >= SV_MIN written without a
# min() function, which the bcftools expression language does not provide.
# bcftools cannot parse a negated compound expression, so SMALL is written out
# as the De Morgan complement of BIG rather than as !(BIG).
BIG="abs(strlen(REF)-strlen(ALT)) >= ${SV_MIN} || (strlen(REF) >= ${SV_MIN} && strlen(ALT) >= ${SV_MIN})"
SMALL="abs(strlen(REF)-strlen(ALT)) < ${SV_MIN} && (strlen(REF) < ${SV_MIN} || strlen(ALT) < ${SV_MIN})"

echo "### SNPs and MNPs"
mtb_bcftools view -i "(TYPE=\"snp\" || TYPE=\"mnp\") && ${SMALL}" "$IN" -Oz -o snps.vcf.gz
mtb_bcftools index -t snps.vcf.gz

echo "### indels < ${SV_MIN} bp"
mtb_bcftools view -i "TYPE=\"indel\" && ${SMALL}" "$IN" -Oz -o indels.vcf.gz
mtb_bcftools index -t indels.vcf.gz

echo "### structural variants >= ${SV_MIN} bp (length change OR block substitution)"
mtb_bcftools view -i "${BIG}" "$IN" -Oz -o svs.vcf.gz
mtb_bcftools index -t svs.vcf.gz

echo "### small variants (the snpEff input): everything not large"
mtb_bcftools view -i "${SMALL}" "$IN" -Oz -o small_variants.vcf.gz
mtb_bcftools index -t small_variants.vcf.gz

printf '\n%-14s %10s\n' "class" "records"
for f in snps indels svs small_variants; do
    printf '%-14s %10s\n' "$f" "$(mtb_bcftools index -n ${f}.vcf.gz)"
done
echo
echo "All three share the coordinate system of ${MTB_REF_PATH}"
echo "and carry genotypes for $(mtb_bcftools query -l "$IN" | wc -l) samples."
