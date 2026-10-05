#!/usr/bin/env bash
# Collapse duplicate allele records left behind by vcfwave, and refill AC/AN/AF.
#
#   bin/vcf_collapse.sh <graph-dir-name> [--in FILE] [--out FILE]
#
# WHY THIS IS NEEDED
#
# vcfwave decomposes each ALT allele of a graph bubble independently and does
# not merge identical decomposed variants across alleles of the same bubble.
# So a SNP carried by several haplotypes of one bubble is emitted once per
# haplotype -- same CHROM/POS/REF/ALT, same ORIGIN, but the carriers split
# across the records.
#
# Example, position 1849 of mtb.complex: six records, all C->A, all
# ORIGIN=...:1592, with AC 24/1/78/1/2/1. The carrier sets are disjoint and sum
# to 107. Reading AC off any single record understates the allele by 4x.
#
# Genome-wide this inflated the SNP file from 66,782 real sites to 168,864
# records. Any tree, GWAS or PastML run taken straight off those records would
# be wrong.
#
# THE FIX
#
#   norm -m +any    merge all records at a position into one multiallelic
#                   record; identical ALTs collapse and genotypes are unioned
#   norm -m -any    split back to biallelic, one record per DISTINCT allele
#   +fill-tags      recompute AC/AN/AF, which the round-trip leaves stale
#
# Verified: at POS 1849 this yields one C->A record with AC=107 AN=416, and at
# POS 3928274 it keeps C->A, C->G and C->N as three records with carriers intact.
#
# `norm -d exact` is NOT a substitute -- it keeps the first record and discards
# the others' genotypes (24 of 107 carriers at POS 1849).
set -euo pipefail

# --- locate config/project_env.sh -----------------------------------------
_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/../config/project_env.sh" \
          "${MTB_WORK:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: cannot locate config/project_env.sh" >&2; exit 1; }
source "$_mtb_env"

[[ $# -ge 1 ]] || { mtb_usage "${BASH_SOURCE[0]}"; exit 2; }
GRAPH_NAME="$1"
GRAPH_DIR="${MTB_GRAPHS}/${GRAPH_NAME}"; shift

IN="${GRAPH_DIR}/all_variants.decomposed.vcf.gz"
OUT="${GRAPH_DIR}/all_variants.collapsed.vcf.gz"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --in)  IN="$2";  shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
mtb_require_file "$IN"

echo "### collapsing duplicate alleles"
echo "    in  : $IN  ($(mtb_bcftools index -n "$IN" 2>/dev/null || echo '?') records)"

mtb_bcftools norm -m +any "$IN" -Ou \
  | mtb_bcftools norm -m -any -Ou \
  | mtb_bcftools +fill-tags -Oz -o "$OUT" -- -t AC,AN,AF
mtb_bcftools index -t "$OUT"

echo "    out : $OUT  ($(mtb_bcftools index -n "$OUT") records)"
echo
echo "Now re-split classes:  bin/vcf_split_classes.sh ${GRAPH_NAME} --in ${OUT}"
