#!/usr/bin/env bash
# Refuse a merged VCF that standard tools cannot read.
#
#   bash bin/vcf_gate.sh <merged.vcf.gz> <H37Rv.fasta> [h37rv_contig]
#
# WHY. Two spec violations sat in every merged VCF this pipeline wrote, and both
# were invisible to the pipeline's own readers: GT=2 on records carrying one ALT
# allele (976 records in scale200; `bcftools +fill-tags` aborts on the first)
# and REF=N on symbolic records (7,291 records that do not match H37Rv). Each
# check below would have stopped the run that produced them.
#
#   1. every GT allele index is at most the record's ALT count
#   2. `bcftools +fill-tags` reads every record (catches 1 and other GT errors)
#   3. `bcftools norm -c e` finds REF equal to H37Rv on the H37Rv contig
#      (node_* contigs carry no sequence here and are not checked)
#   4. records are sorted within each contig
#
# Exit 0 when all pass; otherwise prints what failed and exits 1.
set -euo pipefail

_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: no project_env.sh" >&2; exit 1; }
source "$_mtb_env"

VCF="${1:?usage: vcf_gate.sh <merged.vcf.gz> <H37Rv.fasta> [contig]}"
FA="${2:?usage: vcf_gate.sh <merged.vcf.gz> <H37Rv.fasta> [contig]}"
CTG="${3:-NC_000962.3}"
BCF="${MTB_BCFTOOLS:?MTB_BCFTOOLS is unset}"
[[ -s "$VCF" && -s "${VCF}.tbi" ]] || { echo "FATAL: ${VCF} or its index is missing" >&2; exit 1; }
[[ -s "$FA" ]] || { echo "FATAL: no FASTA ${FA}" >&2; exit 1; }

rc=0
echo "[vcf_gate] ${VCF}"

# 1 -- GT index vs ALT count
bad_gt=$(zcat "$VCF" | awk -F'\t' '!/^#/{
    n = split($5, a, ",")
    for (i = 10; i <= NF; i++) { split($i, g, ":"); if (g[1] != "." && g[1] + 0 > n) { bad++; break } }
} END { print bad + 0 }')
if [[ "$bad_gt" -gt 0 ]]; then
    echo "  FAIL  ${bad_gt} records have a GT index above their ALT count"; rc=1
else
    echo "  ok    every GT index is within its record's ALT alleles"
fi

# 2 -- a standard reader gets through every record
if err=$("$BCF" +fill-tags "$VCF" -Ou -- -t AN,AC 2>&1 >/dev/null) && [[ -z "$err" ]]; then
    echo "  ok    bcftools +fill-tags reads every record"
else
    echo "  FAIL  bcftools +fill-tags: $(head -3 <<< "$err")"; rc=1
fi

# 3 -- REF matches the reference on the H37Rv contig
mm=$("$BCF" norm -c w -f "$FA" -r "$CTG" "$VCF" -Ou 2>&1 >/dev/null \
     | grep -c 'REF_MISMATCH\|does not match' || true)
if [[ "$mm" -gt 0 ]]; then
    echo "  FAIL  ${mm} records on ${CTG} have a REF that does not match ${FA}"; rc=1
else
    echo "  ok    every REF on ${CTG} matches the reference"
fi

# 4 -- sorted within contig
unsorted=$(zcat "$VCF" | awk -F'\t' '!/^#/{ if ($1 == c && $2 + 0 < p) u++; c = $1; p = $2 + 0 } END { print u + 0 }')
if [[ "$unsorted" -gt 0 ]]; then
    echo "  FAIL  ${unsorted} records are out of order within their contig"; rc=1
else
    echo "  ok    records sorted within each contig"
fi

[[ "$rc" -eq 0 ]] && echo "[vcf_gate] passed" || echo "[vcf_gate] FAILED" >&2
exit "$rc"
