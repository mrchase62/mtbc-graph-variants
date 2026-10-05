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
#   norm -m -any    split the few multiallelic records to biallelic
#   COLLAPSE_PY     trim each allele to its minimal REF/ALT (shared suffix, then
#                   shared prefix down to one anchor base), then write ONE
#                   record per (CHROM, POS, REF, ALT) with each sample's GT the
#                   union over the records sharing that key: 1 if any record
#                   says 1, else 0 if any says 0, else missing. Every other
#                   column comes from the first record with that key.
#   sort            trimming moves POS, so the stream is re-sorted
#   +fill-tags      recompute AC/AN/AF from the unioned genotypes
#
# The output has no duplicate keys and no padded SNPs; both are checked, and the
# run fails if either is found.
#
# THE PREVIOUS FIX, AND WHY IT WAS REPLACED (audit GRAPHVCF-5, 2026-10-05)
#
# This used `norm -m +any | norm -m -any`. The merge makes one record per POS
# whose REF is the longest REF there, and the split does not trim it again, so a
# SNP that shared its POS with a longer allele came out padded (CG>TG): 1,399
# records in CX333, 376 of them SNPs, which every reader requiring len(REF)==1
# or an exact (pos,ref,alt) join then lost. And a haploid merged record holds
# one allele per sample, so a genome carrying a SNP plus an indel anchored on
# the same base lost one of them (61 carrier cells read 0).
#
# `norm -d exact` is NOT a substitute -- it keeps the first record and discards
# the others' genotypes (24 of 107 carriers at POS 1849).
#
# THE PRODUCT
#
# all_variants.collapsed.vcf.gz is the graph's variant file for every
# downstream reader. all_variants.decomposed.vcf.gz is vcfwave's per-allele
# intermediate: the same key can appear in several records with the carriers
# split between them, and a reader that does not union them reads wrong
# genotypes (audit Fault A / GRAPHVCF-1). Not left-aligned: POS is vcfwave's,
# trimmed.
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

# Reads a plain-text, biallelic VCF on stdin and writes the collapsed records,
# unsorted, on stdout. tests/test_audit_selection.py runs this text directly.
# --- COLLAPSE_PY begin
read -r -d '' COLLAPSE_PY <<'PY' || true
import sys
def trim(pos, ref, alt):
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt
CODE = {".": 0, "0": 1, "1": 2}     # union order: missing < REF < ALT
GT = ".01"
recs, n_in, n_pad = {}, 0, 0
for line in sys.stdin:
    if line.startswith("##"):
        sys.stdout.write(line); continue
    if line.startswith("#"):
        sys.stdout.write("##MTB_collapse=alleles trimmed to minimal REF/ALT; one "
                         "record per CHROM/POS/REF/ALT; GT union (1 over 0 over .)\n")
        sys.stdout.write(line); continue
    f = line.rstrip("\n").split("\t")
    if "," in f[4]:
        sys.exit(f"FATAL: multiallelic record at {f[0]}:{f[1]}; split first")
    if f[8] != "GT":
        sys.exit(f"FATAL: FORMAT {f[8]} at {f[0]}:{f[1]}; only GT is handled")
    n_in += 1
    pos, ref, alt = trim(int(f[1]), f[3].upper(), f[4].upper())
    if (ref, alt) != (f[3].upper(), f[4].upper()):
        n_pad += 1
    try:
        g = bytearray(CODE[x] for x in f[9:])
    except KeyError as e:
        sys.exit(f"FATAL: genotype {e} at {f[0]}:{f[1]}; expected haploid 0/1/.")
    k = (f[0], pos, ref, alt)
    if k in recs:
        old = recs[k][1]
        for i, v in enumerate(g):
            if v > old[i]:
                old[i] = v
    else:
        recs[k] = ([f[0], str(pos), f[2], ref, alt] + f[5:9], g)
n_noalt = 0
for k, (fix, g) in recs.items():
    n_noalt += 2 not in g
    sys.stdout.write("\t".join(fix + [GT[v] for v in g]) + "\n")
snp_pad = sum(1 for (c, p, r, a) in recs if len(r) > 1 and len(a) > 1
              and (r[0] == a[0] or r[-1] == a[-1]))
print(f"    {n_in} records in, {len(recs)} distinct keys out; {n_pad} trimmed; "
      f"{n_noalt} with no ALT carrier", file=sys.stderr)
if snp_pad:
    sys.exit(f"FATAL: {snp_pad} padded records remain after trimming")
PY
# --- COLLAPSE_PY end

echo "### collapsing duplicate alleles"
echo "    in  : $IN  ($(mtb_bcftools index -n "$IN" 2>/dev/null || echo '?') records)"

SORT_TMP="$(mktemp -d "${OUT}.sort_tmp.XXXXXX")"
trap 'rm -rf "$SORT_TMP"' EXIT
mtb_bcftools norm -m -any "$IN" -Ov \
  | "$MTB_PY" -c "$COLLAPSE_PY" \
  | mtb_bcftools sort -m 2G -T "$SORT_TMP" -Ou - \
  | mtb_bcftools +fill-tags -Oz -o "$OUT" -- -t AC,AN,AF
mtb_bcftools index -f -t "$OUT"

# The product's contract: one record per key. Checked on the written file, so a
# change upstream of the trim (or in the sort) cannot break it silently.
NDUP="$(mtb_bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$OUT" | sort | uniq -d | wc -l)"
[[ "$NDUP" -eq 0 ]] || { echo "FATAL: ${NDUP} duplicate (CHROM,POS,REF,ALT) keys in ${OUT}" >&2; exit 1; }

echo "    out : $OUT  ($(mtb_bcftools index -n "$OUT") records, 0 duplicate keys)"
echo
echo "Now re-split classes:  bin/vcf_split_classes.sh ${GRAPH_NAME} --in ${OUT}"
