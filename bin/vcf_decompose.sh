#!/usr/bin/env bash
#SBATCH --job-name=vcf_decompose
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 0-12:00
#SBATCH -p sapphire
#SBATCH --mem=64G
#SBATCH --output=slurm/decompose_%j.out
#SBATCH --error=slurm/decompose_%j.err
#
# Produce ONE VCF containing SNPs, indels AND structural variants, genotyped
# across every path in the graph.
#
#   sbatch bin/vcf_decompose.sh <graph-dir-name> [--max-allele N] [--chunk N]
#
# This is the step the 2025 pipeline was missing. It ran `vg deconstruct` and
# then filtered to `LV=0 && strlen < 200`, which discards every structural
# variant and every nested allele, and called SVs separately from odgi untangle
# into a TSV with no genotypes. The two tracks could never be joined.
#
# Chain:
#   vg deconstruct -a      nested VCF with LV/PS snarl tags   (already run)
#   vcfbub -l 0 -a MAX     pop nested bubbles; a site with an allele longer
#                          than MAX is REMOVED and its children promoted
#   vcfwave -I 1000        realign each ALT and decompose complex bubbles into
#                          SNPs / indels / SVs, keeping genotypes (records with
#                          an allele over 5 kb always at -I 1000; see below)
#   bcftools sort          -> all_variants.decomposed.vcf.gz (intermediate)
#   bin/vcf_collapse.sh    trim, one record per allele, union genotypes
#                          -> all_variants.collapsed.vcf.gz (THE PRODUCT)
#
# all_variants.collapsed.vcf.gz is what every downstream reader should take.
# The decomposed file holds one record per vcfwave allele: the same
# (pos,ref,alt) recurs with its carriers split across records, and a reader
# that does not union them reads wrong genotypes (audit Fault A / GRAPHVCF-1).
#
# Measured on 20,000 records of the 2025 mtb.complex VCF: vcfwave resolves
# 270 SV records into 2,307, and 17,585 SNPs into 36,687, because complex
# bubbles get unpacked instead of being dropped by a length filter.
#
# vcfwave is the bottleneck per-record (~28 min per 14.5k records at 417 samples,
# and it does not use more than one core), so the input is chunked and run in
# parallel: 47k records finish in under 10 min across 16 cores.
#
# The final `bcftools sort` is then the long pole -- ~1 GB of decomposed VCF with
# 416 samples. It is given 8 GB of sort memory and a temp directory on scratch
# rather than the default 768 MB spilling into /tmp.
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
# resolved now, before the cd into the graph directory below
COLLAPSE_SH="$(cd "$(dirname "$_mtb_env")/.." && pwd)/bin/vcf_collapse.sh"

[[ $# -ge 1 ]] || { mtb_usage "${BASH_SOURCE[0]}"; exit 2; }

GRAPH_NAME="$1"
GRAPH_DIR="${MTB_GRAPHS}/${GRAPH_NAME}"; shift
# vcfbub -a: a site with an allele longer than this is dropped, and its child
# snarls are promoted in its place (keeping their LV). In CX333 that removed one
# 3.4 Mb top-level snarl and promoted 49,120 LV=1 records (audit PGB-11).
MAX_ALLELE=100000
# vcfwave -I is --inv-min: the MINIMUM allele length considered for inverted
# alignment, not a realignment band. There is no band -- vcfwave's length cap is
# -L/--max-length, which defaults to unlimited and is not passed here. Raising
# -I above vcfwave's own default of 64 means inversions shorter than it are never
# inverted-aligned and get decomposed into SNPs and indels instead.
#
# 1000 is what produced the CX333 production VCF (slurm/analyze_an_CX333_45740274,
# 2026-09-11). The default was later changed to 64 here, so a rebuild from this
# script would not have reproduced production (audit PGB-10). Lowering it is a
# choice to make per build with --wave-i, not a silent default.
WAVE_I=1000
# Inverted alignment is quadratic. At a low --inv-min it wedges on a
# multi-kilobase allele (a 23-minute stall on one record at -I 64; see
# bin/vcfwave_rerun.sh). Records whose longest allele exceeds MAX_INV_ALLELE are
# therefore always waved at FALLBACK_I: every record is still decomposed, and
# only the ones that could hang are excluded from inversion detection. At the
# default -I 1000 both arms run alike.
MAX_INV_ALLELE=5000
FALLBACK_I=1000
CHUNK=5000            # records per vcfwave chunk
JOBS="${SLURM_CPUS_PER_TASK:-8}"
IN_VCF=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-allele) MAX_ALLELE="$2"; shift 2 ;;
        --wave-i)     WAVE_I="$2";     shift 2 ;;
        --max-inv-allele) MAX_INV_ALLELE="$2"; shift 2 ;;
        --chunk)      CHUNK="$2";      shift 2 ;;
        --jobs)       JOBS="$2";       shift 2 ;;
        --input)      IN_VCF="$2";     shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ -d "$GRAPH_DIR" ]] || { echo "no such graph dir: $GRAPH_DIR" >&2; exit 1; }
if [[ -z "${IN_VCF:-}" ]]; then
    # prefer the bgzipped form; see pilot_analyze.sh for why both exist
    if   [[ -s "${GRAPH_DIR}/variants.vcf.gz" ]]; then IN_VCF="${GRAPH_DIR}/variants.vcf.gz"
    else IN_VCF="${GRAPH_DIR}/variants.vcf"; fi
fi
mtb_require_file "$MTB_PGGB_SIF" "$IN_VCF"

WORKDIR="${GRAPH_DIR}/decompose"
mkdir -p "$WORKDIR/chunks" "$WORKDIR/sort_tmp"
cd "$GRAPH_DIR"

BUB="${WORKDIR}/bub.vcf"

echo "### [1/4] vcfbub -l 0 -a ${MAX_ALLELE}"
singularity exec $(mtb_bind) "$MTB_PGGB_SIF" \
    vcfbub -l 0 -a "$MAX_ALLELE" -i "$IN_VCF" > "$BUB"
echo "    $(grep -vc '^#' "$BUB") records (from $(grep -vc '^#' "$IN_VCF"))"

echo "### [2/4] splitting into ${CHUNK}-record chunks for vcfwave"
grep '^#' "$BUB" > "${WORKDIR}/header.txt"
rm -f "${WORKDIR}/chunks/"part_* "${WORKDIR}/chunks/"big_* "${WORKDIR}/chunks/"*.body
# Partition by longest allele (as bin/vcfwave_rerun.sh). Each record lands in
# exactly one arm, so the union is the whole input with no duplicates and no gaps.
grep -v '^#' "$BUB" | awk -F'\t' -v cap="$MAX_INV_ALLELE" \
    -v small="${WORKDIR}/chunks/small.body" -v big="${WORKDIR}/chunks/big.body" '
    { m = length($4); n = split($5, A, ",");
      for (i = 1; i <= n; i++) if (length(A[i]) > m) m = length(A[i]);
      if (m > cap) print > big; else print > small }'
touch "${WORKDIR}/chunks/small.body" "${WORKDIR}/chunks/big.body"
echo "    $(wc -l < "${WORKDIR}/chunks/small.body") records <= ${MAX_INV_ALLELE} bp (--inv-min ${WAVE_I})"
echo "    $(wc -l < "${WORKDIR}/chunks/big.body") records >  ${MAX_INV_ALLELE} bp (--inv-min ${FALLBACK_I})"
split -l "$CHUNK" -d -a 4 "${WORKDIR}/chunks/small.body" "${WORKDIR}/chunks/part_"
split -l "$CHUNK" -d -a 4 "${WORKDIR}/chunks/big.body"   "${WORKDIR}/chunks/big_"
rm -f "${WORKDIR}/chunks/"*.body
# nullglob: either arm can be empty, and an unmatched glob must not become a
# literal filename (or, through ls under pipefail, a fatal error)
shopt -s nullglob
for p in "${WORKDIR}/chunks/"part_* "${WORKDIR}/chunks/"big_*; do
    cat "${WORKDIR}/header.txt" "$p" > "${p}.vcf"
    rm -f "$p"
done
CHUNK_VCFS=("${WORKDIR}/chunks/"part_*.vcf "${WORKDIR}/chunks/"big_*.vcf)
shopt -u nullglob
echo "    ${#CHUNK_VCFS[@]} chunks"

echo "### [3/4] vcfwave -I ${WAVE_I} (over ${MAX_INV_ALLELE} bp: -I ${FALLBACK_I}) across ${JOBS} parallel jobs"
BINDS="$(mtb_bind)"
export BINDS WAVE_I FALLBACK_I MTB_PGGB_SIF
wave_one() {
    local v="$1" iv="$WAVE_I"
    case "$(basename "$v")" in big_*) iv="$FALLBACK_I" ;; esac
    singularity exec $BINDS "$MTB_PGGB_SIF" \
        vcfwave -I "$iv" "$v" > "${v%.vcf}.wave.vcf" 2> "${v%.vcf}.wave.err" \
        && echo "    done $(basename "$v") (--inv-min $iv)" \
        || { echo "    FAILED $(basename "$v")" >&2; return 1; }
}
export -f wave_one
printf '%s\n' "${CHUNK_VCFS[@]}" \
  | xargs -P "$JOBS" -I{} bash -c 'wave_one "$@"' _ {}

echo "### [4/4] merging, normalising, indexing"
OUT="${GRAPH_DIR}/all_variants.decomposed.vcf.gz"

# Take the header from a WAVED chunk, not from the vcfbub input: vcfwave adds
# its own ##INFO declarations (ORIGIN, LEN, TYPE, INV, CONFLICT) and reusing the
# pre-wave header leaves them undeclared, which bcftools rejects outright.
FIRST_WAVE=$(ls "${WORKDIR}/chunks/"*.wave.vcf | head -1)
[[ -s "$FIRST_WAVE" ]] || { echo "no vcfwave output produced" >&2; exit 1; }

# The settings go into the header, so a VCF says how it was decomposed.
{
    grep '^##' "$FIRST_WAVE"
    echo "##MTB_decompose=vcfbub -l 0 -a ${MAX_ALLELE}; vcfwave -I ${WAVE_I} (records with an allele over ${MAX_INV_ALLELE} bp: -I ${FALLBACK_I}); input $(basename "$IN_VCF")"
    grep '^#CHROM' "$FIRST_WAVE"
    for w in "${WORKDIR}/chunks/"*.wave.vcf; do grep -v '^#' "$w" || true; done
} | mtb_bcftools sort -m 8G -T "${WORKDIR}/sort_tmp" -Oz -o "$OUT" -
mtb_bcftools index -t "$OUT"


# vcfwave decomposes each ALT allele independently and does not merge identical
# results across alleles of the same bubble, so one SNP can appear once per
# haplotype with the carriers split across records. Collapse them and refill
# AC/AN/AF. One implementation, in bin/vcf_collapse.sh (which explains it), so
# that the step run here and the step run alone cannot drift apart.
echo "### [5/5] collapsing duplicate alleles and refilling AC/AN/AF"
COLLAPSED="${GRAPH_DIR}/all_variants.collapsed.vcf.gz"
bash "$COLLAPSE_SH" "$GRAPH_NAME" \
    --in "$OUT" --out "$COLLAPSED"
echo "    $(mtb_bcftools index -n "$OUT") -> $(mtb_bcftools index -n "$COLLAPSED") records"

echo
echo "### wrote ${OUT}   (intermediate: duplicate keys, carriers split)"
echo "### wrote ${COLLAPSED}   <- THE PRODUCT: every downstream reader takes this one"
mtb_bcftools index -n "$COLLAPSED" | sed 's/^/    records (collapsed): /'
echo "    samples: $(mtb_bcftools query -l "$OUT" | wc -l)"
echo
echo "Split into size classes with:"
echo "  bin/vcf_split_classes.sh ${GRAPH_NAME} --in ${COLLAPSED}"
echo "REMEMBER: bin/sync_back.sh"
