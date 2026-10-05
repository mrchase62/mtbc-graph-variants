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
#   vcfbub -l 0 -a MAX     pop nested bubbles, keep top-level alleles up to MAX
#   vcfwave -I 64          realign each ALT and decompose complex bubbles into
#                          canonical SNPs / indels / SVs, keeping genotypes
#   bcftools norm/sort     canonicalise, index
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

[[ $# -ge 1 ]] || { mtb_usage "${BASH_SOURCE[0]}"; exit 2; }

GRAPH_NAME="$1"
GRAPH_DIR="${MTB_GRAPHS}/${GRAPH_NAME}"; shift
MAX_ALLELE=100000     # vcfbub -a: alleles longer than this stay collapsed
# vcfwave -I is --inv-min: the MINIMUM allele length considered for inverted
# alignment, not a realignment band. There is no band -- vcfwave's length cap is
# -L/--max-length, which defaults to unlimited and is not passed here. Raising
# -I above vcfwave's own default of 64 means inversions shorter than it are never
# inverted-aligned and get decomposed into SNPs and indels instead.
WAVE_I=64
CHUNK=5000            # records per vcfwave chunk
JOBS="${SLURM_CPUS_PER_TASK:-8}"
IN_VCF=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-allele) MAX_ALLELE="$2"; shift 2 ;;
        --wave-i)     WAVE_I="$2";     shift 2 ;;
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
rm -f "${WORKDIR}/chunks/"part_*
grep -v '^#' "$BUB" | split -l "$CHUNK" -d -a 4 - "${WORKDIR}/chunks/part_"
NCHUNK=$(ls "${WORKDIR}/chunks/"part_* | wc -l)
echo "    ${NCHUNK} chunks"

for p in "${WORKDIR}/chunks/"part_*; do
    cat "${WORKDIR}/header.txt" "$p" > "${p}.vcf"
    rm -f "$p"
done

echo "### [3/4] vcfwave -I ${WAVE_I} across ${JOBS} parallel jobs"
BINDS="$(mtb_bind)"
export BINDS WAVE_I MTB_PGGB_SIF
wave_one() {
    local v="$1"
    singularity exec $BINDS "$MTB_PGGB_SIF" \
        vcfwave -I "$WAVE_I" "$v" > "${v%.vcf}.wave.vcf" 2> "${v%.vcf}.wave.err" \
        && echo "    done $(basename "$v")" \
        || { echo "    FAILED $(basename "$v")" >&2; return 1; }
}
export -f wave_one
ls "${WORKDIR}/chunks/"part_*.vcf | xargs -P "$JOBS" -I{} bash -c 'wave_one "$@"' _ {}

echo "### [4/4] merging, normalising, indexing"
OUT="${GRAPH_DIR}/all_variants.decomposed.vcf.gz"

# Take the header from a WAVED chunk, not from the vcfbub input: vcfwave adds
# its own ##INFO declarations (ORIGIN, LEN, TYPE, INV, CONFLICT) and reusing the
# pre-wave header leaves them undeclared, which bcftools rejects outright.
FIRST_WAVE=$(ls "${WORKDIR}/chunks/"part_*.wave.vcf | head -1)
[[ -s "$FIRST_WAVE" ]] || { echo "no vcfwave output produced" >&2; exit 1; }

{
    grep '^#' "$FIRST_WAVE"
    for w in "${WORKDIR}/chunks/"part_*.wave.vcf; do grep -v '^#' "$w" || true; done
} | mtb_bcftools sort -m 8G -T "${WORKDIR}/sort_tmp" -Oz -o "$OUT" -
mtb_bcftools index -t "$OUT"


# vcfwave decomposes each ALT allele independently and does not merge identical
# results across alleles of the same bubble, so one SNP can appear once per
# haplotype with the carriers split across records. Collapse them and refill
# AC/AN/AF. See bin/vcf_collapse.sh for the full explanation.
echo "### [5/5] collapsing duplicate alleles and refilling AC/AN/AF"
COLLAPSED="${GRAPH_DIR}/all_variants.collapsed.vcf.gz"
mtb_bcftools norm -m +any "$OUT" -Ou \
  | mtb_bcftools norm -m -any -Ou \
  | mtb_bcftools +fill-tags -Oz -o "$COLLAPSED" -- -t AC,AN,AF
mtb_bcftools index -t "$COLLAPSED"
echo "    $(mtb_bcftools index -n "$OUT") -> $(mtb_bcftools index -n "$COLLAPSED") records"

echo
echo "### wrote ${OUT}"
echo "### wrote ${COLLAPSED}   <- use this one downstream"
mtb_bcftools index -n "$COLLAPSED" | sed 's/^/    records (collapsed): /'
echo "    samples: $(mtb_bcftools query -l "$OUT" | wc -l)"
echo
echo "Split into size classes with:"
echo "  bin/vcf_split_classes.sh ${GRAPH_NAME} --in ${COLLAPSED}"
echo "REMEMBER: bin/sync_back.sh"
