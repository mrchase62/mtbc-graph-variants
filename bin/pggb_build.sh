#!/usr/bin/env bash
#SBATCH --job-name=pggb
#SBATCH -N 1
#SBATCH -c 112
#SBATCH -t 3-00:00
#SBATCH -p sapphire
#SBATCH --mem=900G
#SBATCH --output=slurm/pggb_%j.out
#SBATCH --error=slurm/pggb_%j.err
#
# RESOURCES: sapphire and intermediate both give 112 cores and ~1 TB per node.
# Earlier runs in this project asked for 48 cores and 128 GB, leaving well over
# half of each node idle -- and wfmash alignment, which is the dominant stage,
# parallelises across threads. Peak memory measured at 205 paths was 42-46 GB, so
# 900 GB is generous headroom rather than a real requirement.
#
# For a build that might exceed 3 days, submit to the intermediate partition
# instead, which allows 14 days (12 nodes rather than 184, so expect to queue):
#
#   sbatch -p intermediate -t 14-00:00 bin/pggb_build.sh <fasta> <name> ...
#
# Build a pangenome graph with pggb.
#
#   sbatch bin/pggb_build.sh <input.fasta.gz> <output-name> [pggb args...]
#
# Reproduces the 2025 mtb.complex settings by default (-s 10000 -p 95 -k 51 -K 21).
# Output lands in $MTB_GRAPHS/<output-name>/.
#
# Examples:
#   sbatch bin/pggb_build.sh mtb.complex.fasta.gz mtb.complex.p95.s10k.k51
#   sbatch bin/pggb_build.sh mtb_lineage4.fasta.gz L4.p99.s10k.k61 -p 99 -k 61
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

if [[ $# -lt 2 ]]; then
    mtb_usage "${BASH_SOURCE[0]}"; exit 2
fi

IN_NAME="$1"; OUT_NAME="$2"; shift 2

IN_HOST="${MTB_FASTAS}/${IN_NAME}"
OUT_HOST="${MTB_GRAPHS}/${OUT_NAME}"

mtb_require_work
mtb_require_file "$MTB_PGGB_SIF" "$IN_HOST"
mkdir -p "$OUT_HOST"

# --- provenance sidecar ------------------------------------------------------
# WHO BUILT THIS GRAPH, AND FROM WHAT CODE. The graph's filename carries the
# pggb component fingerprint, which says what pggb did. It does not say which
# revision of THIS repository drove it, and once panel construction and variant
# calling live in separate repositories that information is not recoverable from
# anywhere else: the downstream repository's history contains no trace of the
# code that built its graph.
#
# So it is written beside the graph, at build time, by the repository that is
# doing the building. p0_prepare.sh reads it into build_info.tsv and
# stamp_build_id.sh puts it in every VCF header, which closes the chain from a
# called variant back to the commit that produced the coordinates it sits in.
#
# Recorded honestly rather than optimistically: a dirty working tree is stamped
# as dirty, because a commit id alone would then be a claim the tree does not
# support.
_prov="${OUT_HOST}/graph_provenance.tsv"
{
    _repo="$(git -C "$(dirname "${BASH_SOURCE[0]}")" config --get remote.origin.url 2>/dev/null || true)"
    _commit="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse HEAD 2>/dev/null || true)"
    _dirty="clean"
    git -C "$(dirname "${BASH_SOURCE[0]}")" diff --quiet HEAD 2>/dev/null || _dirty="dirty"
    [[ -n "$_commit" ]] || { _repo="not-a-git-checkout"; _commit="unknown"; _dirty="unknown"; }
    printf 'panel_repo\t%s\n'   "${_repo:-no-remote-configured}"
    printf 'panel_commit\t%s\n' "$_commit"
    printf 'panel_dirty\t%s\n'  "$_dirty"
    printf 'pggb_container\t%s\n' "$MTB_PGGB_SIF"
    printf 'pggb_args\t%s\n'   "$*"
    printf 'input_fasta\t%s\n' "$(basename "$IN_HOST")"
    printf 'built\t%s\n'       "$(date -Is)"
    printf 'built_by\t%s\n'    "${USER:-unknown}"
} > "$_prov"
echo "[pggb_build] provenance: ${_prov}"

IN_C="$(mtb_in_data "$IN_HOST")"
OUT_C="$(mtb_in_graphs "$OUT_HOST")"

# Number of haplotypes = number of sequences in the input.
if [[ -e "${IN_HOST}.fai" ]]; then
    N_HAP=$(wc -l < "${IN_HOST}.fai")
else
    echo "[pggb_build] no ${IN_HOST}.fai — counting sequences (slower)"
    N_HAP=$(zcat -f "$IN_HOST" | grep -c '^>')
fi
echo "[pggb_build] haplotypes: ${N_HAP}"

# --- parameters, stated explicitly ---------------------------------------
# Every value below is named so that a rerun's settings are visible in this file
# rather than half-inherited from pggb's defaults. Verified against the recovered
# 2025 params.yml files -- see GRAPH_PROVENANCE.md for the audit.
#
#   -s 10000  segment length. pggb default is 1000; 10 kb is 10x that. Long
#             segments cannot anchor inside hypervariable PE_PGRS loci, which is
#             the mechanism behind the 1,632 unresolved block substitutions in
#             the complex graph. Lower this to resolve those regions.
#   -l 30000  minimum mapping block length. pggb default is 5*s = 50000. The 2025
#             builds used 30000, so it MUST be set explicitly: omitting it
#             silently changes the graph.
#   -p 95     mapping identity. Default 90. MTBC is >99.9% identical outside
#             PE/PPE, so 95 is comfortable even including M. canettii.
#   -k 51     seqwish exact-match filter. Default 23. At 51 short shared matches
#             in variable regions are discarded, compounding the -s effect.
#   -K 21     wfmash mash kmer. Default 15. NOTE: only the 2025 complex graph
#             used 21; mtb_lineage4 and pggb_mtb_1_2_4_bovis used the default 15.
#             The three existing graphs are therefore NOT parameter-consistent.
SEGMENT_LENGTH=10000
BLOCK_LENGTH=30000
MAP_PCT_ID=95
MIN_MATCH_LEN=51
MASH_KMER=21

# Anything passed after <output-name> is appended and wins, since pggb takes the
# last occurrence of a flag.
PGGB_ARGS=(-s "$SEGMENT_LENGTH" -l "$BLOCK_LENGTH" -p "$MAP_PCT_ID"
           -k "$MIN_MATCH_LEN" -K "$MASH_KMER" --keep-temp-files)

# pggb -r is --resume: "do not overwrite existing outputs in the given
# directory". All three 2025 builds ran with resume:true, which means their final
# graphs may incorporate intermediates from earlier attempts under other
# parameters -- unverifiable after the fact. Resume is therefore OPT-IN here, via
# MTB_PGGB_RESUME=1, instead of being passed silently on every run.
if [[ "${MTB_PGGB_RESUME:-0}" == "1" ]]; then
    echo "[pggb_build] WARNING: --resume requested; existing intermediates in"
    echo "             ${OUT_HOST} will be reused and NOT rebuilt."
    PGGB_ARGS+=(-r)
elif [[ -n "$(ls -A "$OUT_HOST" 2>/dev/null)" ]]; then
    echo "[pggb_build] FATAL: ${OUT_HOST} is not empty and resume is off." >&2
    echo "             Move it aside, or set MTB_PGGB_RESUME=1 to reuse it." >&2
    exit 1
fi

set -x
singularity exec --cleanenv $(mtb_bind) "$MTB_PGGB_SIF" \
    bash /usr/local/bin/pggb \
        -i "$IN_C" \
        -o "$OUT_C" \
        -n "$N_HAP" \
        -t "${SLURM_CPUS_PER_TASK:-112}" \
        "${PGGB_ARGS[@]}" \
        "$@"
set +x

echo "[pggb_build] finished -> ${OUT_HOST}"
echo "[pggb_build] REMEMBER: bin/sync_back.sh"
