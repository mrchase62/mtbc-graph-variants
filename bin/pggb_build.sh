#!/usr/bin/env bash
#SBATCH --job-name=pggb
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -t 3-00:00
#SBATCH -p sapphire
#SBATCH --mem=900G
#SBATCH --output=slurm/pggb_%j.out
#SBATCH --error=slurm/pggb_%j.err
#
# RESOURCES: 48 cores, as CX333 was built (data/qc/cx333_build.cmd). Going to
# 112 measured 51% parallel efficiency -- 1.18x the speed for 2.33x the cores --
# so the extra cores mostly idle. Peak memory measured at 205 paths was 42-46
# GB, so 900 GB is generous headroom rather than a real requirement.
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
# <input.fasta.gz> is a file name under $MTB_FASTAS, not a path.
# Defaults are CX333's settings: -s 10000 -l 30000 -p 95 -k 23 -K 15, full
# all-against-all mapping, threads from SLURM_CPUS_PER_TASK (48 if unset).
# Sparse mapping is off unless asked for with -x (e.g. -x auto).
# Output lands in $MTB_GRAPHS/<output-name>/, which must be empty or absent.
# Every effective setting is written to <output-name>/graph_provenance.tsv.
#
# Examples:
#   sbatch bin/pggb_build.sh mtb.complex333.fasta.gz CX333.s10k.k23.K15
#   sbatch bin/pggb_build.sh panel.fasta.gz P.s10k.k23.K15.xauto -x auto
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
# The repository the config came from. Under sbatch dirname(BASH_SOURCE) is the
# spool directory, which is not a git checkout, so the provenance below asks
# this one.
_repo_dir="$(cd "$(dirname "$_mtb_env")/.." && pwd)"

if [[ $# -lt 2 ]]; then
    mtb_usage "${BASH_SOURCE[0]}"; exit 2
fi

IN_NAME="$1"; OUT_NAME="$2"; shift 2

IN_HOST="${MTB_FASTAS}/${IN_NAME}"
OUT_HOST="${MTB_GRAPHS}/${OUT_NAME}"

mtb_require_work
mtb_require_file "$MTB_PGGB_SIF" "$IN_HOST"

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
# rather than half-inherited from pggb's defaults. The defaults are the settings
# CX333 was built with (graphs/CX333.s10k.k23.K15/*.params.yml,
# data/qc/cx333_build.cmd) and that the 2026-10 test graphs reused
# (analysis/graph_tests/build_arm.sbatch). They used to be the 2025 mtb.complex
# settings, -k 51 -K 21, so following RUNBOOK -- which passes no arguments --
# gave a different graph under a directory named k23.K15 (audit PGB-2).
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
#   -k 23     seqwish exact-match filter (pggb's default). The 2025 complex
#             graph used 51, which discards short shared matches in variable
#             regions and compounds the -s effect.
#   -K 15     wfmash mash kmer (pggb's default). Only the 2025 complex graph used
#             21.
#   -x        wfmash sparse mapping: keep this fraction of mappings ('auto' for
#             the giant-component heuristic). OFF by default -- every pair of
#             genomes is mapped, as for CX333. On a 52-genome test it gave the
#             same graph within 1% and built 20% faster
#             (analysis/graph_tests/RESULTS.md), but adopting it is a decision
#             not yet taken, so it is only ever an explicit -x.
SEGMENT_LENGTH=10000
BLOCK_LENGTH=30000
MAP_PCT_ID=95
MIN_MATCH_LEN=23
MASH_KMER=15
SPARSE_MAP=""                      # empty = not passed (pggb default 1.0)
THREADS="${SLURM_CPUS_PER_TASK:-48}"

# Arguments after <output-name> override the named settings above. They are
# folded in here rather than appended after them, so that the values pggb runs
# with are the values recorded in the provenance file: appended, pggb took the
# last occurrence and the command line carried both (CX333's log shows
# "-k 51 -K 21 ... -k 23 -K 15"). Anything not named here is passed through
# unchanged, in order, and recorded verbatim.
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--segment-length) SEGMENT_LENGTH="$2"; shift 2 ;;
        -l|--block-length)   BLOCK_LENGTH="$2";   shift 2 ;;
        -p|--map-pct-id)     MAP_PCT_ID="$2";     shift 2 ;;
        -k|--min-match-len)  MIN_MATCH_LEN="$2";  shift 2 ;;
        -K|--mash-kmer)      MASH_KMER="$2";      shift 2 ;;
        -x|--sparse-map)     SPARSE_MAP="$2";     shift 2 ;;
        -t|--threads)        THREADS="$2";        shift 2 ;;
        -n|--n-haplotypes)
            [[ "$2" == "$N_HAP" ]] || echo "[pggb_build] NOTE: -n $2 given; the input has ${N_HAP} sequences"
            N_HAP="$2"; shift 2 ;;
        -r|--resume)
            echo "[pggb_build] FATAL: pass resume as MTB_PGGB_RESUME=1, not $1." >&2; exit 2 ;;
        -i|--input-fasta|-o|--output-dir)
            echo "[pggb_build] FATAL: $1 is set from the positional arguments." >&2; exit 2 ;;
        *) EXTRA_ARGS+=("$1"); shift ;;
    esac
done

# A DIRECTORY NAME THAT STATES SETTINGS MUST STATE THESE ONES. Graph directories
# are named for their parameters (CX333.s10k.k23.K15), and every reader takes
# the name at its word. Each dot-separated token of the form s<N>k, l<N>k, p<N>,
# k<N> or K<N> is checked against the value pggb will run with.
_bad=()
IFS='.' read -r -a _tok <<< "$OUT_NAME"
for _t in "${_tok[@]}"; do
    case "$_t" in
        s[0-9]*k) [[ "${_t:1:-1}000" == "$SEGMENT_LENGTH" ]] || _bad+=("$_t vs -s $SEGMENT_LENGTH") ;;
        l[0-9]*k) [[ "${_t:1:-1}000" == "$BLOCK_LENGTH" ]]   || _bad+=("$_t vs -l $BLOCK_LENGTH") ;;
        p[0-9]*)  [[ "${_t:1}" =~ ^[0-9]+$ ]] && { [[ "${_t:1}" == "$MAP_PCT_ID" ]] || _bad+=("$_t vs -p $MAP_PCT_ID"); } ;;
        k[0-9]*)  [[ "${_t:1}" =~ ^[0-9]+$ ]] && { [[ "${_t:1}" == "$MIN_MATCH_LEN" ]] || _bad+=("$_t vs -k $MIN_MATCH_LEN"); } ;;
        K[0-9]*)  [[ "${_t:1}" =~ ^[0-9]+$ ]] && { [[ "${_t:1}" == "$MASH_KMER" ]] || _bad+=("$_t vs -K $MASH_KMER"); } ;;
    esac
done
if (( ${#_bad[@]} )); then
    echo "[pggb_build] FATAL: the output name ${OUT_NAME} disagrees with the settings:" >&2
    printf '             %s\n' "${_bad[@]}" >&2
    exit 1
fi

PGGB_ARGS=(-s "$SEGMENT_LENGTH" -l "$BLOCK_LENGTH" -p "$MAP_PCT_ID"
           -k "$MIN_MATCH_LEN" -K "$MASH_KMER")
[[ -n "$SPARSE_MAP" ]] && PGGB_ARGS+=(-x "$SPARSE_MAP")
PGGB_ARGS+=(--keep-temp-files)

# pggb -r is --resume: "do not overwrite existing outputs in the given
# directory". All three 2025 builds ran with resume:true, which means their final
# graphs may incorporate intermediates from earlier attempts under other
# parameters -- unverifiable after the fact. Resume is therefore OPT-IN here, via
# MTB_PGGB_RESUME=1, instead of being passed silently on every run.
#
# THIS CHECK COMES BEFORE ANYTHING IS WRITTEN INTO THE DIRECTORY. The provenance
# sidecar used to be written first, so the directory was never empty when the
# check ran and a clean build was impossible without resume (audit PGB-1).
RESUME=0
if [[ "${MTB_PGGB_RESUME:-0}" == "1" ]]; then
    echo "[pggb_build] WARNING: --resume requested; existing intermediates in"
    echo "             ${OUT_HOST} will be reused and NOT rebuilt."
    PGGB_ARGS+=(-r)
    RESUME=1
elif [[ -n "$(ls -A "$OUT_HOST" 2>/dev/null)" ]]; then
    echo "[pggb_build] FATAL: ${OUT_HOST} is not empty and resume is off." >&2
    echo "             Move it aside, or set MTB_PGGB_RESUME=1 to reuse it." >&2
    exit 1
fi
mkdir -p "$OUT_HOST"

IN_C="$(mtb_in_data "$IN_HOST")"
OUT_C="$(mtb_in_graphs "$OUT_HOST")"
PGGB_CMD=(-i "$IN_C" -o "$OUT_C" -n "$N_HAP" -t "$THREADS"
          "${PGGB_ARGS[@]}" "${EXTRA_ARGS[@]}")

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
# support. Every setting pggb receives is a row of its own, the container and
# the input FASTA by sha256 as well as by name, and `finished` is appended only
# when pggb exits 0 -- a sidecar without it describes an attempt, not a graph.
_prov="${OUT_HOST}/graph_provenance.tsv"
{
    _repo="$(git -C "$_repo_dir" config --get remote.origin.url 2>/dev/null || true)"
    _commit="$(git -C "$_repo_dir" rev-parse HEAD 2>/dev/null || true)"
    _dirty="clean"
    git -C "$_repo_dir" diff --quiet HEAD 2>/dev/null || _dirty="dirty"
    [[ -n "$_commit" ]] || { _repo="not-a-git-checkout"; _commit="unknown"; _dirty="unknown"; }
    printf 'panel_repo\t%s\n'   "${_repo:-no-remote-configured}"
    printf 'panel_commit\t%s\n' "$_commit"
    printf 'panel_dirty\t%s\n'  "$_dirty"
    printf 'pggb_container\t%s\n' "$MTB_PGGB_SIF"
    printf 'pggb_container_sha256\t%s\n' "$(sha256sum "$MTB_PGGB_SIF" | cut -d' ' -f1)"
    printf 'pggb_version\t%s\n' "$(singularity exec --cleanenv "$MTB_PGGB_SIF" \
        bash /usr/local/bin/pggb --version 2>/dev/null | head -1 || echo unknown)"
    printf 'pggb_args\t%s\n'   "${PGGB_CMD[*]}"
    printf 'segment_length\t%s\n' "$SEGMENT_LENGTH"
    printf 'block_length\t%s\n'   "$BLOCK_LENGTH"
    printf 'map_pct_id\t%s\n'     "$MAP_PCT_ID"
    printf 'min_match_len\t%s\n'  "$MIN_MATCH_LEN"
    printf 'mash_kmer\t%s\n'      "$MASH_KMER"
    printf 'sparse_map\t%s\n'     "${SPARSE_MAP:-off}"
    printf 'n_haplotypes\t%s\n'   "$N_HAP"
    printf 'threads\t%s\n'        "$THREADS"
    printf 'resume\t%s\n'         "$RESUME"
    printf 'extra_args\t%s\n'     "${EXTRA_ARGS[*]:-}"
    printf 'input_fasta\t%s\n' "$(basename "$IN_HOST")"
    printf 'input_fasta_sha256\t%s\n' "$(sha256sum "$IN_HOST" | cut -d' ' -f1)"
    printf 'built\t%s\n'       "$(date -Is)"
    printf 'built_by\t%s\n'    "${USER:-unknown}"
} > "$_prov"
echo "[pggb_build] provenance: ${_prov}"

set -x
singularity exec --cleanenv $(mtb_bind) "$MTB_PGGB_SIF" \
    bash /usr/local/bin/pggb "${PGGB_CMD[@]}"
set +x

printf 'finished\t%s\n' "$(date -Is)" >> "$_prov"
echo "[pggb_build] finished -> ${OUT_HOST}"
echo "[pggb_build] REMEMBER: bin/sync_back.sh"
