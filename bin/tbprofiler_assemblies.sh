#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 0-01:30
#SBATCH -p shared
#SBATCH --mem=10000
#SBATCH -J tbprof_asm
#SBATCH -o slurm/%x_%A_%a.out
#SBATCH -e slurm/%x_%A_%a.err
#
# Lineage assignment and mixed-lineage screening for ASSEMBLIES.
#
#   sbatch --array=1-N bin/tbprofiler_assemblies.sh          # array over the panel
#   bash bin/tbprofiler_assemblies.sh <accession>            # single genome
#
# Adapted from the lab's read-based tbprofiler runner. Two changes: input is
# `--fasta` rather than `-1/-2`, matching how the 2025 calls were made (their
# sample names are FASTA filenames), and the manifest is the rebuild panel.
#
# Run BEFORE rotation, not after: lineage calling is rotation-independent because
# tb-profiler aligns to H37Rv itself, so screening first avoids rotating genomes
# that are about to be excluded for mixed lineage.
#
# The conda install at envs/mtb_isolates/bin/tb-profiler is BROKEN -- it dies on
# import with ModuleNotFoundError: pydantic_core._pydantic_core. Use the container.
set -euo pipefail

_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh" \
          "${MTB_WORK:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: no project_env.sh" >&2; exit 1; }
source "$_mtb_env"

CONTAINER="${TBPROFILER_CONTAINER:-/n/boslfs02/LABS/sfortune_lab/Lab/containers/tbprofiler_v1.3.sif}"
DB_DIR="${TBPROFILER_DB:-/n/boslfs02/LABS/sfortune_lab/Lab/shared/tbprofiler_db}"
PANEL="${PANEL:-${MTB_DATA}/ncbi/panel.rebuild.tsv}"
ASM_DIR="${ASM_DIR:-${MTB_DATA}/assemblies}"
OUT_DIR="${OUT_DIR:-${MTB_DATA}/tbprofiler}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

[[ -f "$CONTAINER" ]] || { echo "ERROR: container not found: $CONTAINER" >&2; exit 1; }
[[ -d "$DB_DIR" ]]    || { echo "ERROR: db dir not found: $DB_DIR" >&2; exit 1; }
[[ -f "$PANEL" ]]     || { echo "ERROR: panel not found: $PANEL" >&2; exit 1; }

if [[ $# -ge 1 ]]; then
    sample="$1"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "ERROR: no SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    # +1 skips the panel header
    sample="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$PANEL")"
    [[ -n "$sample" ]] || { echo "ERROR: no panel row $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
fi

src="${ASM_DIR}/${sample}.fna.gz"
[[ -f "$src" ]] || { echo "ERROR: assembly not found: $src" >&2; exit 1; }

mkdir -p "$OUT_DIR" slurm
if [[ -s "${OUT_DIR}/${sample}.tbprofiler.done" ]]; then
    echo "already done: $sample"; exit 0
fi

tmp_dir="/tmp/tbp.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$tmp_dir"
trap 'rm -rf "$tmp_dir"' EXIT

# tb-profiler --fasta wants plain FASTA, so decompress into the scratch dir
fa="${tmp_dir}/${sample}.fasta"
zcat "$src" > "$fa"

declare -A _seen
binds=()
_add_bind() {
    local d="$1"; [[ -n "$d" && -d "$d" ]] || return 0
    d="$(cd "$d" && pwd -P)"
    [[ -z "${_seen["$d"]+x}" ]] && { _seen["$d"]=1; binds+=( -B "$d:$d" ); }
    return 0
}
_add_bind "$tmp_dir"; _add_bind "$DB_DIR"; _add_bind "$OUT_DIR"; _add_bind "$ASM_DIR"

echo "sample=$sample  db=$DB_DIR  threads=$THREADS"
singularity exec "${binds[@]}" "$CONTAINER" \
    tb-profiler profile \
        --fasta "$fa" \
        -p "$sample" \
        -d "$tmp_dir" \
        --db_dir "$DB_DIR" \
        --threads "$THREADS" \
        --csv

if [[ -d "${tmp_dir}/results" ]]; then
    mv "${tmp_dir}/results/${sample}".* "$OUT_DIR/" 2>/dev/null || true
else
    echo "ERROR: no results at ${tmp_dir}/results" >&2; exit 1
fi
touch "${OUT_DIR}/${sample}.tbprofiler.done"
echo "done: $sample"
