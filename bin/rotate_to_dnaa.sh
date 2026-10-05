#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=8G
#SBATCH -J rotate_dnaa
#SBATCH -o slurm/rotate_%A_%a.out
#SBATCH -e slurm/rotate_%A_%a.err
#
# Rotate one assembly so that dnaA (Rv0001) starts at position 1, forward strand.
#
#   sbatch --array=1-N bin/rotate_to_dnaa.sh        # array over the panel
#   bash bin/rotate_to_dnaa.sh <accession>          # one genome
#
# WHY THIS IS NOT THE 2025 SCRIPT
# The 2025 rotate_to_dnaa_circlator_patched.sh judged success as "the output file
# is non-empty". That cannot detect the one thing that matters: whether dnaA
# actually ended up at the start, on the forward strand. A circlator run that
# silently did nothing was logged OK. Here every rotation is VERIFIED by mapping
# the dnaA gene back to the rotated sequence, and the achieved offset and strand
# are recorded per genome so the whole panel can be audited afterwards.
#
# It also reads .fna.gz (the 2025 version required plain .fna) and runs as an
# array rather than a serial loop over 491 genomes.
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

CIRCLATOR="${CIRCLATOR_SIF:-/n/boslfs02/LABS/sfortune_lab/Lab/containers/circlator_latest.sif}"
DNAA="${MTB_DATA}/annotation/Mtb_dnaA_coding.fasta"
PANEL="${PANEL:-${MTB_DATA}/ncbi/panel.rebuild.tsv}"
ASM_DIR="${ASM_DIR:-${MTB_DATA}/assemblies}"
OUT_DIR="${OUT_DIR:-${MTB_DATA}/rotated}"
LOG_DIR="${OUT_DIR}/verify"

mtb_require_file "$CIRCLATOR" "$DNAA" "$PANEL" "$MTB_MINIMAP2"

if [[ $# -ge 1 ]]; then
    acc="$1"
else
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "ERROR: no SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    acc="$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$PANEL")"
    [[ -n "$acc" ]] || { echo "ERROR: no panel row $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
fi

src="${ASM_DIR}/${acc}.fna.gz"
mtb_require_file "$src"
mkdir -p "$OUT_DIR" "$LOG_DIR" slurm
final="${OUT_DIR}/${acc}.dnaA_rotated.fasta"
[[ -s "$final" && -s "${LOG_DIR}/${acc}.tsv" ]] && { echo "already done: $acc"; exit 0; }

tmp="$(mktemp -d "/tmp/rot.${SLURM_JOB_ID:-$$}.XXXX")"
trap 'rm -rf "$tmp"' EXIT
zcat "$src" > "${tmp}/in.fasta"

echo "[rotate] $acc"
# circlator runs nucmer in the CURRENT directory and leaves a tmp.run_nucmer.*
# behind, then dies trying to rmdir it. Run with the CWD inside the scratch dir
# so that debris lands there and is cleaned by the trap, rather than piling up in
# the project root.
cp "$DNAA" "${tmp}/dnaa.fasta"
( cd "$tmp" && singularity exec --pwd "$tmp" -B "${tmp}:${tmp}" \
    "$CIRCLATOR" circlator fixstart in.fasta out --genes_fa dnaa.fasta \
    > circlator.log 2>&1 ) || true

if [[ ! -s "${tmp}/out.fasta" ]]; then
    echo -e "${acc}\tNA\tNA\tNA\tNA\tFAILED_no_output" > "${LOG_DIR}/${acc}.tsv"
    echo "FAILED (no output): $acc" >&2
    tail -5 "${tmp}/circlator.log" >&2 || true
    exit 1
fi

# --- verification: where did dnaA actually land, and on which strand?
"$MTB_MINIMAP2" -x asm5 --secondary=no -t "${SLURM_CPUS_PER_TASK:-2}" \
    "${tmp}/out.fasta" "$DNAA" 2>/dev/null > "${tmp}/dnaa.paf" || true

read -r start strand alen tlen <<< "$(awk -F'\t' 'BEGIN{best=0}
    ($10+0)>best { best=$10+0; s=$8; st=$5; al=$10; tl=$7 }
    END{ if(best>0) print s, st, al, tl; else print "NA","NA","0","0" }' "${tmp}/dnaa.paf")"

status="OK"
if [[ "$start" == "NA" ]]; then
    status="FAILED_dnaA_not_found"
elif [[ "$strand" != "+" ]]; then
    status="REVERSE_STRAND"
elif (( start > 100 )); then
    status="OFFSET_${start}"
fi

cp "${tmp}/out.fasta" "$final"
"$MTB_SAMTOOLS" faidx "$final" 2>/dev/null || true
echo -e "${acc}\t${start}\t${strand}\t${alen}\t${tlen}\t${status}" > "${LOG_DIR}/${acc}.tsv"
echo "[rotate] $acc  dnaA at ${start} (${strand})  -> ${status}"
[[ "$status" == OK ]] || echo "NOTE: $acc needs attention: $status" >&2
