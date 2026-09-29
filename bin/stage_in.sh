#!/usr/bin/env bash
# Stage inputs from durable storage onto scratch.
# Safe to re-run: rsync only copies what changed.
#
#   bin/stage_in.sh              # minimal: containers + ref + graph input FASTAs (~6 GB)
#   bin/stage_in.sh rotated      # + rotated_assemblies (~9 GB)
#   bin/stage_in.sh graphs       # + the existing 2025 graphs (~84 GB)
#   bin/stage_in.sh all          # everything above
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

TIER="${1:-minimal}"
RSYNC=(rsync -ah --info=progress2 --no-inc-recursive)

mtb_require_work
mkdir -p "$MTB_CONTAINERS" "$MTB_REFDIR" "$MTB_FASTAS" "$MTB_ROTATED" "$MTB_ASSEMBLIES" "$MTB_GRAPHS"

echo "### staging tier: ${TIER}"
echo "### from: ${MTB_ARCHIVE}"
echo "### to:   ${MTB_WORK}"

echo "--- containers"
"${RSYNC[@]}" "${MTB_ARCHIVE}/pggb_latest.sif"  "$MTB_CONTAINERS/"
"${RSYNC[@]}" "${MTB_ARCHIVE}/vg_v1.69.0.sif"   "$MTB_CONTAINERS/"

echo "--- reference"
"${RSYNC[@]}" "${MTB_ARCHIVE}/data/ref/"        "$MTB_REFDIR/"
"${RSYNC[@]}" "${MTB_ARCHIVE}/H37Rv.fasta"      "$MTB_REFDIR/"

echo "--- graph input FASTAs"
"${RSYNC[@]}" "${MTB_ARCHIVE}/data/fastas/"     "$MTB_FASTAS/"

echo "--- annotation resources (RD intervals, IS6110 loci)"
mkdir -p "$MTB_DATA/annotation"
for f in known_RDs.bed known_RDs_L1L2L4bovis.bed; do
    src="${MTB_ARCHIVE}/data/pggb_mtb_1_2_4_bovis/structural_variants/${f}"
    [[ -e "$src" ]] && "${RSYNC[@]}" "$src" "$MTB_DATA/annotation/"
done
# The GFF ships with the bare accession as seqid; the graph uses the PanSN name.
IS_SRC=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/Databases/H37Rv_IS6110.gff
if [[ -e "$IS_SRC" ]]; then
    awk -v r="$MTB_REF_PATH" 'BEGIN{FS=OFS="\t"} /^#/{print;next} {$1=r; print}' \
        "$IS_SRC" > "$MTB_DATA/annotation/H37Rv_IS6110.pansn.gff"
fi

echo "--- metadata tables"
for f in tbprofiler.txt tbprof_lineages.csv lineage_defs.tsv list.tsv \
         collinearity_synteny_summary.tsv; do
    [[ -e "${MTB_ARCHIVE}/${f}" ]] && "${RSYNC[@]}" "${MTB_ARCHIVE}/${f}" "$MTB_DATA/"
done

if [[ "$TIER" == "rotated" || "$TIER" == "all" ]]; then
    echo "--- rotated assemblies (~9 GB)"
    "${RSYNC[@]}" "${MTB_ARCHIVE}/rotated_assemblies/" "$MTB_ROTATED/"
    "${RSYNC[@]}" "${MTB_ARCHIVE}/data/assemblies/"    "$MTB_ASSEMBLIES/"
fi

if [[ "$TIER" == "graphs" || "$TIER" == "all" ]]; then
    echo "--- existing 2025 graphs (~84 GB)"
    for g in mtb.complex.p95.s10k.k51 pggb_mtb_1_2_4_bovis \
             pggb_mtb_1_2_4_bovis_p99_s10k_k61 mtb_lineage4 pggb_mtb_plcC_glyS; do
        [[ -d "${MTB_ARCHIVE}/data/${g}" ]] || continue
        echo "    ${g}"
        "${RSYNC[@]}" --exclude 'archive/' "${MTB_ARCHIVE}/data/${g}/" "${MTB_GRAPHS}/${g}/"
    done
fi

echo "### done. Disk used on scratch:"
du -sh "$MTB_WORK"
