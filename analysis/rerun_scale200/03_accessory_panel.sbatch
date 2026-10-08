#!/usr/bin/env bash
#SBATCH -J fix1_accpanel
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0-01:00
#SBATCH -o /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/rerun_scale200/logs/%x_%j.out
# Build 7713a8d71d8e-fix1: the accessory panel from the fixed graph VCF, by the
# commands in docs/PANEL_TREE.md section 4, from the run root.
#   sbatch 03_accessory_panel.sbatch extract     (step a)
#   then the blast array (step b), see README.md
#   sbatch 03_accessory_panel.sbatch build       (steps c and d)
set -euo pipefail
R=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
cd "$R/runroot"
source "$R/config/project_env.sh"
I=$MTB_WORK/refbias/build_inputs/7713a8d71d8e-fix1
V=$I/graphvcf/all_variants.collapsed.vcf.gz
P=$I/accessory_panel
mkdir -p "$P/t2"
case "$1" in
extract)
    "$MTB_BCFTOOLS" query -l "$V" > "$P/t2/genomes.txt"
    "$MTB_PY" bin/t2_extract_candidates.py --vcf "$V" --bcftools "$MTB_BCFTOOLS" \
        --out-fasta "$P/t2/candidates.fasta" --out-meta "$P/t2/candidates.tsv"
    echo "genomes $(wc -l < "$P/t2/genomes.txt"), candidates $(grep -c '^>' "$P/t2/candidates.fasta")" ;;
build)
    "$MTB_PY" bin/t2_summary.py --hits "$P/t2/hits" --meta "$P/t2/candidates.tsv" \
        --out "$P/t2/accessory_validation.tsv"
    "$MTB_PY" bin/build_accessory_panel.py --validation "$P/t2/accessory_validation.tsv" \
        --candidates "$P/t2/candidates.tsv" --fasta "$P/t2/candidates.fasta" \
        --genomes "$P/t2/genomes.txt" --h37rv "$MTB_H37RV" --blastn "$MTB_QC_BIN/blastn" \
        --outdir "$P"
    ls -la "$P" ;;
esac
