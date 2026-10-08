#!/usr/bin/env bash
#SBATCH -J fix1_paneltree
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH -t 0-06:00
#SBATCH -o /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/rerun_scale200/logs/%x_%j.out
# Build 7713a8d71d8e-fix1: the panel tree for step ancestral, by
# docs/PANEL_TREE.md section 1 (alignment from the fixed graph VCF with H37Rv
# as an all-REF row, then IQ-TREE GTR+F+ASC+G4 rooted on canettii ET1291).
set -euo pipefail
R=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
cd "$R/runroot"
source "$R/config/project_env.sh"
I=$MTB_WORK/refbias/build_inputs/7713a8d71d8e-fix1
mkdir -p "$I/tree"
"$MTB_PY" bin/vcf_to_alignment.py --vcf "$I/graphvcf/all_variants.collapsed.vcf.gz" \
    --out "$I/tree/cx333.snps.fasta" --sites-out "$I/tree/cx333.sites.tsv" \
    --ref-sample GCF_000195955 --max-missing 0.10
bash bin/build_snp_tree.sh "$I/tree/cx333.snps.fasta" "$MTB_OUTGROUP" "$I/tree/cx333"
