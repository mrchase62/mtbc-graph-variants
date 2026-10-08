#!/usr/bin/env bash
#SBATCH -J fix1_graphvcf
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 0-03:00
#SBATCH -o /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/rerun_scale200/logs/%x_%j.out
# Build 7713a8d71d8e-fix1, step 1: the CX333 graph VCF with the fixed collapse.
# Input: the graph's decomposed VCF (vcfbub -l 0 -a 100000, vcfwave -I 1000,
# 2026-09-11; the same settings bin/vcf_decompose.sh uses now). The fixed
# bin/vcf_collapse.sh (GRAPHVCF-5 trim-and-union, D22 left-alignment) and the
# class split write into the build's input folder; the graph folder is untouched.
set -euo pipefail
R=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
source "$R/config/project_env.sh"
export PATH=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates/bin:$PATH
G=$MTB_GRAPHS/CX333.s10k.k23.K15
O=$MTB_WORK/refbias/build_inputs/7713a8d71d8e-fix1/graphvcf
mkdir -p "$O"
bash "$R/bin/vcf_collapse.sh" CX333.s10k.k23.K15 --in "$G/all_variants.decomposed.vcf.gz" \
    --out "$O/all_variants.collapsed.vcf.gz" --ref "$MTB_REF_FASTA"
bash "$R/bin/vcf_split_classes.sh" CX333.s10k.k23.K15 --in "$O/all_variants.collapsed.vcf.gz" --outdir "$O"
sha256sum "$O"/*.vcf.gz > "$O/sha256.txt"
