#!/usr/bin/env bash
#SBATCH -J fix1_p0
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0-04:00
#SBATCH -o /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/rerun_scale200/logs/%x_%A_%a.out
# Build 7713a8d71d8e-fix1: run P0 steps in order, from the run root.
#   sbatch 02_p0.sbatch stamp paths accessions
#   sbatch --array=1-333 02_p0.sbatch refs        (one reference per task)
# The build is the CX333 graph under a new id, so nothing is written into
# build 7713a8d71d8e. Build inputs made outside P0 are in
# refbias/build_inputs/7713a8d71d8e-fix1/ (see README.md).
set -euo pipefail
R=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
cd "$R/runroot"
source "$R/config/project_env.sh"
export PATH=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates/bin:$PATH
export BUILD_ID=7713a8d71d8e-fix1
export OG=$MTB_GRAPHS/CX333.s10k.k23.K15/mtb.complex333.fasta.gz.f4f5ee2.11fba48.36e68b6.smooth.final.og
I=$MTB_WORK/refbias/build_inputs/7713a8d71d8e-fix1
export GRAPH_VCF=$I/graphvcf/all_variants.collapsed.vcf.gz
export PANEL_SNPS=$I/graphvcf/snps.vcf.gz
export ACCESSORY_DIR=$I/accessory_panel
export INSGT_CLUSTERS=$MTB_WORK/insgt/assets/insertions.tsv
export INSGT_ROUTING=$MTB_WORK/insgt/assets/insertion_routing.tsv
export ANC_TREE=$I/tree/cx333.rooted.nwk ANC_ALN=$I/tree/cx333.snps.fasta ANC_SITES=$I/tree/cx333.sites.tsv
for s in "$@"; do
    echo "=== step $s ($(date -Is))"
    bash bin/p0_prepare.sh --step "$s"
done
