#!/usr/bin/env bash
#SBATCH --job-name=og_test_tree
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 0-03:00
#SBATCH -p sapphire
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Outgroup test: one arm's tree, with the production scale200_cx333 settings
# (data/trees/scale200_cx333.log). Usage: sbatch run_tree.sbatch G|K
set -euo pipefail
ARM="${1:?arm G or K}"
cd "${SLURM_SUBMIT_DIR}"
/n/boslfs02/LABS/sfortune_lab/Lab/software/iqtree-3.0.1-Linux-intel/bin/iqtree3 \
    -s "${ARM}.fasta" -m GTR+F+ASC+G4 -B 1000 -alrt 1000 -o GCF_035581225 \
    -T "${SLURM_CPUS_PER_TASK:-24}" --prefix "${ARM}" -redo
