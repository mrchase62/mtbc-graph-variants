#!/bin/bash
# The association tail for scale200_fix, as a job. Submit from runroot.
# First run: steps 1-3, which submit the tree job and stop. Second run, after
# the tree: steps 3b-5 (the tests). The phenotype is scale200's, unchanged.
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -J assoc_scale200_fix
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export MTB_BUILD_DIR="$(readlink -f refbias/build/7713a8d71d8e-fix1)"
bash /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/assoc/bin/cohort_assoc_tail.sh \
    scale200_fix assoc/scale200_fix/rrdr_carriers.txt
