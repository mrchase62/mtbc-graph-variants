#!/bin/bash
#SBATCH -p shared -c 4 --mem=16G -t 3:00:00 -J sveval_tree
#SBATCH -o /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/sv_eval_1k/out/tree.%j.log
set -euo pipefail
cd /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
source config/project_env.sh
$MTB_PY analysis/sv_eval_1k/build_tree.py
export OMP_NUM_THREADS=4
/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_pangenome_qc/bin/FastTree -nt -gtr -gamma analysis/sv_eval_1k/out/snps.fasta > analysis/sv_eval_1k/out/new1k.fasttree.nwk
