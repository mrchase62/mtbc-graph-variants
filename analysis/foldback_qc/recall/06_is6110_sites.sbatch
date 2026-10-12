#!/bin/bash
# IS6110 junction check (user, 2026-10-12): fold-back junction positions
# (../foldback_v2.py --sites, unchanged) for every part D sample on its
# part D reference (cohort fbrecall, 48 BAMs), so is6110_check.py can ask
# whether reads at new IS6110 insertions are labelled fold-back when the
# labelling is done on a shared reference, as it would be in production.
# Submit from runroot:  sbatch ../analysis/foldback_qc/recall/06_is6110_sites.sbatch
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J is6110_sites
#SBATCH -o ../analysis/foldback_qc/recall/logs/is6110_sites_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
O=../analysis/foldback_qc/recall/out/sites_fbrecall
mkdir -p $O
export O MTB_PY
tail -n +2 refbias/fbrecall/run/p1/refmap.tsv | cut -f1,5 | xargs -P 8 -L 1 bash -c '
  "$MTB_PY" ../analysis/foldback_qc/foldback_v2.py --bam refbias/work/fbrecall_p2/$0.bam \
      --ref refbias/build/7713a8d71d8e-fix1/refs/$1.fasta --sample $0 \
      --out $O/$0.rate.tsv --sites $O/$0.tsv.gz
  echo "$0 done"'
