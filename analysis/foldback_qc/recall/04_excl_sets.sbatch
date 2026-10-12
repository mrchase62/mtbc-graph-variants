#!/bin/bash
# Excluded-isolates check (user, 2026-10-12): could removing fold-back reads
# replace excluding the sample? For the 12 excluded selfcal54 isolates
# (v2 > 4.9), as ../removal/01_sets.sbatch did for the chimeric controls:
#   1. fold-back read names (../foldback_v2.py --names) from the own-genome
#      alignment (selfcal54 P2);
#   2. the QC'd unaligned reads minus every fold-back pair (samtools view
#      -N ^names) -> refbias/fbrecallx100/reads/<isolate>_fb100.bam;
#   3. make_excl100.py writes the runner tables (cohort fbrecallx100: each
#      <isolate>_fb100 on the same reference as in fbrecall).
# Submit from runroot:  sbatch ../analysis/foldback_qc/recall/04_excl_sets.sbatch
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J excl_sets
#SBATCH -o ../analysis/foldback_qc/recall/logs/excl_sets_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
C=refbias/fbrecallx100
X=../analysis/foldback_qc/recall/names
mkdir -p $C/reads $X
export C X MTB_PY MTB_SAMTOOLS
"$MTB_PY" ../analysis/foldback_qc/recall/make_excl100.py --list > $X/excluded.txt
while read s; do awk -F'\t' -v s=$s '$1==s{print s"\t"$5}' refbias/selfcal54/run/p1/refmap.tsv; done < $X/excluded.txt |
  xargs -P 4 -L 1 bash -c '
    "$MTB_PY" ../analysis/foldback_qc/foldback_v2.py --bam refbias/work/selfcal54_p2/$0.bam \
        --ref refbias/build/7713a8d71d8e-fix1/refs/$1.fasta --sample $0 \
        --out $X/$0.rate.tsv --names $X/$0.fb.names.gz
    zcat $X/$0.fb.names.gz | sort -u > $X/$0.fb100.names
    "$MTB_SAMTOOLS" view -b -N "^$X/$0.fb100.names" -o $C/reads/$0_fb100.tmp.bam refbias/selfcal54/reads/$0.bam
    mv -f $C/reads/$0_fb100.tmp.bam $C/reads/$0_fb100.bam
    echo "$0 fold-back pairs $(wc -l < $X/$0.fb100.names)"'
"$MTB_PY" ../analysis/foldback_qc/recall/make_excl100.py
