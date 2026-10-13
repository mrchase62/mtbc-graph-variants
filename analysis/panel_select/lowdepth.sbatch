#!/bin/bash
# Low-depth intervals (depth < 5, all reads counted by samtools depth defaults)
# for each new1k isolate, from P1's kept H37Rv BAM. Used to mark isolate
# positions unknown instead of REF in the panel distance. Submit from this
# directory: sbatch --array=0-9 lowdepth.sbatch
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -J panel_lowdepth
#SBATCH -o logs/lowdepth_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../../config/project_env.sh
W=../../runroot/refbias/work/new1k_p1
mkdir -p lowdepth
export W MTB_SAMTOOLS
tail -n +2 ../../runroot/refbias/cohort.new1k.tsv | cut -f1 | sort \
    | awk -v t=$SLURM_ARRAY_TASK_ID 'NR % 10 == t' \
    | xargs -P 4 -I{} bash -c '
        "$MTB_SAMTOOLS" depth -a "$W/{}.h37rv.bam" \
        | awk -v OFS="\t" "\$3 < 5 { if (\$2 == e + 1 && c == \$1) { e = \$2 } else { if (c != \"\") print c, s - 1, e; c = \$1; s = \$2; e = \$2 } } END { if (c != \"\") print c, s - 1, e }" \
        | gzip > lowdepth/{}.bed.gz.tmp && mv lowdepth/{}.bed.gz.tmp lowdepth/{}.bed.gz'
echo "task $SLURM_ARRAY_TASK_ID done"
