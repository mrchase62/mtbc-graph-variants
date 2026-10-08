#!/usr/bin/env bash
#SBATCH -J resize_test
#SBATCH -p shared
#SBATCH -t 0-03:00
#SBATCH -o logs/%x_%j.out
# Analysis only: does cutting P1/P2 jobs from 8 cores / 16 GB to 4 cores / 8 GB
# change any alignment or call? One scale200 sample, P1's H37Rv step
# (bin/simulate_and_call.sh), run as:
#   prep   extract the reads from the CRAM exactly as P1 does
#   old    the script before the change (git HEAD), 8 threads  (the old P1)
#   new8   the changed script (bwa -K, Java heap 6 GB), 8 threads
#   new4   the changed script, 4 threads, in a 4-core / 8 GB job (the new P1)
#   sbatch -c <cores> --mem=<mem> run_arm.sbatch <arm>
# Writes only under analysis/resize_test/.
set -euo pipefail
R=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
T=$R/analysis/resize_test
SAMPLE=SAMN12126251
source "$R/config/project_env.sh"
W=$MTB_WORK
H37RV=$W/refbias/build/7713a8d71d8e/refs/GCF_000195955.fasta
ARM=$1
cd "$R"
case "$ARM" in
prep)
    REL=$(awk -F'\t' -v s=$SAMPLE '$1==s{print $2}' $W/refbias/scale200.crams.tsv)
    mkdir -p $T/reads
    "$MTB_SAMTOOLS" collate -@ 2 -u -O -T $T/reads/$SAMPLE.collate \
        --reference "$MTB_CRAM_REF" "$MTB_CRAM_ROOT/$REL" \
      | "$MTB_SAMTOOLS" fastq -@ 2 -n -1 $T/reads/${SAMPLE}_1.fq -2 $T/reads/${SAMPLE}_2.fq -0 /dev/null -s /dev/null -
    echo "read pairs: $(( $(wc -l < $T/reads/${SAMPLE}_1.fq) / 4 ))" ;;
old|new8|new4)
    SCRIPT=bin/simulate_and_call.sh
    if [[ $ARM == old ]]; then git show HEAD:bin/simulate_and_call.sh > $T/old_simulate_and_call.sh; SCRIPT=$T/old_simulate_and_call.sh; fi
    NT=8; [[ $ARM == new4 ]] && NT=4
    mkdir -p $T/$ARM
    /usr/bin/time -v env SIM_REQUIRE_FQ=1 SIM_FQ_PREFIX=$T/reads/$SAMPLE SIM_KEEP_FQ=1 SIM_BAM_SUFFIX=.h37rv \
        MTB_THREADS=$NT bash $SCRIPT "$H37RV" "$H37RV" $T/$ARM $SAMPLE 1 150 2> $T/$ARM/time.txt
    grep -E "Maximum resident|Elapsed" $T/$ARM/time.txt ;;
esac
