#!/usr/bin/env bash
# Compare the resize test's arms with each other and with production's P1
# output for the same sample (refbias/work/scale200_p1, 2026-09-23).
# Alignments: every read record (header lines dropped, since they name the
# command line). Calls: every VCF record (header dropped).
set -euo pipefail
R=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
T=$R/analysis/resize_test
source "$R/config/project_env.sh"
S=SAMN12126251
PROD=$MTB_WORK/refbias/work/scale200_p1
BCF=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates/bin/bcftools
bam() { "$MTB_SAMTOOLS" view "$1" | md5sum | cut -c1-12; }
vcf() { "$BCF" view -H "$1" | md5sum | cut -c1-12; }
n()   { "$MTB_SAMTOOLS" view -c "$1"; }
printf "%-12s %-14s %-14s %12s %10s\n" source alignments calls reads calls_n
for a in production old new8 new4; do
    if [[ $a == production ]]; then B=$PROD/$S.h37rv.bam; V=$PROD/$S.h37rv.vcf.gz
    else B=$T/$a/$S.h37rv.bam; V=$T/$a/$S.vcf.gz; fi
    printf "%-12s %-14s %-14s %12s %10s\n" $a "$(bam $B)" "$(vcf $V)" "$(n $B)" "$("$BCF" view -H $V | wc -l)"
done
echo
for a in old new8 new4; do
    echo "$a: $(grep -E 'Maximum resident' $T/$a/time.txt | sed 's/^\s*//'), $(grep -E 'Elapsed' $T/$a/time.txt | sed 's/^\s*//')"
done
