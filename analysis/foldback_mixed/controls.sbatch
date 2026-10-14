#!/bin/bash
# Step 6 of the fold-back vs mixed-call proposal: Mutect2 on the 25 in-panel
# controls' real and synthetic reads (user, 2026-10-13: "run 1-7").
# Real: P1's H37Rv BAM (refbias/work/marinB63_p1/<s>.h37rv.bam). Synthetic:
# refbias/marinSyn25/reads/syn_<s>.bam (wgsim from the isolate's own
# assembly, matched length/insert/depth, 0.2% error), aligned to H37Rv here
# with P1's own commands (simulate_and_call.sh on supplied FASTQs). Both are
# downsampled to <= 100x mean depth, then Mutect2 tumour-only with
# --min-base-quality-score 20 (the collection's settings, per the user) and
# FilterMutectCalls. One task per control. Submit from runroot:
#   sbatch --array=1-25 ../analysis/foldback_mixed/controls.sbatch
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 03:00:00
#SBATCH -J fbmix_ctrl
#SBATCH -o ../analysis/foldback_mixed/logs/ctrl_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../analysis/novel_events/controls/controls.txt)
H37RV=refbias/build/7713a8d71d8e-fix1/refs/GCF_000195955.fasta
O=../analysis/foldback_mixed/controls
W=refbias/work/fbmix_syn
mkdir -p $O $W
G() { singularity exec -B /n/boslfs02 -B /n/netscratch "$MTB_GATK_SIF" gatk --java-options "-Xmx8g" "$@"; }
# synthetic reads to H37Rv, as P1 aligns real reads
if [[ ! -s $W/syn_$S.h37rv.bam ]]; then
    "$MTB_SAMTOOLS" fastq -@ 2 -n -1 $W/syn_${S}_1.fq -2 $W/syn_${S}_2.fq -0 /dev/null -s /dev/null \
        refbias/marinSyn25/reads/syn_$S.bam
    SIM_REQUIRE_FQ=1 SIM_FQ_PREFIX=$W/syn_$S SIM_BAM_SUFFIX=".h37rv" MTB_THREADS=4 \
        bash bin/simulate_and_call.sh $H37RV $H37RV $W syn_$S 1 150
    rm -f $W/syn_${S}_1.fq $W/syn_${S}_2.fq
fi
for arm in real syn; do
    if [[ $arm == real ]]; then B=refbias/work/marinB63_p1/$S.h37rv.bam; else B=$W/syn_$S.h37rv.bam; fi
    D=$("$MTB_SAMTOOLS" depth -a "$B" | awk '{s+=$3} END{print s/NR}')
    F=$(awk -v d="$D" 'BEGIN{f = (d > 100) ? 100 / d : 1; printf "%.4f", f}')
    T=$W/${arm}_$S.ds.bam
    if awk -v f="$F" 'BEGIN{exit !(f < 1)}'; then
        "$MTB_SAMTOOLS" view -@ 4 -b -s "$((SLURM_ARRAY_TASK_ID + 700))${F#0}" -o $T "$B"
    else
        cp "$B" $T
    fi
    "$MTB_SAMTOOLS" index $T
    G AddOrReplaceReadGroups -I $T -O $T.rg.bam -RGID $arm -RGSM ${arm}_$S -RGLB lib -RGPL ILLUMINA -RGPU u \
        > /dev/null 2>&1
    "$MTB_SAMTOOLS" index $T.rg.bam
    G Mutect2 -R $H37RV -I $T.rg.bam --min-base-quality-score 20 -O $O/${arm}_$S.unfiltered.vcf.gz \
        > $O/${arm}_$S.mutect2.log 2>&1
    G FilterMutectCalls -R $H37RV -V $O/${arm}_$S.unfiltered.vcf.gz -O $O/${arm}_$S.vcf.gz \
        > $O/${arm}_$S.filter.log 2>&1
    printf "%s\t%s\t%s\t%s\n" "$S" "$arm" "$D" "$F" > $O/${arm}_$S.depth.tsv
    rm -f $T $T.bai $T.rg.bam $T.rg.bam.bai
done
echo "$S done"
