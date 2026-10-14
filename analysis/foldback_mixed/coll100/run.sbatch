#!/bin/bash
# Collection comparison, 50 clean + 50 chimeric (user, 2026-10-14: "Run the 50
# clean and 50 chimeric too"). One task per line of samples.tsv:
#   1. reads from the CRAM (read in place), as P1 extracts them;
#   2. align to the isolate's matched reference (P1 refmap, build fix1) and
#      call haploid HaplotypeCaller there: bin/simulate_and_call.sh on the
#      real reads, the command P1 and P2 use;
#   3. the user's within-host Mutect2 settings (as controls_user.sbatch) on
#      the matched-reference BAM and on P1's H37Rv BAM; <= 100x, read
#      groups, MarkDuplicates first.
# Submit from runroot after both P1 runs:
#   sbatch --array=1-100%40 ../analysis/foldback_mixed/coll100/run.sbatch
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -J fbmix100
#SBATCH -o ../analysis/foldback_mixed/logs/c100_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
C=../analysis/foldback_mixed/coll100
read -r S GRP LIN FBR REL < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $C/samples.tsv)
COH=$([[ $GRP == clean ]] && echo new1k || echo fbmix50)
BUILD=refbias/build/7713a8d71d8e-fix1
REFID=$(awk -F'\t' -v s=$S '$1==s{print $5}' refbias/$COH/p1/refmap.tsv)
[[ -n "$REFID" ]] || { echo "FATAL: $S has no reference in refbias/$COH/p1/refmap.tsv" >&2; exit 1; }
MREF=$BUILD/refs/$REFID.fasta
H37BAM=refbias/work/${COH}_p1/$S.h37rv.bam
W=refbias/work/fbmix100; mkdir -p $W
O=$C/calls
if [[ ! -s $W/$S.matched.bam ]]; then
    "$MTB_SAMTOOLS" collate -@ 2 -u -O -T $W/$S.collate --reference "$MTB_CRAM_REF" "$MTB_CRAM_ROOT/$REL" \
      | "$MTB_SAMTOOLS" fastq -@ 2 -n -1 $W/${S}_1.fq -2 $W/${S}_2.fq -0 /dev/null -s /dev/null -
    SIM_REQUIRE_FQ=1 SIM_FQ_PREFIX=$W/$S SIM_BAM_SUFFIX=".matched" MTB_THREADS=4 \
        bash bin/simulate_and_call.sh $MREF $MREF $W $S 1 150
    mv -f $W/$S.vcf.gz $O/hc_matched_$S.vcf.gz
    rm -f $W/${S}_1.fq $W/${S}_2.fq $W/$S.vcf.gz.tbi
fi
T=$(mktemp -d /tmp/fb100.$SLURM_ARRAY_TASK_ID.XXXXXX)
trap '[[ -n "$T" && "$T" == /tmp/fb100.* ]] && rm -rf -- "$T"' EXIT
G() { singularity exec -B /n/boslfs02 -B /n/netscratch -B /tmp "$MTB_GATK_SIF" gatk "$@"; }
for FRAME in matched h37rv; do
    if [[ $FRAME == matched ]]; then B=$W/$S.matched.bam; REF=$MREF; else B=$H37BAM; REF=$BUILD/refs/GCF_000195955.fasta; fi
    P=$O/m2_${FRAME}_$S
    D=$("$MTB_SAMTOOLS" depth -a "$B" | awk '{s+=$3} END{print s/NR}')
    F=$(awk -v d="$D" 'BEGIN{printf "%.4f", (d > 100) ? 100 / d : 1}')
    if awk -v f="$F" 'BEGIN{exit !(f < 1)}'; then
        "$MTB_SAMTOOLS" view -@ 4 -b -s "$((SLURM_ARRAY_TASK_ID + 900))${F#0}" -o $T/ds.bam "$B"
    else
        cp "$B" $T/ds.bam
    fi
    G AddOrReplaceReadGroups -I $T/ds.bam -O $T/rg.bam -RGID $S -RGSM $S -RGLB lib -RGPL ILLUMINA -RGPU u > /dev/null 2>&1
    G --java-options "-Xmx10g" MarkDuplicates -I $T/rg.bam -O $T/md.bam --REMOVE_DUPLICATES false \
        --METRICS_FILE $P.dedupMets.txt --VALIDATION_STRINGENCY LENIENT --TMP_DIR $T > /dev/null 2>&1
    "$MTB_SAMTOOLS" index $T/md.bam
    G --java-options "-Xmx12g" Mutect2 -R $REF -I $T/md.bam -O $P.unfiltered.vcf.gz \
        --annotation StrandBiasBySample --read-filter NotSupplementaryAlignmentReadFilter \
        --min-base-quality-score 20 --dont-use-soft-clipped-bases --min-pruning 2 -max-mnp-distance 1 \
        --native-pair-hmm-threads 4 --f1r2-tar-gz $T/f1r2.tar.gz --minimum-allele-fraction 0.05 \
        --kmer-size 21 --kmer-size 31 --active-probability-threshold 0.005 \
        --max-reads-per-alignment-start 0 > $P.mutect2.log 2>&1
    G --java-options "-Xmx4g" LearnReadOrientationModel -I $T/f1r2.tar.gz -O $T/rom.tar.gz > $P.rom.log 2>&1
    G --java-options "-Xmx3g" FilterMutectCalls -V $P.unfiltered.vcf.gz -R $REF -O $P.vcf.gz \
        --microbial-mode --ob-priors $T/rom.tar.gz --min-reads-per-strand 5 --min-allele-fraction 0.05 \
        --max-alt-allele-count 4 --stats $P.unfiltered.vcf.gz.stats > $P.filter.log 2>&1
    rm -f $T/*.bam $T/*.bai $T/*.tar.gz
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$S" "$GRP" "$LIN" "$FRAME" "$(basename $REF .fasta)" "$D" "$F" > $P.info.tsv
done
echo "$S $GRP done"
