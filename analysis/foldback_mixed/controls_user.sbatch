#!/bin/bash
# Step 6, rerun with the user's within-host Mutect2 settings (user,
# 2026-10-14: "yes, run both frames"), in two frames:
#   h37rv   - real: refbias/work/marinB63_p1/<s>.h37rv.bam;
#             syn:  refbias/work/fbmix_syn/syn_<s>.h37rv.bam
#   matched - P2 BAMs on each control's matched panel genome:
#             real: refbias/work/marinB63_p2/<s>.bam;
#             syn:  refbias/work/marinSyn25_p2/syn_<s>.bam
# Each BAM: downsample to <= 100x mean depth (as the default run), read
# groups, MarkDuplicates (not removed), then Mutect2 / LearnReadOrientationModel
# / FilterMutectCalls with the user's arguments, verbatim from their script.
# 100 tasks: task-1 = control index*4 + arm*2 + frame. Submit from runroot:
#   sbatch --array=1-100%40 ../analysis/foldback_mixed/controls_user.sbatch
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -J fbmix_user
#SBATCH -o ../analysis/foldback_mixed/logs/user_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
i=$((SLURM_ARRAY_TASK_ID - 1))
S=$(sed -n "$((i / 4 + 1))p" ../analysis/novel_events/controls/controls.txt)
ARM=$([[ $(( (i / 2) % 2 )) -eq 0 ]] && echo real || echo syn)
FRAME=$([[ $((i % 2)) -eq 0 ]] && echo h37rv || echo matched)
case $ARM/$FRAME in
    real/h37rv)   B=refbias/work/marinB63_p1/$S.h37rv.bam ;;
    syn/h37rv)    B=refbias/work/fbmix_syn/syn_$S.h37rv.bam ;;
    real/matched) B=refbias/work/marinB63_p2/$S.bam ;;
    syn/matched)  B=refbias/work/marinSyn25_p2/syn_$S.bam ;;
esac
if [[ $FRAME == h37rv ]]; then
    REF=refbias/build/7713a8d71d8e-fix1/refs/GCF_000195955.fasta
else
    REF=$("$MTB_SAMTOOLS" view -H "$B" | grep -m1 "^@PG.*bwa mem" | grep -o "refbias/build/[^ ]*\.fasta")
fi
O=../analysis/foldback_mixed/controls_user
P=$O/${FRAME}_${ARM}_$S
T=$(mktemp -d /tmp/fbuser.$SLURM_ARRAY_TASK_ID.XXXXXX)
trap '[[ -n "$T" && "$T" == /tmp/fbuser.* ]] && rm -rf -- "$T"' EXIT
G() { singularity exec -B /n/boslfs02 -B /n/netscratch -B /tmp "$MTB_GATK_SIF" gatk "$@"; }
D=$("$MTB_SAMTOOLS" depth -a "$B" | awk '{s+=$3} END{print s/NR}')
F=$(awk -v d="$D" 'BEGIN{printf "%.4f", (d > 100) ? 100 / d : 1}')
if awk -v f="$F" 'BEGIN{exit !(f < 1)}'; then
    "$MTB_SAMTOOLS" view -@ 4 -b -s "$((SLURM_ARRAY_TASK_ID + 700))${F#0}" -o $T/ds.bam "$B"
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
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$S" "$ARM" "$FRAME" "$(basename $REF .fasta)" "$D" "$F" > $P.info.tsv
echo "$S $ARM $FRAME done"
