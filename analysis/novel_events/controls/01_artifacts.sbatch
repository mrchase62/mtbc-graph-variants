#!/bin/bash
# Library artifacts in the 25 in-panel controls (P2 BAMs, marinB63), one
# array task per control. Submit from runroot:
#   sbatch --array=1-25 ../analysis/novel_events/controls/01_artifacts.sbatch
#
#   1. Split-read orientation, whole genome (samtools view, primary reads):
#      per 1,000 mapped reads, reads whose supplementary alignment (SA tag) is
#      on the same contig, the opposite strand and within 1 kb: inverted local
#      chimeras ("fold-back" reads); and other same-contig split reads.
#      Also the share of reads with the same start, CIGAR, strand and mate
#      start as an earlier read (a duplicate proxy; P2 does not mark them).
#   2. Oxidative damage: GATK 4.1.3 CollectSequencingArtifactMetrics, the
#      standard 8-oxoG test. Damage shows as G>T errors on one read of the
#      pair and C>A on the other ("pre-adapter" errors); the tool reports a
#      Phred-scaled rate per substitution, lower is worse. Reference files
#      are symlinked into a work folder so the .dict and .fai stay out of the
#      build directory.
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 00:40:00
#SBATCH -J ctl_artifacts
#SBATCH -o ../analysis/novel_events/logs/ctl_artifacts_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
O=$N/controls/out
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $N/controls/controls.txt)
R=$(awk -F'\t' -v s="$S" '$1==s{print $5}' refbias/marinB63/run/p1/refmap.tsv)
B=refbias/work/marinB63_p2/$S.bam
mkdir -p $O/ref/$R $O/$S
{ printf "sample\treads\tdup_like\tinv_chimera_per_1k\tother_split_per_1k\n"
  "$MTB_SAMTOOLS" view -@ 1 -F 0x904 $B | awk -v s=$S -f $N/controls/split_reads.awk; } > $O/$S/split_reads.tsv
F=$O/ref/$R/$R.fasta
[[ -e $F ]] || ln -s "$(readlink -f refbias/build/7713a8d71d8e-fix1/refs/$R.fasta)" $F
[[ -e $F.fai ]] || "$MTB_SAMTOOLS" faidx $F
G="singularity exec -B /n/netscratch,/n/boslfs02 /n/boslfs02/LABS/sfortune_lab/Lab/containers/gatk_4.1.3.0.sif gatk"
[[ -e $O/ref/$R/$R.dict ]] || $G CreateSequenceDictionary -R $F -O $O/ref/$R/$R.dict 2>/dev/null || true
$G --java-options "-Xmx4g" CollectSequencingArtifactMetrics -I $B -R $F -O $O/$S/artifacts \
    -EXT .txt 2> $O/$S/gatk.log
cat $O/$S/split_reads.tsv
grep -A20 "^## METRICS" $O/$S/artifacts.pre_adapter_summary_metrics.txt | cut -f3-9 | head -15
