#!/bin/bash
# Real vs synthetic reads, in-panel controls (the user's design, 2026-10-09).
# For each control: synthetic reads from its complete (polished) Marin
# assembly, matched to its real Illumina run (the run Marin used for
# polishing, as aligned in marinB63 P2), so that the only difference between
# the two read sets is that one is real. One array task per line of
# ../controls/controls.txt (1-25). Submit from runroot:
#   sbatch --array=1-25 ../analysis/novel_events/realsyn/01_reads.sbatch
#   1. samtools stats on the real P2 BAM: read count, read length, insert
#      size mean and standard deviation.
#   2. wgsim from the assembly with that read length, insert size and pair
#      count (same depth), 0.2% base error, no mutations (the assembly is the
#      truth), seed fixed per isolate.
#   3. Unaligned BAM refbias/marinSyn25/reads/syn_<sample>.bam.
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -J realsyn_reads
#SBATCH -o ../analysis/novel_events/realsyn/logs/reads_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $N/controls/controls.txt)
FILE=$(awk -F'\t' -v s=$S '$1==s{print $3}' $N/phaseB/truth_genomes.tsv)
O=refbias/marinSyn25
mkdir -p $O/reads $O/params
W=$(mktemp -d /tmp/realsyn.$S.XXXXXX)
trap '[[ -n "$W" && "$W" == /tmp/realsyn.* ]] && rm -rf -- "$W"' EXIT
"$MTB_SAMTOOLS" stats -@ 2 refbias/work/marinB63_p2/$S.bam | grep "^SN" | cut -f2-3 > $W/sn.txt
v() { awk -F'\t' -v k="$1" '$1==k":"{print $2}' $W/sn.txt; }
PAIRS=$(( $(v "raw total sequences") / 2 )); RL=$(v "maximum length")
INS=$(printf "%.0f" "$(v "insert size average")"); SD=$(printf "%.0f" "$(v "insert size standard deviation")")
printf "sample\tpairs\tread_length\tinsert_mean\tinsert_sd\tfile\n%s\t%s\t%s\t%s\t%s\t%s\n" \
    $S $PAIRS $RL $INS $SD $FILE > $O/params/$S.tsv
zcat ../analysis/external_assemblies/genomes/$FILE > $W/g.fa
"$MTB_WGSIM" -N $PAIRS -1 $RL -2 $RL -d $INS -s $SD -e 0.002 -r 0 -R 0 -X 0 \
    -S $((SLURM_ARRAY_TASK_ID + 900)) $W/g.fa $W/r1.fq $W/r2.fq > /dev/null
"$MTB_SAMTOOLS" import -1 $W/r1.fq -2 $W/r2.fq -o $O/reads/syn_$S.tmp.bam
mv -f $O/reads/syn_$S.tmp.bam $O/reads/syn_$S.bam
cat $O/params/$S.tsv
