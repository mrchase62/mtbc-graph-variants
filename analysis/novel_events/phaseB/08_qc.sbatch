#!/bin/bash
# Phase B read QC, matching the user's standard cohort preprocessing
# (mtbc-hybrid-assembly 03a/03b, and the user, 2026-10-09: fastp, 100x,
# minimum 60x, TB-Profiler mixed removed, kraken2 contamination):
#   1. fastp: --detect_adapter_for_pe --correction --qualified_quality_phred 20
#      --length_required 36 (the user's settings). Optional third argument
#      "extra" adds 3' sliding-window trimming (--cut_right, 4-base window,
#      mean Q20) and poly-G trimming (--trim_poly_g). marinQC5 was run with
#      them; the user doubts the production CRAMs used them (production fastp
#      settings are not recorded) and wants them only if they help.
#   2. TB-Profiler 6.7.0 (lab container tbprofiler_v1.3.sif) on the trimmed
#      reads, --no_delly --no_trim, as in 03b. The database (118 MB) is copied
#      to the task's node-local folder first: with 51 tasks reading it on
#      boslfs02 at once, TB-Profiler's freebayes step took a median 310 s
#      against 13 s in a 5-task run (51710211, 2026-10-09)
#   3. rasusa 2.2.2: downsample to 100x (genome 4,411,532 bp, seed 7)
#   4. unaligned BAM for the runner (as 02_prepare.sbatch)
# kraken2 runs separately (09_kraken.sbatch: one job loads the 59 GB database
# once); qc_summary.py combines both into the PASS/FAIL table.
# One array task per line of the sample list ($1). Submit from runroot:
#   sbatch --array=1-N ../analysis/novel_events/phaseB/08_qc.sbatch <list> <cohort> [extra]
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -J marin_qc
#SBATCH -o ../analysis/novel_events/logs/marin_qc_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
LIST="${1:?sample list}"; C="${2:?cohort name}"
EXTRA=()
[[ "${3:-}" == extra ]] && EXTRA=(--cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --trim_poly_g)
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
G=../analysis/novel_events/phaseB/genomes.tsv
col() { awk -F'\t' -v s="$S" -v c="$1" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} $1==s{print $h[c]}' "$G"; }
R1=$(col r1); R2=$(col r2)
IN=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_illumina
CT=/n/boslfs02/LABS/sfortune_lab/Lab/containers
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates/bin
Q=refbias/$C/qc/$S
mkdir -p "$Q" refbias/$C/reads
W=$(mktemp -d /tmp/marinqc.$S.XXXXXX)
trap '[[ -n "$W" && "$W" == /tmp/marinqc.* ]] && rm -rf -- "$W"' EXIT
B="-B /n/netscratch,/n/boslfs02,$W"
singularity exec $B $CT/mtbc-assembly.sif fastp -i "$IN/$R1" -I "$IN/$R2" \
    -o $W/t1.fq.gz -O $W/t2.fq.gz \
    --detect_adapter_for_pe --correction --qualified_quality_phred 20 --length_required 36 \
    "${EXTRA[@]}" --json $Q/fastp.json --html $Q/fastp.html --thread 4
cp -rL /n/boslfs02/LABS/sfortune_lab/Lab/shared/tbprofiler_db $W/tbprofiler_db
singularity exec $B $CT/tbprofiler_v1.3.sif tb-profiler profile -1 $W/t1.fq.gz -2 $W/t2.fq.gz \
    -p "$S" -d $W/tbp --db_dir $W/tbprofiler_db \
    --no_delly --no_trim --threads 4 --csv
cp $W/tbp/results/$S.results.json $Q/tbprofiler.json
$E/rasusa reads -c 100 -g 4411532 -s 7 -o $W/d1.fq.gz -o $W/d2.fq.gz $W/t1.fq.gz $W/t2.fq.gz 2> $Q/rasusa.log
zcat $W/d1.fq.gz $W/d2.fq.gz | awk 'NR%4==2{n+=length($0)} END{print n}' > $Q/final_bases.txt
"$MTB_SAMTOOLS" import -@ 4 -1 $W/d1.fq.gz -2 $W/d2.fq.gz -o refbias/$C/reads/$S.tmp.bam
mv -f refbias/$C/reads/$S.tmp.bam refbias/$C/reads/$S.bam
echo "[$C] $S: $(cat $Q/final_bases.txt) bases after QC and downsampling"
