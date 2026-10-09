#!/bin/bash
# kraken2 on the raw Phase B read pairs, as in mtbc-hybrid-assembly 03a
# (lab standard database). One job for the whole list: the 59 GB database is
# memory-mapped once and stays in the page cache across samples, instead of
# being loaded by every array task. Submit from runroot:
#   sbatch ../analysis/novel_events/phaseB/09_kraken.sbatch <list> <cohort>
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=72G
#SBATCH -t 01:00:00
#SBATCH -J marin_kraken
#SBATCH -o ../analysis/novel_events/logs/marin_kraken_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
LIST="${1:?sample list}"; C="${2:?cohort name}"
G=../analysis/novel_events/phaseB/genomes.tsv
IN=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_illumina
DB=/n/boslfs02/LABS/sfortune_lab/Lab/kraken2standardDB
while read -r S; do
    read -r R1 R2 < <(awk -F'\t' -v s="$S" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} $1==s{print $h["r1"], $h["r2"]}' "$G")
    mkdir -p refbias/$C/qc/$S
    singularity exec -B /n/netscratch,/n/boslfs02 /n/boslfs02/LABS/sfortune_lab/Lab/containers/mtbc-assembly.sif \
        kraken2 --db $DB --memory-mapping --threads 4 --paired \
        --report refbias/$C/qc/$S/kraken2.tsv --output - "$IN/$R1" "$IN/$R2"
    echo "$S $(date +%T)"
done < "$LIST"
