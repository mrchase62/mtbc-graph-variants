#!/bin/bash
# Adapter read-through per sample (adapter_count.py) for one set (user,
# 2026-10-11): phaseB (62 Marin raw FASTQ pairs), production (the 500 random
# production CRAMs with a measured fold-back rate) or selfcal (54 raw FASTQ
# pairs; after the download job). About 1 s per sample. Submit from runroot:
#   sbatch ../analysis/foldback_qc/adapters/01_count.sbatch <phaseB|production|selfcal>
# Output: adapters/<set>.tsv
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -t 00:45:00
#SBATCH -J adapters
#SBATCH -o ../analysis/foldback_qc/adapters/logs/count_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
SET="${1:?set}"
A=../analysis/foldback_qc/adapters
R=$A/rows_$SET
mkdir -p $R
export A R MTB_PY MTB_SAMTOOLS MTB_H37RV MTB_CRAM_ROOT
case "$SET" in
  phaseB)
    D=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_illumina
    awk -F'\t' -v d=$D 'NR>1{print $1"\t"d"/"$7"\t"d"/"$8}' ../analysis/novel_events/phaseB/genomes.tsv ;;
  selfcal)
    D=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/panel_illumina
    awk -F'\t' -v d=$D 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} {print $1"\t"d"/"$h["r1"]"\t"d"/"$h["r2"]}' \
        ../analysis/foldback_qc/selfcal/genomes.tsv ;;
  production)
    "$MTB_PY" -c '
import csv, gzip, json, os
d = json.load(gzip.open(os.environ["MTB_CRAM_ROOT"] + "/metadata/samples.json.gz", "rt"))
for r in csv.DictReader(open("../analysis/novel_events/controls/production_chimera.tsv"), delimiter="\t"):
    print(r["sample"], os.environ["MTB_CRAM_ROOT"] + "/" + d[r["sample"]]["path"], "-", sep="\t")' ;;
  *) echo "unknown set $SET" >&2; exit 1 ;;
esac | xargs -P 4 -L 1 bash -c '
    if [[ $2 == - ]]; then
        "$MTB_PY" $A/adapter_count.py --sample $0 --cram $1 --ref "$MTB_H37RV" --samtools "$MTB_SAMTOOLS" --out $R/$0.tsv
    elif [[ -s $1 && -s $2 ]]; then
        "$MTB_PY" $A/adapter_count.py --sample $0 --fastq $1 $2 --out $R/$0.tsv
    else
        echo "missing FASTQ for $0" >&2
    fi'
{ head -1 $(ls $R/*.tsv | head -1); for f in $R/*.tsv; do tail -n +2 $f; done; } > $A/$SET.tsv
echo "$SET: $(($(wc -l < $A/$SET.tsv) - 1)) samples"
cut -f6 $A/$SET.tsv | tail -n +2 | sort | uniq -c
