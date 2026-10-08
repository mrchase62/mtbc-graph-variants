#!/bin/bash
# Phase B: download one Marin 2024 Illumina run (paired FASTQ) from ENA and
# check its md5. Row = array index + 2 of analysis/novel_events/phaseB/runs.tsv.
# User approved the download on 2026-10-08 (63 runs, 25.6 GB).
# Submit from the repository root:
#   sbatch --array=0-62%10 analysis/novel_events/phaseB/01_download.sbatch
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 0-02:00
#SBATCH -J marin_dl
#SBATCH -o analysis/novel_events/logs/dl_%A_%a.out
set -euo pipefail
D=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_illumina
T=analysis/novel_events/phaseB/runs.tsv
# Fields by column number with awk: bash `read` merges consecutive tabs, so an
# empty column (same_isolate_in_cx333 is often empty) shifted every field
# after it in the first submission (51415076).
n=$((SLURM_ARRAY_TASK_ID + 2))
FTP=$(awk -F'\t' -v n="$n" 'NR==n{print $9}' "$T")
MD5=$(awk -F'\t' -v n="$n" 'NR==n{print $10}' "$T")
[[ "$FTP" == ftp.sra.ebi.ac.uk/* ]] || { echo "FATAL: row $n has no ENA FASTQ link: '$FTP'" >&2; exit 1; }
IFS=';' read -r -a URLS <<< "$FTP"
IFS=';' read -r -a SUMS <<< "$MD5"
for k in "${!URLS[@]}"; do
    f="$D/$(basename "${URLS[$k]}")"
    if [[ -s "$f" ]] && [[ "$(md5sum "$f" | cut -d' ' -f1)" == "${SUMS[$k]}" ]]; then
        echo "$f present, md5 ok"; continue
    fi
    curl -sS --retry 5 --retry-delay 20 -o "$f.part" "https://${URLS[$k]}"
    got=$(md5sum "$f.part" | cut -d' ' -f1)
    [[ "$got" == "${SUMS[$k]}" ]] || { echo "FATAL: md5 mismatch for $f ($got vs ${SUMS[$k]})" >&2; exit 1; }
    mv -f "$f.part" "$f"
    echo "$f downloaded, md5 ok"
done
