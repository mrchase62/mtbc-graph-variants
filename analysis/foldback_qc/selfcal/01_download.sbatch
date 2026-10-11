#!/bin/bash
# Fold-back calibration step 2: download one panel genome's Illumina run
# (paired FASTQ) from ENA and check its md5. As novel_events/phaseB/
# 01_download.sbatch (resume on reset), columns found by header name.
# User approved all 54 runs (35.6 GB), 2026-10-11. Row = array index + 2 of
# selfcal/genomes.tsv. Submit from the repository root:
#   sbatch --array=0-53%8 analysis/foldback_qc/selfcal/01_download.sbatch
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 0-04:00
#SBATCH -J selfcal_dl
#SBATCH -o analysis/foldback_qc/selfcal/logs/dl_%A_%a.out
set -euo pipefail
D=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/panel_illumina
T=analysis/foldback_qc/selfcal/genomes.tsv
mkdir -p "$D"
n=$((SLURM_ARRAY_TASK_ID + 2))
col() { awk -F'\t' -v n="$n" -v c="$1" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} NR==n{print $h[c]}' "$T"; }
FTP=$(col fastq_ftp); MD5=$(col fastq_md5)
[[ "$FTP" == ftp.sra.ebi.ac.uk/* ]] || { echo "FATAL: row $n has no ENA FASTQ link: '$FTP'" >&2; exit 1; }
IFS=';' read -r -a URLS <<< "$FTP"
IFS=';' read -r -a SUMS <<< "$MD5"
for k in "${!URLS[@]}"; do
    f="$D/$(basename "${URLS[$k]}")"
    if [[ -s "$f" ]] && [[ "$(md5sum "$f" | cut -d' ' -f1)" == "${SUMS[$k]}" ]]; then
        echo "$f present, md5 ok"; continue
    fi
    ok=0
    for attempt in 1 2 3 4 5 6 7 8 9 10; do
        if curl -sS -C - --retry 5 --retry-delay 20 -o "$f.part" "https://${URLS[$k]}"; then
            ok=1; break
        fi
        echo "attempt $attempt interrupted at $(stat -c %s "$f.part" 2>/dev/null || echo 0) bytes; resuming" >&2
        sleep 30
    done
    [[ $ok -eq 1 ]] || { echo "FATAL: $f not complete after 10 attempts" >&2; exit 1; }
    got=$(md5sum "$f.part" | cut -d' ' -f1)
    [[ "$got" == "${SUMS[$k]}" ]] || { echo "FATAL: md5 mismatch for $f ($got vs ${SUMS[$k]})" >&2; exit 1; }
    mv -f "$f.part" "$f"
    echo "$f downloaded, md5 ok"
done
