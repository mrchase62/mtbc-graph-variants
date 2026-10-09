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
    # resume a partial file (-C -): ENA resets long transfers (SRR9588558_2,
    # 1.7 GB, twice), and a restart from zero hits the same reset
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
