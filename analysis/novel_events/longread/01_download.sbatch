#!/bin/bash
# setE Nanopore runs (ENA PRJNA720906; same BioSample as the Illumina runs)
# for the tatC-helY element test. One run per array task; row = index + 2 of
# runs.tsv. User approved 2026-10-10 (17 runs, 13.8 GB).
# Submit from the repository root:
#   sbatch --array=0-16%6 analysis/novel_events/longread/01_download.sbatch
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH -t 0-02:00
#SBATCH -J ont_dl
#SBATCH -o analysis/novel_events/logs/ont_dl_%A_%a.out
set -euo pipefail
D=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/data/external_reads/marin_ont
T=analysis/novel_events/longread/runs.tsv
mkdir -p "$D"
n=$((SLURM_ARRAY_TASK_ID + 2))
FTP=$(awk -F'\t' -v n="$n" 'NR==n{print $7}' "$T")
MD5=$(awk -F'\t' -v n="$n" 'NR==n{print $8}' "$T")
[[ "$FTP" == ftp.sra.ebi.ac.uk/* && "$FTP" != *";"* ]] || { echo "FATAL: row $n: '$FTP'" >&2; exit 1; }
f="$D/$(basename "$FTP")"
if [[ -s "$f" ]] && [[ "$(md5sum "$f" | cut -d' ' -f1)" == "$MD5" ]]; then
    echo "$f present, md5 ok"; exit 0
fi
# 2026-10-10: 11 of 17 first downloads were full-size with a bad md5/gzip CRC and
# no interruption logged (silent corruption in transfer), so a file failing md5
# is deleted and fetched again from scratch, up to 3 times.
for full in 1 2 3; do
    ok=0
    for attempt in 1 2 3 4 5 6 7 8 9 10; do  # resume on ENA resets, as in phaseB/01_download.sbatch
        if curl -sS --fail -C - --retry 5 --retry-delay 20 -o "$f.part" "https://$FTP"; then ok=1; break; fi
        echo "attempt $attempt interrupted at $(stat -c %s "$f.part" 2>/dev/null || echo 0) bytes; resuming" >&2
        sleep 30
    done
    [[ $ok -eq 1 ]] || { echo "FATAL: $f not complete after 10 attempts" >&2; exit 1; }
    got=$(md5sum "$f.part" | cut -d' ' -f1)
    if [[ "$got" == "$MD5" ]]; then
        mv -f "$f.part" "$f"
        echo "$f downloaded, md5 ok (pass $full)"
        exit 0
    fi
    echo "pass $full: md5 mismatch for $f ($got vs $MD5); deleting and fetching again" >&2
    rm -f -- "$f.part"
done
echo "FATAL: md5 mismatch for $f after 3 full downloads" >&2
exit 1
