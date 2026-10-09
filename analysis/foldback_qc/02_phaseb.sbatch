#!/bin/bash
# Fold-back chimera rate for the Phase B isolates without one yet (the 37
# novel and setE_polished isolates; the 25 in-panel controls were measured
# by novel_events/controls/01_artifacts.sbatch on the same BAMs). Reads the
# marinB63 P2 BAMs (real reads aligned to each isolate's matched reference),
# with the same count (novel_events/controls/split_reads.awk).
# Submit from runroot:  sbatch ../analysis/foldback_qc/02_phaseb.sbatch
# Output: analysis/foldback_qc/phaseB_rates.tsv (all 62, controls copied in)
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=6G
#SBATCH -t 00:40:00
#SBATCH -J foldback_phaseb
#SBATCH -o ../analysis/foldback_qc/logs/phaseb_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
N=../analysis/novel_events
Q=../analysis/foldback_qc
O=$(mktemp -d /tmp/foldback.phaseb.XXXXXX)
trap '[[ -n "$O" && "$O" == /tmp/foldback.phaseb.* ]] && rm -rf -- "$O"' EXIT
export O N MTB_SAMTOOLS
awk -F'\t' 'NR>1 && $5!="in_panel"{print $1}' $N/phaseB/out/samples.tsv | xargs -P 8 -I{} bash -c '
    "$MTB_SAMTOOLS" view -F 0x904 refbias/work/marinB63_p2/{}.bam \
        | awk -v s={} -f $N/controls/split_reads.awk > $O/{}.tsv'
for s in $(cat $N/controls/controls.txt); do tail -n +2 $N/controls/out/$s/split_reads.tsv > $O/$s.tsv; done
{ printf "sample\treads\tdup_like\tinv_chimera_per_1k\tother_split_per_1k\n"; cat $O/*.tsv; } > $Q/phaseB_rates.tsv
echo "$(($(wc -l < $Q/phaseB_rates.tsv) - 1)) isolates"
