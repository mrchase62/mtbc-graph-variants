#!/bin/bash
# Read-removal experiment, step 2 (../FOLDBACK_SV_EFFECT_PLAN.md, part C).
# After P2 of the <isolate>_all samples (QC'd reads on their own genome):
#   1. fold-back read names (../foldback_v2.py --names, all reads);
#   2. from the QC'd unaligned BAM (reads/<isolate>.bam), three read sets,
#      removing whole pairs by name (samtools view -N ^file):
#        _fb100  every fold-back pair
#        _fb50   half of them (names with an even CRC32)
#        _rand   as many pairs as _fb100, drawn at random (seed 7) from the
#                pairs not counted as fold-back: the depth control
#   3. make_sets.py --stage 2 adds them to the runner tables.
# Submit from runroot:  sbatch --dependency=afterok:<p2 summary> ../analysis/foldback_qc/removal/01_sets.sbatch
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J fbremove_sets
#SBATCH -o ../analysis/foldback_qc/removal/logs/sets_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
C=refbias/fbremove
X=../analysis/foldback_qc/removal/names
mkdir -p $X
export C X MTB_PY MTB_SAMTOOLS
awk -F'\t' 'NR>1{print $1"\t"$5}' $C/run/p1/refmap.tsv | xargs -P 4 -L 1 bash -c '
    s=${0%_all}
    "$MTB_PY" ../analysis/foldback_qc/foldback_v2.py --bam refbias/work/fbremove_p2/$0.bam \
        --ref refbias/build/7713a8d71d8e-fix1/refs/$1.fasta --sample $0 \
        --out $X/$s.rate.tsv --names $X/$s.fb.names.gz
    "$MTB_SAMTOOLS" view $C/reads/$s.bam | cut -f1 | uniq > $X/$s.all.names
    "$MTB_PY" - "$X/$s" <<PY
import gzip, random, sys, zlib
p = sys.argv[1]
fb = set(x.strip() for x in gzip.open(p + ".fb.names.gz", "rt"))
allnames = [x.strip() for x in open(p + ".all.names")]
open(p + ".fb100.names", "w").writelines(x + "\n" for x in sorted(fb))
open(p + ".fb50.names", "w").writelines(x + "\n" for x in sorted(fb) if zlib.crc32(x.encode()) % 2 == 0)
rest = sorted(set(allnames) - fb)
random.seed(7)
open(p + ".rand.names", "w").writelines(x + "\n" for x in sorted(random.sample(rest, len(fb))))
print(p, "pairs", len(set(allnames)), "fold-back", len(fb))
PY
    for k in fb100 fb50 rand; do
        "$MTB_SAMTOOLS" view -b -N "^$X/$s.$k.names" -o $C/reads/${s}_$k.tmp.bam $C/reads/$s.bam
        mv -f $C/reads/${s}_$k.tmp.bam $C/reads/${s}_$k.bam
    done
    rm -f $X/$s.all.names'
"$MTB_PY" ../analysis/foldback_qc/removal/make_sets.py --stage 2
