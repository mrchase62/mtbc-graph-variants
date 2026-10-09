#!/bin/bash
# GRIDSS 2.13.2 timing test (approved 2026-10-09): 2 Phase B isolates,
# mar_QC_6 (107x) and mar_TB1612 (75x), on their P2 BAMs (bwa mem to the
# matched reference R), with the R uncallable mask (callable_mask.py) as the
# exclude list. Report-only comparison; nothing in the pipeline changes.
# GRIDSS writes cache files beside the reference, so it gets symlinks to the
# build's reference files in its own folder, not the build directory itself.
# Submit from runroot:  sbatch ../analysis/novel_events/gridss/01_time.sbatch
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J gridss_time
#SBATCH -o ../analysis/novel_events/logs/gridss_time_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
G=$HOME/.conda/envs/gridss/bin
export PATH=$G:$PATH
N=../analysis/novel_events
O=$N/gridss/out
REFS=refbias/build/7713a8d71d8e-fix1/refs
for s in mar_QC_6 mar_TB1612; do
    r=$(awk -F'\t' -v s=$s '$1==s{print $5}' refbias/marinB63/run/p1/refmap.tsv)
    rd=$O/ref/$r
    mkdir -p $rd $O/$s
    for x in fasta fasta.fai dict fasta.amb fasta.ann fasta.bwt fasta.pac fasta.sa; do
        [[ -e $rd/$r.$x ]] || ln -s "$(readlink -f $REFS/$r.$x)" $rd/$r.$x
    done
    t0=$(date +%s)
    gridss -r $rd/$r.fasta -t 4 -b $N/out/masks/$r.bed \
        -o $O/$s/$s.gridss.vcf -a $O/$s/$s.assembly.bam --workingdir $O/$s/work \
        refbias/work/marinB63_p2/$s.bam > $O/$s/gridss.log 2>&1
    echo "$s $r $(( $(date +%s) - t0 )) s"
done
