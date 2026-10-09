#!/bin/bash
# GRIDSS 2.13.2, one array task per isolate, resources from the timing test
# (seff 51564761: 4 cores at 50% CPU efficiency, 3.4 GB peak of 16 GB; several
# GRIDSS steps are single-threaded Java). 2 cores, 6 GB, 4 GB Java heap.
# Args: <sample list> <cohort> <p1 refmap>  (BAMs: refbias/work/<cohort>_p2/)
# Submit from runroot:
#   sbatch --array=1-N ../analysis/novel_events/gridss/02_array.sbatch <list> <cohort> <refmap>
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 00:15:00
#SBATCH -J gridss
#SBATCH -o ../analysis/novel_events/logs/gridss_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export PATH=$HOME/.conda/envs/gridss/bin:$PATH
LIST="${1:?sample list}"; C="${2:?cohort}"; RM="${3:?refmap}"
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
N=../analysis/novel_events
O=$N/gridss/out/$C
REFS=refbias/build/7713a8d71d8e-fix1/refs
r=$(awk -F'\t' -v s=$s '$1==s{print $5}' "$RM")
rd=$N/gridss/out/ref/$r
mkdir -p $rd $O/$s
for x in fasta fasta.fai dict fasta.amb fasta.ann fasta.bwt fasta.pac fasta.sa; do
    [[ -e $rd/$r.$x ]] || ln -s "$(readlink -f $REFS/$r.$x)" $rd/$r.$x || true
done
t0=$(date +%s)
gridss -r $rd/$r.fasta -t 2 --jvmheap 4g -b $N/out/masks/$r.bed \
    -o $O/$s/$s.gridss.vcf -a $O/$s/$s.assembly.bam --workingdir $O/$s/work \
    refbias/work/${C}_p2/$s.bam > $O/$s/gridss.log 2>&1
echo "$s $r $(( $(date +%s) - t0 )) s"
