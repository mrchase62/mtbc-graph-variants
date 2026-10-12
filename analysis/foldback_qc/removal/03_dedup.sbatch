#!/bin/bash
# Duplicate-marking test (DEDUP_PLAN.md), one isolate per array task. From
# the P2 BAM of all reads (<isolate>_all): mark duplicates (samtools
# markdup -S), remove them, make a random depth control of the same size,
# then run dysgu, our caller v6 + local assembly and the fold-back scan.
# Submit from runroot:
#   sbatch --array=0-8 ../analysis/foldback_qc/removal/03_dedup.sbatch
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -J fbremove_dedup
#SBATCH -o ../analysis/foldback_qc/removal/logs/dedup_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
R=../analysis/foldback_qc/removal
N=../analysis/novel_events
O=$R/out/dedup
W=refbias/work/fbremove_dedup
B=refbias/build/7713a8d71d8e-fix1/refs
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin
DYSGU=refbias/work/svvenv/bin/dysgu
S=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $R/isolates.txt)
[[ -n "$S" ]] || { echo "no isolate for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
REF=$(awk -F'\t' -v s="${S}_all" '$1==s{print $5}' refbias/fbremove/run/p1/refmap.tsv)
[[ -n "$REF" ]] || { echo "no reference for $S" >&2; exit 1; }
FA=$B/$REF.fasta
ALL=refbias/work/fbremove_p2/${S}_all.bam
mkdir -p $O/p2 $O/proto $O/asm $O/asm_work $O/foldback $O/sites $O/dups $W
T=$(mktemp -d "${TMPDIR:-/tmp}/dedup.$S.XXXXXX")

# 1-2. mark, then remove
"$MTB_SAMTOOLS" collate -@ 4 -O -T $T/col $ALL \
  | "$MTB_SAMTOOLS" fixmate -@ 4 -m - - \
  | "$MTB_SAMTOOLS" sort -@ 4 -m 2G -T $T/srt - \
  | "$MTB_SAMTOOLS" markdup -@ 4 -S -T $T/md -f $O/dups/$S.markdup.txt - $W/${S}_dupmark.bam
"$MTB_SAMTOOLS" index $W/${S}_dupmark.bam
"$MTB_SAMTOOLS" view -@ 4 -b -F 1024 -o $W/${S}_dedup.bam $W/${S}_dupmark.bam
"$MTB_SAMTOOLS" index $W/${S}_dedup.bam

# 3. random depth control, 7. duplicate shares
"$MTB_SAMTOOLS" view -f 1024 -F 0x900 $W/${S}_dupmark.bam | cut -f1 | sort -u > $T/dup.names
"$MTB_SAMTOOLS" view -F 0x900 $ALL | cut -f1 | sort -u > $T/all.names
"$MTB_PY" - $T $R/names/$S.fb.names.gz $S $O/dups/$S.tsv <<'PY'
import gzip, random, sys
t, fbf, s, out = sys.argv[1:]
dup = set(x.strip() for x in open(t + "/dup.names"))
allp = [x.strip() for x in open(t + "/all.names")]
fb = set(x.strip() for x in gzip.open(fbf, "rt"))
random.seed(7)
open(t + "/rand.names", "w").writelines(x + "\n" for x in sorted(random.sample(allp, len(dup))))
other = len(allp) - len(fb)
with open(out, "w") as fo:
    fo.write("sample\tpairs\tdup_pairs\tfoldback_pairs\tfoldback_dup\tother_dup\tdup_share_foldback\tdup_share_other\n")
    fd = len(fb & dup)
    fo.write(f"{s}\t{len(allp)}\t{len(dup)}\t{len(fb)}\t{fd}\t{len(dup) - fd}\t"
             f"{fd / max(len(fb), 1):.4f}\t{(len(dup) - fd) / max(other, 1):.4f}\n")
PY
"$MTB_SAMTOOLS" view -@ 4 -b -N "^$T/rand.names" -o $W/${S}_randdup.bam $ALL
"$MTB_SAMTOOLS" index $W/${S}_randdup.bam

# 4. dysgu as in P2 (bin/p2_call.sh)
for k in dupmark dedup randdup; do
    "$DYSGU" run -x -p 4 $FA $T/dy_$k $W/${S}_$k.bam > $O/p2/${S}_$k.dysgu.vcf 2> $O/p2/${S}_$k.dysgu.log
done
# HaplotypeCaller is not rerun; self_noise.py needs a VCF, so link _all's
for k in dupmark dedup randdup; do
    ln -sf "$(readlink -f refbias/fbremove/run/p2/${S}_all.vcf.gz)" $O/p2/${S}_$k.vcf.gz
done

# 5. our caller v6 + local assembly; 6. fold-back scan after removal
for k in dedup randdup; do
    BAM=$W/${S}_$k.bam
    "$MTB_PY" $N/breakpoint_caller.py --bam $BAM --sample ${S}_$k --ref $FA \
        --mask-bed $N/out/masks/$REF.bed --out $O/proto/${S}_$k.events.tsv
    "$MTB_PY" $N/assemble_candidates.py --events $O/proto/${S}_$k.events.tsv --bam $BAM \
        --ref $FA --sample ${S}_$k --spades "$E/python $E/spades.py" --minimap2 $E/minimap2 \
        --bnd-min-frac 0.2 --max-windows 40 \
        --workdir $O/asm_work/${S}_$k --out $O/asm/${S}_$k.events.tsv
done
"$MTB_PY" ../analysis/foldback_qc/foldback_v2.py --bam $W/${S}_dedup.bam --ref $FA --sample ${S}_dedup \
    --out $O/foldback/${S}_dedup.tsv --sites $O/sites/${S}_dedup.tsv.gz
[[ -n "$T" && "$T" == "${TMPDIR:-/tmp}/dedup.$S."* ]] && rm -rf -- "$T"
echo "done $S"
