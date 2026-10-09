#!/bin/bash
# Does polishing with fold-back chimeric short reads work as well as with
# clean ones? One array task per genome x read set (1-9): genomes R27252,
# QC-6, R37765; read sets 0%, 15% inverted-repeat, 15% random (the foldsim
# unaligned BAMs). Submit from runroot:
#   sbatch --array=1-9 ../analysis/foldback_sim/03_polish.sbatch
#   1. make_draft.py: the truth genome with 300 one-base homopolymer indels
#      and 30 substitutions (same draft for every read set of a genome).
#   2. Polypolish 0.6 as in the user's mtbc-hybrid-assembly 05_polish
#      (bwa mem -a per read file, polypolish filter, polypolish polish).
#   3. Pilon 1.24 as in Marin et al. (bwa mem pairs, default --fix all), one
#      round.
#   4. Each polished genome against the truth (minimap2 asm5, paftools call):
#      remaining differences by type, split into draft errors left unfixed
#      and new differences introduced by polishing.
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J foldsim_polish
#SBATCH -o ../analysis/foldback_sim/logs/polish_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
F=../analysis/foldback_sim
i=$((SLURM_ARRAY_TASK_ID - 1))
G=(R27252 QC_6 R37765); FILE=(R27252 QC-6 R37765); SET=(0 ir_150 random_150)
g=${G[$((i / 3))]}; f=${FILE[$((i / 3))]}; rs=${SET[$((i % 3))]}
S=fb_${g}_${rs}
O=$F/out/polish/$S
mkdir -p $O
CT="singularity exec -B /n/netscratch,/n/boslfs02 /n/boslfs02/LABS/sfortune_lab/Lab/containers/mtbc-assembly.sif"
E=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates/bin
W=$(mktemp -d /tmp/foldpol.$S.XXXXXX)
trap '[[ -n "$W" && "$W" == /tmp/foldpol.* ]] && rm -rf -- "$W"' EXIT
zcat ../analysis/external_assemblies/genomes/$f.fna.gz > $W/truth.fa
"$MTB_PY" $F/make_draft.py --truth $W/truth.fa --seed 7 --out $W/draft.fa --errors $O/draft_errors.tsv
"$MTB_SAMTOOLS" fastq -@ 2 -1 $W/r1.fq -2 $W/r2.fq -0 /dev/null -s /dev/null refbias/foldsim/reads/$S.bam
$E/bwa index $W/draft.fa 2>/dev/null
# Polypolish (user's pipeline)
$E/bwa mem -t 4 -a $W/draft.fa $W/r1.fq > $W/a1.sam 2>/dev/null
$E/bwa mem -t 4 -a $W/draft.fa $W/r2.fq > $W/a2.sam 2>/dev/null
$CT polypolish filter --in1 $W/a1.sam --in2 $W/a2.sam --out1 $W/f1.sam --out2 $W/f2.sam 2> $O/polypolish_filter.log
$CT polypolish polish $W/draft.fa $W/f1.sam $W/f2.sam > $W/polypolish.fa 2> $O/polypolish.log
rm -f $W/a1.sam $W/a2.sam $W/f1.sam $W/f2.sam
# Pilon (Marin et al.)
$E/bwa mem -t 4 $W/draft.fa $W/r1.fq $W/r2.fq 2>/dev/null | $E/samtools sort -@ 2 -o $W/pe.bam -
$E/samtools index $W/pe.bam
$E/pilon -Xmx12g --genome $W/draft.fa --frags $W/pe.bam --output pilon --outdir $W --changes > $O/pilon.log 2>&1
cp $W/pilon.changes $O/pilon.changes
# against the truth
for p in draft polypolish pilon; do
    fa=$W/$p.fa; [[ $p == pilon ]] && fa=$W/pilon.fasta
    $E/minimap2 -cx asm5 --cs $W/truth.fa $fa 2>/dev/null | sort -k6,6 -k8,8n \
        | $E/k8 $E/paftools.js call -l 1000 -L 1000 - 2>/dev/null > $O/$p.vs_truth.txt
done
"$MTB_PY" - $O <<'PY'
import sys, csv
o = sys.argv[1]
err = {int(r["truth_pos"]): r["type"] for r in csv.DictReader(open(f"{o}/draft_errors.tsv"), delimiter="\t")}
rows = []
for p in ("draft", "polypolish", "pilon"):
    diffs = []
    for line in open(f"{o}/{p}.vs_truth.txt"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "V" and f[4] == "1":  # covered by exactly one alignment
            diffs.append(int(f[2]) + 1)
    left = sum(1 for d in diffs if any(abs(d - e) <= 10 for e in err))
    rows.append(dict(genome_readset=o.rsplit("/", 1)[1], stage=p, differences=len(diffs),
                     draft_errors_left=left, new_errors=len(diffs) - left))
with open(f"{o}/result.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
for r in rows: print(r)
PY
