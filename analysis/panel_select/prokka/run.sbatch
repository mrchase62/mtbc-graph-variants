#!/bin/bash
# Prokka annotation for the 6 P200 GenBank genomes submitted without one
# (PRJNA555636; user, 2026-10-14: "yes, annotate the 6 with Prokka"). Writes
# <build>/annotation/<acc>.gff.gz (features only, the ##FASTA section dropped)
# so P0 step gff and P6's off-path annotation can use them. Contig names are
# kept as in the assembly, which is what the build's refs/ FASTA uses.
#   sbatch analysis/panel_select/prokka/run.sbatch   (from the repo root)
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J prokka6
#SBATCH -o analysis/panel_select/prokka/prokka_%j.out
set -euo pipefail
source config/project_env.sh
B=runroot/refbias/build/3893a3762ed1/annotation
O=analysis/panel_select/prokka
export PATH="$MTB_QC_BIN:$PATH"
for a in GCA_047782915 GCA_047782945 GCA_047782925 GCA_047782875 GCA_047782895 GCA_047782885; do
    zcat analysis/genbank_scan/genomes/$a.fna.gz > $O/$a.fna
    prokka --outdir $O/$a --prefix $a --locustag ${a/GCA_/P} --kingdom Bacteria --genus Mycobacterium \
        --species tuberculosis --cpus 8 --force --quiet $O/$a.fna
    sed '/^##FASTA/,$d' $O/$a/$a.gff | gzip > $B/$a.gff.gz
    echo "$a: $(zcat $B/$a.gff.gz | awk -F'\t' '$3=="CDS"' | wc -l) CDS"
    rm -f $O/$a.fna
done
