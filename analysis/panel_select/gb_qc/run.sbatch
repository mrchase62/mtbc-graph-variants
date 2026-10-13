#!/bin/bash
# The CX333 selection screens the GenBank scan skipped (collinearity,
# provenance, foreign DNA), run on the 36 GenBank genomes in the top-200 panel,
# after dnaA rotation. Analysis only; repo code (docs/PANEL_BUILD.md 1c-1e).
# Sequencing technology is not queried from NCBI (metadata_unknown.tsv), so
# provenance decides on the IS6110 profile (REFERENCE_STRUCTURE) and reports
# insertions_ge50 / tech_from_sequence for review.
#   cd analysis/panel_select/gb_qc && sbatch run.sbatch
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -J gb_qc
#SBATCH -o logs/gb_qc_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
R=../../..
source $R/config/project_env.sh
G=$(readlink -f ../../genbank_scan/genomes)
mkdir -p rotated rotated_gz
while read -r a; do
    PANEL=$PWD/accessions.txt ASM_DIR=$G OUT_DIR=$PWD/rotated bash $R/bin/rotate_to_dnaa.sh "$a"
    gzip -c rotated/$a.dnaA_rotated.fasta > rotated_gz/$a.fna.gz
done < accessions.txt
cat rotated/verify/*.tsv > rotation_verify.tsv
REF=$MTB_DATA/ref/H37Rv.fasta
$MTB_PY $R/bin/collinearity_qc.py --ref $REF --dir rotated_gz --minimap2 "$MTB_MINIMAP2" \
    --threads 8 --out collinearity_qc.tsv
$MTB_PY $R/bin/assembly_provenance_screen.py --accessions accessions.txt --assembly-dir rotated \
    --ref $REF --is6110 $MTB_DATA/annotation/IS6110.fasta --metadata metadata_unknown.tsv \
    --minimap2 "$MTB_MINIMAP2" --k8 "$MTB_K8" --paftools "$MTB_PAFTOOLS" --out assembly_provenance.tsv
$MTB_PY $R/bin/foreign_insertion_screen.py --accessions accessions.txt --assembly-dir rotated \
    --ref $REF --background ../../external_assemblies/qc/foreign_background.fasta.gz \
    --is6110 $MTB_DATA/annotation/IS6110.fasta --univec "$MTB_UNIVEC" \
    --blastn "$MTB_BLASTN" --makeblastdb "$MTB_MAKEBLASTDB" \
    --minimap2 "$MTB_MINIMAP2" --k8 "$MTB_K8" --paftools "$MTB_PAFTOOLS" \
    --threads 8 --workdir work_foreign --out foreign_insertions.tsv
echo "[gb_qc] done"
