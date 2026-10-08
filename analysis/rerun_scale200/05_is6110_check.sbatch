#!/usr/bin/env bash
#SBATCH -J fix1_is6110chk
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 0-04:00
#SBATCH -o /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/rerun_scale200/logs/%x_%j.out
# Build 7713a8d71d8e-fix1: do the shared IS6110 panel assets
# (is6110/assets/matched_gff, isclean_matched; made 2026-09-18..22) reproduce
# from this build's references with today's code? Both stages are run into a
# separate folder and every product compared with the shared copy. The shared
# folders are only read.
set -euo pipefail
R=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants
cd "$R/runroot"
source "$R/config/project_env.sh"
export MTB_BUILD_DIR=$MTB_WORK/refbias/build/7713a8d71d8e-fix1
C=$MTB_WORK/refbias/build_inputs/7713a8d71d8e-fix1/is6110_check
mkdir -p "$C"
OUTDIR=$C/matched_gff bash is6110/bin/p1i_discover_matched.sh
GFFDIR=$C/matched_gff OUTDIR=$C/isclean_matched bash is6110/bin/p1i_build_matched.sh
S=is6110/assets
n=0; d=0
for f in "$C"/matched_gff/*.is6110.gff "$C"/isclean_matched/*.crossmap.tsv \
         "$C"/isclean_matched/*.isclean.fasta "$C"/isclean_matched/manifest.tsv; do
    rel="${f#$C/}"; n=$((n+1))
    cmp -s "$f" "$S/$rel" || { d=$((d+1)); echo "DIFFERS: $rel"; }
done
echo "compared ${n} files, ${d} differ"
