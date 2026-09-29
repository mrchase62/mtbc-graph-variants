#!/usr/bin/env bash
#SBATCH --job-name=P1i_build
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 0-12:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH --output=slurm/P1ibuild_%A.out
#SBATCH --error=slurm/P1ibuild_%A.err
# Stage 2 of P1i: one IS-clean build per distinct SNP-matched reference.
#
#   sbatch is6110/bin/p1i_build_matched.sh     # or bash, it runs either way
#
# WHY THERE IS AN #SBATCH BLOCK AT ALL. There was none, so sbatch refused the
# script outright and it had to be run with `bash`. That is not just an
# inconvenience: a waiter watching for the job id reported BUILD_DONE although
# nothing had been submitted, because the submission had failed and the wait
# had nothing to wait for. The resource request matches its two siblings,
# p1i_discover_matched.sh and p1i_matched.sh, which ask for four tasks and
# 8 GB on `shared`.
#
# WHY TWELVE HOURS AND NOT ONE. The loop below is SERIAL over every distinct
# reference in the refmap -- this is deliberately not an array, because each
# task of an array would redo the whole set and race the others, which
# happened once and cost 22 cancelled tasks. Per reference the build is fast, a
# median of 2 s and a worst case of 126 s measured across the 97 builds that
# exist, but a cohort large enough to select the whole 333-reference panel
# would be 11.6 h at that worst case. Twelve hours covers the panel; a cohort
# of a hundred selects about 61 references and finishes in minutes.
#
# 18 builds, not 23: five references serve two pilot isolates each. Each build
# excises the intervals is6110_discover_elements.py found (stage 1, graded in
# is6110/docs/P1I_STAGE1.md), appends the canonical element as its own contig,
# writes a crossmap and indexes the result for bwa.
#
# WHY THE INTERVALS ARE DISCOVERED AND NOT READ. Only H37Rv carries annotated
# element spans. The panel assemblies' GFFs carry the transposase as CDS, and
# the reading frame sits 42-51 bp inside each element boundary -- those bases
# are the two 28 bp inverted repeats a junction read clips against, so excising
# the reading frame would leave them in the chromosome.
#
# NOTE ON GATE 1. Stage 1 returns 14 of 16 boundaries base-exact on H37Rv; the
# two misses are copies longer than the 1355 bp query, leaving 22 bp of 21,700.
# That gate was written as no-go and these builds proceed past it by decision,
# not because it passed. See is6110/docs/P1I_STAGE1.md section 4.
set -euo pipefail

for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
source "$_mtb_env"
cd "${SLURM_SUBMIT_DIR:-.}"

REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
REFS="${REFS:-refbias/build/7713a8d71d8e/refs}"
GFFDIR="${GFFDIR:-is6110/assets/matched_gff}"
OUTDIR="${OUTDIR:-is6110/assets/isclean_matched}"
ELEMENT="${ELEMENT:-is6110/assets/IS6110.query.fasta}"
BWA="${MTB_BWA:-${MTB_QC_BIN}/bwa}"
mkdir -p "$OUTDIR"

MAN="${OUTDIR}/manifest.tsv"
printf 'reference\tcontig\tintervals\tbp_removed\tchrom_clean_bp\n' > "$MAN"

for R in $(awk -F'\t' 'NR>1{print $5}' "$REFMAP" | sort -u); do
    REF="${REFS}/${R}.fasta"
    GFF="${GFFDIR}/${R}.is6110.gff"
    OUT="${OUTDIR}/${R}.isclean.fasta"
    XMAP="${OUTDIR}/${R}.crossmap.tsv"
    [[ -s "$REF" ]] || { echo "FATAL: missing ${REF}" >&2; exit 1; }
    [[ -s "$GFF" ]] || { echo "FATAL: missing ${GFF} -- run stage 1 first" >&2; exit 1; }
    CONTIG="$(awk '/^>/{print substr($1,2); exit}' "$REF")"

    if [[ -s "$OUT" && -s "${OUT}.bwt" && -s "$XMAP" ]]; then
        echo "[P1i-build] ${R}: already built"
        # still emit the manifest row. It used to be written only on the build
        # path, so after any incremental run the manifest recorded the last
        # batch rather than the set of builds -- 12 rows for 61 builds.
        N=$(grep -vc '^#' "$GFF" || true)
        BP=$(awk -F'\t' '!/^#/{s+=$5-$4+1} END{print s+0}' "$GFF")
        CL=$("$MTB_SAMTOOLS" faidx "$OUT" "${CONTIG}_isclean" 2>/dev/null \
             | grep -v '^>' | tr -d '\n' | wc -c)
        printf '%s\t%s\t%s\t%s\t%s\n' "$R" "$CONTIG" "$N" "$BP" "$CL" >> "$MAN"
        continue
    fi
    echo "[P1i-build] ${R} (${CONTIG})"
    "$MTB_PY" is6110/bin/is6110_build_isclean.py \
        --reference "$REF" --contig "$CONTIG" --gff "$GFF" \
        --element "$ELEMENT" --out-fasta "${OUT}.tmp" \
        --out-crossmap "${XMAP}.tmp" | sed 's/^/    /'
    mv -f "${OUT}.tmp" "$OUT"; mv -f "${XMAP}.tmp" "$XMAP"
    "$BWA" index "$OUT" 2>/dev/null
    "$MTB_SAMTOOLS" faidx "$OUT"

    # `grep -c` exits 1 when the count is zero, which under `set -e` killed the
    # whole loop on the one reference that carries no IS6110 copies -- after it
    # had been built and indexed correctly. A genome with none is legitimate,
    # so count without letting the count's own exit status end the run.
    N=$(grep -vc '^#' "$GFF" || true)
    BP=$(awk -F'\t' '!/^#/{s+=$5-$4+1} END{print s+0}' "$GFF")
    CL=$("$MTB_SAMTOOLS" faidx "$OUT" "${CONTIG}_isclean" 2>/dev/null | grep -v '^>' | tr -d '\n' | wc -c)
    printf '%s\t%s\t%s\t%s\t%s\n' "$R" "$CONTIG" "$N" "$BP" "$CL" >> "$MAN"
done
echo "written: ${MAN}"
