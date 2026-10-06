#!/usr/bin/env bash
# Refresh the MTBC complete-genome set from NCBI RefSeq.
#
#   bin/refresh_assemblies.sh [--summary-only] [--outdir DIR] [--summary FILE]
#
# --summary FILE re-applies the selection to an archived assembly summary
# (e.g. data/ncbi/assembly_summary_refseq.20260907.txt) instead of fetching.
#
# The 2025 panel was assembled without a download script, so the selection could
# not be reproduced. This records it.
#
# SELECTION CRITERION, stated explicitly:
#   assembly_level  == "Complete Genome"
#   version_status  == "latest"
#   species_taxid   in MTBC_SPECIES_TAXIDS (default 1773 M. tuberculosis,
#                   78331 M. canettii, 1305738 M. orygis)
#     OR organism_name matches ^Mycobacterium (tuberculosis|canetti|orygis)
#
# The name regex alone failed twice: "canettii" with two t's dropped the
# outgroup GCF_035581225 (RefSeq spells it "canetti"), and
# ^Mycobacterium (tuberculosis|canetti) dropped both complete M. orygis genomes,
# GCF_015265495 and GCF_033782915 (audit PGB-13). The taxid alone fails too:
# the six "Mycobacterium tuberculosis complex sp. N0052 ..." genomes each carry
# their own species taxid. So both are used. Every other MTBC member in RefSeq
# (bovis, africanum, microti, caprae, pinnipedii, mungi) is filed under
# species_taxid 1773 as "M. tuberculosis variant ...". After selecting, the
# script lists any complete genome whose name looks like an MTBC member but was
# not selected, so a new species name cannot be dropped silently.
#
# On the 2026-09-07 summary this selects 516: the 514 of panel_manifest.tsv
# plus the two M. orygis genomes.
#
# NOT covered: GenBank-only complete genomes (assembly_summary_genbank.txt);
# this script reads RefSeq only.
#
# The criterion matches 468 assemblies in the 2025-10-22 summary but only 416 were
# built. The additional filter that removed those 52 was never recorded and is NOT
# recoverable -- the collinearity QC table does not explain it (excluded genomes
# carry the same OK/CHECK/WARNING mix as the built ones). So this script defines a
# panel going forward; it does not reproduce the historical one.
set -euo pipefail
_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh" \
          "${MTB_WORK:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: no project_env.sh" >&2; exit 1; }
source "$_mtb_env"

SUMMARY_ONLY=0
OUTDIR="${MTB_DATA}/assemblies"
LOCAL_SUM=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --summary-only) SUMMARY_ONLY=1; shift ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --summary) LOCAL_SUM="$2"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

NCBI="${MTB_DATA}/ncbi"
mkdir -p "$NCBI" "$OUTDIR"
SUM="${NCBI}/assembly_summary_refseq.txt"
STAMP=$(date -u +%Y%m%d)

if [[ -n "$LOCAL_SUM" ]]; then
    [[ -s "$LOCAL_SUM" ]] || { echo "FATAL: no summary $LOCAL_SUM" >&2; exit 1; }
    SUM="$LOCAL_SUM"
    echo "=== using archived RefSeq assembly summary $SUM"
else
    echo "=== fetching RefSeq assembly summary"
    curl -fsS -o "${SUM}.new" \
        https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
    mv "${SUM}.new" "$SUM"
    cp -n "$SUM" "${NCBI}/assembly_summary_refseq.${STAMP}.txt" 2>/dev/null || true
fi
echo "    $(wc -l < "$SUM") rows"

echo "=== applying the selection criterion"
MTBC_SPECIES_TAXIDS="${MTBC_SPECIES_TAXIDS:-1773 78331 1305738}"
awk -F'\t' -v tx="$MTBC_SPECIES_TAXIDS" '
    BEGIN {n = split(tx, t, " "); for (i = 1; i <= n; i++) keep[t[i]] = 1}
    !/^#/ && $12=="Complete Genome" && $11=="latest" &&
    (($7 in keep) || $8 ~ /^Mycobacterium (tuberculosis|canetti|orygis)/) \
    {split($1,a,"."); print a[1]"\t"$1"\t"$8"\t"$20}' "$SUM" \
    | sort -u > "${NCBI}/selected.${STAMP}.tsv"
echo "    $(wc -l < "${NCBI}/selected.${STAMP}.tsv") assemblies selected" \
     "(species taxids: ${MTBC_SPECIES_TAXIDS})"
# names that look like MTBC members but were not selected: review, never drop silently
awk -F'\t' -v tx="$MTBC_SPECIES_TAXIDS" '
    BEGIN {n = split(tx, t, " "); for (i = 1; i <= n; i++) keep[t[i]] = 1}
    !/^#/ && $12=="Complete Genome" && $11=="latest" && !($7 in keep) &&
    $8 !~ /^Mycobacterium (tuberculosis|canetti|orygis)/ &&
    $8 ~ /^Mycobacterium (tubercul|canett|orygis|bovis|africanum|microti|caprae|pinnipedii|mungi|suricattae)/ \
    {print "    NOT SELECTED, check taxid " $7 ": " $1 "  " $8}' "$SUM"

MANIFEST="${NCBI}/panel_manifest.tsv"
if [[ -f "$MANIFEST" ]]; then
    comm -13 <(cut -f1 "$MANIFEST" | sort -u) \
             <(cut -f1 "${NCBI}/selected.${STAMP}.tsv" | sort -u) > "${NCBI}/added.${STAMP}.txt"
    comm -23 <(cut -f1 "$MANIFEST" | sort -u) \
             <(cut -f1 "${NCBI}/selected.${STAMP}.tsv" | sort -u) > "${NCBI}/removed.${STAMP}.txt"
    echo "    vs manifest: $(wc -l < "${NCBI}/added.${STAMP}.txt") added, "\
"$(wc -l < "${NCBI}/removed.${STAMP}.txt") gone"
else
    cut -f1 "${NCBI}/selected.${STAMP}.tsv" > "${NCBI}/added.${STAMP}.txt"
    echo "    no manifest yet: treating all as new"
fi

[[ $SUMMARY_ONLY -eq 1 ]] && { echo "=== --summary-only, stopping"; exit 0; }

echo "=== downloading genomes to ${OUTDIR}"
n=0; skip=0; fail=0
while read -r acc; do
    [[ -n "$acc" ]] || continue
    line=$(awk -F'\t' -v a="$acc" '$1==a' "${NCBI}/selected.${STAMP}.tsv" | head -1)
    ftp=$(cut -f4 <<<"$line"); full=$(cut -f2 <<<"$line")
    [[ -n "$ftp" && "$ftp" != "na" ]] || { fail=$((fail+1)); continue; }
    out="${OUTDIR}/${acc}.fna.gz"
    [[ -s "$out" ]] && { skip=$((skip+1)); continue; }
    url="${ftp}/$(basename "$ftp")_genomic.fna.gz"
    if curl -fsS -o "${out}.part" "$url"; then mv "${out}.part" "$out"; n=$((n+1))
    else rm -f "${out}.part"; echo "    FAILED $acc" >&2; fail=$((fail+1)); fi
done < "${NCBI}/added.${STAMP}.txt"
echo "    downloaded $n, already present $skip, failed $fail"

cp "${NCBI}/selected.${STAMP}.tsv" "$MANIFEST"
echo "=== manifest updated -> $MANIFEST"
echo
echo "NOT done by this script, and each needs a decision:"
echo "  - rotation to a common origin (2025 used rotated_assemblies/)"
echo "  - collinearity/synteny QC and what to exclude on it"
echo "  - whether to keep lab-derived strains (BCG, H37Ra, mc2 6030, Erdman)"
echo "  - concatenation into a PanSN fasta.gz for pggb"
