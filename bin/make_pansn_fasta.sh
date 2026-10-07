#!/usr/bin/env bash
# Concatenate the rotated panel into one PanSN-named FASTA for pggb.
#
#   bin/make_pansn_fasta.sh --panel FILE --out NAME [--rotated DIR]
#
# --panel is build_panel.py's --out list (header, accession in column 1); --out
# is the bare name, written to ${MTB_FASTAS}/NAME.fasta.gz. Both are required:
# the old defaults (panel.rebuild.tsv, the 484-genome first round, and
# mtb.complex490) built a panel no decision ever chose. An existing
# NAME.fasta.gz is never overwritten (audit PGB-3): the production panel
# mtb.complex333.fasta.gz is the graph's input and must stay as built.
#
# PanSN naming is sample#haplotype#contig, e.g. GCF_000195955#1#NC_000962.3.
# pggb splits on the first '#' (via --exclude-delim '#') to group haplotypes, and
# vg deconstruct later reduces a path name to its SAMPLE field -- which is why
# every downstream join in this project is on the GCF accession.
#
# Haplotype is always 1: these are single-contig complete bacterial genomes, one
# chromosome each, verified before concatenation.
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

PANEL=""
ROT="${MTB_DATA}/rotated"
NAME=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --panel) PANEL="$2"; shift 2 ;;
        --rotated) ROT="$2"; shift 2 ;;
        --out) NAME="$2"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$PANEL" && -n "$NAME" ]] || { echo "usage: $0 --panel FILE --out NAME [--rotated DIR]" >&2; exit 2; }
mtb_require_file "$PANEL"
BG="${MTB_QC_BIN}/bgzip"
mtb_require_file "$BG"

OUT="${MTB_FASTAS}/${NAME}.fasta"
mkdir -p "$MTB_FASTAS"
if [[ -e "${OUT}.gz" ]]; then
    echo "FATAL: ${OUT}.gz exists; refusing to overwrite a built panel. Choose a new --out name." >&2
    exit 1
fi
: > "$OUT"
n=0; missing=0; multi=0
while read -r acc _rest; do
    f="${ROT}/${acc}.dnaA_rotated.fasta"
    if [[ ! -s "$f" ]]; then echo "  MISSING rotated: $acc" >&2; missing=$((missing+1)); continue; fi
    c=$(grep -c '^>' "$f")
    if [[ "$c" -ne 1 ]]; then echo "  SKIP multi-contig ($c): $acc" >&2; multi=$((multi+1)); continue; fi
    contig=$(head -1 "$f" | sed 's/^>//' | awk '{print $1}')
    # header becomes sample#1#contig; sequence is passed through untouched
    awk -v h=">${acc}#1#${contig}" 'NR==1{print h; next} {print}' "$f" >> "$OUT"
    n=$((n+1))
done < <(tail -n+2 "$PANEL")

echo "[pansn] wrote ${n} sequences (${missing} missing, ${multi} multi-contig skipped)"
"$BG" -f -@ "${SLURM_CPUS_PER_TASK:-8}" "$OUT"
"$MTB_SAMTOOLS" faidx "${OUT}.gz"

dup=$(cut -f1 "${OUT}.gz.fai" | sort | uniq -d | wc -l)
echo "[pansn] paths: $(wc -l < "${OUT}.gz.fai")   duplicate names: ${dup}"
awk '{s+=$2} END{printf "[pansn] total %.1f Mb, mean %.0f bp, min %d, max %d\n", s/1e6, s/NR, m?m:$2, x}
     {if(!m||$2<m)m=$2; if($2>x)x=$2}' "${OUT}.gz.fai"
if grep -qx "$(printf '%s\t[0-9]*' "$MTB_REF_PATH" | cut -f1)" <(cut -f1 "${OUT}.gz.fai"); then
    echo "[pansn] reference path ${MTB_REF_PATH}: PRESENT"
else
    echo "[pansn] WARNING: reference path ${MTB_REF_PATH} NOT FOUND -- pggb and every"
    echo "        downstream step reference this exact name" >&2
fi
echo "[pansn] -> ${OUT}.gz"
