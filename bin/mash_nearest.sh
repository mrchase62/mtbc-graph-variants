#!/usr/bin/env bash
# Nearest-neighbour distance for every assembly, by mash.
#
# WHAT IT IS FOR. Two genomes that are nearly identical are either a genuine
# clonal pair or the same submission twice, and the QC master table carries the
# distance so that a panel decision can take account of it. It is a screening
# distance, not a phylogeny: mash compares k-mer sketches, so it is fast and
# approximate and should never be used where a SNP distance is meant.
#
# WHY THIS SCRIPT EXISTS. data/qc/mash_nearest.tsv was produced ad hoc and had
# no script in the repository, which was found while inventorying segment 1 of
# PIPELINE_SEGMENTS.md. The table was therefore not reproducible and its
# parameters were unrecorded. They are now: sketch size 10000, k-mer 21.
#
# SCOPE. It runs over every assembly present, not over the current panel, so
# the output covers the whole candidate set and stays valid as the panel
# changes. The earlier table covered only the 490-genome round.
#
#   bash bin/mash_nearest.sh [assembly-dir] [out.tsv]
set -euo pipefail
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { source "$_c"; break; }
done

ASM="${1:-data/assemblies}"
OUT="${2:-data/qc/mash_nearest.tsv}"
MASH="${MTB_MASH:-$(command -v mash || echo "$HOME/bin/mash")}"
SKETCH_SIZE="${SKETCH_SIZE:-10000}"
KMER="${KMER:-21}"

[[ -x "$MASH" ]] || { echo "FATAL: mash not found at ${MASH}" >&2; exit 1; }
[[ -d "$ASM" ]]  || { echo "FATAL: no assembly directory ${ASM}" >&2; exit 1; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/mashnn.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

ls "$ASM"/*.fna.gz > "${WORK}/list" 2>/dev/null || \
    { echo "FATAL: no *.fna.gz under ${ASM}" >&2; exit 1; }
N=$(wc -l < "${WORK}/list")
echo "[mash] sketching ${N} assemblies, s=${SKETCH_SIZE} k=${KMER}"
"$MASH" sketch -s "$SKETCH_SIZE" -k "$KMER" -o "${WORK}/all" -l "${WORK}/list" 2>&1 \
    | tail -2

echo "[mash] all-against-all distances"
"$MASH" dist "${WORK}/all.msh" "${WORK}/all.msh" > "${WORK}/dist.tsv"

# mash prints reference<TAB>query<TAB>distance<TAB>p<TAB>shared; take the
# smallest non-self distance for each query.  Accession is the basename up to
# the first dot, which is how every other QC table keys.
mkdir -p "$(dirname "$OUT")"
awk -F'\t' '
  function acc(p,  n,b) { n=split(p,b,"/"); split(b[n],c,"."); return c[1] }
  { q=acc($2); r=acc($1); if (q==r) next
    if (!(q in best) || $3+0 < best[q]+0) { best[q]=$3; who[q]=r } }
  END { printf "accession\tnearest_mash_dist\tnearest_neighbour\n"
        for (q in best) printf "%s\t%.8f\t%s\n", q, best[q], who[q] }
' "${WORK}/dist.tsv" | { read -r h; echo "$h"; sort; } > "$OUT"

echo "[mash] wrote $(( $(wc -l < "$OUT") - 1 )) rows to ${OUT}"
