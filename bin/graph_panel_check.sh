#!/usr/bin/env bash
# Does the graph contain exactly the panel it was built from?
#
# WHY THIS EXISTS. Nothing checked it. Every coordinate this project produces
# downstream of the graph is a position on one of its paths, and graph node ids
# are build-scoped, so a graph that silently lost or gained a path would
# invalidate work without any error appearing. Found as a gap while
# inventorying segment 2 of PIPELINE_SEGMENTS.md.
#
# It compares path names from `odgi paths -L` against the sequence names in the
# panel FASTA, as SETS and not as counts -- a count check would pass a graph
# that dropped one genome and duplicated another.
#
#   bash bin/graph_panel_check.sh [graph.og] [panel.fasta.gz]
set -euo pipefail
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { source "$_c"; break; }
done

GRAPH="${1:-graphs/CX333.s10k.k23.K15/mtb.complex333.fasta.gz.f4f5ee2.11fba48.36e68b6.smooth.final.og}"
PANEL="${2:-data/fastas/mtb.complex333.fasta.gz}"
ODGI="${MTB_ODGI:-odgi}"

[[ -s "$GRAPH" ]] || { echo "FATAL: no graph at ${GRAPH}" >&2; exit 1; }
[[ -s "$PANEL" ]] || { echo "FATAL: no panel at ${PANEL}" >&2; exit 1; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/gpc.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

"$ODGI" paths -i "$GRAPH" -L | sort > "${WORK}/graph.paths"
zcat -f "$PANEL" | awk '/^>/{print substr($1,2)}' | sort > "${WORK}/panel.paths"

NG=$(wc -l < "${WORK}/graph.paths"); NP=$(wc -l < "${WORK}/panel.paths")
echo "graph : ${GRAPH}"
echo "panel : ${PANEL}"
echo "paths : ${NG} in the graph, ${NP} in the panel"

# duplicates would make a count check pass while the sets differ
DG=$(sort "${WORK}/graph.paths" | uniq -d | wc -l)
DP=$(sort "${WORK}/panel.paths" | uniq -d | wc -l)
[[ "$DG" -gt 0 ]] && echo "  WARNING: ${DG} duplicated path names in the graph"
[[ "$DP" -gt 0 ]] && echo "  WARNING: ${DP} duplicated sequence names in the panel"

MISSING=$(comm -13 "${WORK}/graph.paths" "${WORK}/panel.paths" | head -20)
EXTRA=$(comm -23 "${WORK}/graph.paths" "${WORK}/panel.paths" | head -20)
NM=$(comm -13 "${WORK}/graph.paths" "${WORK}/panel.paths" | wc -l)
NE=$(comm -23 "${WORK}/graph.paths" "${WORK}/panel.paths" | wc -l)

if [[ "$NM" -eq 0 && "$NE" -eq 0 && "$DG" -eq 0 && "$DP" -eq 0 ]]; then
    echo "PASS: the two sets are identical"
    exit 0
fi
echo "FAIL"
[[ "$NM" -gt 0 ]] && { echo "  ${NM} in the panel but NOT in the graph:"; echo "$MISSING" | sed 's/^/    /'; }
[[ "$NE" -gt 0 ]] && { echo "  ${NE} in the graph but NOT in the panel:"; echo "$EXTRA" | sed 's/^/    /'; }
exit 1
