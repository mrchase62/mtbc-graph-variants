#!/usr/bin/env bash
# Fetch NCBI's UniVec_Core for the foreign screen's vector check (D29).
#
#   bash bin/fetch_univec.sh [DEST_DIR]     default: <repo>/data/univec
#
# Writes UniVec_Core, NCBI's README.uv and UniVec_Core.sha256 (with the build
# named in the README), and refuses to overwrite a copy of another build, so
# a run's inputs.key can be traced to one release. Fetched 2026-10-07:
# build 10.0, 3,155 sequences, sha256 ea5d1524...92c5c.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${1:-${ROOT}/data/univec}"
URL="https://ftp.ncbi.nlm.nih.gov/pub/UniVec"
mkdir -p "$DEST"
TMP="$(mktemp -d "${DEST}/.fetch.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
curl -sS -f -o "${TMP}/UniVec_Core" "${URL}/UniVec_Core"
curl -sS -f -o "${TMP}/README.uv" "${URL}/README.uv"
BUILD="$(awk '$1=="UniVec_Core"{print $2, $3; exit}' "${TMP}/README.uv")"
[[ -n "$BUILD" ]] || { echo "FATAL: no UniVec_Core build in README.uv" >&2; exit 1; }
N="$(grep -c '^>' "${TMP}/UniVec_Core")"
SHA="$(sha256sum "${TMP}/UniVec_Core" | cut -d' ' -f1)"
if [[ -s "${DEST}/UniVec_Core.sha256" ]]; then
    OLD="$(cut -f1 "${DEST}/UniVec_Core.sha256")"
    if [[ "$OLD" != "$SHA" ]]; then
        echo "FATAL: ${DEST} holds another UniVec_Core ($(cut -f2 "${DEST}/UniVec_Core.sha256"));" \
             "NCBI now serves ${BUILD} (${SHA}). Fetch into a new directory." >&2
        exit 1
    fi
fi
mv -f "${TMP}/UniVec_Core" "${TMP}/README.uv" "$DEST/"
printf '%s\t%s\t%s sequences\t%s\n' "$SHA" "$BUILD" "$N" "$(date +%F)" > "${DEST}/UniVec_Core.sha256"
echo "UniVec_Core ${BUILD}: ${N} sequences, sha256 ${SHA} -> ${DEST}"
