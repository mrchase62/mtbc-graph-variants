#!/usr/bin/env bash
#
# Stamp the P0 graph-build provenance into a VCF header.
#
#   bash bin/stamp_build_id.sh out.vcf.gz [more.vcf.gz ...]
#   MTB_BUILD_DIR=refbias/build/<id> bash bin/stamp_build_id.sh out.vcf.gz
#
# P0 names its asset directory after the graph's checksum, so ASSETS from
# different graph builds cannot mix. That is only half of the guarantee: a VCF
# produced against those assets does not itself say which build made it, and the
# graph build changes because bugs are found rather than by choice. Without the
# stamp, output from a superseded graph is indistinguishable from output from the
# fixed one, and the distinction is unrecoverable after the fact -- which is why
# this is worth doing before any output is produced at scale rather than after.
#
# Idempotent: a VCF already carrying ##MTB_graph_build is left alone, so
# re-running a pipeline stage does not accumulate duplicate header lines.
#
# Resolution order for the build directory: $MTB_BUILD_DIR, else the single
# directory under $BUILD_ROOT if there is exactly one. If no build directory can
# be resolved the file is left UNSTAMPED and this exits 0 with a notice, so that
# existing harnesses which predate P0 keep working. If a build directory IS
# resolved and stamping then fails, that is an error and exits non-zero: a
# provenance step that fails quietly is worse than one that is absent.
set -euo pipefail

_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: no project_env.sh" >&2; exit 1; }
source "$_mtb_env"

# Real binaries, not the container fallbacks: the pggb container path in
# project_env.sh resolves under the project tree and is not present here, so
# mtb_bcftools would die on a missing image rather than reheader anything.
BCFTOOLS="${MTB_BCFTOOLS:-bcftools}"
TABIX="${MTB_TABIX:-${MTB_QC_BIN}/tabix}"
[[ -x "$BCFTOOLS" ]] || { echo "FATAL: bcftools not found at ${BCFTOOLS}" >&2; exit 1; }

BUILD_ROOT="${BUILD_ROOT:-refbias/build}"
BUILD="${MTB_BUILD_DIR:-}"
if [[ -z "$BUILD" ]]; then
    mapfile -t _cands < <(find "$BUILD_ROOT" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
    if [[ "${#_cands[@]}" -eq 1 ]]; then
        BUILD="${_cands[0]}"
    elif [[ "${#_cands[@]}" -gt 1 ]]; then
        echo "[stamp] ${#_cands[@]} builds under ${BUILD_ROOT}; set MTB_BUILD_DIR" >&2
        exit 1
    fi
fi
if [[ -z "$BUILD" || ! -s "${BUILD}/build_info.tsv" ]]; then
    echo "[stamp] no P0 build directory resolved; leaving $# file(s) unstamped" >&2
    exit 0
fi

_info() { awk -F'\t' -v k="$1" '$1==k{print $2; exit}' "${BUILD}/build_info.tsv"; }
BUILD_ID="$(_info build_id)"
GRAPH="$(_info graph)"
GRAPH_SHA="$(_info graph_sha256)"
[[ -n "$BUILD_ID" ]] || { echo "[stamp] ${BUILD}/build_info.tsv has no build_id" >&2; exit 1; }
MANIFEST_SHA="-"
[[ -s "${BUILD}/manifest.tsv" ]] && \
    MANIFEST_SHA="$(sha256sum "${BUILD}/manifest.tsv" | cut -d' ' -f1)"

for vcf in "$@"; do
    [[ -s "$vcf" ]] || { echo "[stamp] missing: $vcf" >&2; exit 1; }
    if "$BCFTOOLS" view -h "$vcf" 2>/dev/null | grep -q '^##MTB_graph_build='; then
        echo "[stamp] already stamped: $vcf"
        continue
    fi
    hdr="${vcf}.hdr.$$"
    new="${vcf}.stamped.$$"
    "$BCFTOOLS" view -h "$vcf" > "$hdr"
    # Insert before #CHROM. Header lines must precede it, and appending after it
    # would produce a VCF that parsers silently mis-read rather than reject.
    awk -v id="$BUILD_ID" -v g="$GRAPH" -v gs="$GRAPH_SHA" \
        -v ms="$MANIFEST_SHA" -v bd="$BUILD" '
        /^#CHROM/ && !done {
            printf "##MTB_graph_build=%s\n", id
            printf "##MTB_graph=%s\n", g
            printf "##MTB_graph_sha256=%s\n", gs
            printf "##MTB_p0_build_dir=%s\n", bd
            printf "##MTB_p0_manifest_sha256=%s\n", ms
            done = 1
        }
        { print }
    ' "$hdr" > "${hdr}.new"
    grep -q '^##MTB_graph_build=' "${hdr}.new" \
        || { rm -f "$hdr" "${hdr}.new"; echo "[stamp] no #CHROM line in $vcf" >&2; exit 1; }
    "$BCFTOOLS" reheader -h "${hdr}.new" -o "$new" "$vcf"
    mv -f "$new" "$vcf"
    rm -f "$hdr" "${hdr}.new"
    if [[ -s "${vcf}.tbi" && -x "$TABIX" ]]; then "$TABIX" -f -p vcf "$vcf"; fi
    echo "[stamp] ${vcf}: build ${BUILD_ID}"
done
