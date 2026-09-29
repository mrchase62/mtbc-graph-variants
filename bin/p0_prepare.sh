#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH -J p0_prepare
#SBATCH -o slurm/%x_%A_%a.out
#SBATCH -e slurm/%x_%A_%a.err
#
# P0: per-graph preprocessing. Build every asset the per-sample stages need,
# once per graph build, and stamp a build id so a rebuild is detectable.
#
#   bash bin/p0_prepare.sh                     # all cheap steps, in order
#   bash bin/p0_prepare.sh --step gff          # one step
#   bash bin/p0_prepare.sh --list              # what exists and what does not
#   sbatch --array=1-333 bin/p0_prepare.sh --step refs   # the expensive step
#
# Why a separate stage. Assets belong to a graph BUILD, not to a sample. The
# graph build changes because bugs are found, not because anyone chose to version
# it, so what matters is invalidation and provenance: knowing which build a result
# came from, and being able to regenerate everything for a fixed build without
# touching per-sample logic. Two failures this is meant to prevent:
#
#   - a per-sample task re-deriving a per-graph asset. T18 ran a 24-task array to
#     get sub-lineage that one column of an existing table already held. Asset
#     derivation belongs here, once, where it is visible.
#   - output from a superseded graph being pooled with output from the fixed one.
#     Every asset directory is named by the build id, so they cannot mix silently.
#
# Steps are idempotent: each writes a .done marker and is skipped if present.
# Deleting a marker re-runs that step alone.
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

GRAPH_DIR="${GRAPH_DIR:-graphs/CX333.s10k.k23.K15}"
OG="${OG:-$(ls "${GRAPH_DIR}"/*.smooth.final.og 2>/dev/null | head -1)}"
PANEL_SNPS="${PANEL_SNPS:-${GRAPH_DIR}/snps.vcf.gz}"
H37RV_PATH="${H37RV_PATH:-GCF_000195955#1#NC_000962.3}"
H37RV_ACC="${H37RV_ACC:-GCF_000195955}"
ASM_DIR="${ASM_DIR:-data/assemblies}"
NCBI_SUM="${NCBI_SUM:-data/ncbi/assembly_summary_refseq.txt}"
MASK="${MASK:-data/annotation/H37Rv_repeat_mask.bed}"
# The FASTA the graph was built from. Its sequences are dnaA-rotated, so for a
# minority of accessions the graph's coordinate frame differs from refs/; the
# frames step measures that difference once per build.
PANEL_FASTA="${PANEL_FASTA:-data/fastas/mtb.complex333.fasta.gz}"
BWA="${MTB_BWA:-${MTB_QC_BIN}/bwa}"
GATK_SIF="${MTB_GATK_SIF:?MTB_GATK_SIF is unset; see config/project_env.sh}"
gatk_run() { singularity exec -B "${MTB_WORK}:${MTB_WORK}" "$GATK_SIF" gatk "$@"; }
ACCESSORY_DIR="${ACCESSORY_DIR:-refbias/panel}"
BUILD_ROOT="${BUILD_ROOT:-refbias/build}"
STEP=""
LIST=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --step) STEP="$2"; shift 2 ;;
        --list) LIST=1; shift ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

[[ -s "$OG" ]] || { echo "FATAL: no graph .og found under ${GRAPH_DIR}" >&2; exit 1; }

# --- build id -----------------------------------------------------------------
# Derived from the graph's own bytes, so the same graph always yields the same id
# and a rebuilt graph always yields a different one. Not a timestamp: two runs of
# P0 on the same graph must land in the same directory, or idempotence is a lie.
BUILD_ID="${BUILD_ID:-}"
if [[ -z "$BUILD_ID" ]]; then
    BUILD_ID="$(sha256sum "$OG" | cut -c1-12)"
fi
BUILD="${BUILD_ROOT}/${BUILD_ID}"
mkdir -p "$BUILD"/{annotation,refs,assets,logs} slurm

_done()  { [[ -s "${BUILD}/logs/$1.done" ]]; }
_mark()  { date -Is > "${BUILD}/logs/$1.done"; }
_say()   { echo "[P0:${BUILD_ID}] $*"; }

if [[ "$LIST" == 1 ]]; then
    echo "graph    : $OG"
    echo "build id : $BUILD_ID"
    echo "build dir: $BUILD"
    for s in stamp paths accessions gff refs frames assets manifest; do
        if _done "$s"; then printf "  %-12s done   %s\n" "$s" "$(cat "${BUILD}/logs/$s.done")"
        else printf "  %-12s MISSING\n" "$s"; fi
    done
    exit 0
fi

# --- step: stamp --------------------------------------------------------------
step_stamp() {
    _done stamp && { _say "stamp: already done"; return 0; }
    {
        echo "build_id\t${BUILD_ID}"
        echo "graph\t$(cd "$(dirname "$OG")" && pwd -P)/$(basename "$OG")"
        echo "graph_sha256\t$(sha256sum "$OG" | cut -d' ' -f1)"
        echo "graph_bytes\t$(stat -c%s "$OG")"
        echo "created\t$(date -Is)"
        echo "created_by\t${USER:-unknown}"
        # CARRY THE PANEL REPOSITORY'S IDENTITY FORWARD. pggb_build.sh writes
        # graph_provenance.tsv beside the graph recording which repository and
        # commit built it. That is the only place that information exists once
        # panel construction is a separate repository, so it is copied into the
        # build stamp here and from there into every VCF header.
        #
        # When the sidecar is absent -- any graph built before this existed --
        # the keys are NOT invented. An explicit marker is written instead, so a
        # reader can tell "nobody recorded this" from "nobody asked".
        _gprov="$(dirname "$OG")/graph_provenance.tsv"
        if [[ -s "$_gprov" ]]; then
            awk -F'\t' '$1=="panel_repo"||$1=="panel_commit"||$1=="panel_dirty"||$1=="pggb_args"{print $1"\t"$2}' "$_gprov"
        else
            echo "panel_provenance\tabsent:graph_predates_provenance_sidecar"
        fi
    } | sed 's/\\t/\t/g' > "${BUILD}/build_info.tsv"
    _say "stamp: ${BUILD}/build_info.tsv"
    _mark stamp
}

# --- step: paths -------------------------------------------------------------
step_paths() {
    _done paths && { _say "paths: already done"; return 0; }
    "$MTB_ODGI" paths -i "$OG" -L > "${BUILD}/assets/paths.txt"
    local n; n=$(wc -l < "${BUILD}/assets/paths.txt")
    [[ "$n" -gt 0 ]] || { echo "FATAL: odgi paths returned nothing" >&2; exit 1; }
    grep -qxF "$H37RV_PATH" "${BUILD}/assets/paths.txt" \
        || { echo "FATAL: H37Rv path ${H37RV_PATH} absent from the graph" >&2; exit 1; }
    _say "paths: ${n} paths, H37Rv path present"
    _mark paths
}

# --- step: accessions --------------------------------------------------------
step_accessions() {
    _done accessions && { _say "accessions: already done"; return 0; }
    # PanSN path names are <accession>#<haplotype>#<contig>
    cut -d'#' -f1 "${BUILD}/assets/paths.txt" | sort -u > "${BUILD}/assets/accessions.txt"
    local n; n=$(wc -l < "${BUILD}/assets/accessions.txt")
    _say "accessions: ${n} distinct panel accessions"
    _mark accessions
}

# --- step: gff ---------------------------------------------------------------
# The off-path annotation source. A record that does not map to the H37Rv path is
# annotated in the MATCHED REFERENCE's own coordinates using that reference's
# RefSeq GFF -- no projection, because the off-path sequence IS the reference's
# sequence. data/assemblies currently holds 333 .fna.gz and zero GFFs, so this is
# the one genuinely missing asset.
step_gff() {
    _done gff && { _say "gff: already done"; return 0; }
    [[ -s "$NCBI_SUM" ]] || { echo "FATAL: no assembly summary at ${NCBI_SUM}" >&2; exit 1; }
    local n=0 skip=0 fail=0 acc ftp url out
    while read -r acc; do
        [[ -n "$acc" ]] || continue
        out="${BUILD}/annotation/${acc}.gff.gz"
        if [[ -s "$out" ]]; then skip=$((skip+1)); continue; fi
        # column 20 is ftp_path; accessions in the summary carry a .N version
        ftp=$(awk -F'\t' -v a="$acc" '$1 ~ "^"a"\\." {print $20; exit}' "$NCBI_SUM")
        if [[ -z "$ftp" || "$ftp" == "na" ]]; then
            echo "${acc}\tno_ftp_path" >> "${BUILD}/logs/gff.failures.tsv"
            fail=$((fail+1)); continue
        fi
        ftp="${ftp%/}"
        url="${ftp}/$(basename "$ftp")_genomic.gff.gz"
        if curl -fsS --retry 3 --retry-delay 2 -o "${out}.part" "$url"; then
            mv "${out}.part" "$out"; n=$((n+1))
        else
            rm -f "${out}.part"
            echo "${acc}\t${url}" >> "${BUILD}/logs/gff.failures.tsv"
            fail=$((fail+1))
        fi
    done < "${BUILD}/assets/accessions.txt"
    _say "gff: downloaded ${n}, already present ${skip}, failed ${fail}"
    # A partial GFF set is not a usable asset: off-path records for a reference
    # with no GFF cannot be annotated at all, so refuse to mark the step done.
    if [[ "$fail" -gt 0 ]]; then
        echo "[P0] ${fail} accessions have no GFF; see ${BUILD}/logs/gff.failures.tsv" >&2
        echo "[P0] step 'gff' NOT marked done -- rerun to retry only the missing ones" >&2
        return 1
    fi
    _mark gff
}

# --- step: refs --------------------------------------------------------------
# Stage each panel genome as plain FASTA with faidx and a bwa index, so P2 never
# builds an index inside the sample loop. Expensive and embarrassingly parallel:
#   sbatch --array=1-$(wc -l < <build>/assets/accessions.txt) bin/p0_prepare.sh --step refs
step_refs() {
    local list="${BUILD}/assets/accessions.txt"
    local accs=()
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        accs=( "$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n' "$list")" )
    else
        mapfile -t accs < "$list"
    fi
    local acc src fa dict lock ext
    for acc in "${accs[@]}"; do
        [[ -n "$acc" ]] || continue
        fa="${BUILD}/refs/${acc}.fasta"
        dict="${fa%.*}.dict"
        src="${ASM_DIR}/${acc}.fna.gz"
        if [[ ! -s "$src" ]]; then
            echo "${acc}	no_assembly" >> "${BUILD}/logs/refs.failures.tsv"
            echo "[P0] ${acc}: no assembly at ${src}" >&2; continue
        fi
        [[ -s "$fa" ]] || zcat "$src" > "$fa"
        [[ -s "${fa}.fai" ]] || "$MTB_SAMTOOLS" faidx "$fa"
        # Race-safe index build, copied from bin/simulate_and_call.sh: bwa index
        # writes fixed sibling filenames, so two tasks on the same FASTA collide.
        # Index a private copy and move the products into place, since mv is
        # atomic and a loser simply overwrites with an identical file.
        if [[ ! -s "${fa}.bwt" ]]; then
            lock="${fa}.bwaindex.$$"
            cp -f "$fa" "$lock"
            "$BWA" index "$lock" > "${BUILD}/logs/bwa_index.${acc}.log" 2>&1
            # bwt LAST: it is what `[[ -s ${fa}.bwt ]]` treats as "indexed"
            for ext in amb ann pac sa bwt; do
                [[ -f "${lock}.${ext}" ]] && mv -f "${lock}.${ext}" "${fa}.${ext}"
            done
            rm -f "$lock"
        fi
        # GATK needs a sequence dictionary; building it here keeps it out of the
        # per-sample loop, which is the whole point of P0.
        if [[ ! -s "$dict" ]]; then
            if gatk_run CreateSequenceDictionary -R "$fa" -O "${dict}.$$" \
                    > "${BUILD}/logs/dict.${acc}.log" 2>&1; then
                mv -f "${dict}.$$" "$dict"
            else
                rm -f "${dict}.$$"
                echo "${acc}	dict_failed" >> "${BUILD}/logs/refs.failures.tsv"
            fi
        fi
    done
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        local want built
        want=$(wc -l < "$list")
        built=$(ls "${BUILD}/refs"/*.bwt 2>/dev/null | wc -l)
        _say "refs: ${built} of ${want} indexed"
        if [[ "$built" -eq "$want" ]]; then _mark refs
        else echo "[P0] step 'refs' NOT marked done -- ${want} expected, ${built} built" >&2; fi
    fi
}

# --- step: frames ------------------------------------------------------------
# The panel-vs-refs coordinate frame of every accession. graph_frame.py reads it
# to convert a refs-frame coordinate before it is handed to odgi, which
# interprets coordinates in the PANEL frame. It is a property of this build's
# graph and refs, so it is measured here and nowhere else; refbias_run.sh points
# MTB_GRAPH_FRAMES at it. Needs the refs step to have finished.
step_frames() {
    _done frames && { _say "frames: already done"; return 0; }
    _done refs || { echo "FATAL: step 'frames' needs step 'refs' done first" >&2; return 1; }
    mtb_require_file "$PANEL_FASTA"
    local out="${BUILD}/assets/graph_frame_offsets.tsv"
    "$MTB_PY" graphframe/bin/graph_frame_offsets.py \
        --panel "$PANEL_FASTA" --refs "${BUILD}/refs" \
        --samtools "$MTB_SAMTOOLS" --out "${out}.tmp"
    mv -f "${out}.tmp" "$out"
    # Refuse to mark done if any accession has no strand: graph_frame.py treats
    # such a row as unknown, which silently drops that accession's projections.
    # A note other than "ok" with a strand is fine -- some probes hit repeats and
    # the unique ones agreed.
    local bad
    bad=$(awk -F'\t' 'NR>1 && $4==""' "$out" | wc -l)
    if [[ "$bad" -gt 0 ]]; then
        echo "[P0] frames: ${bad} accessions unresolved (empty strand) -- see ${out}" >&2
        echo "[P0] step 'frames' NOT marked done" >&2
        return 1
    fi
    _say "frames: $(( $(wc -l < "$out") - 1 )) accessions -> ${out}"
    _mark frames
}

# --- step: assets ------------------------------------------------------------
step_assets() {
    _done assets && { _say "assets: already done"; return 0; }
    # the panel SNP matrix that reference selection reads, and the masks and
    # accessory sequence the scorers read -- linked, not copied, with their
    # checksums recorded in the manifest so a changed input is detectable
    ln -sfn "$(cd "$(dirname "$PANEL_SNPS")" && pwd -P)/$(basename "$PANEL_SNPS")" \
        "${BUILD}/assets/panel_snps.vcf.gz"
    [[ -s "${PANEL_SNPS}.tbi" ]] && ln -sfn \
        "$(cd "$(dirname "$PANEL_SNPS")" && pwd -P)/$(basename "$PANEL_SNPS").tbi" \
        "${BUILD}/assets/panel_snps.vcf.gz.tbi"
    [[ -s "$MASK" ]] && ln -sfn "$(cd "$(dirname "$MASK")" && pwd -P)/$(basename "$MASK")" \
        "${BUILD}/assets/repeat_mask.bed"
    # The accessory LOCUS TABLE, not just the sequence. Its locus ids are keyed
    # on the H37Rv anchor of each graph bubble (ACC_0001594 carries pos=1594), so
    # it is the naming layer for off-path records and is shared across samples
    # regardless of which reference they were called against. P0 previously
    # linked only the FASTAs, which carry sequence but no anchors.
    [[ -s "${ACCESSORY_DIR}/panel_manifest.tsv" ]] && ln -sfn \
        "$(cd "$ACCESSORY_DIR" && pwd -P)/panel_manifest.tsv" \
        "${BUILD}/assets/accessory_loci.tsv"
    # IS6110 catalogue. These currently live under refbias/t11/, a test
    # directory: the extraction audit flagged the pipeline depending on test
    # artefacts, and the IS6110 locus list and anchor sets are per-graph inputs
    # like any other, so they belong in the build's assets.
    [[ -s "${IS6110_LOCI:-refbias/t11/loci_clean.tsv}" ]] && ln -sfn \
        "$(cd "$(dirname "${IS6110_LOCI:-refbias/t11/loci_clean.tsv}")" && pwd -P)/$(basename "${IS6110_LOCI:-refbias/t11/loci_clean.tsv}")" \
        "${BUILD}/assets/is6110_loci.tsv"
    [[ -s "${IS6110_ANCHORS:-refbias/t11/anchor_sets.tsv}" ]] && ln -sfn \
        "$(cd "$(dirname "${IS6110_ANCHORS:-refbias/t11/anchor_sets.tsv}")" && pwd -P)/$(basename "${IS6110_ANCHORS:-refbias/t11/anchor_sets.tsv}")" \
        "${BUILD}/assets/is6110_anchors.tsv"
    [[ -s "${IS6110_GFF:-data/annotation/H37Rv_IS6110.pansn.gff}" ]] && ln -sfn \
        "$(cd "$(dirname "${IS6110_GFF:-data/annotation/H37Rv_IS6110.pansn.gff}")" && pwd -P)/$(basename "${IS6110_GFF:-data/annotation/H37Rv_IS6110.pansn.gff}")" \
        "${BUILD}/assets/is6110_elements.gff"
    local f
    for f in "${ACCESSORY_DIR}"/accessory_*.fasta; do
        [[ -s "$f" ]] && ln -sfn "$(cd "$(dirname "$f")" && pwd -P)/$(basename "$f")" \
            "${BUILD}/assets/$(basename "$f")"
    done
    _say "assets: $(ls -1 "${BUILD}/assets" | wc -l) entries"
    _mark assets
}

# --- step: manifest ----------------------------------------------------------
step_manifest() {
    local out="${BUILD}/manifest.tsv"
    { printf 'asset\tpath\tbytes\tsha256\n'
      local p
      while IFS= read -r p; do
          [[ -f "$p" ]] || continue
          printf '%s\t%s\t%s\t%s\n' "$(basename "$p")" "$p" \
              "$(stat -Lc%s "$p")" "$(sha256sum "$p" | cut -d' ' -f1)"
      # The GFFs are external NCBI downloads and can change upstream without
      # notice, so they are the single most important thing to pin. The staged
      # reference FASTAs are omitted deliberately: they are decompressed
      # byte-for-byte from data/assemblies, which carries its own provenance, and
      # hashing 333 x 4.4 MB on every manifest run buys nothing.
      done < <(find "$BUILD" -maxdepth 2 \
                 \( -name '*.tsv' -o -name '*.txt' -o -name '*.bed' \
                    -o -name '*.vcf.gz' -o -name '*.fasta' -o -name '*.gff.gz' \) \
                 -not -path "${BUILD}/refs/*" -not -name manifest.tsv | sort)
    } > "$out"
    _say "manifest: $(( $(wc -l < "$out") - 1 )) assets recorded"
    _say "annotation GFFs: $(ls -1 "${BUILD}/annotation"/*.gff.gz 2>/dev/null | wc -l)"
    _say "indexed refs   : $(ls -1 "${BUILD}/refs"/*.bwt 2>/dev/null | wc -l)"
    _say "sequence dicts : $(ls -1 "${BUILD}/refs"/*.dict 2>/dev/null | wc -l)"
    # A manifest is a provenance record, so it must not claim completeness while
    # assets are still being produced. Write it either way -- it is useful as a
    # snapshot -- but only mark the step done when every prior step is done, so a
    # partial build cannot be mistaken for a finished one later.
    local incomplete=""
    local s2
    for s2 in stamp paths accessions gff refs frames assets; do
        _done "$s2" || incomplete="${incomplete} ${s2}"
    done
    if [[ -n "$incomplete" ]]; then
        echo "[P0] manifest written as a SNAPSHOT only; steps not done:${incomplete}" >&2
        echo "[P0] step 'manifest' NOT marked done -- rerun when those complete" >&2
        return 1
    fi
    _mark manifest
}

if [[ -n "$STEP" ]]; then
    "step_${STEP}"
else
    step_stamp; step_paths; step_accessions; step_assets
    _say "cheap steps complete. Remaining, run explicitly:"
    _say "  bash bin/p0_prepare.sh --step gff"
    _say "  sbatch --array=1-\$(wc -l < ${BUILD}/assets/accessions.txt) bin/p0_prepare.sh --step refs"
    _say "  bash bin/p0_prepare.sh --step frames"
    _say "  bash bin/p0_prepare.sh --step manifest"
fi
