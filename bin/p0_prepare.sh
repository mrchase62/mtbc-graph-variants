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
#   OG=<graph.og> ACCESSORY_DIR=<panel dir> bash bin/p0_prepare.sh   # all cheap steps
#   OG=<graph.og> bash bin/p0_prepare.sh --step gff          # one step
#   OG=<graph.og> bash bin/p0_prepare.sh --list              # what exists and what does not
#   sbatch --array=1-333 bin/p0_prepare.sh --step refs   # the expensive step (OG exported)
#
# Every input is the graph's or named: OG (or GRAPH_DIR) has no default, nor
# has ACCESSORY_DIR (docs/PANEL_TREE.md section 4); the outgroup is
# config/project_env.sh's MTB_OUTGROUP, stamped into build_info.tsv.
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
# Steps are idempotent: each writes a .done marker recording the build id and
# the checksums of its inputs and code, and is skipped when they still match.
# A marker that no longer matches is refused, not silently rerun or reused;
# deleting a marker re-runs that step alone. Assets are copied into the build,
# never linked, so a finished build cannot change underneath its manifest.
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

# THE GRAPH AND ITS DIRECTORY ARE ONE CHOICE (audit P0P2-3). OG alone used to
# leave GRAPH_DIR at CX333, so a new graph got a new build id and the old
# graph's SNP matrix. Given OG, the directory is OG's; given both, they must
# agree; given GRAPH_DIR alone, the one graph in it. Given neither, P0
# refuses: it used to fall back to graphs/CX333..., so a run that forgot OG
# re-prepared the old graph under its old build id without a word.
if [[ -z "${OG:-}" && -z "${GRAPH_DIR:-}" ]]; then
    echo "FATAL: set OG (the graph's .og) or GRAPH_DIR (the directory holding it)" >&2
    exit 1
fi
if [[ -n "${OG:-}" ]]; then
    _og_dir="$(cd "$(dirname "$OG")" 2>/dev/null && pwd -P)" \
        || { echo "FATAL: OG=${OG} is not readable" >&2; exit 1; }
    if [[ -n "${GRAPH_DIR:-}" && "$(cd "$GRAPH_DIR" 2>/dev/null && pwd -P)" != "$_og_dir" ]]; then
        echo "FATAL: OG=${OG} is not in GRAPH_DIR=${GRAPH_DIR}; set one, or both consistently" >&2
        exit 1
    fi
    GRAPH_DIR="${GRAPH_DIR:-$(dirname "$OG")}"
else
    mapfile -t _ogs < <(ls "${GRAPH_DIR}"/*.smooth.final.og 2>/dev/null)
    [[ "${#_ogs[@]}" -le 1 ]] \
        || { echo "FATAL: ${#_ogs[@]} graphs in ${GRAPH_DIR}; set OG" >&2; exit 1; }
    OG="${_ogs[0]:-}"
fi
PANEL_SNPS="${PANEL_SNPS:-${GRAPH_DIR}/snps.vcf.gz}"
# The graph's collapsed variant file (bin/vcf_decompose.sh's product), the
# source of panel allele frequencies (P5) and the outgroup's alleles
# (panel_polarity). Copied into the build like every other asset.
GRAPH_VCF="${GRAPH_VCF:-${GRAPH_DIR}/all_variants.collapsed.vcf.gz}"
# The outgroup: config/project_env.sh's MTB_OUTGROUP unless OUTGROUP is set
# (empty: the graph has none). Recorded in build_info.tsv by step stamp, and
# step panel_polarity refuses a different one, so the polarity table and the
# association tail's tree root are always the same genome.
OUTGROUP="${OUTGROUP-${MTB_OUTGROUP-}}"
# The IS6110-clean builds whose crossmaps give P1's tie-break its interval
# counts. (The accessory catalogue is made by step 'catalogue'.)
ISCLEAN_DIR="${ISCLEAN_DIR:-is6110/assets/isclean_matched}"
H37RV_PATH="${H37RV_PATH:-GCF_000195955#1#NC_000962.3}"
H37RV_ACC="${H37RV_ACC:-GCF_000195955}"
ASM_DIR="${ASM_DIR:-data/assemblies}"
NCBI_SUM="${NCBI_SUM:-data/ncbi/assembly_summary_refseq.txt}"
MASK="${MASK:-data/annotation/H37Rv_repeat_mask.bed}"
# The FASTA the graph was built from. Its sequences are dnaA-rotated, so for a
# minority of accessions the graph's coordinate frame differs from refs/; the
# frames step measures that difference once per build.
# pggb names its outputs after the input FASTA (<input>.<hashes>.smooth.final
# .og), so the default is derived from the graph, not fixed to CX333; the
# frames step also checks its sequence names against the graph's paths.
PANEL_FASTA="${PANEL_FASTA:-data/fastas/$(basename "${OG:-x}" | sed -n 's/^\(.*\.fa\(sta\)\{0,1\}\.gz\)\..*/\1/p')}"
BWA="${MTB_BWA:-${MTB_QC_BIN}/bwa}"
GATK_SIF="${MTB_GATK_SIF:?MTB_GATK_SIF is unset; see config/project_env.sh}"
gatk_run() { singularity exec -B "${MTB_WORK}:${MTB_WORK}" "$GATK_SIF" gatk "$@"; }
# The accessory panel directory: panel_manifest.tsv (the locus table) and
# accessory_{novel,mosaic}.fasta, made from THIS graph by the commands in
# docs/PANEL_TREE.md section 4. No default: it was refbias/panel, CX333's, which a
# new graph's build would have copied in. Step assets requires it and checks
# its genome record against the build's.
ACCESSORY_DIR="${ACCESSORY_DIR:-}"
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
GRAPH_SHA="$(sha256sum "$OG" | cut -d' ' -f1)"
BUILD_ID="${BUILD_ID:-}"
if [[ -z "$BUILD_ID" ]]; then
    BUILD_ID="${GRAPH_SHA:0:12}"
fi
BUILD="${BUILD_ROOT}/${BUILD_ID}"
# An existing build directory belongs to the graph recorded in it. BUILD_ID
# set by hand must not let a second graph write into it.
if [[ -s "${BUILD}/build_info.tsv" ]]; then
    _rec="$(awk -F'\t' '$1=="graph_sha256"{print $2}' "${BUILD}/build_info.tsv")"
    [[ "$_rec" == "$GRAPH_SHA" ]] || {
        echo "FATAL: ${BUILD} was built from graph ${_rec}, not ${OG} (${GRAPH_SHA})" >&2; exit 1; }
fi
mkdir -p "$BUILD"/{annotation,refs,assets,logs} slurm

# A STEP IS DONE ONLY FOR THE INPUTS IT RAN ON (audit TP-4). The marker used
# to hold a date, so after a fix to ancestral_alleles.py, or a new tree,
# `--step ancestral` printed "already done" and kept the faulty table. Each
# marker now records the build id and, for steps that read anything outside
# the graph, the sha256 of those inputs and of the code. A marker whose record
# differs from the current inputs is REFUSED, not rerun over and not reused:
# delete the marker (and the step's products) to redo the step deliberately.
_STEP_KEYS=()           # "key<TAB>value" lines for the step being checked
_inputs() {             # _inputs key file [key file ...]: record sha256s
    _STEP_KEYS=("build_id	${BUILD_ID}")
    while [[ $# -ge 2 ]]; do
        if [[ -e "$2" ]]; then
            _STEP_KEYS+=("$1	$(sha256sum "$2" | cut -d' ' -f1)")
        else
            _STEP_KEYS+=("$1	absent")
        fi
        shift 2
    done
}
_done() {
    local m="${BUILD}/logs/$1.done" kv k v have
    [[ -s "$m" ]] || return 1
    [[ "${#_STEP_KEYS[@]}" -gt 0 ]] || _STEP_KEYS=("build_id	${BUILD_ID}")
    for kv in "${_STEP_KEYS[@]}"; do
        k="${kv%%	*}"; v="${kv#*	}"
        have="$(mtb_kv "$m" "$k")"
        # a date-only marker inside the directory named by this graph's
        # checksum is this build's; only recorded INPUTS can be stale
        [[ "$k" == build_id && -z "$have" ]] && continue
        if [[ "$have" != "$v" ]]; then
            echo "FATAL: step '$1' was done with ${k}=${have:-<not recorded>}, now ${v}." >&2
            echo "       Its products are from other inputs or code. Delete ${m} and" >&2
            echo "       the step's products to redo it, or use a new build." >&2
            exit 1
        fi
    done
}
_mark() {
    local kv
    # A MANIFEST PREDATES ANY STEP DONE AFTER IT (review 2, R2-INT-3): a step
    # re-run once the build was marked complete changes what the manifest
    # vouches for, so the build is incomplete again until it is re-made.
    [[ "$1" == manifest ]] || rm -f "${BUILD}/logs/manifest.done"
    [[ "${#_STEP_KEYS[@]}" -gt 0 ]] || _STEP_KEYS=("build_id	${BUILD_ID}")
    { printf 'completed\t%s\n' "$(date -Is)"
      for kv in "${_STEP_KEYS[@]}"; do printf '%s\n' "$kv"; done
    } > "${BUILD}/logs/$1.done.tmp"
    mv -f "${BUILD}/logs/$1.done.tmp" "${BUILD}/logs/$1.done"
    _STEP_KEYS=()
}
_say()   { echo "[P0:${BUILD_ID}] $*"; }
# Copy, never link: a linked asset changes when its source does, and the
# build's record then describes a file nobody can recover (audit P0P2-2).
_copy()  { cp -fL "$1" "${2}.tmp.$$" && mv -f "${2}.tmp.$$" "$2"; }

if [[ "$LIST" == 1 ]]; then
    echo "graph    : $OG"
    echo "build id : $BUILD_ID"
    echo "build dir: $BUILD"
    for s in stamp paths accessions gff refs frames assets catalogue nodes is6110_intervals \
             panel_polarity ancestral manifest; do
        if [[ -s "${BUILD}/logs/$s.done" ]]; then
            printf "  %-16s done   %s\n" "$s" "$(head -1 "${BUILD}/logs/$s.done")"
        else printf "  %-16s MISSING\n" "$s"; fi
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
        # WHICH SEQUENCES refs/ HOLDS (audit PGB-15, kept by decision). They
        # are the deposited assemblies, not the graph's dnaA-rotated paths:
        # the RefSeq GFFs in annotation/ and the IS6110 crossmaps are in the
        # deposited frame, and assets/graph_frame_offsets.tsv converts.
        echo "refs_source\t${ASM_DIR} (deposited frame; graph frame via assets/graph_frame_offsets.tsv)"
        # the panel genome that roots the trees and polarises the variants
        # (config MTB_OUTGROUP); 'none' when the graph has no outgroup
        echo "outgroup\t${OUTGROUP:-none}"
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
        # to a temporary name: a task killed mid-write left a truncated FASTA
        # that the next run's -s test accepted and indexed (audit P0P2-9)
        [[ -s "$fa" ]] || { zcat "$src" > "${fa}.tmp.$$" && mv -f "${fa}.tmp.$$" "$fa"; }
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
        local dicts; dicts=$(ls "${BUILD}/refs"/*.dict 2>/dev/null | wc -l)
        _say "refs: ${built} of ${want} indexed, ${dicts} with a dictionary"
        # a logged dict_failed is a missing asset too (audit P0P2-9)
        if [[ "$built" -eq "$want" && "$dicts" -eq "$want" ]]; then _mark refs
        else echo "[P0] step 'refs' NOT marked done -- ${want} expected, ${built} built" >&2; fi
    fi
}

# A dependency is met when the step's marker exists; whether that step's own
# inputs are current is its own check, run when it is run.
_have()  { [[ -s "${BUILD}/logs/$1.done" ]]; }

# --- step: frames ------------------------------------------------------------
# The panel-vs-refs coordinate frame of every accession. graph_frame.py reads it
# to convert a refs-frame coordinate before it is handed to odgi, which
# interprets coordinates in the PANEL frame. It is a property of this build's
# graph and refs, so it is measured here and nowhere else; refbias_run.sh points
# MTB_GRAPH_FRAMES at it. Needs the refs step to have finished.
step_frames() {
    mtb_require_file "$PANEL_FASTA"
    _inputs panel_fasta_sha256 "$PANEL_FASTA" \
            code_sha256 graphframe/bin/graph_frame_offsets.py
    _done frames && { _say "frames: already done"; return 0; }
    _have refs || { echo "FATAL: step 'frames' needs step 'refs' done first" >&2; return 1; }
    # The FASTA must be the one this graph was built from: its sequence names
    # are the graph's path names, all of them and nothing else.
    "$MTB_PY" bin/p0_check.py fasta-names --fasta "$PANEL_FASTA" \
        --paths "${BUILD}/assets/paths.txt" || return 1
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
# COPIED INTO THE BUILD, NOT LINKED (audit P0P2-2). The links resolved to
# files that kept changing: the repeat mask went from 181 to 395 intervals on
# Sep 28 under a manifest that recorded the old checksum, and every VCF stamp
# since then names a mask that was not used. A build is now immutable: its
# assets are its own bytes, the manifest hashes them, and refbias_run.sh
# re-hashes them before every chain (bin/p0_check.py verify).
#
# And each graph-derived asset is CHECKED AGAINST THE GRAPH (audit P0P2-3):
# the panel SNP matrix and the collapsed graph VCF must have exactly the
# build's genomes as samples, so another graph's files are refused.
step_assets() {
    local acc="${BUILD}/assets/accessions.txt"
    _inputs panel_snps_sha256 "$PANEL_SNPS" graph_vcf_sha256 "$GRAPH_VCF" \
            mask_sha256 "$MASK" \
            accessory_loci_sha256 "${ACCESSORY_DIR}/panel_manifest.tsv" \
            accessory_genomes_sha256 "${ACCESSORY_DIR}/panel_manifest.genomes.txt" \
            is6110_loci_sha256 "${IS6110_LOCI:-}" \
            is6110_anchors_sha256 "${IS6110_ANCHORS:-}" \
            is6110_gff_sha256 "${IS6110_GFF:-data/annotation/H37Rv_IS6110.pansn.gff}"
    _done assets && { _say "assets: already done"; return 0; }
    _have accessions || { echo "FATAL: step 'assets' needs step 'accessions' done first" >&2; return 1; }
    [[ -n "$ACCESSORY_DIR" ]] || {
        echo "FATAL: step 'assets' needs ACCESSORY_DIR: this graph's accessory panel" >&2
        echo "       (panel_manifest.tsv and its genome record; docs/PANEL_TREE.md section 4)" >&2
        return 1; }
    mtb_require_file "$PANEL_SNPS" "$GRAPH_VCF" "$MASK" \
        "${ACCESSORY_DIR}/panel_manifest.tsv" || return 1
    "$MTB_PY" bin/p0_check.py vcf-samples --vcf "$PANEL_SNPS" \
        --accessions "$acc" --ref-acc "$H37RV_ACC" || return 1
    "$MTB_PY" bin/p0_check.py vcf-samples --vcf "$GRAPH_VCF" \
        --accessions "$acc" --ref-acc "$H37RV_ACC" || return 1
    # The accessory panel is made outside P0 (an array of blasts), so it must
    # say which genomes it was made from: build_accessory_panel.py writes
    # panel_manifest.genomes.txt, the graph VCF's samples, and they must be
    # this build's. A table with no record (CX333's refbias/panel) is refused.
    "$MTB_PY" bin/p0_check.py accessory-panel --dir "$ACCESSORY_DIR" \
        --accessions "$acc" --ref-acc "$H37RV_ACC" || return 1
    # replace links left by an earlier version of this step, never write
    # through them into the file they point at
    find "${BUILD}/assets" -maxdepth 1 -type l -delete
    # the panel SNP matrix that reference selection reads, the graph's variant
    # file, and the masks and accessory sequence the scorers read
    _copy "$PANEL_SNPS" "${BUILD}/assets/panel_snps.vcf.gz"
    [[ -s "${PANEL_SNPS}.tbi" ]] && _copy "${PANEL_SNPS}.tbi" "${BUILD}/assets/panel_snps.vcf.gz.tbi"
    _copy "$GRAPH_VCF" "${BUILD}/assets/graph_collapsed.vcf.gz"
    [[ -s "${GRAPH_VCF}.tbi" ]] && _copy "${GRAPH_VCF}.tbi" "${BUILD}/assets/graph_collapsed.vcf.gz.tbi"
    _copy "$MASK" "${BUILD}/assets/repeat_mask.bed"
    # The accessory LOCUS TABLE, not just the sequence. Its locus ids are keyed
    # on the H37Rv anchor of each graph bubble (ACC_0001594 carries pos=1594), so
    # it is the naming layer for off-path records and is shared across samples
    # regardless of which reference they were called against. P0 previously
    # linked only the FASTAs, which carry sequence but no anchors.
    _copy "${ACCESSORY_DIR}/panel_manifest.tsv" "${BUILD}/assets/accessory_loci.tsv"
    _copy "${ACCESSORY_DIR}/panel_manifest.genomes.txt" "${BUILD}/assets/accessory_loci.genomes.txt"
    # The accessory catalogue is no longer copied from accessory/assets/: it
    # is made from this build's own files by step 'catalogue' below.
    # IS6110 locus list and anchor sets: copied ONLY WHEN NAMED. They were
    # taken by default from refbias/t11/, a pilot test directory, though no
    # pass of the chain reads them (only the working tree's p1b/p1c
    # benchmarks, which are not production passes) and nothing in this
    # repository makes them. A build of a new graph must not carry CX333's.
    if [[ -n "${IS6110_LOCI:-}" ]]; then
        mtb_require_file "$IS6110_LOCI" || return 1
        _copy "$IS6110_LOCI" "${BUILD}/assets/is6110_loci.tsv"
    fi
    if [[ -n "${IS6110_ANCHORS:-}" ]]; then
        mtb_require_file "$IS6110_ANCHORS" || return 1
        _copy "$IS6110_ANCHORS" "${BUILD}/assets/is6110_anchors.tsv"
    fi
    [[ -s "${IS6110_GFF:-data/annotation/H37Rv_IS6110.pansn.gff}" ]] && \
        _copy "${IS6110_GFF:-data/annotation/H37Rv_IS6110.pansn.gff}" "${BUILD}/assets/is6110_elements.gff"
    local f
    for f in "${ACCESSORY_DIR}"/accessory_*.fasta; do
        [[ -s "$f" ]] && _copy "$f" "${BUILD}/assets/$(basename "$f")"
    done
    _say "assets: $(ls -1 "${BUILD}/assets" | wc -l) entries, none of them links"
    _mark assets
}

# --- step: catalogue ---------------------------------------------------------
# The accessory catalogue locus_presence.py genotypes against (level 1) and
# the merge reads for level-1 records. It had no producer in the repository:
# accessory/assets/accessory_catalogue.{tsv,fasta} were copied into the build
# from a file made by hand. Its producer is accessory/bin/merge_catalogues.py
# -- run on CX333's decomposed VCF with insgt's clusters and routing and
# gwas1000's census, it reproduces both files byte for byte -- so it is made
# here from the build's own accessory locus table and collapsed graph VCF
# (the collapsed file, audit PGB-9: on CX333, 157 of 802 rows differ from
# the decomposed file's catalogue, 23 of them in the representative allele's
# length, the rest in carrier sets).
#
# The census (one cohort's call counts) is not a property of the build and
# is not passed; its three columns are empty. The insgt clusters and routing
# (graph-specific, hand-made) are passed only when INSGT_CLUSTERS and
# INSGT_ROUTING name this graph's; otherwise the route columns are empty and
# the merge writes no ACCROUTE.
step_catalogue() {
    # h37rv_cov95 (review 2, R2-IS-1) is measured against the build's H37Rv
    local h37fa="${BUILD}/refs/${H37RV_ACC}.fasta"
    local blastn="${MTB_BLASTN:-${MTB_QC_BIN}/blastn}"
    _inputs loci_sha256 "${BUILD}/assets/accessory_loci.tsv" \
            graph_vcf_sha256 "${BUILD}/assets/graph_collapsed.vcf.gz" \
            clusters_sha256 "${INSGT_CLUSTERS:-}" routing_sha256 "${INSGT_ROUTING:-}" \
            h37rv_sha256 "$h37fa" \
            code_sha256 accessory/bin/merge_catalogues.py
    _done catalogue && { _say "catalogue: already done"; return 0; }
    _have assets || { echo "FATAL: step 'catalogue' needs step 'assets' done first" >&2; return 1; }
    _have refs || { echo "FATAL: step 'catalogue' needs step 'refs' done first (it measures h37rv_cov95 against ${h37fa})" >&2; return 1; }
    [[ -x "$blastn" ]] || { echo "FATAL: no blastn at ${blastn} (set MTB_BLASTN)" >&2; return 1; }
    mtb_require_file "${BUILD}/assets/accessory_loci.tsv" \
        "${BUILD}/assets/graph_collapsed.vcf.gz" || return 1
    local out="${BUILD}/assets/accessory_catalogue"
    rm -f "${out}".fasta.*          # an index of an earlier catalogue
    "$MTB_PY" accessory/bin/merge_catalogues.py \
        --loci "${BUILD}/assets/accessory_loci.tsv" \
        --graph-vcf "${BUILD}/assets/graph_collapsed.vcf.gz" \
        --clusters "${INSGT_CLUSTERS:-}" --routing "${INSGT_ROUTING:-}" \
        --census "" --out "${out}.tsv.tmp" --out-fasta "${out}.fasta.tmp" \
        --h37rv "$h37fa" --blastn "$blastn" \
        || { rm -f "${out}".*.tmp; return 1; }
    "$MTB_PY" bin/p0_check.py catalogue --tsv "${out}.tsv.tmp" \
        --accessions "${BUILD}/assets/accessions.txt" || { rm -f "${out}".*.tmp; return 1; }
    mv -f "${out}.tsv.tmp" "${out}.tsv"
    mv -f "${out}.fasta.tmp" "${out}.fasta"
    # indexed once here: locus_presence.py indexes a FASTA without one, and
    # an array of tasks would each write the index into the build at once
    "$BWA" index "${out}.fasta" 2>/dev/null
    _say "catalogue: $(( $(wc -l < "${out}.tsv") - 1 )) loci -> ${out}.{tsv,fasta}"
    _mark catalogue
}

# --- step: nodes -------------------------------------------------------------
# Node-to-path membership and node offsets, which p5_states.py needs for every
# node-frame key. They were made by hand from one GFA into accessory/assets/,
# outside any build, and P5 read them by a relative default; node ids are
# graph-specific, so on a new graph P5 would have read the old graph's
# membership with no error (audit P4P5-7). Made here from the build's own
# .og (odgi view -g), over every node the H37Rv path does not traverse --
# which is where every node-frame key lies -- and passed explicitly by
# p5_merge.sh. About five minutes on CX333.
step_nodes() {
    _inputs code_sha256 accessory/bin/node_path_membership.py
    _done nodes && { _say "nodes: already done"; return 0; }
    local gfa="${BUILD}/logs/graph.$$.gfa"
    "$MTB_ODGI" view -i "$OG" -g > "${gfa}.tmp"
    mv -f "${gfa}.tmp" "$gfa"
    if ! "$MTB_PY" accessory/bin/node_path_membership.py --gfa "$gfa" \
            --exclude-path "$H37RV_PATH" \
            --out "${BUILD}/assets/node_paths.tsv" \
            --out-positions "${BUILD}/assets/node_positions.tsv"; then
        rm -f "$gfa"; return 1
    fi
    rm -f "$gfa"
    _mark nodes
}

# --- step: is6110_intervals --------------------------------------------------
# The P1 tie-break's per-reference IS6110 interval counts. p1_summary.py read
# them from is6110/assets/isclean_matched, outside the build, and a missing
# crossmap silently skipped the rule (audit P0P2-4). They are a property of
# the panel's assemblies, so they are a build asset: one row per panel genome
# but H37Rv (never a candidate), each checked against this build's refs/.
step_is6110_intervals() {
    _inputs isclean_manifest_sha256 "${ISCLEAN_DIR}/manifest.tsv" \
            code_sha256 bin/p0_check.py
    _done is6110_intervals && { _say "is6110_intervals: already done"; return 0; }
    _have refs || { echo "FATAL: step 'is6110_intervals' needs step 'refs' done first" >&2; return 1; }
    local out="${BUILD}/assets/is6110_intervals.tsv"
    "$MTB_PY" bin/p0_check.py is6110-intervals --isclean-dir "$ISCLEAN_DIR" \
        --refs "${BUILD}/refs" --accessions "${BUILD}/assets/accessions.txt" \
        --skip "$H37RV_ACC" --out "${out}.tmp" || { rm -f "${out}.tmp"; return 1; }
    mv -f "${out}.tmp" "$out"
    _mark is6110_intervals
}

# --- step: panel_polarity ----------------------------------------------------
# The outgroup's allele at every graph variant, which the association chain's
# event writer reads for polarity. It was refbias/assets/panel_polarity.tsv,
# made by hand from whichever VCF; it is a property of this graph, so it is
# made here from the build's own collapsed VCF. OUTGROUP= (empty) records that
# the graph has no outgroup, and the step writes nothing.
step_panel_polarity() {
    # keys are left-aligned against the build's own H37Rv, exactly as the
    # cohort VCF's are (review 2, R2-INT-1), so that sequence is an input too
    local h37fa="${BUILD}/refs/${H37RV_ACC}.fasta"
    _inputs graph_vcf_sha256 "${BUILD}/assets/graph_collapsed.vcf.gz" \
            h37rv_sha256 "$h37fa" \
            code_sha256 bin/panel_polarity.py \
            norm_code_sha256 bin/mtb_norm.py
    _STEP_KEYS+=("outgroup	${OUTGROUP:-none}")
    _done panel_polarity && { _say "panel_polarity: already done"; return 0; }
    _have assets || { echo "FATAL: step 'panel_polarity' needs step 'assets' done first" >&2; return 1; }
    # the outgroup the build was stamped with is the one the association tail
    # roots its trees on; a table made with another would disagree with it
    local rec; rec="$(mtb_kv "${BUILD}/build_info.tsv" outgroup)"
    if [[ -n "$rec" && "$rec" != "${OUTGROUP:-none}" ]]; then
        echo "FATAL: build ${BUILD_ID} was stamped with outgroup ${rec}, not '${OUTGROUP:-none}'" >&2
        return 1
    fi
    if [[ -z "$OUTGROUP" ]]; then
        _say "panel_polarity: OUTGROUP is empty; no table (the event writer needs --panel-polarity '')"
        _mark panel_polarity; return 0
    fi
    _have refs || { echo "FATAL: step 'panel_polarity' needs step 'refs' done first (it left-aligns against ${h37fa})" >&2; return 1; }
    local out="${BUILD}/assets/panel_polarity.tsv"
    "$MTB_PY" bin/panel_polarity.py --panel-vcf "${BUILD}/assets/graph_collapsed.vcf.gz" \
        --h37rv "$h37fa" --h37rv-contig "${H37RV_PATH##*#}" \
        --outgroup "$OUTGROUP" --bcftools "$MTB_BCFTOOLS" --out "${out}.tmp" \
        || { rm -f "${out}.tmp"; return 1; }
    mv -f "${out}.tmp" "$out"
    _mark panel_polarity
}

# --- step: ancestral ---------------------------------------------------------
# The ancestral allele at every panel SNP site, for the merged VCF's AA tag. A
# property of the panel and its tree, not of any cohort, so it is a build asset
# and every cohort's VCF joins the same table. It used to be a file in
# data/trees that nothing regenerated; p5_finish.sh now reads it from here and
# refuses a build without it. Seconds to compute.
#
# Its inputs are a tree and alignment made outside P0 by the panel-tree step
# (docs/PANEL_TREE.md: vcf_to_alignment.py on this build's
# assets/graph_collapsed.vcf.gz, then build_snp_tree.sh), so they are checked
# against the panel -- the alignment's taxa must be exactly the build's
# genomes -- copied into the build, and recorded in the marker with the
# code's checksum (audit TP-4: after the polarity fixes the old marker kept
# the faulty table).
#
# THERE ARE NO DEFAULT INPUTS. They were data/trees/cx333.*, CX333's tree,
# which the taxa check refuses on any other panel but which a rebuild of a
# panel with the same genomes would have taken silently. ANC_TREE, ANC_ALN
# and ANC_SITES must name this build's panel-tree products.
step_ancestral() {
    local t="${ANC_TREE:-}" al="${ANC_ALN:-}" st="${ANC_SITES:-}"
    if [[ -z "$t" || -z "$al" || -z "$st" ]]; then
        echo "FATAL: step 'ancestral' needs ANC_TREE, ANC_ALN and ANC_SITES: this build's" >&2
        echo "       panel tree, alignment and sites table (docs/PANEL_TREE.md)" >&2
        return 1
    fi
    mtb_require_file "$t"; mtb_require_file "$al"; mtb_require_file "$st"
    _inputs tree_sha256 "$t" alignment_sha256 "$al" sites_sha256 "$st" \
            code_sha256 bin/ancestral_alleles.py
    # the tree's outgroup leaves; unset, ancestral_alleles.py's default (the
    # two CX333 canettii), which it refuses if they are not leaves of the tree
    local -a og_args=()
    if [[ -n "${ANC_OUTGROUPS:-}" ]]; then
        og_args=(--outgroups "$ANC_OUTGROUPS")
        _STEP_KEYS+=("outgroups	${ANC_OUTGROUPS}")
    fi
    _done ancestral && { _say "ancestral: already done"; return 0; }
    _have accessions || { echo "FATAL: step 'ancestral' needs step 'accessions' done first" >&2; return 1; }
    "$MTB_PY" bin/p0_check.py taxa --fasta "$al" \
        --accessions "${BUILD}/assets/accessions.txt" || return 1
    mkdir -p "${BUILD}/assets/ancestral_inputs"
    _copy "$t"  "${BUILD}/assets/ancestral_inputs/tree.nwk"
    _copy "$al" "${BUILD}/assets/ancestral_inputs/alignment.fasta"
    _copy "$st" "${BUILD}/assets/ancestral_inputs/sites.tsv"
    local out="${BUILD}/assets/ancestral.tsv"
    "$MTB_PY" bin/ancestral_alleles.py --tree "${BUILD}/assets/ancestral_inputs/tree.nwk" \
        --alignment "${BUILD}/assets/ancestral_inputs/alignment.fasta" \
        --sites "${BUILD}/assets/ancestral_inputs/sites.tsv" "${og_args[@]}" \
        --out "${out}.tmp"
    mv -f "${out}.tmp" "$out"
    _say "ancestral: $(( $(wc -l < "$out") - 1 )) sites -> ${out}"
    _mark ancestral
}

# --- step: manifest ----------------------------------------------------------
# Every file under assets/ and annotation/, plus build_info.tsv. The proj
# stores and refs_* trees are caches that grow by design, so they are not in
# it; bin/p0_check.py verify re-hashes exactly these rows at every use.
step_manifest() {
    local out="${BUILD}/manifest.tsv"
    # an earlier complete manifest must not vouch for this one, which may be
    # a snapshot (review 2, R2-INT-3)
    rm -f "${BUILD}/logs/manifest.done"
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
      done < <( { echo "${BUILD}/build_info.tsv"
                  find "${BUILD}/assets" "${BUILD}/annotation" -type f; } | sort)
    } > "${out}.tmp"
    mv -f "${out}.tmp" "$out"
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
    # ancestral too (review 2, R2-BUILD-2): without it the build was marked
    # complete and the chain failed only at its last pass, p5vcf
    for s2 in stamp paths accessions gff refs frames assets catalogue nodes \
              is6110_intervals panel_polarity ancestral; do
        _have "$s2" || incomplete="${incomplete} ${s2}"
    done
    if [[ -n "$incomplete" ]]; then
        echo "[P0] manifest written as a SNAPSHOT only; steps not done:${incomplete}" >&2
        echo "[P0] step 'manifest' NOT marked done -- rerun when those complete" >&2
        return 1
    fi
    # an asset left as a link would make the build mutable again
    if [[ -n "$(find "${BUILD}/assets" "${BUILD}/annotation" -type l | head -1)" ]]; then
        echo "[P0] assets/ still holds links; rerun --step assets (delete logs/assets.done first)" >&2
        return 1
    fi
    _inputs; _mark manifest
}

_run() { _STEP_KEYS=(); "step_$1"; }
if [[ -n "$STEP" ]]; then
    _run "$STEP"
else
    _run stamp; _run paths; _run accessions; _run assets; _run catalogue; _run nodes; _run panel_polarity
    _say "cheap steps complete. Remaining, run explicitly:"
    _say "  bash bin/p0_prepare.sh --step gff"
    _say "  sbatch --array=1-\$(wc -l < ${BUILD}/assets/accessions.txt) bin/p0_prepare.sh --step refs"
    _say "  bash bin/p0_prepare.sh --step frames"
    _say "  bash bin/p0_prepare.sh --step is6110_intervals   # crossmaps for every panel genome first:"
    _say "      MTB_BUILD_DIR=${BUILD} bash is6110/bin/p1i_discover_matched.sh"
    _say "      MTB_BUILD_DIR=${BUILD} bash is6110/bin/p1i_build_matched.sh"
    _say "  ANC_TREE=.. ANC_ALN=.. ANC_SITES=.. bash bin/p0_prepare.sh --step ancestral   # docs/PANEL_TREE.md"
    _say "  bash bin/p0_prepare.sh --step manifest"
fi
