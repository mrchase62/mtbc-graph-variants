#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Mtb pangenome project — single source of truth for all paths.
#
# Every script in bin/ sources this file. Nothing anywhere else should contain
# a literal /n/netscratch or /n/boslfs02 path.
#
#   source "$(dirname "$0")/../config/project_env.sh"
#
# Override any variable from the environment, e.g. to run somewhere else:
#   MTB_WORK=/some/other/place sbatch bin/pggb_build.sh ...
# ---------------------------------------------------------------------------

# --- Site file --------------------------------------------------------------
# A SITE FILE, GITIGNORED, holds the values for whichever cluster this is. That
# is what keeps the published configuration free of one site's paths while a
# working checkout still runs.
#
# It must be sourced FIRST, before any default below. It used to be sourced
# after the Roots and working subdirectories were assigned, so a site file that
# set MTB_WORK moved MTB_WORK alone while MTB_DATA, MTB_GRAPHS, MTB_LOGS and the
# rest still pointed into the old scratch tree.
#
# It is found relative to THIS file, not to the current directory. The old
# `config/site.local.sh` resolved against $PWD, so sourcing the config from
# anywhere but the repository root skipped the site file without a word and
# left MTB_CRAM_ROOT empty.
_mtb_cfg_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
_mtb_site="${MTB_SITE_FILE:-${_mtb_cfg_dir}/site.local.sh}"
if [[ -r "$_mtb_site" ]]; then
    source "$_mtb_site"
elif [[ -n "${MTB_SITE_FILE:-}" ]]; then
    echo "[project_env] WARNING: MTB_SITE_FILE=${MTB_SITE_FILE} is not readable" >&2
fi
unset _mtb_site

# --- Roots ------------------------------------------------------------------

# Working root: fast scratch. PURGED PERIODICALLY. Nothing here is safe.
: "${MTB_WORK:=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome}"

# Persistent root: durable lab storage. Results are synced here by bin/sync_back.sh.
: "${MTB_PERSIST:=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/MtbPangenome}"

# Archive: the Oct-Dec 2025 project, read-only source for staging inputs.
: "${MTB_ARCHIVE:=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/bkups260111/Assemblies}"

# --- Working subdirectories -------------------------------------------------

: "${MTB_BIN:=${MTB_WORK}/bin}"
: "${MTB_CONFIG:=${MTB_WORK}/config}"
: "${MTB_DATA:=${MTB_WORK}/data}"
: "${MTB_FASTAS:=${MTB_DATA}/fastas}"
: "${MTB_REFDIR:=${MTB_DATA}/ref}"
: "${MTB_ASSEMBLIES:=${MTB_DATA}/assemblies}"
: "${MTB_ROTATED:=${MTB_DATA}/rotated_assemblies}"
: "${MTB_GRAPHS:=${MTB_WORK}/graphs}"
: "${MTB_LOGS:=${MTB_WORK}/logs}"
: "${MTB_SLURM:=${MTB_WORK}/slurm}"

# --- Containers -------------------------------------------------------------

# --- Read source ------------------------------------------------------------
#
# THE PIPELINE CONSUMES ALIGNED READS; WHERE THEY COME FROM IS SITE CONFIG.
# Every cohort so far has been run off a shared CRAM collection staged by a
# colleague on netscratch. That was an opportunistic reuse that saved a great
# deal of time, and it is not the pipeline's input contract:
#
#   it is transient -- netscratch is purged without warning, and it took the
#     pggb and vg containers on 2026-09-28;
#   it is not ours -- the path carries another person's directory layout, which
#     should not be baked into scripts, least of all published ones;
#   it is one source among several -- FASTQ from new sequencing runs will not
#     look like it, and neither will the next collection.
#
# So this is declared once, here, with NO default. A script that needs it reads
# MTB_CRAM_ROOT and fails with a clear message when it is unset, rather than
# falling back to a path that happens to exist on one cluster today. Set it in
# your environment or in a site file; the 2026 runs used the collection under
# ${MTB_SCRATCH_ROOT:-/n/netscratch/sfortune_lab/Lab}/<collection>.
: "${MTB_CRAM_ROOT:=}"
# The reference the CRAMs were encoded against, needed to decode them.
: "${MTB_CRAM_REF:=}"
# The collection's parent. The 2026 collection carried a per-isolate results
# table, pipeline_inputs/ and trees/ there. Another read source need not have
# any of that, so anything depending on it must degrade rather than assume.
: "${MTB_CRAM_META_ROOT:=}"

# Fail loudly and early rather than letting a stale default resolve to nothing.
mtb_require_cram_root() {
    if [[ -z "${MTB_CRAM_ROOT:-}" ]]; then
        echo "FATAL: MTB_CRAM_ROOT is unset. The pipeline needs a collection of" >&2
        echo "  aligned reads and does not assume one. Set it to the root of your" >&2
        echo "  CRAM (or BAM) collection, e.g." >&2
        echo "    export MTB_CRAM_ROOT=/path/to/collection" >&2
        echo "  and MTB_CRAM_REF to the reference those files were encoded" >&2
        echo "  against. Both are declared in config/project_env.sh." >&2
        return 1
    fi
    [[ -d "$MTB_CRAM_ROOT" ]] || {
        echo "FATAL: MTB_CRAM_ROOT=${MTB_CRAM_ROOT} does not exist." >&2
        echo "  If this was a scratch collection it may have been purged." >&2
        return 1; }
}

: "${MTB_CONTAINERS:=${MTB_WORK}/containers}"
: "${MTB_PGGB_SIF:=${MTB_CONTAINERS}/pggb_latest.sif}"
: "${MTB_VG_SIF:=${MTB_CONTAINERS}/vg_v1.69.0.sif}"
# GraphAligner 1.0.19, for aligning long sequences (complete assemblies, de novo
# contigs) to a pggb GFA. vg's own aligners are short-read oriented; GraphAligner
# handles general cyclic GFA, which is what pggb produces.
: "${MTB_GRAPHALIGNER_SIF:=${MTB_CONTAINERS}/graphaligner.sif}"
: "${MTB_LAB_CONTAINERS:=/n/boslfs02/LABS/sfortune_lab/Lab/containers}"

# --- External tools (lab installs, not containerised) -----------------------

# Declared here because scripts read them: each used to be defaulted inside the
# script that needed it, so show_config could not report them and a site file
# was the only way to learn they existed.
: "${MTB_GATK_SIF:=${MTB_LAB_CONTAINERS}/mtb_gatk.sif}"
: "${MTB_DELLY_ENV:=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/tb-profiler}"

# Python with pandas/numpy/scikit-learn, for the sv_*.py tools.
: "${MTB_PY:=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_pangenome_qc/bin/python}"
# Python that can import mtbvartools and dendropy, for the association arm's
# bytestream writer. MTB_PY cannot: mtbvartools is an editable install present
# only in this environment, so a script needing it must name this interpreter
# explicitly rather than assume the pipeline python.
: "${MTB_PY_VT:=/n/home11/mchase/.conda/envs/mtb_isolates_cluster/bin/python}"
# odgi binary (untangle/pav/sort/flip). The pggb container also carries one.
: "${MTB_ODGI:=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/odgi/bin/odgi}"

# samtools / bedtools are NOT on the login-node PATH. They live in the
# mtb_pangenome_qc conda env (and in the pggb container). odgi_pav.sh needs both.
: "${MTB_QC_BIN:=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_pangenome_qc/bin}"
: "${MTB_SAMTOOLS:=${MTB_QC_BIN}/samtools}"
# bgzip lives in the same QC env. It was used by name in an ad-hoc command and
# expanded to the empty string, which silently produced a 0-byte FASTA that then
# failed to index -- the same failure mode MTB_MINIMAP2 had before it was defined.
: "${MTB_BGZIP:=${MTB_QC_BIN}/bgzip}"
# tabix ships with bgzip and was never named here, so anything needing it fell
# back to bare `tabix`, which is not on PATH on this cluster.
: "${MTB_TABIX:=${MTB_QC_BIN}/tabix}"
: "${MTB_BEDTOOLS:=${MTB_QC_BIN}/bedtools}"
# bcftools is NOT in the QC env and NOT on the login-node PATH; the mtb_bcftools
# shell function below falls back to the pggb container. Python helpers cannot
# call a shell function, so they need a real path -- passing "$MTB_BCFTOOLS"
# undefined produced `PermissionError: [Errno 13] Permission denied: ''`, the
# same failure mode MTB_MINIMAP2, MTB_BGZIP and MTB_K8 each had before.
: "${MTB_BCFTOOLS:=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates/bin/bcftools}"
# minimap2 lives in the same QC conda env. It was referenced by name in the
# project notes but never actually defined here, so anything using it got an
# empty string and failed with a bare PermissionError.
: "${MTB_MINIMAP2:=${MTB_QC_BIN}/minimap2}"
: "${MTB_BWA:=${MTB_QC_BIN}/bwa}"
: "${MTB_WGSIM:=${MTB_QC_BIN}/wgsim}"
: "${MTB_TMPBASE:=${TMPDIR:-/tmp}}"
# paftools.js needs the k8 javascript shell; neither is in the QC env, both ship
# with minimap2 and live in ~/bin here. Same failure mode again: they were passed
# by bare name to is6110_pairwise_call.py, which requires them as arguments, so
# an undefined variable became an empty argument.
: "${MTB_K8:=${HOME}/bin/k8}"
: "${MTB_PAFTOOLS:=${HOME}/bin/paftools.js}"
# BLAST+ for the foreign screen's UniVec check (D29); not in the QC env
: "${MTB_BLASTN:=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin/blastn}"
: "${MTB_MAKEBLASTDB:=/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin/makeblastdb}"
# NCBI UniVec_Core, fetched into the repository's data/ by bin/fetch_univec.sh
: "${MTB_UNIVEC:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/univec/UniVec_Core}"

: "${MTB_SNPEFF_JAR:=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/snpEff/snpEff.jar}"
: "${MTB_JAVA:=/n/boslfs02/LABS/sfortune_lab/Lab/software/jdk-20.0.1/bin/java}"
: "${MTB_SNPEFF_DB:=Mycobacterium_tuberculosis_h37rv}"
# Contig name the snpEff database expects (PanSN path name is stripped to this).
: "${MTB_SNPEFF_CONTIG:=NC_000962}"

# --- Reference --------------------------------------------------------------

# PanSN-named H37Rv, used as the reference PATH inside every graph.
: "${MTB_REF_FASTA:=${MTB_REFDIR}/GCF_000195955.pansn.fasta}"
# The exact path name as it appears in the GFA / .og. vg deconstruct -P wants this.
: "${MTB_REF_PATH:=GCF_000195955#1#NC_000962.3}"
# Plain H37Rv (non-PanSN header) for nucmer/MUMmer QC.
: "${MTB_H37RV:=${MTB_REFDIR}/H37Rv.fasta}"

# --- Outgroup ---------------------------------------------------------------
# THE PANEL GENOME THAT ROOTS EVERY TREE AND POLARISES EVERY VARIANT, named in
# one place. It was GCF_035581225 (M. canettii) written into P0, the
# association tail's tree and event steps, add_outgroup.py, panel_polarity.py
# and build_snp_tree.sh separately, so a graph without that genome could not
# be run without editing five files. P0 records the value it used in the
# build's build_info.tsv (key `outgroup`), and the association tail reads it
# from there, so a cohort is always rooted on its own build's outgroup.
# Set it EMPTY for a graph with no outgroup (P0 then makes no polarity table).
# The panel tree's own outgroup leaves for the AA tag are separate: P0 step
# ancestral, ANC_OUTGROUPS (decision D44).
#
# Leaves of a cohort tree that are outside the MTBC besides the outgroup: the
# panel's second canettii and the read-based canettii isolate. The
# association tail pins the MRCA of every other leaf -- the MTBC node -- to
# the panel's AA (R2-TREES-6); names not in a tree are ignored.
: "${MTB_NON_MTBC_TIPS=GCF_000253375,canettii}"
: "${MTB_OUTGROUP=GCF_035581225}"

# --- Container bind convention ---------------------------------------------
#
# All containerised steps bind  ${MTB_DATA} -> /data  and  ${MTB_GRAPHS} -> /graphs.
# Use mtb_bind() to build the flags, and mtb_in_data()/mtb_in_graphs() to convert
# a host path to its in-container equivalent.

: "${MTB_BIND_DATA:=/data}"
: "${MTB_BIND_GRAPHS:=/graphs}"

mtb_bind() {
    # /data and /graphs are the conventional mounts that pggb's -i/-o use.
    # MTB_WORK is additionally bound at its own path, so a host path passed to a
    # containerised tool resolves verbatim and no translation is needed.
    printf -- '-B %s:%s -B %s:%s -B %s:%s' \
        "${MTB_DATA}" "${MTB_BIND_DATA}" \
        "${MTB_GRAPHS}" "${MTB_BIND_GRAPHS}" \
        "${MTB_WORK}" "${MTB_WORK}"
}

# bcftools is not installed on the login nodes; the pggb container ships it.
# Prefer a host binary when one exists, else run the container's.
mtb_bcftools() {
    if command -v bcftools >/dev/null 2>&1; then
        bcftools "$@"
    else
        singularity exec $(mtb_bind) "${MTB_PGGB_SIF}" bcftools "$@"
    fi
}

# Same for tabix/bgzip, which travel with bcftools.
mtb_tabix() {
    if command -v tabix >/dev/null 2>&1; then tabix "$@"
    else singularity exec $(mtb_bind) "${MTB_PGGB_SIF}" tabix "$@"; fi
}

# host path under $MTB_DATA -> /data/...
mtb_in_data() {
    local p="$1"
    case "$p" in
        "${MTB_DATA}"/*) printf '%s' "${MTB_BIND_DATA}/${p#"${MTB_DATA}"/}" ;;
        *) echo "[project_env] ERROR: '$p' is not under MTB_DATA (${MTB_DATA})" >&2; return 1 ;;
    esac
}

# host path under $MTB_GRAPHS -> /graphs/...
mtb_in_graphs() {
    local p="$1"
    case "$p" in
        "${MTB_GRAPHS}"/*) printf '%s' "${MTB_BIND_GRAPHS}/${p#"${MTB_GRAPHS}"/}" ;;
        *) echo "[project_env] ERROR: '$p' is not under MTB_GRAPHS (${MTB_GRAPHS})" >&2; return 1 ;;
    esac
}

# --- Compute defaults -------------------------------------------------------

: "${MTB_THREADS:=${SLURM_CPUS_PER_TASK:-8}}"
: "${MTB_PARTITION:=sapphire}"

# --- Export ------------------------------------------------------------------
# Every MTB_* variable above is set with `: "${VAR:=default}"`, which assigns
# but does NOT export. Sourcing this file therefore gave a shell script `$MTB_PY`
# by ordinary expansion while leaving a CHILD process nothing: 39 variables were
# set and none reached a subprocess.
#
# That stayed invisible because the Python entry points read them as
# `os.environ.get("MTB_ODGI", "odgi")` and silently fall back to a bare command
# name. Where the fallback happens to be on PATH the run succeeds against a
# binary nobody chose; where it is not, the run dies with a bare
# FileNotFoundError. On 2026-09-22 is6110_project_sites.py failed exactly that
# way -- `odgi` is not on PATH -- while is6110_breaks_without_calls.py had been
# quietly using /n/home11/mchase/bin/minimap2 instead of the configured one.
# Both are 2.28-r1209, checked, so no earlier result is affected, but that was
# luck rather than design.
#
# Indirect prefix expansion rather than a hand-written list, because a list is
# the thing that goes stale: the display list in mtb_show_config() below names
# 23 of the 39 and was mistaken for an export more than once.
export ${!MTB_*}

# --- Guard rails ------------------------------------------------------------

# Refuse to run if scratch has been purged out from under us.
mtb_require_work() {
    if [[ ! -d "${MTB_WORK}" ]]; then
        echo "[project_env] FATAL: MTB_WORK does not exist: ${MTB_WORK}" >&2
        echo "               netscratch may have been purged. Re-create it and run:" >&2
        echo "               ${MTB_PERSIST}/bin/stage_in.sh" >&2
        return 1
    fi
}

mtb_require_file() {
    local f
    for f in "$@"; do
        if [[ ! -e "$f" ]]; then
            echo "[project_env] FATAL: required file missing: $f" >&2
            return 1
        fi
    done
}

# --- Build identity of products (audit section B, rerun safety) --------------
# A skip-if-exists guard that tests only existence keeps a previous build's
# output when the chain is rerun on a new graph (audit P0P2-1, TP-4, P3IS-8).
# These helpers let a guard ask WHICH build a product came from.

# The single P0 build to use: $MTB_BUILD_DIR, else the one build under
# $BUILD_ROOT whose manifest step completed. A half-built directory is never
# picked by default (audit P0P2-13), and two completed builds are refused.
mtb_resolve_build() {
    local root="${BUILD_ROOT:-refbias/build}" b
    if [[ -n "${MTB_BUILD_DIR:-}" ]]; then
        printf '%s\n' "${MTB_BUILD_DIR%/}"; return 0
    fi
    local -a cands=()
    for b in "${root}"/*/; do
        [[ -s "${b}logs/manifest.done" ]] && cands+=("${b%/}")
    done
    if [[ "${#cands[@]}" -ne 1 ]]; then
        echo "[project_env] FATAL: ${#cands[@]} completed builds (logs/manifest.done) under ${root}; set MTB_BUILD_DIR" >&2
        return 1
    fi
    printf '%s\n' "${cands[0]}"
}

# The ##MTB_graph_build stamp of a VCF (plain or bgzipped); empty if none.
mtb_vcf_build_id() {
    local f="$1"
    [[ -s "$f" ]] || return 0
    { if [[ "$f" == *.gz ]]; then gzip -cd "$f"; else cat "$f"; fi; } 2>/dev/null \
        | awk '/^#CHROM/{exit} /^##MTB_graph_build=/{sub(/^##MTB_graph_build=/,""); print; exit}' \
        || true
}

# A value from a key<TAB>value file (a .done marker or provenance sidecar);
# empty when the file or the key is absent.
mtb_kv() {
    [[ -s "$1" ]] || return 0
    awk -F'\t' -v k="$2" '$1==k{print $2; exit}' "$1"
}

# Print a script's leading comment block as its usage message.
# Immune to line-number drift, unlike `sed -n '11,17p' "$0"`.
mtb_usage() {
    local script="${1:-${BASH_SOURCE[1]}}"
    awk '
        NR==1 && /^#!/  { next }               # skip shebang
        /^#SBATCH/      { next }               # skip sbatch directives
        /^#/            { sub(/^# ?/, ""); print; seen=1; next }
        seen            { exit }               # stop at first non-comment after the block
        { next }
    ' "$script" >&2
}

# Print the resolved configuration (bin/show_config.sh calls this).
mtb_show_config() {
    local v
    for v in MTB_WORK MTB_PERSIST MTB_ARCHIVE MTB_DATA MTB_GRAPHS \
             MTB_CONTAINERS MTB_PGGB_SIF MTB_VG_SIF MTB_GRAPHALIGNER_SIF \
             MTB_REF_FASTA MTB_REF_PATH MTB_H37RV MTB_OUTGROUP MTB_NON_MTBC_TIPS \
             MTB_SNPEFF_JAR MTB_JAVA MTB_SNPEFF_DB MTB_PY MTB_PY_VT MTB_ODGI MTB_MINIMAP2 \
             MTB_SAMTOOLS MTB_BGZIP MTB_TABIX MTB_BEDTOOLS MTB_BCFTOOLS MTB_K8 MTB_PAFTOOLS \
             MTB_BLASTN MTB_MAKEBLASTDB MTB_UNIVEC \
             MTB_BWA MTB_WGSIM MTB_GATK_SIF MTB_DELLY_ENV \
             MTB_CRAM_ROOT MTB_CRAM_REF MTB_BUILD_DIR MTB_GRAPH_FRAMES \
             MTB_THREADS MTB_PARTITION; do
        printf '%-18s %s\n' "$v" "${!v}"
    done
}
