#!/usr/bin/env bash
#SBATCH --job-name=P5_finish
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH -t 0-08:00
# Eight hours, not one. At 997 isolates this step writes a 783 MB matrix
# and then validates, sanity-checks, clusters the SV arm and annotates on
# top of it; gwas1000 was killed at the one-hour limit having finished the
# matrix and nothing after it, which cost the whole tail of the chain. The
# earlier OUT_OF_MEMORY on the same job was a different fault, fixed in
# p5_matrix.py -- peak RSS is now 9.2 GB against the 32 GB requested.
#SBATCH -p shared
#SBATCH --mem=32000
# 32 GB, not 8. gwas1000's matrix step died OUT_OF_MEMORY at the 8 GB
# ceiling 98 seconds in: p5_matrix.py held every sample's whole state
# table as a dict of strings, which at 997 isolates and 173,685 keys is
# 173 million Python dict entries. That is fixed -- the matrix is now one
# uint8 array of 173 MB -- but the keys table and the per-row write still
# want room, and a job that dies here costs the whole tail of the chain.
#SBATCH --output=slurm/P5fin_%A.out
#SBATCH --error=slurm/P5fin_%A.err
#
# Assemble the matrix and run every downstream check, in one job, so the whole
# chain derives from ONE generation of the upstream stages.
#
# This exists because of a stale-input mistake worth not repeating: the matrix
# was once built from one generation of P4's output while P6 was annotated from a
# later one. The graph build id cannot catch that -- both generations carry the
# same build id -- so the protection has to be procedural: rebuild keys, states,
# matrix and every check together, never piecemeal.
set -euo pipefail

_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
source "$_mtb_env"

# Exactly one build, or MTB_BUILD_DIR. `find | head -1` picked whichever build
# the filesystem listed first when there were several, while p4_place.sh
# refused -- so two passes of one chain could run against different builds.
BUILD="${MTB_BUILD_DIR:-}"
if [[ -z "$BUILD" ]]; then
    mapfile -t _c < <(find "${BUILD_ROOT:-refbias/build}" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
    [[ "${#_c[@]}" -eq 1 ]] || { echo "FATAL: set MTB_BUILD_DIR (${#_c[@]} builds)" >&2; exit 1; }
    BUILD="${_c[0]}"
fi
export MTB_BUILD_DIR="$BUILD"

# Paths were hard-coded to the pilot layout, so this script silently checked
# the pilot's freshness while finishing another cohort. They are variables now,
# with the pilot as the default, and bin/refbias_run.sh sets them per cohort.
P1DIR="${P1DIR:-refbias/p1}"
P2DIR="${P2DIR:-refbias/p2}"
P4DIR="${P4DIR:-refbias/p4}"
P4BDIR="${P4BDIR:-refbias/p4b}"
OUTDIR="${OUTDIR:-refbias/p5}"
REFMAP="${REFMAP:-${P1DIR}/refmap.tsv}"
COHORT="${COHORT:-refbias/cohort.pilot.tsv}"
export P1DIR P2DIR P4DIR P4BDIR OUTDIR REFMAP COHORT

# the step, read here as well as below: the freshness guard depends on it
P5FSTEP="${1:-${P5FSTEP:---all}}"

# --- freshness guard --------------------------------------------------------
# EVERY find HERE IS -maxdepth 1. Without it they recurse, and refbias/p5
# contains with_is6110/, a subdirectory of state files from an earlier
# IS6110-augmented run. Those are legitimately old, so the guard counted 23 of
# them as stale states of the CURRENT cohort and refused to build a matrix
# that was in fact correct -- keys at 12:43:33 against states from 12:55
# onward. A guard that blocks a good run is as costly as one that passes a bad
# one, and this one was silent about which directory the offending files were
# in, which is what made it slow to diagnose.
# This job once ran after its upstream array was CANCELLED, because the
# dependency was afterany, and rebuilt the matrix from keys and states files that
# were a day old. The output was byte-identical to the pre-fix matrix, which is
# exactly what a stale rebuild looks like: nothing errors, and the numbers simply
# do not move. Dependencies are now afterok, and this checks anyway.
newest_placed=$(find "$P4DIR" -maxdepth 1 -name '*.placed.tsv' -newer "${OUTDIR}/keys.tsv" 2>/dev/null | wc -l)
if [[ "$newest_placed" -gt 0 ]]; then
    echo "FATAL: ${newest_placed} P4 outputs are NEWER than ${OUTDIR}/keys.tsv;" >&2
    echo "       the key set is stale -- rerun --keys before finishing" >&2
    exit 1
fi
# P4b was not in this guard, and that is how both cohorts' sv_placed.tsv files
# sat two days behind bin/p4b_place_sv.py -- missing the entire `inherited`
# component, 93 of 159 rows per isolate -- with nothing reporting it. It was
# not caught downstream either, because the SV matrix was built by hand and
# this script did not run it. It does now, so a stale P4b is FATAL here rather
# than a warning: the matrix would otherwise be built on 41% of the input.
# Compared against the script's own mtime, because P4b has no later stage to
# be older than.
# Only the steps that READ P4b need this: --pre builds sv_matrix.tsv from it.
# The merge reads sv_matrix.tsv, not P4b, so refusing it here blocked a VCF
# rebuild for every cohort whenever p4b_place_sv.py changed, even when the SV
# matrix was deliberately being kept.
if [[ -d "$P4BDIR" && ( "$P5FSTEP" == "--all" || "$P5FSTEP" == "--pre" ) ]]; then
    stale_p4b=$(find "$P4BDIR" -maxdepth 1 -name '*.sv_placed.tsv' \
        ! -newer bin/p4b_place_sv.py 2>/dev/null | wc -l)
    if [[ "$stale_p4b" -gt 0 ]]; then
        echo "FATAL: ${stale_p4b} P4b outputs in ${P4BDIR} are OLDER than" >&2
        echo "       bin/p4b_place_sv.py, so they predate the current code." >&2
        echo "       Rerun P4b: bash bin/refbias_run.sh <cohort> --only p4b" >&2
        exit 1
    fi
fi

stale_states=$(find "$OUTDIR" -maxdepth 1 -name '*.states.tsv' ! -newer "${OUTDIR}/keys.tsv" 2>/dev/null | wc -l)
if [[ "$stale_states" -gt 0 ]]; then
    echo "FATAL: ${stale_states} state files are OLDER than keys.tsv;" >&2
    echo "       they were computed against a different key set" >&2
    exit 1
fi
n_states=$(ls "${OUTDIR}"/*.states.tsv 2>/dev/null | wc -l)
n_expect=$(( $(wc -l < "$REFMAP") - 1 ))
if [[ "$n_states" -ne "$n_expect" ]]; then
    echo "FATAL: ${n_states} state files but ${n_expect} cohort samples" >&2
    exit 1
fi
echo "freshness: keys newer than all P4 output, ${n_states} state files all newer than keys"
echo

# TWO STEPS NOW, AND THE SPLIT IS LOAD-BEARING. Everything up to and including
# P6 annotation builds the tables; the merged VCF consumes them. Between the two
# sit a per-sample array (pass p5svgt, which measures deletion absence) and the
# IS6110 arm's own P5 stages, neither of which can run inside this job. Before
# the split, p5finish depended only on p5states, so a fresh cohort built its
# merged VCF BEFORE the IS6110 arm had produced anything and the insertion
# sites were simply missing from it.
#
#   --pre     matrix, validation, sanity, SV matrix, P6 annotation
#   --merge   the merged cohort VCF; with VCF_SHARDS > 1, every shard in turn
#             and then the assembly, in this one job (for running by hand)
#   --merge-shard     one shard of the merged VCF: shard SLURM_ARRAY_TASK_ID-1
#                     of VCF_SHARDS (or VCF_SHARD), written under vcf_parts/
#   --merge-assemble  join the VCF_SHARDS shards into merged.vcf.gz, then
#                     stamp and gate it
#   --all     both, the default, for running it by hand
P5FSTEP="${1:-${P5FSTEP:---all}}"
case "$P5FSTEP" in
  --all|--pre|--merge|--merge-shard|--merge-assemble) ;;
  *) echo "usage: $0 [--all|--pre|--merge|--merge-shard|--merge-assemble]" >&2; exit 2 ;;
esac

if [[ "$P5FSTEP" == "--all" || "$P5FSTEP" == "--pre" ]]; then
echo "=== matrix ==="
bash bin/p5_merge.sh --matrix
echo
echo "=== lineage recovery (the pilot's endpoint check) ==="
# Every cohort consumer below reads the memory-mapped states array and the
# per-site table that --matrix just wrote, not a dense keys x samples text
# matrix: at 10,000 isolates that matrix is about 23 GB, and the sanity check
# and the old merge each held all of it in memory.
ARR=(--states-array "$OUTDIR" --keys "${OUTDIR}/keys.tsv" --sites "${OUTDIR}/sites.tsv")
"$MTB_PY" bin/p5_validate.py "${ARR[@]}" \
    --refmap "$REFMAP" --out "${OUTDIR}/validation.tsv"
echo
echo "=== sanity checks ==="
"$MTB_PY" bin/p5_sanity.py "${ARR[@]}" \
    --refmap "$REFMAP" --cohort "$COHORT" \
    --p2-summary "${P2DIR}/p2_summary.tsv" --out "${OUTDIR}/sanity.tsv"
echo
echo "=== SV matrix (the structural-variant arm) ==="
# Kept a separate table from the small-variant matrix on purpose: SV rows key
# on tolerance clustering rather than an exact allele, they never assert REF
# because no-call is mostly not-detected, and they carry a per-caller
# confidence band. bin/p5_sv_matrix.py says why at length.
"$MTB_PY" bin/p5_sv_matrix.py --refmap "$REFMAP" --dir "$P4BDIR" \
    --calibration "${SV_CALIBRATION:-refbias/STAGE11.calibration.tsv}" \
    --out "${OUTDIR}/sv_matrix.tsv"
echo
echo "=== P6 annotation ==="
# --out is passed explicitly. p6_annotate.py defaults it to refbias/p6, which
# is the pilot's path, so every other cohort silently overwrote the pilot's
# table -- or, as happened on 2026-09-23, the pilot overwrote scale100's.
P6DIR="${P6DIR:-refbias/p6}"
mkdir -p "$P6DIR"
# --p4-dir is passed too. It was not, so P6 read its default, the PILOT's
# refbias/p4: for any other cohort every off-path site was annotated from the
# pilot's carriers, or -- where that directory was absent -- from none at all
# (11,533 scale200 sites came out offpath_no_carrier).
"$MTB_PY" bin/p6_annotate.py --build "$BUILD" \
    --sites "${OUTDIR}/sites.tsv" --refmap "$REFMAP" --p4-dir "$P4DIR" \
    --out "${P6DIR}/annotated.tsv"
fi   # end --pre

if [[ "$P5FSTEP" == "--pre" ]]; then
    echo "--pre complete; run --merge after p5svgt and the IS6110 P5 stages"
    exit 0
fi

echo
echo "=== merged cohort VCF ==="
# The cohort-level deliverable. Built here rather than as its own pass because
# it is a transformation of files this step has just written, and a separate
# pass would be one more dependency to wire wrong. The IS6110 key table comes
# from pass p1iv, which may not have run yet -- the merge simply omits that
# class if the table is absent, and says so in its own output.
# ONE tag, derived once, used by everything below. It used to be derived
# separately per lookup, and COHORT_NAME was derived differently again as
# basename(OUTDIR) -- which is the literal "p5" for every cohort, since every
# outroot ends in /p5. That stamped ##cohort=p5 on the VCF and sent
# merge_cohort_vcf.py looking for a "p5" cohort's IS6110 stage-2 table.
_tag="$(basename "$(dirname "$OUTDIR")")"
[[ "$_tag" == "refbias" ]] && _tag="pilot"
COHORT_NAME="${COHORT_NAME:-$_tag}"
echo "  cohort name: ${COHORT_NAME}"

# NOT resolved here either, for the same reason as IS6110_STATES below. This
# loop used to fall back to the unprefixed is6110/results/p1i_cohort_keys.tsv
# -- the PILOT's table -- for any cohort without its own, and because it then
# passed --is6110-keys explicitly it bypassed the guard in merge_cohort_vcf.py
# written to stop exactly that. merge_cohort_vcf.py maps the cohort name to its
# own table and emits no IS6110 block, loudly, when there is none.
IS6110_KEYS="${IS6110_KEYS:-}"
# The two state tables that turn presence-only blocks into genotypes. Both are
# optional and the merge says so loudly when one is absent, because without
# them the SV and IS6110 blocks carry no REF at all and nothing downstream that
# has to contrast presence with absence can use them.
# p5_svgt.sh writes svgt_iv_states.tsv in interval mode and svgt_states.tsv
# in the original matrix mode; sv_states.tsv is a name nothing has ever
# written, so this lookup silently found nothing for every cohort. Interval
# mode is preferred where both exist, because it genotypes every sample at
# every interval rather than trusting caller breakpoints.
# THE INTERVAL CATALOGUE, which supersedes the caller deletions when present.
# Until 2026-09-27 the merge read only sv_states.tsv, keyed on caller clusters,
# while the rebuilt arm wrote svgt_iv_states.tsv keyed on catalogued intervals
# -- two key spaces with nothing in common, so the interval arm never reached
# the VCF at all. Both are passed now: the catalogue supplies the deletions and
# the caller matrix still supplies the insertions it cannot genotype.
IVTAB="${IVTAB:-refbias/assets/sv_intervals.tsv}"
IVSTATES="${IVSTATES:-${OUTDIR}/svgt_iv_states.tsv}"
IV_ARGS=()
if [[ -s "$IVTAB" && -s "$IVSTATES" ]]; then
    IV_ARGS=(--sv-intervals "$IVTAB" --sv-interval-states "$IVSTATES")
    echo "  deletion block from ${IVTAB} genotyped in ${IVSTATES}"
else
    echo "  no interval states at ${IVSTATES}; the deletion block falls back to"
    echo "  caller clusters, which are 80.4% singletons at a 10 bp tolerance"
fi

# SV_STATES OVERLAYS THE CALLER MATRIX and is therefore keyed on `key`, the
# caller cluster id. svgt_iv_states.tsv is keyed on `interval` and belongs to
# --sv-interval-states above; passing it here raised KeyError: 'key', because
# the earlier version of this resolver preferred it. The two key spaces have
# nothing in common and must not be crossed.
if [[ -z "${SV_STATES:-}" ]]; then
    for _c in "${OUTDIR}/svgt_states.tsv" "${OUTDIR}/sv_states.tsv"; do
        [[ -s "$_c" ]] && { SV_STATES="$_c"; break; }
    done
fi
[[ -s "${SV_STATES:-}" ]] || SV_STATES=""
[[ -n "$SV_STATES" ]] && echo "  SV states: ${SV_STATES}"
# NOT resolved here. merge_cohort_vcf.py maps a cohort name to its stage-2
# table itself, and deliberately refuses to fall back to the unprefixed pilot
# file for a cohort that has none -- a shell loop that tried the pilot's table
# second would reintroduce exactly the bug that comment guards against, which
# is joining the 23-isolate pilot's states onto a different cohort.
IS6110_STATES="${IS6110_STATES:-}"

# LEVEL 1, the accessory presence block. One biallelic record per accessory
# locus, every sample stated, from accessory/bin/locus_presence.py. It is
# separate from the node-frame small variants inside those loci because those
# are a different character: for a sample that does not carry the insert, a
# variant inside it is inapplicable, not reference and not unknown.
ACCPRES="${ACCPRES:-accessory/${COHORT_NAME}}"
ACC_ARGS=()
if compgen -G "${ACCPRES}/*.presence.tsv" >/dev/null; then
    ACC_ARGS=(--accessory-presence "$ACCPRES")
    echo "  level-1 accessory presence from ${ACCPRES}" \
         "($(ls "${ACCPRES}"/*.presence.tsv | wc -l) tables)"
else
    echo "  no level-1 presence tables at ${ACCPRES}; the accessory arm will"
    echo "  carry within-insert variation only, with no presence character"
fi

# REF bases and the build id come from the build, not from the merge's
# defaults: symbolic records used to carry REF=N, and --build-id defaulted to a
# hard-coded id that any rebuild would have silently kept.
H37RV_FASTA="${BUILD}/refs/GCF_000195955.fasta"
[[ -s "$H37RV_FASTA" ]] || { echo "FATAL: no H37Rv FASTA at ${H37RV_FASTA}" >&2; exit 1; }
BUILD_ID="$(awk -F'\t' '$1=="build_id"{print $2}' "${BUILD}/build_info.tsv")"
[[ -n "$BUILD_ID" ]] || { echo "FATAL: no build_id in ${BUILD}/build_info.tsv" >&2; exit 1; }

# The small-variant block comes from the states array --matrix wrote. A cohort
# finished before the array existed still has matrix.tsv, which is used then.
if [[ -s "${OUTDIR}/states.u8.npy" && -s "${OUTDIR}/states.meta.tsv" ]]; then
    SMALL=(--states-array "$OUTDIR" --keys "${OUTDIR}/keys.tsv")
elif [[ -s "${OUTDIR}/matrix.tsv" ]]; then
    echo "  no states array in ${OUTDIR}; reading the legacy matrix.tsv"
    SMALL=(--matrix "${OUTDIR}/matrix.tsv" --keys "${OUTDIR}/keys.tsv")
else
    echo "FATAL: neither states.u8.npy nor matrix.tsv in ${OUTDIR}; run --pre" >&2
    exit 1
fi
MERGE=("$MTB_PY" bin/merge_cohort_vcf.py "${SMALL[@]}" \
    --sv-matrix "${OUTDIR}/sv_matrix.tsv" \
    --h37rv-fasta "$H37RV_FASTA" --build-id "$BUILD_ID" \
    ${IS6110_KEYS:+--is6110-keys "$IS6110_KEYS"} \
    ${IS6110_STATES:+--is6110-states "$IS6110_STATES"} \
    ${SV_STATES:+--sv-states "$SV_STATES"} \
    "${IV_ARGS[@]}" \
    "${ACC_ARGS[@]}" \
    --cohort-name "$COHORT_NAME")

# SHARDED BY GENOME REGION. One process holds every record's cells before it
# writes: 8.5 GB at 997 isolates, about 260 GB projected at 10,000. Each shard
# owns a contiguous H37Rv range and a contiguous range of node contigs, loads
# only its own rows of every input, and writes record-only parts; the assembly
# streams them after one header. The assembled file is byte-identical to the
# single-process one (tests/run_tests.py checks it).
VCF_SHARDS="${VCF_SHARDS:-1}"
PARTS="${OUTDIR}/vcf_parts/s"
mkdir -p "${OUTDIR}/vcf_parts"
case "$P5FSTEP" in
  --merge-shard)
    I="${VCF_SHARD:-$(( ${SLURM_ARRAY_TASK_ID:?--merge-shard needs an array task or VCF_SHARD} - 1 ))}"
    "${MERGE[@]}" --n-shards "$VCF_SHARDS" --shard "$I" --part-prefix "$PARTS" \
        --out /dev/null
    exit 0 ;;
  --merge-assemble|--merge|--all)
    if [[ "$VCF_SHARDS" -gt 1 ]]; then
        if [[ "$P5FSTEP" != "--merge-assemble" ]]; then
            for ((I = 0; I < VCF_SHARDS; I++)); do
                "${MERGE[@]}" --n-shards "$VCF_SHARDS" --shard "$I" \
                    --part-prefix "$PARTS" --out /dev/null
            done
        fi
        # Parts from an earlier generation must not be assembled with this one.
        # Only testable against the states array: a cohort still on the legacy
        # matrix has no states.meta.tsv, and `find -newer <missing file>` fails,
        # which under set -e and pipefail ended this job silently.
        if [[ -s "${OUTDIR}/states.meta.tsv" ]]; then
            stale=$(find "${OUTDIR}/vcf_parts" -maxdepth 1 -name 's.*.meta.json' \
                    ! -newer "${OUTDIR}/states.meta.tsv" | wc -l)
            if [[ "$stale" -gt 0 ]]; then
                echo "FATAL: ${stale} VCF shard parts predate the states array; rerun the shards" >&2
                exit 1
            fi
        fi
        "${MERGE[@]}" --assemble --n-shards "$VCF_SHARDS" --part-prefix "$PARTS" \
            --out "${OUTDIR}/merged.vcf.gz"
    else
        "${MERGE[@]}" --out "${OUTDIR}/merged.vcf.gz"
    fi ;;
esac

# STAMP THE GRAPH PROVENANCE. bin/stamp_build_id.sh existed and was called by
# p2_call.sh for the per-sample caller VCFs, but never here -- so the cohort
# deliverable, the one file anybody actually uses, carried no record of which
# graph produced its coordinates. Both merged VCFs in the tree on 2026-09-29 had
# zero ##MTB headers.
#
# This matters more now that the graph is built by a separate repository: its
# identity is no longer anywhere in this repository's history, and the path
# recorded in build_info.tsv points into scratch, which is purged. The stamp
# carries the pggb fingerprint parsed from the graph filename, which is the part
# that actually identifies the graph.
bash bin/stamp_build_id.sh "${OUTDIR}/merged.vcf.gz"

# THE DELIVERABLE MUST BE READABLE BY STANDARD TOOLS, checked here and not by a
# consumer months later. bin/vcf_gate.sh fails on a GT index above the ALT
# count, a record bcftools cannot parse, a REF that does not match H37Rv, or an
# unsorted contig -- each of which the pipeline's own readers tolerated.
bash bin/vcf_gate.sh "${OUTDIR}/merged.vcf.gz" "$H37RV_FASTA"

echo
echo "=== chain rebuilt from one generation ==="
