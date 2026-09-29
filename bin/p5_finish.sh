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

BUILD="${MTB_BUILD_DIR:-$(find refbias/build -mindepth 1 -maxdepth 1 -type d | head -1)}"
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
if [[ -d "$P4BDIR" ]]; then
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
#   --merge   the merged cohort VCF
#   --all     both, the default, for running it by hand
P5FSTEP="${1:-${P5FSTEP:---all}}"
case "$P5FSTEP" in
  --all|--pre|--merge) ;;
  *) echo "usage: $0 [--all|--pre|--merge]" >&2; exit 2 ;;
esac

if [[ "$P5FSTEP" != "--merge" ]]; then
echo "=== matrix ==="
bash bin/p5_merge.sh --matrix
echo
echo "=== lineage recovery (the pilot's endpoint check) ==="
"$MTB_PY" bin/p5_validate.py --matrix "${OUTDIR}/matrix.tsv" \
    --refmap "$REFMAP" --out "${OUTDIR}/validation.tsv"
echo
echo "=== sanity checks ==="
"$MTB_PY" bin/p5_sanity.py --matrix "${OUTDIR}/matrix.tsv" \
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
"$MTB_PY" bin/p6_annotate.py --build "$BUILD" \
    --matrix "${OUTDIR}/matrix.tsv" --refmap "$REFMAP" \
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

IS6110_KEYS="${IS6110_KEYS:-}"
if [[ -z "$IS6110_KEYS" ]]; then
    _kt="$_tag"; [[ "$_kt" == "pilot" ]] && _kt=""
    for _c in "is6110/results/${_kt}_p1i_cohort_keys.tsv" \
              "is6110/results/p1i_cohort_keys.tsv"; do
        [[ -s "$_c" ]] && { IS6110_KEYS="$_c"; break; }
    done
fi
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

"$MTB_PY" bin/merge_cohort_vcf.py \
    --matrix "${OUTDIR}/matrix.tsv" --sv-matrix "${OUTDIR}/sv_matrix.tsv" \
    ${IS6110_KEYS:+--is6110-keys "$IS6110_KEYS"} \
    ${IS6110_STATES:+--is6110-states "$IS6110_STATES"} \
    ${SV_STATES:+--sv-states "$SV_STATES"} \
    "${IV_ARGS[@]}" \
    "${ACC_ARGS[@]}" \
    --cohort-name "$COHORT_NAME" \
    --out "${OUTDIR}/merged.vcf.gz"

echo
echo "=== chain rebuilt from one generation ==="
