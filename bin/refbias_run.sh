#!/usr/bin/env bash
# One entry point for the short-read passes P1 to P5.
#
# WHY THIS EXISTS. Every pass was launched by hand with five or six environment
# variables that have to agree with each other, and nothing checked that they
# did. Two bugs of exactly that class were hit on 2026-09-21 in the IS6110
# work: a glob picked up a second reference frame's files, and the mismatch
# showed up only as agreement collapsing from 0.87 to 0.016, never as an error.
# A wrong path in this pipeline does not fail, it produces a plausible number.
# See SEGMENT5_SHORT_READ_PASSES.md section 6.
#
# WHAT IT GUARANTEES.
#   * every pass in a chain gets the SAME cohort, refmap, build and graph;
#   * the graph comes from the build directory's own build_info.tsv rather than
#     from a glob over a hard-coded directory name, which is what six of the
#     pass scripts still fall back to;
#   * dependencies are afterok, never afterany -- p5_finish.sh exists because a
#     matrix was once rebuilt from a day-old key set after its upstream array
#     was cancelled, and nothing errored;
#   * the structural-variant arm is part of the chain: p5_finish.sh builds
#     refbias/<root>/p5/sv_matrix.tsv from P4b, which for two days nothing did
#     -- see refbias/archive/README.md;
#   * every pass's --summary step runs after its array, because those steps are
#     not reporting extras: P1's writes refmap.tsv, which every later pass
#     reads, and P2's writes p2_summary.tsv, which p5_finish.sh reads. The
#     first version of this script omitted them and P2 stopped the chain with
#     "no P1 refmap"; afterok then cancelled the other six jobs, which is the
#     behaviour wanted but was found the hard way;
#   * --dry-run prints the exact sbatch lines and submits nothing;
#   * --throttle N caps how many tasks of one pass run at once (default 50,
#     0 for unthrottled), so a large cohort does not take a shared partition
#     for itself.
#
#   bash bin/refbias_run.sh --list
#   bash bin/refbias_run.sh pilot --dry-run
#   bash bin/refbias_run.sh pilot
#   bash bin/refbias_run.sh scale100 --from p3
#   bash bin/refbias_run.sh scale1000 --throttle 100
#   bash bin/refbias_run.sh pilot --only p4b
#   bash bin/refbias_run.sh pilot --status
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
for _c in "${MTB_ENV_FILE:-}" "${HERE}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { source "$_c"; break; }
done

REGISTRY="${REGISTRY:-refbias/cohorts.tsv}"
BUILD_ROOT="${BUILD_ROOT:-refbias/build}"
# p1g, p1i and p1iv are the IS6110 insertion-site arm. They were run by hand
# for months and the inventory shows exactly what that costs: the fixed arm
# exists for three cohorts, the matched arm for two, the VCFs for two, and no
# cohort has all three. A step that must be remembered is a step that will be
# forgotten, so they are passes now. They sit after p1 because they need its
# refmap, and they are independent of p2 to p5, which is why the chain does not
# serialise them behind it.
# THE ORDER CHANGED, AND THE REASON IS THE MERGED VCF. p5 used to end with
# p5_finish.sh building merged.vcf.gz, depending on nothing but p5states, so a
# fresh cohort wrote its cohort-level deliverable BEFORE the IS6110 arm had run
# and before any deletion had been genotyped for absence. p5 now stops after
# the tables (`--pre`); p1is earns the IS6110 REFs, p5svgt measures deletion
# absence, and p5vcf builds the VCF last, depending on all three.
ORDER=(p1 p2 p3 p4 p4b p5 p1g p1i p1iv p1is p5svgt p5vcf)

die() { echo "FATAL: $*" >&2; exit 1; }

# --- registry ---------------------------------------------------------------
list_cohorts() {
    printf "%-10s %-34s %-18s %s\n" cohort cohort_table outroot note
    awk -F'\t' '!/^#/ && $1!="cohort" && NF>=6 {printf "%-10s %-34s %-18s %s\n",$1,$2,$4,$6}' "$REGISTRY"
}
lookup() {  # cohort field -> value
    awk -F'\t' -v c="$1" -v f="$2" '!/^#/ && $1!="cohort" && $1==c {print $f; found=1}
        END{ if(!found) exit 3 }' "$REGISTRY"
}

COHORT_NAME="" FROM="" TO="" ONLY="" DRY=0 STATUS=0
# How many array tasks of one pass may run at once, appended to the array spec
# as sbatch's %N. Every pass here is embarrassingly parallel and short, so a
# thousand-isolate cohort would otherwise queue a thousand tasks at once in
# each of six passes -- which does not break anything, but lands the whole
# fair-share cost of one cohort in one instant on a shared partition and
# leaves no room for anyone else's work, including this project's own reruns.
# Set to 0 to submit unthrottled.
THROTTLE="${THROTTLE:-50}"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --list)     [[ -s "$REGISTRY" ]] || die "no registry at $REGISTRY"
                    list_cohorts; exit 0 ;;
        --from)     FROM="$2"; shift 2 ;;
        --to)       TO="$2"; shift 2 ;;
        --only)     ONLY="$2"; shift 2 ;;
        --dry-run)  DRY=1; shift ;;
        --throttle) THROTTLE="$2"; shift 2 ;;
        --status)   STATUS=1; shift ;;
        # README and INPUTS.md document `--cohort <name>`; the positional form
        # is kept as well so existing invocations do not break.
        --cohort)   [[ -z "$COHORT_NAME" ]] || die "one cohort at a time"
                    COHORT_NAME="$2"; shift 2 ;;
        -h|--help)  sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)         die "unknown option $1" ;;
        *)          [[ -z "$COHORT_NAME" ]] || die "one cohort at a time"
                    COHORT_NAME="$1"; shift ;;
    esac
done
[[ -n "$COHORT_NAME" ]] || { list_cohorts; die "name a cohort"; }
[[ -s "$REGISTRY" ]] || die "no registry at $REGISTRY"

COHORT="$(lookup "$COHORT_NAME" 2)" || die "cohort '${COHORT_NAME}' is not in ${REGISTRY}"
CRAMS="$(lookup "$COHORT_NAME" 3)"
OUTROOT="$(lookup "$COHORT_NAME" 4)"
PASSES="$(lookup "$COHORT_NAME" 5)"
WORKPFX="$(lookup "$COHORT_NAME" 7)"
[[ -n "$WORKPFX" ]] || die "cohort '${COHORT_NAME}' has no workprefix in ${REGISTRY}"
[[ -s "$COHORT" ]] || die "cohort table ${COHORT} does not exist"
[[ -s "$CRAMS" ]]  || die "CRAM table ${CRAMS} does not exist"
# Every per-sample task reads row TASK_ID+1 of the raw file, so N must count
# exactly the rows after the header, and nothing but trailing blank lines may
# sit among them. `grep -vc '^#'` counted blank lines, so a trailing newline
# gave N+1 tasks and the last one failed with "no cohort row", which afterok
# turned into a cancelled chain. A blank or comment line BETWEEN rows is worse:
# every later task would read its neighbour's sample, so it is refused.
N="$(awk -F'\t' 'NR==1{next} $1!="" && $1!~/^#/{last=NR; n++} END{print n+0}' "$COHORT")"
_gap="$(awk -F'\t' 'NR==1{next} $1!="" && $1!~/^#/{if(bad){print bad; exit} next}
                    {if(!bad) bad=NR}' "$COHORT")"
[[ -z "$_gap" ]] || die "cohort table ${COHORT} has a blank or comment line at line ${_gap} between sample rows; task indices would shift"
[[ "$N" -gt 0 ]] || die "cohort table ${COHORT} has no rows"

# --- build and graph, from the build record and not from a glob -------------
# Not the newest by mtime: re-running any P0 step on an old build rewrites a
# file in it and made that build "newest", so a chain would run against the
# old graph. Use MTB_BUILD_DIR when set; otherwise require exactly one build
# whose manifest step completed, and refuse to guess between several.
if [[ -n "${MTB_BUILD_DIR:-}" ]]; then
    BUILD="$MTB_BUILD_DIR"
else
    mapfile -t _builds < <(for _b in "${BUILD_ROOT}"/*/; do
        [[ -s "${_b}logs/manifest.done" ]] && printf '%s\n' "${_b%/}"; done)
    case "${#_builds[@]}" in
        1) BUILD="${_builds[0]}" ;;
        0) die "no completed build (logs/manifest.done) under ${BUILD_ROOT}; run bin/p0_prepare.sh or set MTB_BUILD_DIR" ;;
        *) die "${#_builds[@]} completed builds under ${BUILD_ROOT} (${_builds[*]}); set MTB_BUILD_DIR to choose one" ;;
    esac
fi
BUILD="${BUILD%/}"
[[ -d "$BUILD" ]] || die "no build directory ${BUILD}"
INFO="${BUILD}/build_info.tsv"
[[ -s "$INFO" ]] || die "no build_info.tsv in ${BUILD}; run bin/p0_prepare.sh --step stamp"
OG="$(awk -F'\t' '$1=="graph"{print $2}' "$INFO")"
BUILD_ID="$(awk -F'\t' '$1=="build_id"{print $2}' "$INFO")"
[[ -s "$OG" ]] || die "the graph recorded in ${INFO} is not readable: ${OG}"

REFMAP="${OUTROOT}/p1/refmap.tsv"

# --- which passes ------------------------------------------------------------
declare -a RUN=()
if [[ -n "$ONLY" ]]; then
    # Accept a comma-separated list, not just one pass. The tail of a chain
    # is three passes -- p1is, p5svgt, p5vcf -- and passing them as one string
    # produced `no rule for pass 'p1is,p5svgt,p5vcf'` after the submission
    # banner had already printed, which reads like the passes do not exist.
    IFS=',' read -r -a RUN <<< "$ONLY"
else
    started=0
    for p in "${ORDER[@]}"; do
        [[ -n "$FROM" && "$p" != "$FROM" && "$started" -eq 0 ]] && continue
        started=1
        RUN+=("$p")
        [[ -n "$TO" && "$p" == "$TO" ]] && break
    done
    [[ ${#RUN[@]} -gt 0 ]] || die "--from '${FROM}' is not one of: ${ORDER[*]}"
fi
# The registry's `passes` column now FILTERS rather than merely warns. It used
# to print a note and run the pass anyway, which made the column decorative:
# `refbias_run.sh l49` would have submitted p1 through p5 for a cohort
# registered for the IS6110 arm alone, and the only sign would have been a note
# on stderr. An explicit --only or --from still overrides, because that is a
# deliberate instruction from whoever typed it.
# --from and --to select a RANGE within the registered set; only --only
# overrides the registry. Letting --from disable the filter meant
# `refbias_run.sh scale100_h37rv --from p3` submitted the three IS6110 passes
# for a cohort registered p2..p5, which then had to be cancelled by hand.
# "Start later in the chain" is not the same instruction as "ignore what this
# cohort is for".
if [[ -z "$ONLY" ]]; then
    declare -a SEL=()
    for p in "${RUN[@]}"; do
        if [[ ",${PASSES}," == *",${p},"* ]]; then
            SEL+=("$p")
        else
            echo "  skipping '${p}': not in the registry's passes for ${COHORT_NAME}" >&2
        fi
    done
    [[ ${#SEL[@]} -gt 0 ]] || die "none of ORDER is listed in passes for ${COHORT_NAME} (${PASSES})"
    RUN=("${SEL[@]}")
else
    for p in "${RUN[@]}"; do
        [[ ",${PASSES}," == *",${p},"* ]] || \
            echo "  note: '${p}' is not in the registry's passes for ${COHORT_NAME} (${PASSES}); running it anyway because it was named explicitly" >&2
    done
fi

cat <<EOS
cohort   : ${COHORT_NAME}   ${N} isolates
table    : ${COHORT}
crams    : ${CRAMS}
outroot  : ${OUTROOT}
workpfx  : ${WORKPFX}
refmap   : ${REFMAP}
build    : ${BUILD}   id ${BUILD_ID}
graph    : ${OG}
passes   : ${RUN[*]}
EOS

# --- status ------------------------------------------------------------------
if [[ "$STATUS" -eq 1 ]]; then
    echo
    printf "%-5s %-26s %8s %8s\n" pass outdir present expected
    for p in "${ORDER[@]}"; do
        d="${OUTROOT}/${p}"
        # One pattern per pass, and each must match exactly one file per
        # sample. P2 writes both <sample>.vcf.gz and <sample>.g.vcf.gz, so a
        # bare *.vcf.gz counts every isolate twice -- the same glob mistake
        # that cost two debugging rounds in the IS6110 work on 2026-09-21.
        # Globs, not `ls | wc -l`: under pipefail an empty directory made ls
        # fail and --status exited at the first pass with no outputs. And
        # every pass needs its own case -- p1g onward used to fall through and
        # reuse p5's pattern.
        skip=""
        case "$p" in
            p1)   pat="*.candidates.tsv" ;;
            p2)   pat="*.vcf.gz" ; skip=".g.vcf.gz" ;;
            p3)   pat="*.copynumber.tsv" ;;
            p4)   pat="*.placed.tsv" ;;
            p4b)  pat="*.sv_placed.tsv" ;;
            p5)   pat="*.states.tsv" ;;
            p1g|p1i) pat="*.junctions.tsv" ;;
            *)    printf "%-5s %-26s %8s %8s\n" "$p" "(cohort-level pass)" "-" "-"
                  continue ;;
        esac
        n=0
        shopt -s nullglob
        for _f in "$d"/$pat; do
            [[ -n "$skip" && "$_f" == *"$skip" ]] && continue
            n=$((n + 1))
        done
        shopt -u nullglob
        if [[ "$n" -eq "$N" ]]; then note=""
        elif [[ "$n" -lt "$N" ]]; then note="   <-- $((N-n)) missing"
        else note="   <-- $((n-N)) MORE than the cohort; check the pattern"; fi
        printf "%-5s %-26s %8s %8s%s\n" "$p" "$d" "$n" "$N" "$note"
    done
    exit 0
fi

# --- submit ------------------------------------------------------------------
# EVERY cross-pass path, exported to every pass. Passing only OUTDIR and WORK
# was a real bug: p3 reads P2WORK, p4 reads P1WORK and P2DIR, p4b reads P2DIR,
# and all of those default to the PILOT's directories. So a non-default cohort
# silently read the pilot's inputs. It surfaced on 2026-09-21 when a scale100
# P4b re-run produced header-only files for the 77 isolates that are not in the
# pilot, and real output for the 23 that are -- no error, just a smaller
# number. Set them all here so no pass can fall back to a default.
DIRS="P1DIR=${OUTROOT}/p1,P2DIR=${OUTROOT}/p2,P3DIR=${OUTROOT}/p3"
DIRS="${DIRS},P4DIR=${OUTROOT}/p4,P4BDIR=${OUTROOT}/p4b,P5DIR=${OUTROOT}/p5"
# P6DIR is exported although no pass in ORDER writes it: p5_finish.sh runs the
# P6 annotation as its last step, and bin/p6_annotate.py defaults its --out to
# a hard-coded refbias/p6/annotated.tsv. Every cohort therefore wrote the same
# file, and the pilot re-run of 2026-09-23 replaced the scale100 generation --
# 61,049 rows became 40,444 -- with nothing failing and sync_back committing
# the result. Scoping it to the cohort's own outroot is the fix.
DIRS="${DIRS},P6DIR=${OUTROOT}/p6"
DIRS="${DIRS},P1GDIR=${OUTROOT}/p1g,P1IDIR=${OUTROOT}/p1i"
DIRS="${DIRS},P1IVCF=${OUTROOT}/p1i/vcf"
DIRS="${DIRS},P1WORK=${WORKPFX}p1,P2WORK=${WORKPFX}p2"
# The panel-vs-refs frame table is a build asset (p0_prepare.sh --step frames).
# Older builds predate that step; they fall back to graph_frame.py's default
# table, with a warning, rather than failing a chain that ran before.
FRAMES="${MTB_GRAPH_FRAMES:-${BUILD}/assets/graph_frame_offsets.tsv}"
if [[ ! -s "$FRAMES" ]]; then
    echo "WARNING: no frame table at ${FRAMES}; run bin/p0_prepare.sh --step frames." >&2
    echo "         Falling back to graphframe/results/graph_frame_offsets.tsv." >&2
    FRAMES="${HERE}/graphframe/results/graph_frame_offsets.tsv"
    [[ -s "$FRAMES" ]] || die "no frame table there either; run bin/p0_prepare.sh --step frames"
fi
EXPORT="ALL,MTB_BUILD_DIR=${BUILD},BUILD_ROOT=${BUILD_ROOT},OG=${OG},MTB_GRAPH_FRAMES=${FRAMES}"
EXPORT="${EXPORT},COHORT=${COHORT},CRAMMAP=${CRAMS},CRAMS=${CRAMS},REFMAP=${REFMAP}"
EXPORT="${EXPORT},${DIRS}"
mkdir -p slurm

deps() {  # join the job ids that exist, for --dependency=afterok:a:b:c
    local out="" j
    for j in "$@"; do
        [[ -n "$j" ]] && out="${out:+$out:}$j"
    done
    printf '%s' "$out"
}

submit() {  # name script extra_export array dep
    local name="$1" script="$2" extra="$3" array="$4" dep="$5"
    local -a args=(-J "rb_${COHORT_NAME}_${name}")
    # %N goes on the array spec, never on a non-array job: sbatch rejects
    # `--array=%50` outright, and the summary steps are submitted with an
    # empty array argument.
    if [[ -n "$array" ]]; then
        if [[ "$THROTTLE" -gt 0 ]]; then
            args+=(--array="${array}%${THROTTLE}")
        else
            args+=(--array="$array")
        fi
    fi
    [[ -n "$dep" ]] && args+=(--dependency="afterok:${dep}")
    args+=(--export="${EXPORT}${extra:+,$extra}")
    # stdout of this function is captured as the job id by the caller, so the
    # human-readable line must go to stderr or it becomes part of the id and
    # cascades into every downstream --dependency.
    if [[ "$DRY" -eq 1 ]]; then
        printf "  sbatch %s %s\n" "${args[*]}" "$script" >&2
        echo "DRYRUN_${name}"
        return
    fi
    # A rejected submission (QOS or submit limit, bad partition) must stop the
    # chain. Inside the caller's $(...) set -e does not apply, so an sbatch
    # failure used to return an empty id with status 0; deps() then dropped it
    # and the next pass went in with NO dependency, running on stale inputs.
    # Die here, in the subshell, AND print nothing, so the caller's check fails.
    local id
    if ! id="$(sbatch --parsable "${args[@]}" "$script")"; then
        echo "FATAL: sbatch rejected ${name} (${script}); stopping the chain" >&2
        [[ ${#JOB[@]} -gt 0 ]] && \
            echo "       already submitted, cancel if unwanted: scancel ${JOB[*]}" >&2
        return 1
    fi
    id="${id%%;*}"                       # --parsable may append ";cluster"
    if [[ ! "$id" =~ ^[0-9]+$ ]]; then
        echo "FATAL: sbatch returned no job id for ${name}: '${id}'" >&2
        return 1
    fi
    echo "$id"
}

# Check the I/O contract before submitting anything. A pass validates its own
# inputs, but inside the job -- after the array is submitted and after afterok
# has committed the rest of the chain to waiting on it. That cost 104 cancelled
# jobs on 2026-09-24 when the pinned-reference arm lacked P1's H37Rv-frame
# calls. IO_CHECK=0 skips it; a failure is advisory on --dry-run and fatal
# otherwise, because a dry run is how you find out what is missing.
if [[ "${IO_CHECK:-1}" != "0" && -s refbias/io_contract.tsv ]]; then
    # Tell the checker which passes are about to run, so an input that an
    # earlier pass in THIS chain produces is not reported as missing.
    _willrun="$(IFS=,; echo "${RUN[*]}")"
    if ! "$MTB_PY" refbias/bin/io_contract.py "$COHORT_NAME" \
            --will-run "$_willrun" > "${TMPDIR:-/tmp}/io_$$.txt" 2>&1; then
        sed 's/^/  /' "${TMPDIR:-/tmp}/io_$$.txt" >&2
        rm -f "${TMPDIR:-/tmp}/io_$$.txt"
        if [[ "$DRY" -eq 1 ]]; then
            echo "  (--dry-run: reporting only)" >&2
        else
            die "I/O contract check failed; fix the inputs or set IO_CHECK=0"
        fi
    fi
    rm -f "${TMPDIR:-/tmp}/io_$$.txt"
fi

echo
[[ "$DRY" -eq 1 ]] && echo "--dry-run: nothing will be submitted"
declare -A JOB=()
for p in "${RUN[@]}"; do
    case "$p" in
      p1)
        EX="OUTDIR=${OUTROOT}/p1,WORK=${WORKPFX}p1"
        JOB[p1a]=$(submit p1 bin/p1_select_reference.sh "$EX" "1-${N}" "")
        # P1's summary IS refmap.tsv, not a report, and every later pass
        # reads it. Omitting it is what broke the first run of this script.
        JOB[p1]=$(submit p1sum bin/p1_select_reference.sh \
                "${EX},P1STEP=--summary" "" "${JOB[p1a]}") ;;
      p2)
        EX="OUTDIR=${OUTROOT}/p2,WORK=${WORKPFX}p2"
        JOB[p2a]=$(submit p2 bin/p2_call.sh "$EX" "1-${N}" "${JOB[p1]:-}")
        # p2_summary.tsv is read by p5_sanity.py through p5_finish.sh
        JOB[p2]=$(submit p2sum bin/p2_call.sh \
                "${EX},P2STEP=--summary" "" "${JOB[p2a]}") ;;
      p3)
        EX="OUTDIR=${OUTROOT}/p3,WORK=${WORKPFX}p3"
        JOB[p3a]=$(submit p3 bin/p3_accessory.sh "$EX" "1-${N}" "${JOB[p2]:-}")
        JOB[p3]=$(submit p3sum bin/p3_accessory.sh \
                "${EX},P3STEP=--summary" "" "${JOB[p3a]}") ;;
      p4)
        EX="OUTDIR=${OUTROOT}/p4,WORK=${WORKPFX}p4"
        JOB[p4a]=$(submit p4 bin/p4_place.sh "$EX" "1-${N}" "${JOB[p2]:-}")
        JOB[p4]=$(submit p4sum bin/p4_place.sh \
                "${EX},P4STEP=--summary" "" "${JOB[p4a]}") ;;
      p4b)
        EX="OUTDIR=${OUTROOT}/p4b,WORK=${WORKPFX}p4b"
        JOB[p4ba]=$(submit p4b bin/p4b_place_sv.sh "$EX" "1-${N}" "${JOB[p2]:-}")
        JOB[p4b]=$(submit p4bsum bin/p4b_place_sv.sh \
                "${EX},P4BSTEP=--summary" "" "${JOB[p4ba]}") ;;
      p5)
        # three steps, and the order is load-bearing: keys once over all of
        # P4's output, then states per sample against those keys, then the
        # matrix. p5_finish.sh refuses to run if the generations disagree.
        P5EX="OUTDIR=${OUTROOT}/p5,WORK=${WORKPFX}p5,COHORT_NAME=${COHORT_NAME}"
        JOB[p5keys]=$(submit p5keys bin/p5_merge.sh "${P5EX},P5STEP=--keys" "" "${JOB[p4]:-}")
        JOB[p5states]=$(submit p5states bin/p5_merge.sh "${P5EX},P5STEP=--states" "1-${N}" "${JOB[p5keys]:-}")
        # --pre only: matrix, validation, sanity, SV matrix, P6. The merged VCF
        # is pass p5vcf, after the two genotyping passes below.
        # --pre also builds sv_matrix.tsv from P4b, and p4 and p4b run in
        # parallel off p2, so it must wait for BOTH or it reads a partial P4b.
        JOB[p5]=$(submit p5pre bin/p5_finish.sh "${P5EX},P5FSTEP=--pre" "" \
                  "$(deps "${JOB[p5states]:-}" "${JOB[p4b]:-}")")
        ;;
      p1g)
        # fixed-reference arm: every isolate against the same cut-down H37Rv.
        EX="OUTDIR=${OUTROOT}/p1g,WORK=${WORKPFX}p1g"
        JOB[p1ga]=$(submit p1g is6110/bin/p1g_isclean.sh "$EX" "1-${N}" "${JOB[p1]:-}")
        JOB[p1g]=$(submit p1gsum is6110/bin/p1g_isclean.sh \
                "${EX},P1GSTEP=--summary" "" "${JOB[p1ga]}") ;;
      p1i)
        # matched arm: each isolate against its own reference, cut down. The
        # cut-down builds are panel assets and already exist for all 333, so
        # nothing here builds them; if the panel changes, run
        # p1i_discover_matched.sh then p1i_build_matched.sh, both serial.
        # Array only: p1i_matched.sh has no --summary mode, because the
        # matched arm's per-cohort aggregation IS pass p1iv. Submitting a
        # summary step anyway is what made the first back-filled l7 run fail
        # with "FATAL: no task id" -- a non-array job has no
        # SLURM_ARRAY_TASK_ID, so the script fell through to its per-sample
        # branch and refused.
        EX="OUTDIR=${OUTROOT}/p1i,WORK=${WORKPFX}p1i"
        JOB[p1i]=$(submit p1i is6110/bin/p1i_matched.sh "$EX" "1-${N}" "${JOB[p1]:-}") ;;
      p1iv)
        # label, reconcile, then write the VCFs. One job, three steps in order,
        # because each reads the previous one's whole-cohort table and there is
        # nothing per-sample left to parallelise.
        EX="OUTDIR=${OUTROOT}/p1i,WORK=${WORKPFX}p1iv,P1IVCF=${OUTROOT}/p1i/vcf"
        JOB[p1iv]=$(submit p1iv bin/p1i_vcf.sh "$EX" "" "${JOB[p1i]:-}") ;;
      p1is)
        # Put the IS6110 arm into P5's key space and earn a REF for every
        # isolate that did not report a site. Both scripts existed and were
        # validated long before this pass; what was missing was anything that
        # ran them, so the merged VCF's IS6110 block was 98.7% NOCALL.
        EX="P1IDIR=${OUTROOT}/p1i,COHORT_NAME=${COHORT_NAME}"
        JOB[p1is]=$(submit p1is bin/p1i_p5states.sh "$EX" "" "${JOB[p1iv]:-}") ;;
      p5svgt)
        # Measure deletion absence, so the SV block is not presence-only.
        # Needs sv_matrix.tsv, which p5 --pre writes.
        EX="OUTDIR=${OUTROOT}/p5,WORK=${WORKPFX}p5svgt,SVDIR=${OUTROOT}/p5/svgt"
        JOB[svgtp]=$(submit p5svgtprobes bin/p5_svgt.sh \
                "${EX},SVGTSTEP=--probes" "" "${JOB[p5]:-}")
        JOB[svgta]=$(submit p5svgt bin/p5_svgt.sh \
                "${EX},SVGTSTEP=--states" "1-${N}" "${JOB[svgtp]}")
        JOB[p5svgt]=$(submit p5svgtmerge bin/p5_svgt.sh \
                "${EX},SVGTSTEP=--merge" "" "${JOB[svgta]}") ;;
      p5vcf)
        # The cohort-level deliverable, last, after everything that contributes
        # a genotype to it. Depending on only some of these is how a merged VCF
        # came to be written with no insertion sites in it at all.
        P5EX="OUTDIR=${OUTROOT}/p5,WORK=${WORKPFX}p5,COHORT_NAME=${COHORT_NAME}"
        _up="$(deps "${JOB[p5]:-}" "${JOB[p5svgt]:-}" "${JOB[p1is]:-}")"
        # Sharded by genome region above a few hundred isolates: one shard per
        # 500 isolates by default (VCF_SHARDS overrides), so a 10,000-isolate
        # cohort runs 20 shards each holding a twentieth of the records.
        SH="${VCF_SHARDS:-$(( (N + 499) / 500 ))}"
        if [[ "$SH" -gt 1 ]]; then
            P5EX="${P5EX},VCF_SHARDS=${SH}"
            JOB[p5vcfsh]=$(submit p5vcfshard bin/p5_finish.sh \
                    "${P5EX},P5FSTEP=--merge-shard" "1-${SH}" "$_up")
            JOB[p5vcf]=$(submit p5vcf bin/p5_finish.sh \
                    "${P5EX},P5FSTEP=--merge-assemble" "" "${JOB[p5vcfsh]}")
        else
            JOB[p5vcf]=$(submit p5vcf bin/p5_finish.sh "${P5EX},P5FSTEP=--merge" "" "$_up")
        fi ;;
      *) die "no rule for pass '${p}'" ;;
    esac
    # P5 submits three jobs rather than an array plus a summary, so the
    # generic line printed only p5finish twice and left keys and states
    # invisible -- which matters when one of them is the job that failed.
    if [[ "$p" == "p5" ]]; then
        echo "  -> p5 keys ${JOB[p5keys]:-?}, states ${JOB[p5states]:-?}, pre ${JOB[p5]:-?}"
    elif [[ "$p" == "p5svgt" ]]; then
        echo "  -> p5svgt probes ${JOB[svgtp]:-?}, states ${JOB[svgta]:-?}, merge ${JOB[p5svgt]:-?}"
    else
        echo "  -> ${p} array ${JOB[${p}a]:-${JOB[$p]:-?}}, summary ${JOB[$p]:-?}"
    fi
done

if [[ "$DRY" -eq 1 ]]; then
    echo
    echo "NOTE: each pass takes its step as \$1, and sbatch --export cannot"
    echo "      supply a positional argument, so it is passed as P1STEP,"
    echo "      P2STEP, ... P5STEP. Each script reads its own, with \$1"
    echo "      winning if both are given."
fi
