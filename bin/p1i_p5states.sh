#!/usr/bin/env bash
#SBATCH --job-name=P1i_p5states
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 0-08:00
#SBATCH -p shared
#SBATCH --mem=16000
#SBATCH --output=slurm/P1ip5_%j.out
#SBATCH --error=slurm/P1ip5_%j.err
#
# Pass p1is: put the IS6110 arm into P5's key space, then EARN a REF for every
# isolate that did not report a site.
#
#   sbatch bin/p1i_p5states.sh          # or bash, it runs either way
#
# WHY THIS IS ITS OWN PASS. Both scripts already existed, were validated in
# is6110/docs/P5_MERGE_PLAN.md, and were never wired into bin/refbias_run.sh --
# so stage 2 had been run by hand for the pilot and scale100 and for nobody
# else, and merge_cohort_vcf.py read the stage-1 key table regardless. The
# effect was measured on scale200: the IS6110 block of the merged VCF was 0.32%
# ALT, 0.95% REF and 98.73% NOCALL, and with almost no leaf in the reference
# state a parsimony reconstruction could say nothing -- of 828 insertion sites,
# 425 had exactly one derived leaf and only 30 had two or more, so the
# convergence question the arm exists to ask could not be asked at all.
#
# Stage 1 (is6110_p5_merge.py) translates the arm's keys into P5's key space
# and leaves every non-reporting isolate NOCALL, asserting nothing.
# Stage 2 (is6110_p5_stage2.py) projects a carrier's coordinate onto each
# non-carrier's own path, converts it into that isolate's element-free frame,
# and reads depth there: depth and no junction candidate earns REF, a position
# that does not project is ABSENT, anything else stays NOCALL. The depth floor
# matches p5_states.py's, so an IS6110 REF means the same strength of evidence
# as every other REF in the matrix.
set -euo pipefail

_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: no project_env.sh" >&2; exit 1; }
source "$_mtb_env"
cd "${SLURM_SUBMIT_DIR:-.}"

REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
P1IDIR="${P1IDIR:-refbias/p1i}"
COHORT_NAME="${COHORT_NAME:-}"
RES="${RES:-is6110/results}"
# The pilot's tables are unprefixed and scale100's carry the `scale` tag; every
# other cohort uses its own name. merge_cohort_vcf.py maps the same three cases.
case "$COHORT_NAME" in
  ""|pilot|pilot_rerun) TAG="" ;;
  scale100)             TAG="scale_" ;;
  *)                    TAG="${COHORT_NAME}_" ;;
esac
KEYTAB="${KEYTAB:-${RES}/${TAG}p1i_cohort_keys.tsv}"
[[ -s "$KEYTAB" ]] || { echo "FATAL: no IS6110 cohort key table at ${KEYTAB} -- run pass p1iv first" >&2; exit 1; }

OUT1K="${RES}/${TAG}p5_is6110_keys.tsv"
OUT1S="${RES}/${TAG}p5_is6110_states.tsv"
OUT2="${RES}/${TAG}p5_is6110_states.stage2.tsv"

# STEPS. Stage 2 used to be one serial job: an odgi load per distinct reference
# and a samtools pass per sample, all in one process -- 29 minutes at 200
# isolates, 99 at 997, and on course for a day at 10,000. It is now three
# steps the runner submits as two arrays and a join; --all still runs every
# step in this one job.
#
#   --stage1    stage 1, cohort-wide and cheap; also one table per sample
#   --project   stage 2, one reference per array task (task i = i-th distinct
#               reference of the refmap; tasks past the last exit 0)
#   --sample    stage 2, one sample per array task
#   --merge     stage 2, the cohort table
#   --all       all of the above, in order (default)
P1ISSTEP="${1:-${P1ISSTEP:---all}}"
case "$P1ISSTEP" in
  --all|--stage1|--project|--sample|--merge) ;;
  *) echo "usage: $0 [--all|--stage1|--project|--sample|--merge]" >&2; exit 2 ;;
esac
STAGE1DIR="${P1IDIR}/p5stage1"
STAGE2=("$MTB_PY" is6110/bin/is6110_p5_stage2.py
    --cohort-keys "$KEYTAB" --stage1-dir "$STAGE1DIR" --refmap "$REFMAP"
    --isclean-dir "$P1IDIR" --workdir "${P1IDIR}/p5stage2" --out "$OUT2"
    --threads "${SLURM_CPUS_PER_TASK:-1}")

if [[ "$P1ISSTEP" == "--all" || "$P1ISSTEP" == "--stage1" ]]; then
    echo "=== stage 1: into P5's key space ==="
    "$MTB_PY" is6110/bin/is6110_p5_merge.py \
        --cohort-keys "$KEYTAB" --refmap "$REFMAP" \
        --out-keys "$OUT1K" --out-states "$OUT1S" --out-states-dir "$STAGE1DIR"
    [[ "$P1ISSTEP" == "--stage1" ]] && exit 0
fi

# The graph and its path list come from the build, never from the script's
# defaults, which are a relative glob and a hard-coded build id.
BUILD="${MTB_BUILD_DIR:?MTB_BUILD_DIR is unset; run through bin/refbias_run.sh}"
OG="${OG:-$(awk -F'\t' '$1=="graph"{print $2}' "${BUILD}/build_info.tsv")}"
[[ -s "$OG" ]] || { echo "FATAL: no graph at '${OG}'" >&2; exit 1; }
[[ -s "${BUILD}/assets/paths.txt" ]] || { echo "FATAL: no ${BUILD}/assets/paths.txt" >&2; exit 1; }
STAGE2+=(--graph "$OG" --paths "${BUILD}/assets/paths.txt"
         --store "${BUILD}/proj_is6110")

case "$P1ISSTEP" in
  --project)
    "${STAGE2[@]}" --mode project --index "${SLURM_ARRAY_TASK_ID:?--project is an array step}" ;;
  --sample)
    "${STAGE2[@]}" --mode sample --index "${SLURM_ARRAY_TASK_ID:?--sample is an array step}" ;;
  --merge)
    "${STAGE2[@]}" --mode merge
    echo "written: ${OUT2}" ;;
  --all)
    echo
    echo "=== stage 2: earn REF for the non-carriers ==="
    "${STAGE2[@]}" --mode all
    echo "written: ${OUT2}" ;;
esac
