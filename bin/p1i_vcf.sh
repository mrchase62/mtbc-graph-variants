#!/usr/bin/env bash
#SBATCH --job-name=P1i_vcf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=16000
#SBATCH --output=slurm/P1ivcf_%A.out
#SBATCH --error=slurm/P1ivcf_%A.err
#
# The IS6110 arm's cohort-level tail: label the stacks, reconcile the two
# sides, write the VCFs. Submitted by bin/refbias_run.sh as pass p1iv.
#
#   OUTDIR=refbias/scale/p1i bash bin/p1i_vcf.sh
#
# FIVE STEPS, NOT THREE. The first version of this script ran promote,
# reconcile and write_vcf, and the writer died on a missing
# <tag>_p1i_sites_flank_all.tsv. The cohort tail is longer than it looks:
# promote labels the stacks, reconcile joins the element side to the junction
# side, project_sites carries each call into the H37Rv frame through the graph,
# place_by_flank decides which calls survive independent flank placement, and
# only then can the VCF be written. Missing a step produced a FileNotFoundError
# rather than a wrong number, which is the good case; it is recorded here so the
# order is not rediscovered.
#
# ONE JOB, FIVE STEPS IN ORDER. Each reads the previous step's whole-cohort
# table, so there is nothing per-sample left to parallelise and splitting them
# into three Slurm jobs would only add three chances for a dependency to be
# wired wrong. That is not hypothetical: this pipeline has already lost a run
# to a pass whose --summary step was omitted, and another to an afterok chain
# reading a stale generation.
set -euo pipefail

for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
source "$_mtb_env"
cd "${SLURM_SUBMIT_DIR:-.}"

P1IDIR="${P1IDIR:-${OUTDIR:-refbias/p1i}}"
REFMAP="${REFMAP:-refbias/p1/refmap.tsv}"
VCFDIR="${P1IVCF:-${P1IDIR}/vcf}"
TAG="${TAG:-$(basename "$(dirname "$P1IDIR")")}"
[[ "$TAG" == "refbias" ]] && TAG="pilot"
RES="${RES:-is6110/results}"
OG="${OG:-$(awk -F'\t' '$1=="graph"{print $2}' "${MTB_BUILD_DIR:-refbias/build/7713a8d71d8e}/build_info.tsv")}"
mkdir -p "$VCFDIR" "$RES"

# The PILOT's tables are unprefixed (p1i_cohort_keys.tsv and so on), and that
# is what p1i_p5states.sh and merge_cohort_vcf.py read for it. Writing
# pilot_p1i_cohort_keys.tsv here meant a pilot rerun left the old table in
# place and stage 2 silently used it. Every other cohort is <tag>_.
PFX="${TAG}_"; [[ "$TAG" == "pilot" ]] && PFX=""
SITES="${RES}/sites_${TAG}_p1i.tsv"
CALLS="${RES}/calls_${TAG}_p1i.tsv"
RECON="${RES}/${PFX}p1i_reconcile.tsv"
FLANK="${RES}/${PFX}p1i_sites_flank_all.tsv"
H37RV="${RES}/${PFX}p1i_sites_h37rv_all.tsv"
KEYS="${RES}/${PFX}p1i_cohort_keys.tsv"

echo "=== p1iv: cohort ${TAG}, stacks in ${P1IDIR}"
n=$(ls "${P1IDIR}"/*.elstacks.tsv 2>/dev/null | wc -l)
[[ "$n" -gt 0 ]] || { echo "FATAL: no elstacks in ${P1IDIR}; run pass p1i first" >&2; exit 1; }
echo "    ${n} per-sample stack files"

echo "=== 0/5 DR-array rescue"
# Copies inside the CRISPR DR array have junction reads that align equally well
# to many identical repeats, so no stack reaches the mapping-quality threshold
# and the copy is never called. This merges such stacks into one locus-level
# site; samples whose DR copy was already called are left exactly as they were.
# See is6110/bin/is6110_repeat_rescue.py.
"$MTB_PY" is6110/bin/is6110_repeat_rescue.py --dir "$P1IDIR" --refmap "$REFMAP" \
    --clean-dir "${CLEANDIR:-is6110/assets/isclean_matched}"

echo "=== 1/5 labelling"
"$MTB_PY" is6110/bin/is6110_promote_sites.py \
    --dir "$P1IDIR" --refmap "$REFMAP" --out "$SITES" --calls-out "$CALLS"

echo "=== 2/5 reconcile"
# --check-against is deliberately empty. Its default points at the PILOT's
# isclean_summary.tsv, which has 23 rows, so any other cohort is refused -- a
# real refusal that cost a run before it was understood.
"$MTB_PY" is6110/bin/is6110_reconcile.py \
    --elside-dir "$P1IDIR" --junc-dir "$P1IDIR" --refmap "$REFMAP" \
    --check-against "" --out "$RECON" \
    --summary-out "${RECON%.tsv}_summary.tsv"

echo "=== 3/5 project into the H37Rv frame"
# --all-stacks: project EVERY reconciled stack, not only the A/B tiers. The
# output is named _all and the writer needs a flank row for every site; without
# the flag one-sided and unconfirmed two-sided-wide sites got no row, were keyed
# "node:" and were dropped by the P5 merge (1,800 gwas1000 rows). The pilot and
# scale100 were run by hand WITH the flag, so this matches what they did.
"$MTB_PY" is6110/bin/is6110_project_sites.py --all-stacks \
    --reconcile "$RECON" --refmap "$REFMAP" --graph "$OG" \
    --workdir "${WORK:-refbias/work/p1iv}" --out "$H37RV"

echo "=== 4/5 flank placement"
"$MTB_PY" is6110/bin/is6110_place_by_flank.py \
    --sites "$H37RV" --minimap2 "$MTB_MINIMAP2" \
    --workdir "${WORK:-refbias/work/p1iv}" --out "$FLANK"

echo "=== 5/5 VCFs"
"$MTB_PY" is6110/bin/is6110_write_vcf.py \
    --reconcile "$RECON" --flank "$FLANK" --refmap "$REFMAP" \
    --outdir "$VCFDIR" --keys-out "$KEYS"

echo "=== p1iv done"
printf '    %-38s %s\n' "$SITES" "$(wc -l < "$SITES") lines"
printf '    %-38s %s\n' "$CALLS" "$(wc -l < "$CALLS") lines"
printf '    %-38s %s\n' "$RECON" "$(wc -l < "$RECON") lines"
printf '    %-38s %s\n' "$KEYS" "$(wc -l < "$KEYS") lines"
echo "    VCFs: $(ls "${VCFDIR}"/*.is6110.vcf 2>/dev/null | wc -l) matched-frame, $(ls "${VCFDIR}"/*.h37rv.vcf 2>/dev/null | wc -l) H37Rv-frame"
