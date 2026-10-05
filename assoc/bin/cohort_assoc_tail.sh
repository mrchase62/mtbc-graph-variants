#!/usr/bin/env bash
# From a finished cohort chain to the association tests, in one command.
#
#   cd /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/runroot
#   bash /n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/assoc/bin/cohort_assoc_tail.sh \
#       gwas1000 [phenotype.txt]
#
# RUN FROM RUNROOT, BY THE REPOSITORY PATH. The working tree's root also holds
# an assoc/bin/ and a bin/, with stale copies of this script, the scan and the
# audit; a relative `bash assoc/bin/...` from there runs those. runroot links
# assoc/bin and bin to the repository and its data directories to the working
# tree, so the code is the repository's and the data are production's. This
# script refuses to run anywhere whose assoc/bin and bin are not its own.
#
# The five steps between p5vcf and an answer, which were run by hand for
# scale200 and are collected here so a second cohort does not have to
# rediscover them. Each step skips if its product already exists, so the
# script is safe to re-run after a failure.
#
#   1  alignment from the merged VCF, with the reference as a tip
#   2  the outgroup read off the panel VCF -- without it the cohort root is
#      inferred from the cohort alone and deep branches go undetermined
#   3  the tree, rooted on that outgroup
#   4  the event matrix: per-branch gains, the input every test reads
#   5  the tests -- variant-level scan, then IS6110 collapsed to genes
set -euo pipefail
_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
source "$_mtb_env"

# THE SIBLINGS MUST BE THIS SCRIPT'S OWN. Every assoc/bin and bin script below
# is called through these two absolute paths, and the run directory's relative
# assoc/bin and bin must resolve to the same places -- true in runroot, false
# at the working-tree root, whose copies are stale.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO="$(cd "${HERE}/../.." && pwd -P)"
_here_rel="$(cd assoc/bin 2>/dev/null && pwd -P || true)"
_bin_rel="$(cd bin 2>/dev/null && pwd -P || true)"
if [[ "$_here_rel" != "$HERE" || "$_bin_rel" != "${REPO}/bin" ]]; then
    echo "FATAL: run this from runroot. Here ($(pwd)) assoc/bin is" \
         "'${_here_rel:-absent}' and bin is '${_bin_rel:-absent}', but this" \
         "script's siblings are ${HERE} and ${REPO}/bin" >&2
    exit 1
fi
for _s in "${HERE}/add_outgroup.py" "${HERE}/write_event_matrix.py" \
          "${HERE}/assoc_scan.py" "${HERE}/is6110_gene_burden.py" \
          "${REPO}/bin/vcf_to_alignment.py" "${REPO}/bin/build_snp_tree.sh" \
          "${REPO}/bin/retier_intervals.py" "${REPO}/bin/audit_chain.py"; do
    [[ -s "$_s" ]] || { echo "FATAL: missing ${_s}" >&2; exit 1; }
done

# The cohort's own table, as the registry records it, for the lineage null.
lookup_cohort_table() {
    awk -F'\t' -v c="$C" '!/^#/ && $1!="cohort" && $1==c {print $2}' \
        refbias/cohorts.tsv 2>/dev/null | head -1
}

C="${1:?usage: $0 <cohort> [phenotype.txt]}"
# LEVEL 1 presence tables for this cohort, if they exist. They make the
# accessory arm's within-insert variation conditional on carrying the insert.
ACCPRES="${ACCPRES:-accessory/${C}}"
compgen -G "${ACCPRES}/*.presence.tsv" >/dev/null || ACCPRES=""
PHENO="${2:-}"
OUT="refbias/${C}"
[[ "$C" == "pilot" ]] && OUT="refbias"
VCF="${OUT}/p5/merged.vcf.gz"
[[ -s "$VCF" ]] || { echo "FATAL: no ${VCF}; the cohort chain is not finished" >&2; exit 1; }
mkdir -p "assoc/${C}"

echo "=== 1. alignment"
if [[ ! -s "data/trees/${C}.snps.fasta" ]]; then
    "$MTB_PY" "${REPO}/bin/vcf_to_alignment.py" --vcf "$VCF" \
        --out "data/trees/${C}.snps.fasta" \
        --sites-out "data/trees/${C}.sites.tsv" \
        --ref-sample H37Rv --max-missing 0.10
else echo "  already built"; fi

echo "=== 2. outgroup"
if [[ ! -s "data/trees/${C}.og.fasta" ]]; then
    "$MTB_PY" "${HERE}/add_outgroup.py" \
        --alignment "data/trees/${C}.snps.fasta" \
        --sites "data/trees/${C}.sites.tsv" \
        --out "data/trees/${C}.og.fasta"
else echo "  already built"; fi

echo "=== 3. tree"
if [[ ! -s "data/trees/${C}.rooted.nwk" ]]; then
    J=$(sbatch --parsable "${REPO}/bin/build_snp_tree.sh" "data/trees/${C}.og.fasta" \
            GCF_035581225 "data/trees/${C}")
    echo "  submitted ${J}; re-run this script when it finishes"
    exit 0
else echo "  already built"; fi

echo "=== 4. event matrix"
# GUARD ON THE LAST FILE WRITTEN, NOT THE FIRST. labelled.nwk is written
# before the reconstruction starts and summary.txt after everything else, so a
# run that dies in between leaves a directory this step used to treat as
# finished. gwas1000 did exactly that -- the writer refused a duplicate key,
# left labelled.nwk behind, and the next attempt skipped step 4 and failed in
# step 5 on a missing variants.tsv.
if [[ ! -s "assoc/${C}/events/summary.txt" ]]; then
    "$MTB_PY_VT" "${HERE}/write_event_matrix.py" --vcf "$VCF" \
        --tree "data/trees/${C}.rooted.nwk" --out "assoc/${C}/events" \
        --ref-sample H37Rv --outgroup-name GCF_035581225 \
        --outgroup-fasta "data/trees/${C}.og.fasta" \
        --outgroup-sites "data/trees/${C}.sites.tsv" \
        ${ACCPRES:+--accessory-presence "$ACCPRES"} \
        --dedupe suffix
else echo "  already built"; fi

[[ -n "$PHENO" ]] || { echo "no phenotype given; stopping before the tests"; exit 0; }
# LEVEL 2 is conditional on level 1, so the presence tables are passed to the
# event matrix as well as to the VCF. Without them a variant inside an insert is
# reconstructed over every sample, including the ones that do not carry the
# sequence it lives in.
echo "=== 5a. variant-level scan"
# THE LINEAGE NULL NEEDS A LINEAGE TABLE, and without --lineages the scan
# silently writes an empty p_lineage column rather than failing. gwas1000's
# first scan came back with all 14,033 rows blank in that column, which reads
# as "nothing survived the lineage null" when in fact it was never run. The
# cohort phenotype table carries the lineage call, so it is passed here.
LINTAB="${LINTAB:-refbias/${C}.phenotype.tsv}"
[[ -s "$LINTAB" ]] || LINTAB="$(lookup_cohort_table)"
# No table is fatal, not a blank column: the survivor rule needs the lineage
# null, and a run without it reports nothing surviving rather than an error.
if [[ -z "$LINTAB" || ! -s "$LINTAB" ]]; then
    echo "FATAL: no lineage table for ${C} (tried refbias/${C}.phenotype.tsv" \
         "and refbias/cohorts.tsv); set LINTAB" >&2
    exit 1
fi
echo "  lineage table: ${LINTAB}"
# THE SV EVIDENCE TIERS, REBUILT FROM THIS COHORT'S OWN CATALOGUE EVERY RUN.
# The scan keys each catalogued deletion by its evidence tier. The chain
# stopped passing the table, so from 2026-10-02 every svi: deletion was nulled
# in one `other` pool -- and the only tiered table on disk
# (refbias/assets/sv_intervals.retiered.tsv) was built from the previous
# catalogue and held 292 of scale200's 467 current IDs. Rebuilding it here
# takes seconds and cannot go stale; the scan refuses any svi: ID it cannot
# tier. The evidence is this cohort's own interval genotypes.
SVTIERS="assoc/${C}/sv_intervals.retiered.tsv"
for _f in "${OUT}/p5/sv_intervals.tsv" "${OUT}/p5/svgt_iv_states.tsv" \
          "${OUT}/p1/refmap.tsv"; do
    [[ -s "$_f" ]] || { echo "FATAL: no ${_f} for the SV evidence tiers" >&2; exit 1; }
done
"$MTB_PY" "${REPO}/bin/retier_intervals.py" \
    --intervals "${OUT}/p5/sv_intervals.tsv" \
    --iv-states "${OUT}/p5/svgt_iv_states.tsv" \
    --refmap "${OUT}/p1/refmap.tsv" --cohort-tag "$C" \
    --out "$SVTIERS"
# THE SCAN NEEDS THE PRESENCE TABLES TOO. Without them it measures the
# callability floor against the whole tree and runs no level-2 null, so every
# variant inside an accessory locus is judged on branches it can never occupy.
# gwas1000's 2026-10-02 scan ran that way: cond_carriers and p_cond were blank
# for every accessory variant, and the audit flagged 6 records with two or more
# origins that the floor had dropped.
"$MTB_PY_VT" "${HERE}/assoc_scan.py" --events "assoc/${C}/events" \
    --phenotype "$PHENO" --lineages "$LINTAB" \
    ${ACCPRES:+--accessory-presence "$ACCPRES"} \
    --sv-intervals "$SVTIERS" \
    --out "assoc/${C}/scan.tsv"
echo "=== 5b. IS6110, collapsed to genes"
# The burden takes the same lineage table as the scan. Without it the lineage
# null is skipped and p_lineage is written blank, which gwas1000 and scale200
# both did through 2026-10-02.
"$MTB_PY_VT" "${HERE}/is6110_gene_burden.py" --events "assoc/${C}/events" \
    --phenotype "$PHENO" --lineages "$LINTAB" \
    --out "assoc/${C}/is6110_gene.tsv"
# THE SAME BURDEN FOR SMALL VARIANTS AND SVs. These were run by hand from
# 2026-09-25 and never added here, so the 2026-10-02 rerun archived them and
# produced neither. The small-variant burden is where loss-of-function
# resistance shows up (pncA, ethA, gid, the fabG1 and embA promoters): many
# distinct rare alleles, none testable alone, many origins per gene.
for _cls in small sv; do
    echo "=== 5c. ${_cls} variants, collapsed to genes"
    "$MTB_PY_VT" "${HERE}/is6110_gene_burden.py" --events "assoc/${C}/events" \
        --phenotype "$PHENO" --lineages "$LINTAB" --cls "$_cls" \
        --out "assoc/${C}/${_cls}_gene.tsv"
done

echo
echo "=== 6. chain audit: is every class and key space still accounted for?"
# THE STANDING CHECK. Three times in one week a downstream stage silently
# ignored something upstream had added -- the svi: key space, the
# accessory_presence class, and the callability denominator on conditional
# variants -- and none raised an error. This reconciles every (class, region)
# stratum and every ID key space across both boundaries, and fails when a loss
# is not fully explained by the scan's own two filters. It runs last, after
# every output is written, so a failure reports a problem rather than
# destroying the run's results.
AUDIT_RC=0
# The audit is the repo's copy, found next to this script; it reads the event
# matrix's tree, through assoc/bin's tree code, to count applicable branches.
"$MTB_PY_VT" "${REPO}/bin/audit_chain.py" --cohort "$C" \
    --assoc-bin "$HERE" \
    ${ACCPRES:+--accessory-presence "$ACCPRES"} \
    --out "assoc/${C}/chain_audit.txt" || AUDIT_RC=$?
if [[ "$AUDIT_RC" -ne 0 ]]; then
    echo
    echo "############################################################"
    echo "# CHAIN AUDIT FAILED. A class, region or key space present  #"
    echo "# upstream has no representative downstream, and the loss   #"
    echo "# is not explained by the scan's documented filters. The    #"
    echo "# outputs above were still written -- read them knowing an  #"
    echo "# arm may be missing from them.                             #"
    echo "#   assoc/${C}/chain_audit.txt                              #"
    echo "############################################################"
fi
exit "$AUDIT_RC"
