#!/usr/bin/env bash
#SBATCH --job-name=archive_aln
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -t 0-01:00
#SBATCH -p shared
#SBATCH --mem=4000
#SBATCH --output=slurm/archive_%A_%a.out
#SBATCH --error=slurm/archive_%A_%a.err
#
# Keep the alignments later passes read, as lossless CRAM on durable storage,
# and optionally drop the ones nothing reads.
#
#   sbatch --array=1-N bin/archive_alignments.sh --archive   # one sample per task
#   sbatch --array=1-N bin/archive_alignments.sh --restore   # CRAM back to BAM
#   bash bin/archive_alignments.sh --status                  # what exists where
#   bash bin/archive_alignments.sh --archive <sample>        # one sample, by hand
#
# WHY. sync_back.sh excludes alignments, and the ones later passes read are
# lost when scratch is purged: the matched-reference BAM (P3, SV genotyping)
# and the element-free matched BAM (IS6110 stage 2). Recomputing them means
# realigning -- about 0.7 CPU-hours per sample for P2 alone -- while as CRAM
# against the build's own references they are 37% of the BAM size, about 90 MB
# per sample for both. The graph is stable, so they stay valid until it changes.
#
# WHAT IS KEPT, per sample                  reference it is encoded against
#   P2WORK/<s>.bam          -> .archive.cram   ${BUILD}/refs/<R>.fasta
#   P1IDIR/<s>.isclean.bam  -> .archive.cram   ${CLEANDIR}/<R>.isclean.fasta
#
# WHAT NOTHING READS, and is deleted only with ARCHIVE_DROP_UNUSED=1
#   P1WORK/<s>.h37rv.bam    P1's own alignment; P1 keeps its VCF, nothing
#                           after P1 opens the BAM
#   P1GDIR/<s>.isclean.bam  the fixed-reference IS6110 arm; read only by p1g
#                           itself, and the merged VCF does not use that arm
#
# VERIFIED BEFORE ANYTHING IS REMOVED. Each CRAM is decoded and compared with
# its BAM: the record count, every record with its auxiliary tags sorted (CRAM
# regenerates NM and MD on decoding, so they move within a record but not in
# value), and the header apart from the M5 and UR fields CRAM adds to each @SQ.
# A mismatch fails the task and leaves the BAM alone. With ARCHIVE_DELETE_BAM=1
# a verified BAM and its index are then removed; by default they are kept.
#
# The CRAMs are named *.archive.cram so sync_back.sh can include them while it
# still excludes other alignments, including the read collection's own CRAMs.
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

BUILD="${MTB_BUILD_DIR:?MTB_BUILD_DIR is unset; run through bin/refbias_run.sh}"
REFMAP="${REFMAP:?REFMAP is unset}"
P1WORK="${P1WORK:-refbias/work/p1}"
P2WORK="${P2WORK:-refbias/work/p2}"
P1IDIR="${P1IDIR:-refbias/p1i}"
P1GDIR="${P1GDIR:-refbias/p1g}"
CLEANDIR="${CLEANDIR:-is6110/assets/isclean_matched}"
ST="${MTB_SAMTOOLS:?MTB_SAMTOOLS is unset}"
NT="${SLURM_CPUS_PER_TASK:-2}"

STEP="${1:-${ARCHIVESTEP:-}}"
case "$STEP" in
  --archive|--restore|--status) ;;
  *) echo "usage: $0 --archive|--restore|--status [sample]" >&2; exit 2 ;;
esac

# (bam, reference) pairs kept for one sample
kept() {
    local s="$1" r
    r="$(awk -F'\t' -v s="$s" '$1==s{print $5; exit}' "$REFMAP")"
    [[ -n "$r" ]] || { echo "FATAL: ${s} not in ${REFMAP}" >&2; return 1; }
    printf '%s\t%s\n' "${P2WORK}/${s}.bam" "${BUILD}/refs/${r}.fasta"
    printf '%s\t%s\n' "${P1IDIR}/${s}.isclean.bam" "${CLEANDIR}/${r}.isclean.fasta"
}

cram_of() { printf '%s' "${1%.bam}.archive.cram"; }

# records with auxiliary tags sorted, so tag order does not count as a change
records_md5() {
    awk -F'\t' '{n = 0; delete t; for (i = 12; i <= NF; i++) t[++n] = $i
                 asort(t); s = $1; for (i = 2; i <= 11; i++) s = s "\t" $i
                 for (i = 1; i <= n; i++) s = s "\t" t[i]; print s}' \
        | md5sum | cut -d' ' -f1
}
header_md5() {
    grep -v '^@PG' | sed -E 's/\t(M5|UR):[^\t]*//g' | md5sum | cut -d' ' -f1
}

archive_one() {
    local bam="$1" fa="$2" cram tmp
    cram="$(cram_of "$bam")"
    if [[ ! -s "$bam" ]]; then
        if [[ -s "$cram" ]]; then echo "  ${cram}: archived, BAM already removed"
        else echo "  WARNING: neither ${bam} nor its archive exists"; fi
        return 0
    fi
    [[ -s "$fa" ]] || { echo "FATAL: no reference ${fa} for ${bam}" >&2; return 1; }
    tmp="${cram}.tmp.$$"
    if [[ -s "$cram" && "$cram" -nt "$bam" ]]; then
        echo "  ${cram}: already archived and newer than the BAM"
    else
        "$ST" view -@ "$NT" -C -T "$fa" -o "$tmp" "$bam"
        mv -f "$tmp" "$cram"
        "$ST" index "$cram"
    fi
    # verify, every time: an archive is only as good as its last check
    local c1 c2 r1 r2 h1 h2
    c1="$("$ST" view -c "$bam")"
    c2="$("$ST" view -c -T "$fa" "$cram")"
    r1="$("$ST" view "$bam" | records_md5)"
    r2="$("$ST" view -T "$fa" "$cram" | records_md5)"
    h1="$("$ST" view -H "$bam" | header_md5)"
    h2="$("$ST" view -H -T "$fa" "$cram" | header_md5)"
    if [[ "$c1" != "$c2" || "$r1" != "$r2" || "$h1" != "$h2" ]]; then
        echo "FATAL: ${cram} does not reproduce ${bam}" \
             "(records ${c1}/${c2}, record md5 ${r1:0:12}/${r2:0:12}," \
             "header md5 ${h1:0:12}/${h2:0:12}); BAM kept" >&2
        return 1
    fi
    echo "  ${cram}: verified, ${c1} records," \
         "$(( $(stat -L -c%s "$cram") * 100 / $(stat -L -c%s "$bam") ))% of the BAM"
    if [[ "${ARCHIVE_DELETE_BAM:-0}" == 1 ]]; then
        rm -f "$bam" "${bam}.bai"
        echo "  ${bam}: removed (ARCHIVE_DELETE_BAM=1)"
    fi
}

restore_one() {
    local bam="$1" fa="$2" cram tmp
    cram="$(cram_of "$bam")"
    if [[ -s "$bam" ]]; then echo "  ${bam}: present"; return 0; fi
    [[ -s "$cram" ]] || { echo "FATAL: no BAM and no archive for ${bam}" >&2; return 1; }
    tmp="${bam}.tmp.$$"
    "$ST" view -@ "$NT" -b -T "$fa" -o "$tmp" "$cram"
    mv -f "$tmp" "$bam"
    "$ST" index "$bam"
    echo "  ${bam}: restored from ${cram}"
}

samples() {
    if [[ -n "${2:-}" ]]; then echo "$2"; return; fi
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n{print $1}' "$REFMAP"
    else
        awk -F'\t' 'NR>1{print $1}' "$REFMAP"
    fi
}

if [[ "$STEP" == "--status" ]]; then
    printf '%-14s %8s %8s %8s\n' kind bam archive neither
    for kind in matched isclean; do
        nb=0; na=0; nn=0
        while read -r s; do
            while IFS=$'\t' read -r bam fa; do
                case "$kind" in matched) [[ "$bam" == *isclean* ]] && continue ;;
                                isclean) [[ "$bam" == *isclean* ]] || continue ;; esac
                if [[ -s "$bam" ]]; then nb=$((nb + 1))
                elif [[ -s "$(cram_of "$bam")" ]]; then na=$((na + 1))
                else nn=$((nn + 1)); fi
            done < <(kept "$s")
        done < <(awk -F'\t' 'NR>1{print $1}' "$REFMAP")
        printf '%-14s %8d %8d %8d\n' "$kind" "$nb" "$na" "$nn"
    done
    exit 0
fi

rc=0
while read -r s; do
    [[ -n "$s" ]] || { echo "FATAL: no refmap row for this task" >&2; exit 1; }
    echo "[archive] ${s} ${STEP}"
    while IFS=$'\t' read -r bam fa; do
        if [[ "$STEP" == "--archive" ]]; then archive_one "$bam" "$fa" || rc=1
        else restore_one "$bam" "$fa" || rc=1; fi
    done < <(kept "$s")
    if [[ "$STEP" == "--archive" && "${ARCHIVE_DROP_UNUSED:-0}" == 1 && "$rc" -eq 0 ]]; then
        for f in "${P1WORK}/${s}.h37rv.bam" "${P1GDIR}/${s}.isclean.bam"; do
            if [[ -e "$f" ]]; then
                rm -f "$f" "${f}.bai"
                echo "  ${f}: removed (ARCHIVE_DROP_UNUSED=1; nothing after its own pass reads it)"
            fi
        done
    fi
done < <(samples "$STEP" "${2:-}")
exit "$rc"
