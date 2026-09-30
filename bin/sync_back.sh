#!/usr/bin/env bash
# Copy results from scratch back to durable storage.
#
# RUN THIS AFTER EVERY SUBSTANTIVE RUN. netscratch is purged without warning;
# the 2025 incarnation of this project was lost exactly this way.
#
#   bin/sync_back.sh             # graphs + logs + scripts, excluding bulky intermediates
#   bin/sync_back.sh --full      # everything, including seqwish/smoothxg intermediates
#   bin/sync_back.sh --dry-run   # show what would move
#
# Also commits changes to TRACKED files and pushes git history to a bare mirror
# at ${MTB_PERSIST}.git (override with MTB_GIT_MIRROR). rsync propagates the
# current state; only the git history can recover an overwritten script.
#
# Untracked files are rsynced but NOT committed -- they are listed instead, so
# that `git add` stays a human decision. This matters because several sessions
# may share one working tree, and an `add -A` here commits someone else's
# in-progress work under this script's identity.
set -euo pipefail
# --- locate config/project_env.sh -----------------------------------------
# Under sbatch, BASH_SOURCE[0] is Slurm's spool copy of this script rather than
# the file in bin/, so the relative lookup alone is not enough. Fall back to the
# submit directory and then to an already-exported MTB_WORK.
_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/../config/project_env.sh" \
          "${MTB_WORK:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || {
    echo "FATAL: cannot locate config/project_env.sh." >&2
    echo "       Submit from the project root, or export MTB_ENV_FILE." >&2
    exit 1
}
source "$_mtb_env"

FULL=0; DRY=()
for a in "$@"; do
    case "$a" in
        --full)    FULL=1 ;;
        --dry-run) DRY=(--dry-run) ;;
        *) echo "unknown option: $a" >&2; exit 2 ;;
    esac
done

mtb_require_work
mkdir -p "$MTB_PERSIST"/{results,logs,bin,config}

# rsync exit 24 means "some files vanished before they could be transferred",
# which is normal here: this tree is synced while Slurm jobs are still writing
# and deleting scratch files, and a caller's temp directory disappearing
# mid-copy is not a failure of the backup. Under `set -e` that nonzero status
# aborted the whole script, and because the run was invoked through a pipe the
# pipeline's status came from the last command, so a sync that never reached
# its git push reported success. That happened on 2026-09-23 and left the
# mirror four commits behind while the run looked clean. rs() tolerates 24 and
# nothing else.
RSYNC=(rsync -ah --info=progress2 --no-inc-recursive "${DRY[@]}")
# SYNC_BACKUP=1 keeps the durable copy of any file an rsync overwrites, under
# .sync_backup/<timestamp>/ in the destination. Without it a file that
# REGRESSED on scratch (truncated, or rewritten by a bad run) silently replaced
# the good durable copy. Off by default only because it costs space.
if [[ "${SYNC_BACKUP:-0}" == 1 ]]; then
    RSYNC+=(--backup --backup-dir=".sync_backup/$(date -u +%Y%m%dT%H%M%SZ)")
fi
# Any failure below sets this, and the script then exits nonzero and does not
# append to .last_sync -- a failed push or copy used to be reported as success.
SYNC_FAILED=0
rs() {
    local st=0
    "${RSYNC[@]}" "$@" || st=$?
    if [[ $st -ne 0 && $st -ne 24 ]]; then
        echo "FATAL: rsync exited ${st} for: $*" >&2
        return "$st"
    fi
    [[ $st -eq 24 ]] && echo "    (rsync: some files vanished mid-copy; continuing)"
    return 0
}

# Intermediates pggb can always regenerate from the final graph + params.
EXCLUDES=()
if [[ $FULL -eq 0 ]]; then
    EXCLUDES=(
        --exclude '*.seqwish.gfa'
        --exclude '*.seqwish.gfa.prep.*.gfa'
        --exclude '*.seqwish.gfa.smooth.*.gfa'
        --exclude '*.smooth.gfa'
        --exclude '*.smooth.fix.gfa'
        --exclude '*.alignments.wfmash.paf'
        --exclude '*.mappings.wfmash.paf'
        --exclude 'seqwish-*/'
        --exclude 'tmp.norm.vcf'
        --exclude 'decompose/'
        # staged inputs already live in MTB_ARCHIVE; don't duplicate them
        # only the UNCOMPRESSED form is excluded; variants.vcf.gz is synced
        # because it is cheap (0.56 GB total) and saves re-running deconstruct.
        # It is a convenience, not a last copy: the source graphs all live in
        # $MTB_ARCHIVE/data/<graph>/ -- see GRAPH_PROVENANCE.md section 1.
        --exclude 'variants.vcf'
        # superseded by all_variants.collapsed.vcf.gz, and regenerable from
        # variants.vcf via bin/vcf_decompose.sh -- see VARIANT_STRATEGY.md sec 8
        --exclude 'all_variants.decomposed.vcf.gz'
        --exclude 'all_variants.decomposed.vcf.gz.tbi'
    )
    echo "### excluding regenerable intermediates (use --full to keep them)"
fi

echo "--- scripts and config (always)"
rs "$MTB_BIN/"    "$MTB_PERSIST/bin/"
rs "$MTB_CONFIG/" "$MTB_PERSIST/config/"

# Documents, at their repository-relative paths. This step was missing, and the
# omission was quiet in the worst way: the root-level .md files on durable
# storage all carried one date, 2026-09-07, from a one-time copy of the whole
# root, and nothing written afterwards ever landed there as a file. Of 130
# tracked documents only 99 were present and only 83 were byte-identical, so a
# reader opening the persist copy of a pipeline index was reading a version two
# weeks old with no indication of it. The git push below does make them
# recoverable, but recoverable-from-history is not the same as readable, and the
# point of the persist tree is that someone can read it without a checkout.
#
# The tracked set is the definition of "a document this project stands behind",
# which is why the list comes from git rather than a glob: an untracked draft in
# the working tree is deliberately left out, and a document that was renamed or
# deleted in git stops being copied. --delete is NOT used, because the same
# persist tree also holds documents from eras before this repository's history.
echo "--- documents"
if [[ -d "$MTB_WORK/.git" ]]; then
    git -C "$MTB_WORK" ls-files -z '*.md' \
        | rs --from0 --files-from=- "$MTB_WORK/" "$MTB_PERSIST/"
else
    echo "    (no git checkout -- skipping)"
fi

# TOP-LEVEL TRACKED FILES. The documents block above takes every tracked *.md,
# and the directory blocks take bin, config, data, logs and the output trees --
# which between them miss a tracked file sitting in the project root with any
# other extension. Two were missing from the mirror when this was checked on
# 2026-09-27: .gitignore and accessory.txt. Both are recoverable from the git
# mirror, so nothing was lost, but the file tree is supposed to stand on its own.
# The "NOT COPIED" warning below could not catch them either, because it tests
# directories and these are files.
echo "--- top-level tracked files"
if git -C "$MTB_WORK" rev-parse --git-dir >/dev/null 2>&1; then
    # tr -d DELETES the separators and joins every path into one string; the
    # first version of this did that and copied nothing while printing no error.
    # The list is built first, tolerating ONLY grep's no-match status, and the
    # copies then run in this shell so a failed rsync is seen. The old
    # `... | while ...; done || true` also swallowed every rsync failure.
    mapfile -t _top < <(git -C "$MTB_WORK" ls-files --directory -z \
        | tr '\000' '\n' | { grep -v '/' || true; } | { grep -v '\.md$' || true; })
    for _f in "${_top[@]}"; do
        [[ -f "$MTB_WORK/$_f" ]] || continue
        rs "$MTB_WORK/$_f" "$MTB_PERSIST/$_f" || SYNC_FAILED=1
    done
fi
# `|| true` because this runs under `set -e` with pipefail: a grep that matches
# nothing exits 1 and took the whole script down with it, silently, after the
# section header had already printed. The first run of this block copied nothing
# AND skipped every later section including the git push, while still exiting 0
# from the caller's point of view.

echo "--- graphs and results"
rs "${EXCLUDES[@]}" "$MTB_GRAPHS/" "$MTB_PERSIST/results/"

# loci/ holds every per-locus analysis output -- the PPE38, IS6110, CRISPR,
# prophage and W148 results that the project's documents cite. It was missing
# from this script entirely, so after a session of analysis the documents were on
# durable storage via git and the numbers behind them were not. data/annotation
# and data/qc are curated inputs that cannot be regenerated by rerunning anything:
# the 43 canonical spoligotype spacers, the tiered RD table, the repeat mask, the
# provenance flags and the panel definitions.
# --- pipeline output trees ---------------------------------------------------
# refbias/, is6110/ and graphframe/ did not exist when this script was written,
# so for a long time the entire pipeline output tree was outside its scope while
# the script reported success. Their TRACKED contents ride the git mirror, but
# every per-sample table is gitignored by project convention and so lived in
# exactly one place on scratch.
#
# Alignments are excluded by default and the difference is not marginal:
# refbias/ is 73 GB with them and 14 GB without. They are regenerable from P1
# and P2, which is hours of array time but no lost information, whereas an
# output table produced by code that has since changed is not regenerable at
# all. --full includes them.
# Covered on request 2026-09-20 after the coverage audit named them:
#   containers    178 MB, graphaligner.sif and its recipes. A container image
#                 is not regenerable from anything in this repository -- the
#                 recipe may rebuild, but not to the same bits -- so it is the
#                 kind of thing durable storage is for.
#   dist          36 KB, the packaged IS6110 scripts and their README.
#   dr_elements   337 MB, the DR catalogue and its blast outputs.
#   assoc         the association arm: the phyoverlap2 event-matrix writer and
#                 the CallBytestreams it produces. The bytestreams are
#                 regenerable from a merged VCF and a tree, but the writer
#                 itself lives nowhere else -- assoc/bin is not under bin/, so
#                 the scripts section above does not reach it.
OUTPUT_TREES=(refbias is6110 graphframe truth7 containers dist dr_elements assoc sv2frame insgt accessory giraffe)
echo "--- pipeline outputs (${OUTPUT_TREES[*]})"
ALIGN_EXCLUDES=()
if [[ $FULL -eq 0 ]]; then
    ALIGN_EXCLUDES=(
        # The archived alignments ARE copied: bin/archive_alignments.sh
        # writes lossless CRAMs of the two alignments later passes read, at
        # 37% of the BAM size, precisely so they survive a scratch purge.
        # rsync takes the first rule that matches, so these come first.
        --include '*.archive.cram' --include '*.archive.cram.crai'
        --exclude '*.bam'   --exclude '*.bam.bai'
        --exclude '*.cram'  --exclude '*.cram.crai'
        --exclude '*.sam'
        --exclude '*.fastq' --exclude '*.fastq.gz'
        --exclude '*.fq'    --exclude '*.fq.gz'
        --exclude '*.depth'
    )
fi
for _t in "${OUTPUT_TREES[@]}"; do
    [[ -d "$MTB_WORK/$_t" ]] && rs "${ALIGN_EXCLUDES[@]}" \
        "$MTB_WORK/$_t/" "$MTB_PERSIST/$_t/"
done

# A HARDCODED LIST GOES STALE SILENTLY, WHICH IS HOW refbias/ WENT UNSYNCED AND
# HOW truth7/ WOULD HAVE. Anything at the top level that no section of this
# script copies is reported, so a new directory announces itself instead of
# being quietly dropped. This warns rather than guessing, because what belongs
# in durable storage is a decision, not a default.
# `loci` is copied by the analysis-outputs section below, not here
_covered=" bin config graphs data logs slurm loci ${OUTPUT_TREES[*]} "
_unsynced=()
for _d in "$MTB_WORK"/*/; do
    _b="$(basename "$_d")"
    [[ "$_covered" == *" ${_b} "* ]] && continue
    case "$_b" in .*) continue ;; esac
    _unsynced+=("$_b")
done
if [[ ${#_unsynced[@]} -gt 0 ]]; then
    echo "    NOT COPIED by any section of this script -- add to OUTPUT_TREES"
    echo "    if they should be durable:"
    printf '      %s\n' "${_unsynced[@]}"
fi

echo "--- analysis outputs and curated inputs"
[[ -d "$MTB_WORK/loci" ]] && rs "$MTB_WORK/loci/" "$MTB_PERSIST/loci/"
[[ -d "$MTB_DATA/annotation" ]] && rs "$MTB_DATA/annotation/" "$MTB_PERSIST/data/annotation/"
[[ -d "$MTB_DATA/qc" ]] && rs "$MTB_DATA/qc/" "$MTB_PERSIST/data/qc/"
# data/ref holds the H37Rv reference the whole project is coordinate-anchored
# to, and the Actinomycetota comparison genomes. ~123M, and it was on
# purgeable netscratch only -- losing it would break every locus script.
[[ -d "$MTB_DATA/ref" ]] && rs "$MTB_DATA/ref/" "$MTB_PERSIST/data/ref/"
# data/trees was not in this list, so the CX333 phylogeny reached the mirror
# only as git objects and was absent as a file -- and its 24 MB SNP alignment,
# which is untracked, reached it not at all. The tree cost twelve minutes on
# twenty-four cores and is an input to the ancestral-allele work, so it is
# worth the 25 MB. The same omission as the documents step, found the same way:
# by checking for the file on persist rather than trusting a clean exit.
[[ -d "$MTB_DATA/trees" ]] && rs "$MTB_DATA/trees/" "$MTB_PERSIST/data/trees/"

# Assemblies are deliberately NOT synced by default: data/rotated is 2.1 GB
# regenerable from data/assemblies by rotate_to_dnaa.sh, and data/assemblies is
# 626 MB re-downloadable from NCBI. --full takes them.
if [[ $FULL -eq 1 ]]; then
    echo "--- assemblies (--full)"
    [[ -d "$MTB_DATA/assemblies" ]] && rs "$MTB_DATA/assemblies/" "$MTB_PERSIST/data/assemblies/"
    [[ -d "$MTB_DATA/rotated"    ]] && rs "$MTB_DATA/rotated/"    "$MTB_PERSIST/data/rotated/"
fi

echo "--- logs"
[[ -d "$MTB_LOGS"  ]] && rs "$MTB_LOGS/"  "$MTB_PERSIST/logs/"
[[ -d "$MTB_SLURM" ]] && rs "$MTB_SLURM/" "$MTB_PERSIST/logs/slurm/"

# --- git history to the durable bare mirror --------------------------------
# rsync copies the CURRENT state of bin/; it cannot recover a script that was
# overwritten. That happened once already: the graph-derived IS6110 matrix
# builder was lost by reusing its filename, and rsync had faithfully propagated
# the replacement. Pushing the history is what makes that recoverable, so it
# belongs here rather than in someone's memory.
GIT_MIRROR="${MTB_GIT_MIRROR:-${MTB_PERSIST}.git}"
if [[ -d "$MTB_WORK/.git" ]]; then
    echo "--- git history -> $GIT_MIRROR"
    # Only TRACKED modifications are committed, via `add -u`. This used to be
    # `add -A`, which swept untracked files in too, and that is actively harmful
    # when more than one session shares this working tree: another person's
    # half-finished output gets committed under the sync_back identity with an
    # opaque timestamp message, attributed to nobody and described as nothing.
    # It has already happened twice -- a set of hely_tatc_rnaseq_* scripts landed
    # inside an unrelated IS6110 commit, and a refbias/t11/scan_anchored/ run was
    # committed mid-flight. Untracked files are reported and left alone; adding a
    # new file to the repository is a decision for a human, whereas protecting an
    # edit to a file already under version control is what this script is for.
    if [[ ${#DRY[@]} -gt 0 ]]; then
        echo "    (dry run) would commit tracked changes and push main"
    else
        if [[ -n "$(git -C "$MTB_WORK" status --porcelain --untracked-files=no)" ]]; then
            echo "    tracked changes present -- committing them"
            git -C "$MTB_WORK" add -u
            git -C "$MTB_WORK" -c user.name="${MTB_GIT_NAME:-sync_back}" \
                -c user.email="${MTB_GIT_EMAIL:-sync_back@localhost}" \
                commit -q -m "sync_back: $(date -u '+%Y-%m-%dT%H:%M:%SZ')" || true
        fi
        # Untracked files are still rsynced to durable storage above; they are
        # simply not committed. Name them so nothing new is silently left out of
        # version control.
        _untracked=$(git -C "$MTB_WORK" ls-files --others --exclude-standard)
        if [[ -n "$_untracked" ]]; then
            echo "    NOT committed (untracked -- 'git add' them yourself if they belong in git):"
            # Read the first 20 from a here-string rather than piping printf into
            # head.  Under `set -o pipefail` a producer killed by SIGPIPE makes the
            # pipeline return 141, and `set -e` then terminated this script here --
            # before the push below -- while still reporting success to the caller.
            # head reading a here-string has no producer to signal.
            _n_untracked=$(wc -l <<<"$_untracked")
            sed -n '1,20p' <<<"$_untracked" | sed 's/^/      /'
            if (( _n_untracked > 20 )); then
                echo "      ... and $(( _n_untracked - 20 )) more"
            fi
        fi
        [[ -d "$GIT_MIRROR" ]] || {
            echo "    creating bare mirror"
            git init -q --bare "$GIT_MIRROR"
            git -C "$GIT_MIRROR" symbolic-ref HEAD refs/heads/main
        }
        git -C "$MTB_WORK" remote get-url persist >/dev/null 2>&1 \
            || git -C "$MTB_WORK" remote add persist "$GIT_MIRROR"
        if git -C "$MTB_WORK" push -q persist main; then
            echo "    pushed $(git -C "$MTB_WORK" rev-parse --short HEAD)"
        else
            echo "    ERROR: git push failed; history is NOT mirrored" >&2
            SYNC_FAILED=1
        fi
    fi
else
    echo "--- git history: no repository at $MTB_WORK/.git, skipping" >&2
fi

if [[ "$SYNC_FAILED" -ne 0 ]]; then
    echo "### FAILED: see the errors above; .last_sync NOT updated" >&2
    exit 1
fi
echo "### done."
if [[ ${#DRY[@]} -eq 0 ]]; then
    date -u '+%Y-%m-%dT%H:%M:%SZ synced from '"$MTB_WORK" >> "$MTB_PERSIST/.last_sync"
    du -sh "$MTB_PERSIST"
fi
