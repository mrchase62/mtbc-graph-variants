#!/usr/bin/env bash
# Copy the key scale200 outputs, old run and fixed rerun, from netscratch
# (90-day purge) to the durable lab mirror, under names that say what they
# are. The user asked for this on 2026-10-08 ("copy the key outputs to the
# mirror. directories will need to be renamed, results or out. refbias doesn't
# make sense"). BAMs and per-sample working files are not copied.
#
#   results/cohorts/scale200/
#     README.md
#     COMPARISON.md                   old vs new report (analysis/rerun_scale200)
#     2026-10-01_original/            build 7713a8d71d8e
#     2026-10-08_audit_fixes/         build 7713a8d71d8e-fix1
#       calls/                 merged VCF (+ index), site and key tables, P5 checks
#       structural_variants/   SV interval catalogue and SV matrix
#       reference_choice/      P1 matched reference per sample, P2 summary
#       association/           scan, gene burdens, event matrix, chain audit, phenotype
#       accessory_presence/    per-sample accessory presence tables
#       trees/                 cohort and cohort+panel trees, alignments
#
# rsync -a, never --delete: re-running adds or refreshes, it removes nothing.
# Run from runroot.
set -euo pipefail
M=/n/boslfs02/LABS/sfortune_lab/Lab/mchase/MtbPangenome/results/cohorts/scale200
[[ -d refbias && -d assoc && -d accessory ]] || { echo "run from runroot" >&2; exit 1; }

copy_run() {  # <pipeline cohort name> <destination folder name> <tree prefixes...>
    local c="$1" d="$M/$2"; shift 2
    mkdir -p "$d"/{calls,structural_variants,reference_choice,association,accessory_presence,trees}
    rsync -a refbias/"$c"/p5/{merged.vcf.gz,merged.vcf.gz.tbi,sites.tsv,keys.tsv,validation.tsv,sanity.tsv,states.meta.tsv} "$d/calls/"
    [[ -e refbias/"$c"/.mtb_build ]] && rsync -a refbias/"$c"/.mtb_build "$d/calls/build_stamp.txt"
    rsync -a refbias/"$c"/p5/{sv_intervals.tsv,sv_matrix.tsv} "$d/structural_variants/"
    rsync -a refbias/"$c"/p1/refmap.tsv* "$d/reference_choice/"
    rsync -a refbias/"$c"/p2/p2_summary.tsv "$d/reference_choice/"
    rsync -aL assoc/"$c"/ "$d/association/"
    rsync -aL accessory/"$c"/ "$d/accessory_presence/"
    for p in "$@"; do
        compgen -G "data/trees/${p}.*" >/dev/null && rsync -a data/trees/"$p".* "$d/trees/"
    done
    echo "$c -> $d: $(du -sh "$d" | cut -f1)"
}

copy_run scale200     2026-10-01_original    scale200 scale200_cx333
copy_run scale200_fix 2026-10-08_audit_fixes scale200_fix
rsync -a ../analysis/rerun_scale200/COMPARISON.md "$M/COMPARISON.md"
rsync -a ../analysis/rerun_scale200/archive_README.md "$M/README.md"
