#!/bin/bash
# Real-data check of our caller wired into production (user, 2026-10-12).
# Cohort fbrecall (48 isolates, P2 finished, BAMs kept):
#   1. bin/p2_call.sh per sample, with the runner's folders: a finished
#      sample runs only the new step (mask, caller, assembly, element
#      genotype). 27 of 29 masks were copied from the evaluation
#      (novel_events/out/masks); GCF_002357935 and GCF_002357955 are built
#      here, to compare.
#   2. bin/p4b_place_sv.sh per sample into a test folder (run/p4b_test).
# compare.py then checks the events against recall/out/asm_af/keep_na_0.1
# and scores the P4b rows. Submit from runroot:
#   sbatch ../analysis/sv_wiring/01_fbrecall_check.sbatch
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 03:00:00
#SBATCH -J sv_wiring
#SBATCH -o ../analysis/sv_wiring/logs/fbrecall_check_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ../config/project_env.sh
export MTB_BUILD_DIR=refbias/build/7713a8d71d8e-fix1
export REFMAP=refbias/fbrecall/run/p1/refmap.tsv CRAMMAP=refbias/fbrecall/crams.tsv
L=../analysis/sv_wiring/logs
tail -n +2 $REFMAP | cut -f1 | \
  OUTDIR=refbias/fbrecall/run/p2 WORK=refbias/work/fbrecall_p2 SLURM_CPUS_PER_TASK=1 \
  xargs -P 8 -I{} bash -c 'bash bin/p2_call.sh {} > '$L'/p2.{}.log 2>&1 || echo "P2 FAILED {}"'
echo "P2 step finished: $(ls refbias/fbrecall/run/p2/*.ours.events.tsv | wc -l) of 48 ours tables"
tail -n +2 $REFMAP | cut -f1 | \
  P2DIR=refbias/fbrecall/run/p2 OUTDIR=refbias/fbrecall/run/p4b_test WORK=refbias/work/fbrecall_p4b_test \
  SLURM_CPUS_PER_TASK=2 xargs -P 4 -I{} bash -c 'bash bin/p4b_place_sv.sh {} > '$L'/p4b.{}.log 2>&1 || echo "P4b FAILED {}"'
echo "P4b step finished: $(ls refbias/fbrecall/run/p4b_test/*.sv_placed.tsv | wc -l) of 48 tables"
