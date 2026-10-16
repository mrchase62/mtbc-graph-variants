#!/usr/bin/env bash
#SBATCH -J chain_watch
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=200M
#SBATCH -t 1-00:00
#SBATCH -o slurm/chain_watch_%j.out
#
# Watch one cohort's chain while it runs, so lost run time is seen in minutes,
# not the next morning (user, 2026-10-15: "put some checks into the scripts to
# prevent all this lost run time"). What new1k_p200 lost, and what this does:
#
#   FAILED / TIMEOUT / OUT_OF_MEMORY tasks sat unnoticed for hours (a P5
#     states array hit its 6 h limit 32 times; p1iv failed and every pass
#     after it waited on it). Each one is reported the cycle it appears, and
#     with --stop-on-fail the cohort's remaining jobs are cancelled, because
#     afterok dependents of a failed task can never start anyway.
#   An array sat PENDING (JobArrayTaskLimit) with none of its tasks running:
#     Slurm's per-array running count had wedged (new1k_p200 p5svgtproj and
#     p5sv2frame, 2026-10-15). Detected and cleared by re-setting the array's
#     throttle (scontrol update ArrayTaskThrottle), which is what freed them.
#   A task running far past its pass's typical length is reported (more than
#     --slow-factor times the median of that pass's completed tasks, and over
#     --slow-min minutes), so a slow pass is seen before it times out.
#
#   sbatch bin/chain_watch.sh <cohort> [--every MIN] [--stop-on-fail]
#   bash   bin/chain_watch.sh <cohort> --once        # one check, then exit
#
# Jobs are found by name, rb_<cohort>_*, the names bin/refbias_run.sh gives
# them; refbias_run.sh --watch submits this alongside a chain. It exits when
# no job of the cohort is left in the queue, after printing each pass's wall
# time and CPU-hours.
set -uo pipefail
COHORT="${1:?usage: $0 <cohort> [--every MIN] [--stop-on-fail] [--once] [--slow-factor X] [--slow-min MIN]}"; shift
EVERY=5; STOP=0; ONCE=0; SLOWF=5; SLOWMIN=30
while [[ $# -gt 0 ]]; do
    case "$1" in
        --every) EVERY="$2"; shift 2 ;;
        --stop-on-fail) STOP=1; shift ;;
        --once) ONCE=1; shift ;;
        --slow-factor) SLOWF="$2"; shift 2 ;;
        --slow-min) SLOWMIN="$2"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
PFX="rb_${COHORT}_"
START="$(date +%Y-%m-%dT%H:%M:%S)"
START="${CHAIN_WATCH_SINCE:-$START}"
declare -A REPORTED=()
say() { echo "[watch ${COHORT} $(date +%H:%M)] $*"; }

check() {
    # 1. failures since the watch began
    while read -r id name state; do
        [[ -n "$id" && -z "${REPORTED[$id]:-}" ]] || continue
        REPORTED[$id]=1
        say "FAILURE: ${name} task ${id} ${state}"
        if [[ "$STOP" -eq 1 ]]; then
            say "--stop-on-fail: cancelling every queued job of ${COHORT}"
            squeue -u "$USER" -h -o "%A %j" | awk -v p="$PFX" 'index($2,p)==1{print $1}' \
                | sort -u | xargs -r scancel
        fi
    done < <(sacct -u "$USER" -S "$START" -n -X -P -o JobID,JobName,State 2>/dev/null \
             | awk -F'|' -v p="$PFX" 'index($2,p)==1 && $3 ~ /^(FAILED|TIMEOUT|OUT_OF_ME|NODE_FAIL)/ {print $1, $2, $3}')
    # 2. an array waiting on its own task limit with nothing of it running
    while read -r aid name; do
        [[ -n "$aid" ]] || continue
        running="$(squeue -u "$USER" -h -r -j "$aid" -t RUNNING 2>/dev/null | wc -l)"
        if [[ "$running" -eq 0 ]]; then
            thr="$(scontrol show job "$aid" 2>/dev/null | grep -o 'ArrayTaskThrottle=[0-9]*' | head -1 | cut -d= -f2)"
            thr="${thr:-50}"
            say "STALL: ${name} (${aid}) pending on JobArrayTaskLimit with no task running; resetting its throttle (${thr})"
            scontrol update jobid="$aid" ArrayTaskThrottle="$(( thr + 1 ))" >/dev/null 2>&1 \
                && scontrol update jobid="$aid" ArrayTaskThrottle="$thr" >/dev/null 2>&1
        fi
    done < <(squeue -u "$USER" -h -o "%F %j %T %r" 2>/dev/null \
             | awk -v p="$PFX" 'index($2,p)==1 && $3=="PENDING" && $4=="JobArrayTaskLimit" {print $1, $2}' | sort -u)
    # 3. tasks running far longer than their pass's completed median
    while read -r id name el; do
        [[ -n "$id" && -z "${REPORTED[slow$id]:-}" ]] || continue
        REPORTED[slow$id]=1
        say "SLOW: ${name} task ${id} running $(( el / 60 )) min, over ${SLOWF}x its pass's median"
    done < <(
        sacct -u "$USER" -S "$START" -n -X -P -o JobID,JobName,State,ElapsedRaw 2>/dev/null \
          | awk -F'|' -v p="$PFX" -v f="$SLOWF" -v m="$SLOWMIN" '
              index($2,p)!=1 {next}
              $3=="COMPLETED" {n[$2]++; v[$2,n[$2]]=$4; next}
              $3=="RUNNING"   {r[$1]=$2; e[$1]=$4}
              END {
                for (j in n) { k=n[j]; for(a=1;a<=k;a++) s[a]=v[j,a];
                  for(a=1;a<=k;a++) for(b=a+1;b<=k;b++) if(s[b]<s[a]){t=s[a];s[a]=s[b];s[b]=t}
                  med[j]=s[int((k+1)/2)] }
                for (i in r) if ((r[i] in med) && e[i] > f*med[r[i]] && e[i] > m*60) print i, r[i], e[i]
              }')
}

summary() {
    say "no job of ${COHORT} left in the queue; per pass since ${START}:"
    sacct -u "$USER" -S "$START" -n -X -P -o JobName,State,Start,End,ElapsedRaw,AllocCPUS 2>/dev/null \
      | awk -F'|' -v p="$PFX" 'index($1,p)==1 {
            pass=substr($1,length(p)+1); n[pass]++; st[pass" "$2]++
            if (!(pass in s) || $3<s[pass]) s[pass]=$3; if ($4>e[pass]) e[pass]=$4; c[pass]+=$5*$6/3600; tot+=$5*$6/3600 }
          END { for (q in n) printf "  %-14s tasks %5d  %s .. %s  CPU-h %7.1f\n", q, n[q], s[q], e[q], c[q]
                printf "  total CPU-h %.1f\n", tot }' | sort
}

while true; do
    check
    [[ "$ONCE" -eq 1 ]] && exit 0
    if ! squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -q "^${PFX}"; then
        summary; exit 0
    fi
    sleep "$(( EVERY * 60 ))"
done
