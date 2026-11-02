#!/usr/bin/env bash
# cramit run — convert and validate the run's BAMs on this machine.
#
# Usage: cramit run [-j PARALLEL] [-t THREADS] [--retry-failed] RUN
#   -j PARALLEL     BAMs converted at the same time (default 1)
#   -t THREADS      threads per BAM (default 4)
#   --retry-failed  also retry BAMs that failed before (once you have fixed them)
#
# Safe to stop and start again: finished BAMs are skipped, and an interrupted
# BAM starts over cleanly.
set -euo pipefail
source "$CRAMIT_HOME/lib/common.sh"

PAR=1; THREADS=4; RETRY=0; ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -j) PAR="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        --retry-failed) RETRY=1; shift ;;
        -h|--help) show_usage "$0" ;;
        -*) die "unknown option: $1" ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} -eq 1 ] || show_usage "$0" 1
[[ "$PAR" =~ ^[1-9][0-9]*$ && "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "-j and -t must be positive integers"
load_run "${ARGS[0]}"
need_samtools

TODO=$(pending_samples "$RETRY")
N=$(grep -c . <<< "$TODO" || true)
if [ "$N" -eq 0 ]; then
    echo "nothing to do — every BAM is done, running, or failed (see: cramit status -v $RUN)"
    exit 0
fi
echo "converting $N BAM(s), $PAR at a time, $THREADS thread(s) each"

run_one() {
    local key bam st
    IFS=$'\t' read -r key bam _ <<< "$1"
    bash "$CRAMIT_HOME/lib/convert.sh" -r "$REF" -l "$RUN/samples/$key" -t "$THREADS" "$bam" > /dev/null 2>&1 || true
    st=$(cut -f1,2 "$RUN/samples/$key/status" 2>/dev/null | tr '\t' ' ')
    case "$st" in ok*) st=done ;; "failed "*) st="failed: ${st#failed }" ;; *) st="stopped" ;; esac
    printf '%s  %-24s %s\n' "$(date +%T)" "$st" "$bam"
}
export -f run_one
export RUN REF THREADS

printf '%s\n' "$TODO" | xargs -d '\n' -P "$PAR" -I{} bash -c 'run_one "$1"' _ {}
echo ""
bash "$CRAMIT_HOME/lib/status.sh" "$RUN"
