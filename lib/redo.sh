#!/usr/bin/env bash
# cramit redo — convert some BAMs of a run again.
#
# Usage: cramit redo RUN SAMPLE...
#   SAMPLE  the BAM's path, its file name (with or without .bam), or its key
#           (the folder name in RUN/samples/)
#
# For a BAM that is done, redo deletes its CRAM (and .crai) and resets it, so
# the next `cramit run RUN` or `cramit submit RUN` converts and validates it
# again. Use it when a CRAM has to be made again: it changed after validation,
# or it was made while the reference sat at another path and you want the
# header to record the new one.
# For a BAM that failed or stopped, redo only resets it, so it is retried
# without --retry-failed. A CRAM that cramit did not validate is left alone.
#
# redo never touches a BAM that is still running, or one whose BAM is gone
# (its CRAM is then the only copy of the reads).
set -euo pipefail
source "$CRAMIT_HOME/lib/common.sh"

ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) show_usage "$0" ;;
        -*) die "unknown option: $1" ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} -ge 2 ] || show_usage "$0" 1
load_run "${ARGS[0]}"

STATES=$(sample_states)
reset=0; refused=0
for q in "${ARGS[@]:1}"; do
    row=$(find_sample "$q") || { refused=$((refused + 1)); continue; }
    IFS=$'\t' read -r key bam cram <<< "$row"
    state=$(awk -F'\t' -v k="$key" '$1 == k { print $2 }' <<< "$STATES")
    d="$RUN/samples/$key"
    case "$state" in
        running)     echo "  kept   $key: still running"; refused=$((refused + 1)); continue ;;
        not_started) echo "  -      $key: not converted yet, nothing to redo"; continue ;;
    esac
    if [ ! -e "$bam" ]; then
        echo "  kept   $key: its BAM is gone, so the CRAM is the only copy"; refused=$((refused + 1)); continue
    fi
    if [ "$state" = done ]; then
        rm -f -- "$cram" "$cram.crai"
        what="CRAM deleted, will be converted again"
    else
        what="reset from '$state', will be retried"
        [ ! -e "$cram" ] || what="$what (the CRAM next to it was not validated by cramit and is left alone)"
    fi
    rm -f -- "$d/status"
    printf '%s  cramit redo: %s\n' "$(date '+%F %T')" "$what" >> "$d/convert.log"
    echo "  redo   $key: $what"
    reset=$((reset + 1))
done

echo ""
if [ "$reset" -gt 0 ]; then
    echo "$reset BAM(s) to convert again. Next:  cramit run $RUN   or   cramit submit $RUN"
fi
[ "$refused" -eq 0 ]
