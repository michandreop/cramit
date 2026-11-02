#!/usr/bin/env bash
# cramit status — how far a run has got, and what needs attention.
#
# Usage: cramit status [-v] RUN
#   -v  list every BAM that is not done, with the reason and what to do
#
# Also writes RUN/status.tsv (key, state, reason, detail, bam, cram) for your
# own scripts. Reasons are explained in docs/SPECIAL_CASES.md.
set -euo pipefail
source "$CRAMIT_HOME/lib/common.sh"

VERBOSE=0; ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -v) VERBOSE=1; shift ;;
        -h|--help) show_usage "$0" ;;
        -*) die "unknown option: $1" ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} -eq 1 ] || show_usage "$0" 1
load_run "${ARGS[0]}"

TSV="$RUN/status.tsv"
{ printf 'key\tstate\treason\tdetail\tbam\tcram\n'; sample_states; } > "$TSV.tmp" && mv -f "$TSV.tmp" "$TSV"

# BAMs already deleted, among the done ones (one stat call for all of them)
GONE=$(awk -F'\t' 'NR>1 && $2=="done"{print $5}' "$TSV" > "$RUN/.done_bams"
       n=$(wc -l < "$RUN/.done_bams")
       e=$(xargs -r -d '\n' stat -c '%n' < "$RUN/.done_bams" 2>/dev/null | wc -l || true)
       rm -f "$RUN/.done_bams"; echo $((n - e)))

hint() {
    case "$1" in
        bam_missing)           echo "the BAM path no longer exists" ;;
        bam_corrupt)           echo "re-download or regenerate the BAM, then: cramit submit --retry-failed (or run --retry-failed)" ;;
        wrong_reference)       echo "aligned to another reference — convert it in a separate run with the right -r" ;;
        not_sorted)            echo "sort it (samtools sort) or fix the header's SO tag, then retry with --retry-failed" ;;
        existing_cram_differs) echo "a different CRAM is already at that path — move or delete it, then retry with --retry-failed" ;;
        cram_differs)          echo "the new CRAM did not match and was discarded — retry with --retry-failed; if it happens again, report it" ;;
        timeout)               echo "resubmit with more time: cramit submit -T 48:00:00 RUN" ;;
        out_of_memory)         echo "resubmit with more memory: cramit submit -m 16G RUN" ;;
        cancelled|incomplete)  echo "cramit submit / cramit run picks it up again" ;;
        *)                     echo "see the log in RUN/samples/<key>/convert.log" ;;
    esac
}

awk -F'\t' -v gone="$GONE" -v run="$RUN" -v attn="$(wc -l < "$RUN/attention.txt" 2>/dev/null || echo 0)" '
NR == 1 { next }
{ n++; st[$2]++; if ($2 == "failed") why[$3]++ }
function row(label, k, note) { if (st[k]) printf "  %-16s %6d   %s\n", label, st[k], note }
END {
    printf "cramit status: %s\n", run
    printf "  %-16s %6d\n", "BAMs in run", n
    row("done", "done", gone ? "(" gone " BAM(s) already deleted)" : "")
    row("running", "running", "")
    row("not started", "not_started", "")
    row("failed", "failed", "needs a fix — cramit status -v")
    for (k in why) printf "      %-22s %d\n", k, why[k]
    row("timeout", "timeout", "cramit submit -T <longer> RUN")
    row("out of memory", "out_of_memory", "cramit submit -m <more> RUN")
    row("cancelled", "cancelled", "picked up again by the next submit/run")
    row("incomplete", "incomplete", "picked up again by the next submit/run")
    if (attn > 0) printf "  %-16s %6d   %s\n", "set aside", attn, "by scan — see " run "/attention.txt"
    left = st["done"] - gone
    if (left > 0) printf "\n%d BAM(s) can be deleted now: cramit delete-bams %s\n", left, run
}' "$TSV"

if [ "$VERBOSE" -eq 1 ]; then
    echo ""
    awk -F'\t' -v OFS='\t' 'NR>1 && $2!="done" && $2!="running" && $2!="not_started" {print ($2=="failed" ? $3 : $2), $1, $4}' "$TSV" \
    | sort | while IFS=$'\t' read -r why key detail; do
        printf '%s  %s\n    %s\n    → %s\n' "$why" "$key" "$detail" "$(hint "$why" | sed "s#RUN#$RUN#")"
    done
fi
