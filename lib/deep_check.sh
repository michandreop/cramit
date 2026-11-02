#!/usr/bin/env bash
# cramit deep-check — compare a BAM and its CRAM record by record, aux tags included.
#
# Usage: cramit deep-check -r REF [-n N] BAM CRAM
#   -r REF  reference the CRAM was made with
#   -n N    compare only the first N records (default: all)
#
# The standard validation compares read names, sequences, qualities and flags,
# not aux tags (UMI tags, MM/ML methylation, MD/NM, ...). Use this where tags
# matter. CRAM may store tags in another order, so each record's tags are
# sorted before comparing.
#
# Exit status: 0 identical, 1 different (first difference printed), 2 error.
set -uo pipefail
source "$CRAMIT_HOME/lib/common.sh"

REF=""; N=0
while getopts "r:n:h" opt; do
    case $opt in r) REF="$OPTARG" ;; n) N="$OPTARG" ;; h) show_usage "$0" ;; *) show_usage "$0" 2 ;; esac
done
shift $((OPTIND - 1))
[ $# -eq 2 ] && [ -n "$REF" ] || show_usage "$0" 2
BAM="$1"; CRAM="$2"
for f in "$BAM" "$CRAM" "$REF"; do [ -f "$f" ] || { echo "cramit: not found: $f" >&2; exit 2; }; done
[[ "$N" =~ ^[0-9]+$ ]] || { echo "cramit: -n must be a number" >&2; exit 2; }
need_samtools

# records, with each record's aux tags (fields 12+) sorted
records() {
    local lim=cat; [ "$N" -gt 0 ] && lim="head -n $N"
    samtools view "$@" 2>/dev/null | $lim | awk -F'\t' 'BEGIN{OFS="\t"} {
        n = 0
        for (i = 12; i <= NF; i++) { v = $i; j = n; while (j > 0 && t[j] > v) { t[j+1] = t[j]; j-- } t[j+1] = v; n++ }
        s = $1; for (i = 2; i <= 11 && i <= NF; i++) s = s OFS $i
        for (i = 1; i <= n; i++) s = s OFS t[i]
        print s }'
}

out=$(cmp <(records "$BAM") <(records -T "$REF" "$CRAM") 2>&1)
rc=$?
scope=$([ "$N" -gt 0 ] && echo "first $N records" || echo "all records")
if [ "$rc" -eq 0 ]; then
    echo "IDENTICAL ($scope, aux tags included): $CRAM"
    exit 0
fi
line=$(sed -n 's/.*line \([0-9]*\).*/\1/p' <<< "$out")
if [ -n "$line" ]; then
    echo "DIFFERENT at record $line ($scope): $CRAM"
    echo "  BAM : $(records "$BAM" | sed -n "${line}p" | cut -c1-300)"
    echo "  CRAM: $(records -T "$REF" "$CRAM" | sed -n "${line}p" | cut -c1-300)"
else
    echo "DIFFERENT ($scope): $out"
fi
exit 1
