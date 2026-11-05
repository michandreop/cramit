#!/usr/bin/env bash
# cramit deep-check — compare a BAM and its CRAM record by record, every field
# and aux tag included.
#
# Usage: cramit deep-check [-n N] [-e K] -r REF BAM CRAM
#        cramit deep-check [-n N] [-e K] RUN SAMPLE
#   -r REF  reference the CRAM was made with
#   RUN SAMPLE  take the BAM, CRAM and reference from a run; SAMPLE is the
#           BAM's path, its file name (with or without .bam) or its key
#   -n N    compare only the first N records (default: all)
#   -e K    show the first K records that differ (default 3)
#
# The standard validation compares read names, sequences, qualities and flags.
# deep-check also compares positions, CIGAR, MAPQ, mate fields and aux tags
# (UMI tags, MM/ML methylation, MD/NM, ...). It reads both files to the end,
# counts the differences by field and shows the first few. Tags are sorted
# before comparing, since CRAM may store them in another order.
#
# Two differences are expected and counted on their own lines, not as errors
# (details in docs/SPECIAL_CASES.md):
#   - CRAM does not store the CIGAR or MAPQ of unmapped reads (flag 0x4), so
#     they come back as * and 0. Some aligners (older BWA) leave values there;
#     the SAM spec says they carry no meaning.
#   - samtools adds MD and NM tags, computed from the reference, to reads that
#     had none when it reads a CRAM.
#
# Takes 1 to 1.5 hours per 100 million reads. On SLURM:
#   sbatch -c 3 --mem 4G -t 8:00:00 --wrap "cramit deep-check RUN SAMPLE"
#
# Exit status: 0 identical (apart from that expected difference),
# 1 different, 2 error.
set -uo pipefail
source "$CRAMIT_HOME/lib/common.sh"

REF=""; N=0; EX=3
while getopts "r:n:e:h" opt; do
    case $opt in
        r) REF="$OPTARG" ;; n) N="$OPTARG" ;; e) EX="$OPTARG" ;;
        h) show_usage "$0" ;; *) show_usage "$0" 2 ;;
    esac
done
shift $((OPTIND - 1))
[ $# -eq 2 ] || show_usage "$0" 2
err() { echo "cramit: $*" >&2; exit 2; }
[[ "$N" =~ ^[0-9]+$ ]]  || err "-n must be a number"
[[ "$EX" =~ ^[0-9]+$ ]] || err "-e must be a number"

if [ -z "$REF" ] && [ -f "$1/cramit.conf" ]; then
    load_run "$1"
    row=$(find_sample "$2") || exit 2
    IFS=$'\t' read -r _ BAM CRAM <<< "$row"
else
    [ -n "$REF" ] || err "give -r REF BAM CRAM, or RUN SAMPLE"
    BAM="$1"; CRAM="$2"
fi
for f in "$BAM" "$CRAM" "$REF"; do [ -f "$f" ] || err "not found: $f"; done
( need_samtools ) >/dev/null 2>&1 || err "samtools not found — put it on PATH, or set CRAMIT_SETUP"
need_samtools

T=$(mktemp -d "${TMPDIR:-/tmp}/cramit-deep.XXXXXX") || err "cannot make a temporary folder"
pa=""; pb=""
trap '[ -z "$pa" ] || kill $pa $pb 2>/dev/null; rm -rf "$T"' EXIT
mkfifo "$T/bam" "$T/cram" || err "cannot make named pipes in $T"
samtools view "$BAM" > "$T/bam" 2> "$T/bam.err" & pa=$!
samtools view -T "$REF" "$CRAM" > "$T/cram" 2> "$T/cram.err" & pb=$!

# Records side by side: BAM record \001 CRAM record. When one file runs out
# first, its side is empty.
paste -d $'\001' "$T/bam" "$T/cram" | awk -F'\001' -v limit="$N" -v ex="$EX" -v early="$T/early" '
function sorted(r, from, to,   i, j, t, s, v, k) {
    k = 0
    for (i = from; i <= to; i++) { v = r[i]; j = k; while (j > 0 && t[j] > v) { t[j+1] = t[j]; j-- } t[j+1] = v; k++ }
    s = ""; for (i = 1; i <= k; i++) s = s "\t" t[i]; return s
}
function cut_(s) { return length(s) > 120 ? substr(s, 1, 117) "..." : s }
function tagdiff(   i, x, y, s) {
    delete x; delete y; s = ""
    for (i = 12; i <= na; i++) x[a[i]] = 1
    for (i = 1; i <= nc; i++) y[c[i]] = 1
    for (i in x) if (!(i in y)) s = s "\n        only in BAM:  " cut_(i)
    for (i in y) if (!(i in x)) s = s "\n        only in CRAM: " cut_(i)
    return s
}
function note(what) { bad = bad (bad == "" ? "" : ", ") what }
BEGIN {
    split("QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL", fname, " ")
}
{
    if ($1 == "" || $2 == "" || NF != 2) {
        stop = $2 == "" ? "the CRAM has fewer records than the BAM" : "the CRAM has more records than the BAM"
        stopat = recs + 1; print "" > early; exit
    }
    na = split($1, a, "\t"); nb = split($2, b, "\t")
    if (a[1] != b[1]) {
        stop = "the read names stop lining up (a read is missing, extra or out of order)"
        stopat = recs + 1; sa = $1; sb = $2; print "" > early; exit
    }
    recs++
    unm = int(a[2] / 4) % 2
    if (unm) unmapped++
    bad = ""; dropped = 0
    for (f = 2; f <= 11; f++) {
        if (a[f] == b[f]) continue
        # expected: CRAM does not store the CIGAR/MAPQ of unmapped reads
        if (unm && f == 6 && b[6] == "*") { cig++; dropped = 1; continue }
        if (unm && f == 5 && b[5] == 0)   { mq++;  dropped = 1; continue }
        d[fname[f]]++; note(fname[f])
    }
    if (dropped) norm++
    # expected: the CRAM reader adds MD/NM tags to reads that had none
    md = 0; nm = 0
    for (i = 12; i <= na; i++) { p = substr(a[i], 1, 3); if (p == "MD:") md = 1; else if (p == "NM:") nm = 1 }
    nc = 0; add = 0
    for (i = 12; i <= nb; i++) {
        p = substr(b[i], 1, 3)
        if ((p == "MD:" && !md) || (p == "NM:" && !nm)) { add = 1; continue }
        c[++nc] = b[i]
    }
    if (add) added++
    if (sorted(a, 12, na) != sorted(c, 1, nc)) { d["tags"]++; note("tags") }
    if (bad != "") {
        other++
        if (other <= ex) {
            shown = shown sprintf("\n    record %d (%s): %s differ", recs, a[1], bad)
            for (f = 2; f <= 11; f++) if (index(", " bad ",", ", " fname[f] ","))
                shown = shown sprintf("\n        %-6s BAM: %s   CRAM: %s", fname[f], cut_(a[f]), cut_(b[f]))
            if (index(bad, "tags")) shown = shown tagdiff()
        }
    }
    if (limit && recs >= limit) { print "" > early; exit }
}
END {
    printf "  records compared        %d%s\n", recs, limit && recs >= limit ? " (the first " limit ")" : (stop ? "" : " (all)")
    printf "  unmapped reads          %d\n", unmapped
    if (norm || added) {
        print "  expected differences (see docs/SPECIAL_CASES.md):"
        if (norm)  printf "    unmapped reads whose CIGAR/MAPQ CRAM does not store  %d (CIGAR %d, MAPQ %d)\n", norm, cig, mq
        if (added) printf "    reads given MD/NM tags by the CRAM reader            %d\n", added
    }
    printf "  records that differ     %d\n", other
    for (k in d) printf "    %-8s %d\n", k, d[k]
    if (shown != "") printf "  first differences:%s\n", shown
    if (stop != "") {
        printf "  stopped at record %d: %s; nothing after it was compared\n", stopat, stop
        if (sa != "") printf "        BAM:  %s\n        CRAM: %s\n", cut_(sa), cut_(sb)
    }
    if (other || stop != "") { print "DIFFERENT"; exit 1 }
    print (norm || added) ? "IDENTICAL apart from the expected differences above" : "IDENTICAL"
    exit 0
}' > "$T/report"
rc=${PIPESTATUS[1]}
wait "$pa"; ra=$?; wait "$pb"; rb=$?; pa=""
# samtools must end cleanly; when the comparison stopped early it is stopped
# by a broken pipe (exit 141), which is fine.
for side in bam:$ra cram:$rb; do
    r=${side#*:}
    if [ "$r" -ne 0 ] && ! { [ "$r" -eq 141 ] && [ -e "$T/early" ]; }; then
        echo "cramit deep-check: samtools could not read the ${side%%:*} (exit $r):" >&2
        sed -n '1,5p' "$T/${side%%:*}.err" >&2
        exit 2
    fi
done
[ "$rc" -le 1 ] || err "comparison failed (awk exit $rc)"
echo "cramit deep-check"
echo "  BAM   $BAM"
echo "  CRAM  $CRAM"
echo "  REF   $REF"
cat "$T/report"
for side in bam cram; do
    [ -s "$T/$side.err" ] && { echo "  samtools said about the $side:"; sed -n '1,3p' "$T/$side.err" | sed 's/^/    /'; }
done
exit "$rc"
