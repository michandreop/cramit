#!/usr/bin/env bash
# cramit convert — convert one BAM to CRAM and prove the CRAM holds the same reads.
#
# Usage: cramit convert -r REF [-o OUTDIR] [-l LOGDIR] [-t THREADS] BAM
#   -r REF      reference FASTA the BAM was aligned to (.fa, or bgzipped .fa.gz)
#   -o OUTDIR   write the CRAM here (default: next to the BAM, as NAME.cram)
#   -l LOGDIR   log, stats and status file (default: ./cramit_logs/<name>__<hash>)
#   -t THREADS  threads (default: $SLURM_CPUS_PER_TASK, else 4)
#
# The BAM is read once: the same stream feeds the CRAM writer, samtools stats
# and samtools flagstat. The new CRAM is then decoded once for its own stats
# and flagstat. It is published (renamed to NAME.cram) only if both match, so
# a NAME.cram written by cramit is always complete and validated.
#
# If NAME.cram already exists, it is not overwritten: it is checked against
# the BAM and kept if it matches.
#
# Exit status: 0 validated CRAM in place, 1 not converted (input or conversion
# problem), 2 the CRAM does not match the BAM.
set -euo pipefail
source "$CRAMIT_HOME/lib/common.sh"

REF=""; OUTDIR=""; LOGDIR=""; THREADS="${SLURM_CPUS_PER_TASK:-4}"
while getopts "r:o:l:t:h" opt; do
    case $opt in
        r) REF="$OPTARG" ;; o) OUTDIR="$OPTARG" ;; l) LOGDIR="$OPTARG" ;; t) THREADS="$OPTARG" ;;
        h) show_usage "$0" ;; *) show_usage "$0" 1 ;;
    esac
done
shift $((OPTIND - 1))
[ $# -eq 1 ] || show_usage "$0" 1
[ -n "$REF" ] || die "-r REF is required"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "-t must be a positive integer"
[[ "$1" == *.bam ]] || die "not a .bam file: $1"

BAM=$(abspath "$1"); REF=$(abspath "$REF"); NAME=$(basename "$BAM" .bam)
if [ -n "$OUTDIR" ]; then CRAM="$(abspath "$OUTDIR")/$NAME.cram"; else CRAM="${BAM%.bam}.cram"; fi
LOGDIR=$(abspath "${LOGDIR:-cramit_logs/$(sample_key "$BAM")}")
LOG="$LOGDIR/convert.log"
mkdir -p "$LOGDIR"

need_samtools
index_reference "$REF"

log() { printf '%s  %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"; }

# Record why this BAM was not converted, and stop.
fail() {   # exit_code reason detail
    write_status "$LOGDIR" failed "$2" "$3"
    log "FAILED ($2): $3"
    echo "cramit: $NAME: $2 — $3 (log: $LOG)" >&2
    exit "$1"
}

first_error() { grep -m1 -E '\[E::|rror|fail' "$1" 2>/dev/null || tail -n1 "$1" 2>/dev/null || true; }

# ── Already done? (checked before anything else, so it holds after the BAM is deleted)
if [ -f "$CRAM" ] && [ -f "$LOGDIR/status" ] \
   && [ "$(cut -f1 "$LOGDIR/status")" = ok ] && [ "$(cut -f4 "$LOGDIR/status")" = "$(stat -c%s "$CRAM")" ]; then
    echo "already done: $CRAM (validated earlier, unchanged since)"
    exit 0
fi

write_status "$LOGDIR" running
START=$(date +%s)
TMP="${CRAM%.cram}.cramit-tmp.$$.cram"
FIFO=""
cleanup() { [ -z "$FIFO" ] || rm -rf "$FIFO"; rm -f "$TMP" "$TMP.crai" "$CRAM.crai.cramit-tmp.$$"; }
on_exit() {
    local rc=$?
    cleanup
    # an unexpected error (not a signal, not fail()) must not leave "running" behind
    if [ "$rc" -ne 0 ] && [ "$rc" -lt 128 ] && [ "$(cut -f1 "$LOGDIR/status" 2>/dev/null)" = running ]; then
        write_status "$LOGDIR" failed internal_error "exit $rc — see $LOG"
    fi
}
trap on_exit EXIT
trap 'exit 143' TERM
trap 'exit 130' INT

{ echo ""; echo "=== cramit convert  $(date '+%F %T')"; } >> "$LOG"
log "BAM:  $BAM"
log "CRAM: $CRAM"
log "REF:  $REF"
log "$(samtools --version | head -1), $THREADS threads, ${CRAMIT_WHERE:-local:$(hostname):$$}"

# ── Inspect the BAM header ────────────────────────────────────
[ -e "$BAM" ] || fail 1 bam_missing "BAM not found: $BAM"
HDR=$(samtools view -H "$BAM" 2> "$LOGDIR/header.err") \
    || fail 1 bam_corrupt "cannot read the BAM header: $(first_error "$LOGDIR/header.err")"
SO=$(awk -F'\t' '$1=="@HD"{for(i=2;i<=NF;i++) if($i~/^SO:/) print substr($i,4)}' <<< "$HDR")
NSQ=$(awk '$1=="@SQ"{n++} END{print n+0}' <<< "$HDR")
if   [ "$NSQ" -eq 0 ];                       then INDEX=no;    WHY="unaligned BAM"
elif [ "$SO" = coordinate ];                 then INDEX=yes;   WHY="coordinate-sorted"
elif [[ "$SO" =~ ^(queryname|unsorted)$ ]];  then INDEX=no;    WHY="$SO-sorted BAM"
else                                              INDEX=maybe; WHY="sort order not stated"
fi
log "index: $INDEX ($WHY)"
QC_OPT=""; [ "$NSQ" -gt 0 ] || QC_OPT="-u"    # quickcheck needs -u for files without @SQ

# Every BAM contig must be in the reference at the same length. Otherwise
# samtools does not fail: it quietly stores reference sequence in the CRAM.
if [ "$NSQ" -gt 0 ]; then
    BAD=$(awk -F'\t' -v fai="$REF.fai" '
        BEGIN { while ((getline l < fai) > 0) { split(l, a, "\t"); len[a[1]] = a[2] } }
        $1 == "@SQ" { sn = ""; ln = ""
            for (i = 2; i <= NF; i++) { if ($i ~ /^SN:/) sn = substr($i, 4); else if ($i ~ /^LN:/) ln = substr($i, 4) }
            if (!(sn in len)) print sn " not in reference"
            else if (len[sn] != ln) print sn " has length " ln " in the BAM, " len[sn] " in the reference" }' <<< "$HDR")
    [ -z "$BAD" ] || fail 1 wrong_reference \
        "$(wc -l <<< "$BAD") contig(s) do not match, e.g. $(head -3 <<< "$BAD" | paste -sd ';' -)"
fi

# ── One decode pass over a file → stats + flagstat (+ optionally a CRAM) ──
# Sets RC_READ RC_PIPE RC_WRITE RC_STATS RC_FLAGSTAT.
FIFO=$(mktemp -d "${TMPDIR:-/tmp}/cramit.XXXXXX")
one_pass() {   # tag input [cram_out]
    local tag="$1" in="$2" out="${3:-}" ps pf p
    rm -f "$FIFO/s" "$FIFO/f"; mkfifo "$FIFO/s" "$FIFO/f"
    samtools stats - < "$FIFO/s" > "$LOGDIR/$tag.stats" 2> "$LOGDIR/$tag.stats.err" & ps=$!
    samtools flagstat - < "$FIFO/f" > "$LOGDIR/$tag.flagstat" 2> "$LOGDIR/$tag.flagstat.err" & pf=$!
    set +e
    if [ -n "$out" ]; then
        # 1 decompression thread is plenty; CRAM encoding gets the rest
        samtools view -h -u -@ 1 "$in" 2> "$LOGDIR/$tag.read.err" \
            | tee "$FIFO/s" "$FIFO/f" \
            | samtools view -C -T "$REF" --output-fmt-option version=3.0 -@ "$THREADS" -o "$out" - 2> "$LOGDIR/cram.write.err"
        p=("${PIPESTATUS[@]}"); RC_WRITE=${p[2]}
    else
        samtools view -h -u -@ "$THREADS" -T "$REF" "$in" 2> "$LOGDIR/$tag.read.err" | tee "$FIFO/s" > "$FIFO/f"
        p=("${PIPESTATUS[@]}"); RC_WRITE=0
    fi
    RC_READ=${p[0]}; RC_PIPE=${p[1]}
    wait "$ps"; RC_STATS=$?
    wait "$pf"; RC_FLAGSTAT=$?
    set -e
}

# Did the BAM pass read every record? samtools can exit 0 after a read error.
check_bam_pass() {
    if [ "$RC_READ" -ne 0 ] || grep -qE 'Read block operation failed|error reading file' "$LOGDIR/bam.read.err"; then
        fail 1 bam_corrupt "the BAM could not be read to the end (truncated or damaged): $(first_error "$LOGDIR/bam.read.err")"
    fi
    if [ "$RC_STATS" -ne 0 ] || [ "$RC_FLAGSTAT" -ne 0 ]; then
        fail 1 validate_error "samtools stats/flagstat failed on the BAM: $(first_error "$LOGDIR/bam.stats.err")"
    fi
}

# Same reads? flagstat must be identical, and so must the samtools stats CHK
# checksums over read names, sequences and base qualities.
same_reads() {
    local a b
    a=$(grep '^CHK' "$LOGDIR/bam.stats" | cut -f2-); b=$(grep '^CHK' "$LOGDIR/cram.stats" | cut -f2-)
    if [ -z "$a" ] || [ -z "$b" ]; then DIFF="no CHK line in samtools stats output"; return 1; fi
    if [ "$a" != "$b" ]; then DIFF="read checksums differ (BAM CHK $a, CRAM CHK $b)"; return 1; fi
    if ! cmp -s "$LOGDIR/bam.flagstat" "$LOGDIR/cram.flagstat"; then DIFF="flagstat differs (see $LOGDIR/*.flagstat)"; return 1; fi
}

# Build an index via a temp name, so a killed run never leaves a partial .crai.
# samtools index fails when the reads are not in position order.
make_index() {   # cram final_crai
    local t="$2.cramit-tmp.$$"
    if samtools index -@ "$THREADS" -o "$t" "$1" 2>> "$LOGDIR/index.err"; then mv -f "$t" "$2"; else rm -f "$t"; return 1; fi
}
NOT_SORTED="the header says coordinate-sorted, but the reads are not in position order (samtools index failed)"

# ── A CRAM is already there: check it, never overwrite it ─────
if [ -f "$CRAM" ]; then
    log "a CRAM already exists — checking it against the BAM (it will not be modified)"
    one_pass bam "$BAM"
    check_bam_pass
    samtools quickcheck $QC_OPT "$CRAM" 2>/dev/null \
        || fail 2 existing_cram_differs "the existing CRAM is incomplete (fails samtools quickcheck): $CRAM"
    one_pass cram "$CRAM"
    if [ "$RC_READ" -ne 0 ] || [ "$RC_STATS" -ne 0 ] || [ "$RC_FLAGSTAT" -ne 0 ]; then
        fail 2 existing_cram_differs "the existing CRAM cannot be decoded with this reference: $(first_error "$LOGDIR/cram.read.err")"
    fi
    same_reads || fail 2 existing_cram_differs "the existing CRAM does not match the BAM: $DIFF"
    if [ "$INDEX" != no ] && [ ! -s "$CRAM.crai" ]; then
        if make_index "$CRAM" "$CRAM.crai"; then log "indexed the existing CRAM"
        elif [ "$INDEX" = yes ]; then fail 1 not_sorted "$NOT_SORTED"; fi
    fi
    write_status "$LOGDIR" ok - - "$(stat -c%s "$CRAM")"
    log "OK: the existing CRAM matches the BAM — kept as is"
    exit 0
fi

# ── Convert (single pass over the BAM) ────────────────────────
rm -f "${CRAM%.cram}".cramit-tmp.*        # left by an attempt that was killed
mkdir -p "$(dirname "$CRAM")"
log "converting (BAM read once: CRAM + stats + flagstat) → $(basename "$TMP")"
one_pass bam "$BAM" "$TMP"
check_bam_pass
if grep -q 'Enabling embed_ref' "$LOGDIR/cram.write.err"; then
    fail 1 wrong_reference "samtools had to embed reference sequence in the CRAM — the reference does not match the BAM"
fi
if [ "$RC_WRITE" -ne 0 ] || [ "$RC_PIPE" -ne 0 ]; then
    fail 1 convert_error "writing the CRAM failed: $(first_error "$LOGDIR/cram.write.err")"
fi
log "converted"

# ── Index, as far as the read order allows ────────────────────
case "$INDEX" in
    yes)
        make_index "$TMP" "$TMP.crai" || fail 1 not_sorted "$NOT_SORTED"
        log "indexed" ;;
    maybe)
        if make_index "$TMP" "$TMP.crai"; then log "reads are in position order — indexed"
        else log "reads are not in position order — no index (expected)"; fi ;;
    no)
        log "no index: a $WHY cannot be indexed (expected)" ;;
esac

# ── Validate (single pass over the CRAM) ──────────────────────
samtools quickcheck $QC_OPT "$TMP" 2>/dev/null || fail 1 validate_error "the new CRAM fails samtools quickcheck"
log "validating (CRAM decoded once: stats + flagstat)"
one_pass cram "$TMP"
if [ "$RC_READ" -ne 0 ] || [ "$RC_STATS" -ne 0 ] || [ "$RC_FLAGSTAT" -ne 0 ]; then
    fail 1 validate_error "decoding the new CRAM failed: $(first_error "$LOGDIR/cram.read.err")"
fi
same_reads || fail 2 cram_differs "$DIFF — the CRAM was discarded, the BAM is untouched"

# ── Publish: index first, CRAM last (the CRAM's presence marks completion) ──
if [ -f "$TMP.crai" ]; then mv -f "$TMP.crai" "$CRAM.crai"; else rm -f "$CRAM.crai"; fi
mv -f "$TMP" "$CRAM"
BAM_B=$(stat -c%s "$BAM"); CRAM_B=$(stat -c%s "$CRAM")
write_status "$LOGDIR" ok - - "$CRAM_B"
SECS=$(( $(date +%s) - START ))
log "OK: $(human_bytes "$BAM_B") → $(human_bytes "$CRAM_B") ($(awk -v a="$BAM_B" -v b="$CRAM_B" 'BEGIN{printf "%.0f", 100*(1-b/a)}')% smaller) in $((SECS/3600))h$(printf '%02d' $((SECS%3600/60)))m"
