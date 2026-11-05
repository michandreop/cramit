#!/usr/bin/env bash
# cramit scan — find the BAMs, check each one, and set up a run directory.
#
# Usage: cramit scan -r REF -o RUN [-j JOBS] BAM_DIR | BAM_LIST
#   -r REF   reference FASTA the BAMs were aligned to (one reference per run)
#   -o RUN   run directory to create: lists, logs and progress live here
#   -j JOBS  BAMs checked in parallel (default 8)
#
# BAM_DIR is searched recursively for *.bam (symlinks are skipped). BAM_LIST
# is a text file with one BAM path per line.
#
# Each BAM gets a quick check (header + end of file only, no full read):
#   OK               will be converted (listed in RUN/samples.tsv)
#   MISSING          path does not exist
#   CORRUPT          header unreadable, or no end-of-file marker (truncated,
#                    or still being written)
#   WRONG_REFERENCE  a contig is missing from REF or has another length or checksum
#   NOT_WRITABLE     the CRAM cannot be written next to the BAM
# Everything but OK goes to RUN/attention.txt and is left alone.
# Re-running scan on the same RUN refreshes the lists and keeps all progress.
set -euo pipefail
source "$CRAMIT_HOME/lib/common.sh"

REF=""; RUN_ARG=""; JOBS=8
while getopts "r:o:j:h" opt; do
    case $opt in
        r) REF="$OPTARG" ;; o) RUN_ARG="$OPTARG" ;; j) JOBS="$OPTARG" ;;
        h) show_usage "$0" ;; *) show_usage "$0" 1 ;;
    esac
done
shift $((OPTIND - 1))
[ $# -eq 1 ] || show_usage "$0" 1
[ -n "$REF" ] && [ -n "$RUN_ARG" ] || die "-r REF and -o RUN are required"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || die "-j must be a positive integer"
IN="$1"; REF=$(abspath "$REF"); RUN=$(abspath "${RUN_ARG%/}")

need_samtools
index_reference "$REF"

if [ -f "$RUN/cramit.conf" ]; then
    OLD_REF=$(source "$RUN/cramit.conf"; echo "$REF")
    [ "$OLD_REF" = "$REF" ] || die "$RUN was set up with reference $OLD_REF — use a new run directory for another reference"
fi
mkdir -p "$RUN/samples" "$RUN/slurm"

# ── Collect BAM paths (absolute, unique) ──────────────────────
LINKS=0
if [ -d "$IN" ]; then
    IN=$(cd "$IN" && pwd)
    find "$IN" -type f -name '*.bam' | sort > "$RUN/.found"
    LINKS=$(find "$IN" -type l -name '*.bam' | wc -l)
    WHERE="in $IN"
elif [ -f "$IN" ]; then
    list_lines "$IN" | while IFS= read -r p; do abspath "$p"; done | sort -u > "$RUN/.found"
    WHERE="listed in $IN"
else
    die "not a directory or a file: $IN"
fi
N=$(wc -l < "$RUN/.found")
[ "$N" -gt 0 ] || { rm -f "$RUN/.found"; die "no BAMs $WHERE"; }

# Reference checksums (M5), to catch a reference with the right contig names
# and lengths but different sequence. Uses an existing .dict, else builds one.
DICT=""
for d in "$REF.dict" "${REF%.gz}.dict" "${REF%.fa*}.dict" "${REF%.fasta*}.dict" "$RUN/reference.dict"; do
    [ -f "$d" ] && { DICT="$d"; break; }
done
if [ -z "$DICT" ]; then
    echo "computing reference checksums (once per run)..."
    samtools dict -o "$RUN/reference.dict.tmp" "$REF" && mv -f "$RUN/reference.dict.tmp" "$RUN/reference.dict"
    DICT="$RUN/reference.dict"
fi

# ── Check one BAM → status reason sort_order index size cram_exists key bam ──
check_bam() {
    local bam="$1" fai="$2" dict="$3" hdr so nsq index="-" size=0 refm cram status=OK reason="-" key
    key=$(sample_key "$bam"); cram="${bam%.bam}.cram"
    out() { printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$status" "$reason" "${so:--}" "$index" "$size" \
                "$([ -e "$cram" ] && echo yes || echo no)" "$key" "$bam"; }
    if [ ! -e "$bam" ]; then status=MISSING; reason="path does not exist"; out; return; fi
    size=$(stat -L -c%s "$bam" 2>/dev/null || echo 0)
    if ! hdr=$(samtools view -H "$bam" 2>/dev/null); then status=CORRUPT; reason="header unreadable"; out; return; fi
    so=$(awk -F'\t' '$1=="@HD"{for(i=2;i<=NF;i++) if($i~/^SO:/) print substr($i,4)}' <<< "$hdr")
    nsq=$(awk '$1=="@SQ"{n++} END{print n+0}' <<< "$hdr")
    if [ "$nsq" -eq 0 ]; then index=no; so="${so:-none} (unaligned)"
    else case "$so" in coordinate) index=yes ;; queryname|unsorted) index=no ;; *) index=maybe ;; esac
    fi
    if ! samtools quickcheck $([ "$nsq" -eq 0 ] && echo -u) "$bam" 2>/dev/null; then
        status=CORRUPT; reason="no end-of-file marker: truncated, or still being written"; out; return
    fi
    if [ "$nsq" -gt 0 ]; then
        refm=$(awk -F'\t' -v fai="$fai" -v dict="$dict" '
            BEGIN {
                while ((getline l < fai) > 0) { split(l, a, "\t"); len[a[1]] = a[2] }
                while ((getline l < dict) > 0) if (l ~ /^@SQ/) {
                    n = split(l, a, "\t"); sn = ""; m = ""
                    for (i = 2; i <= n; i++) { if (a[i] ~ /^SN:/) sn = substr(a[i], 4); else if (a[i] ~ /^M5:/) m = substr(a[i], 4) }
                    if (sn != "") md5[sn] = tolower(m)
                }
            }
            function add(list, n, item) { return n <= 3 ? list (n > 1 ? ", " : "") item : list }
            function more(n) { return n > 3 ? ", +" (n - 3) " more" : "" }
            $1 == "@SQ" { sn = ""; ln = ""; m = ""; as = ""; ur = ""
                for (i = 2; i <= NF; i++) {
                    if ($i ~ /^SN:/) sn = substr($i, 4); else if ($i ~ /^LN:/) ln = substr($i, 4)
                    else if ($i ~ /^M5:/) m = tolower(substr($i, 4))
                    else if ($i ~ /^AS:/) as = $i; else if ($i ~ /^UR:/) ur = $i }
                if (names == "") names = (as ur == "") ? "-" : as (as != "" && ur != "" ? " " : "") ur
                if (!(sn in len)) mlist = add(mlist, ++miss, sn)
                else if (len[sn] != ln) llist = add(llist, ++lenbad, sn " is " ln " in the BAM, " len[sn] " in the reference")
                else if (m != "" && md5[sn] != "" && md5[sn] != m) clist = add(clist, ++m5bad, sn) }
            END {
                if (miss)        printf "%d contig(s) not in the reference (%s%s)", miss, mlist, more(miss)
                else if (lenbad) printf "%d contig(s) have another length (%s%s)", lenbad, llist, more(lenbad)
                else if (m5bad)  printf "%d contig(s) have another sequence checksum, M5 (%s%s)", m5bad, clist, more(m5bad)
                if ((miss || lenbad || m5bad) && names != "-") printf "; the BAM names its reference %s", names
            }' <<< "$hdr")
        if [ -n "$refm" ]; then status=WRONG_REFERENCE; reason="$refm"; out; return; fi
    fi
    if [ ! -w "$(dirname "$bam")" ]; then status=NOT_WRITABLE; reason="no write permission in $(dirname "$bam")"; fi
    out
}
export -f check_bam sample_key

echo "checking $N BAMs $WHERE ($JOBS at a time)..."
xargs -d '\n' -P "$JOBS" -I{} bash -c 'check_bam "$1" "$2" "$3"' _ {} "$REF.fai" "$DICT" \
    < "$RUN/.found" | sort -t$'\t' -k8,8 > "$RUN/.checked"
rm -f "$RUN/.found"

{ printf 'status\treason\tsort_order\tindex\tsize_bytes\tcram_exists\tkey\tbam\n'; cat "$RUN/.checked"; } > "$RUN/preflight.tsv"
{ printf 'key\tbam\tcram\n'
  awk -F'\t' -v OFS='\t' '$1=="OK"{c=$8; sub(/\.bam$/, ".cram", c); print $7, $8, c}' "$RUN/.checked"; } > "$RUN/samples.tsv"
awk -F'\t' -v OFS='\t' '$1!="OK"{print $1, $2, $8}' "$RUN/.checked" > "$RUN/attention.txt"
rm -f "$RUN/.checked"
printf 'REF=%q\n' "$REF" > "$RUN/cramit.conf"

# ── Summary ───────────────────────────────────────────────────
NOK=$(( $(wc -l < "$RUN/samples.tsv") - 1 ))
{
    echo ""
    echo "cramit scan  $(date '+%F %T')"
    echo "  BAMs found $WHERE: $N ($(human_bytes "$(awk -F'\t' 'NR>1{s+=$5} END{printf "%d", s}' "$RUN/preflight.tsv")"))"
    echo "  reference: $REF"
    echo ""
    echo "  ready to convert: $NOK    → $RUN/samples.tsv"
    echo "  set aside:        $(wc -l < "$RUN/attention.txt")    → $RUN/attention.txt"
    awk -F'\t' '$1!="OK"{n[$1]++} END{for (k in n) printf "      %-16s %d\n", k, n[k]}' "$RUN/attention.txt"
    noidx=$(awk -F'\t' 'NR>1 && $1=="OK" && $4=="no"' "$RUN/preflight.tsv" | wc -l)
    maybe=$(awk -F'\t' 'NR>1 && $1=="OK" && $4=="maybe"' "$RUN/preflight.tsv" | wc -l)
    have=$(awk -F'\t' 'NR>1 && $1=="OK" && $6=="yes"' "$RUN/preflight.tsv" | wc -l)
    [ "$noidx" -eq 0 ] || echo "  note: $noidx BAM(s) are queryname-sorted, unsorted or unaligned — their CRAMs will have no index"
    [ "$maybe" -eq 0 ] || echo "  note: $maybe BAM(s) do not state a sort order — indexed only if the reads turn out to be sorted"
    [ "$have" -eq 0 ]  || echo "  note: $have BAM(s) already have a .cram next to them — it will be checked against the BAM, never overwritten"
    [ "$LINKS" -eq 0 ] || echo "  note: $LINKS symlinked .bam path(s) skipped — scan the folder they point to instead"
    echo ""
    echo "Next:  cramit run $RUN        (this machine)"
    echo "   or  cramit submit $RUN     (SLURM)"
} | tee "$RUN/scan_summary.txt"
