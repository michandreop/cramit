#!/usr/bin/env bash
# cramit delete-bams — delete the BAMs whose CRAM has been validated.
#
# Usage: cramit delete-bams [--force] RUN
#   (no option)  dry run: list what would be deleted, delete nothing
#   --force      delete
#
# Every check is redone on disk at the moment of deletion. A BAM (and its
# .bai) is deleted only if ALL of these hold:
#   1. cramit validated its CRAM, and the CRAM still has the size it had then
#   2. the BAM is a real file, not a symlink (deleting a link frees nothing)
#   3. if the BAM is coordinate-sorted, its CRAM has a non-empty .crai index
#   4. the CRAM passes samtools quickcheck
#   5. the CRAM lists the same contigs (names and lengths) as the BAM
# Each decision is written to RUN/delete-bams.<time>.<pid>.log.
set -uo pipefail
source "$CRAMIT_HOME/lib/common.sh"

FORCE=0; ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --force) FORCE=1; shift ;;
        -h|--help) show_usage "$0" ;;
        -*) die "unknown option: $1" ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} -eq 1 ] || show_usage "$0" 1
load_run "${ARGS[0]}"
need_samtools

LOG="$RUN/delete-bams.$(date +%Y%m%d_%H%M%S).$$.log"
{ echo "# cramit delete-bams  $(date '+%F %T')  mode: $([ $FORCE -eq 1 ] && echo DELETE || echo DRY-RUN)"
  echo "# action <TAB> reason <TAB> bam"; } > "$LOG" || die "cannot write $LOG"

contigs() { awk -F'\t' '$1=="@SQ"{sn="";ln=""; for(i=2;i<=NF;i++){if($i~/^SN:/)sn=$i; else if($i~/^LN:/)ln=$i} print sn"\t"ln}'; }

ndel=0; nskip=0; ngone=0; nopen=0; bytes=0
skip() { printf 'SKIP\t%s\t%s\n' "$1" "$bam" >> "$LOG"; nskip=$((nskip + 1)); }

while IFS=$'\t' read -r key bam cram; do
    [ "$key" = key ] && continue
    IFS=$'\t' read -r state _ _ okbytes _ < "$RUN/samples/$key/status" 2>/dev/null || state=""
    if [ "$state" != ok ]; then nopen=$((nopen + 1)); continue; fi
    if [ ! -e "$bam" ] && [ ! -L "$bam" ]; then ngone=$((ngone + 1)); continue; fi
    [ -L "$bam" ]  && { skip bam_is_symlink; continue; }
    [ -f "$cram" ] || { skip cram_missing; continue; }
    [ "$(stat -c%s "$cram")" = "$okbytes" ] || { skip cram_changed_since_validation; continue; }
    bam_hdr=$(samtools view -H "$bam" 2>/dev/null)   || { skip bam_unreadable; continue; }
    cram_hdr=$(samtools view -H "$cram" 2>/dev/null) || { skip cram_unreadable; continue; }
    so=$(awk -F'\t' '$1=="@HD"{for(i=2;i<=NF;i++) if($i~/^SO:/) print substr($i,4)}' <<< "$bam_hdr")
    if [ "$so" = coordinate ] && [ ! -s "$cram.crai" ]; then skip crai_missing; continue; fi
    if [ -e "$cram.crai" ] && [ ! -s "$cram.crai" ]; then skip crai_empty; continue; fi
    qc=""; [ -n "$(contigs <<< "$cram_hdr")" ] || qc="-u"
    samtools quickcheck $qc "$cram" 2>/dev/null || { skip cram_quickcheck_failed; continue; }
    [ "$(contigs <<< "$bam_hdr")" = "$(contigs <<< "$cram_hdr")" ] || { skip contigs_differ; continue; }

    size=$(stat -c%s "$bam"); bais=()
    for b in "$bam.bai" "${bam%.bam}.bai"; do [ -f "$b" ] && { bais+=("$b"); size=$((size + $(stat -c%s "$b"))); }; done
    if [ "$FORCE" -eq 1 ]; then
        rm -f "$bam" "${bais[@]}"
        [ -e "$bam" ] && { skip delete_failed; continue; }
        printf 'DELETED\t-\t%s\n' "$bam" >> "$LOG"
    else
        printf 'WOULD_DELETE\t-\t%s\n' "$bam" >> "$LOG"
    fi
    ndel=$((ndel + 1)); bytes=$((bytes + size))
done < "$RUN/samples.tsv"

echo ""
if [ "$FORCE" -eq 1 ]; then
    echo "deleted:         $ndel BAM(s), $(human_bytes "$bytes") freed"
else
    echo "DRY RUN — nothing deleted. Run again with --force to delete."
    echo "would delete:    $ndel BAM(s), $(human_bytes "$bytes")"
fi
echo "kept (a check failed): $nskip$([ $nskip -gt 0 ] && echo "   ← reasons in the log")"
echo "not validated yet:     $nopen"
echo "deleted earlier:       $ngone"
echo "log: $LOG"
