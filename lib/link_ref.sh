#!/usr/bin/env bash
# cramit link-ref — put a link to the reference next to every validated CRAM.
#
# Usage: cramit link-ref [--name NAME] RUN
#   --name NAME   link name, without extension (default ref_for_bam_reconstruction)
#
# A CRAM can only be read with the reference it was made with. Each CRAM's
# header records that reference's path and checksums, but nothing in the
# folder shows it. This adds, in every folder that holds a validated CRAM:
#   NAME.fa.gz      -> the run's reference (extension as the reference's)
#   NAME.fa.gz.fai  -> its index, and NAME.fa.gz.gzi if it is bgzipped
# so anyone can rebuild a BAM with
#   samtools view -b -T NAME.fa.gz -o sample.bam sample.cram
#
# Links that already point to this reference are kept. A NAME that exists and
# is anything else is reported and left alone. Safe to run again at any time.
set -euo pipefail
source "$CRAMIT_HOME/lib/common.sh"

NAME=ref_for_bam_reconstruction; ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --name) NAME="${2:?--name needs a value}"; shift 2 ;;
        -h|--help) show_usage "$0" ;;
        -*) die "unknown option: $1" ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} -eq 1 ] || show_usage "$0" 1
[[ "$NAME" =~ ^[A-Za-z0-9._-]+$ ]] || die "--name must be a plain file name"
load_run "${ARGS[0]}"
if [ ! -f "$REF.fai" ] || { [[ "$REF" == *.gz ]] && [ ! -f "$REF.gzi" ]; }; then
    need_samtools; index_reference "$REF"
fi

# The link keeps the reference's extension: hs37d5.fa.gz -> NAME.fa.gz
base=$(basename "$REF"); gz=""
case "$base" in *.gz) gz=.gz; base=${base%.gz} ;; esac
case "$base" in *.*) ext=".${base##*.}$gz" ;; *) ext="$gz" ;; esac
LINK="$NAME$ext"
names=("$LINK" "$LINK.fai"); targets=("$REF" "$REF.fai")
if [ -n "$gz" ]; then names+=("$LINK.gzi"); targets+=("$REF.gzi"); fi

linked=0; already=0; taken=0; failed=0; problems=""
while IFS= read -r dir; do
    made=0; bad=""
    for i in "${!names[@]}"; do
        p="$dir/${names[$i]}"; t="${targets[$i]}"
        if [ -L "$p" ] && [ "$(readlink "$p")" = "$t" ]; then
            continue
        elif [ -e "$p" ] || [ -L "$p" ]; then
            bad="name taken: $p"; break
        elif ln -s "$t" "$p" 2>/dev/null; then
            made=1
        else
            bad="cannot write: $p"; break
        fi
    done
    if [ -n "$bad" ]; then
        case "$bad" in name*) taken=$((taken + 1)) ;; *) failed=$((failed + 1)) ;; esac
        problems+="  $bad"$'\n'
    elif [ "$made" -eq 1 ]; then linked=$((linked + 1))
    else already=$((already + 1)); fi
done < <(sample_states | awk -F'\t' '$2 == "done" {print $6}' \
            | while IFS= read -r c; do [ -f "$c" ] && dirname "$c"; done | sort -u)

total=$((linked + already + taken + failed))
echo "cramit link-ref: $RUN"
printf '  %-22s %s\n' "reference" "$REF" "link" "$LINK (+ index)"
printf '  %-22s %d\n' "folders with CRAMs" "$total" "  linked now" "$linked" \
    "  already linked" "$already" "  kept, name taken" "$taken" "  could not write" "$failed"
[ -z "$problems" ] || { echo ""; printf '%s' "$problems"; }
[ "$total" -eq 0 ] || { echo ""; echo "rebuild a BAM:  samtools view -b -T $LINK -o sample.bam sample.cram   (in the CRAM's folder)"; }
[ $((taken + failed)) -eq 0 ]
