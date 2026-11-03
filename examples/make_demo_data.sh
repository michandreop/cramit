#!/usr/bin/env bash
# Builds the tiny demo data used in examples/README.md: a small reference and
# a folder of six small BAMs. Three convert cleanly; three have a problem that
# cramit catches.
#
# Usage: bash examples/make_demo_data.sh DIR      (needs samtools and awk)
set -euo pipefail
DIR="${1:?usage: make_demo_data.sh DIR}"
[ ! -e "$DIR" ] || { echo "$DIR already exists — give a new folder" >&2; exit 1; }
if [ -n "${CRAMIT_SETUP:-}" ]; then eval "$CRAMIT_SETUP" >/dev/null 2>&1 || true; fi
command -v samtools >/dev/null || { echo "samtools not found" >&2; exit 1; }

mkdir -p "$DIR/ref" "$DIR/bams/batch2" "$DIR/.build"
T="$DIR/.build"

# Reference: two random 5 kb contigs. Reads: 400 x 100 bp taken from them.
awk -v ref="$DIR/ref/ref.fa" -v sam="$T/reads.sam" 'BEGIN {
    srand(7); split("A C G T", b, " "); L = 5000
    for (c = 1; c <= 2; c++) {
        s = ""; for (i = 0; i < L; i++) s = s b[int(rand() * 4) + 1]; seq[c] = s
        print ">chr" c > ref; for (i = 1; i <= L; i += 60) print substr(s, i, 60) > ref
    }
    q = sprintf("%100s", ""); gsub(/ /, "I", q)
    print "@HD\tVN:1.6" > sam
    for (c = 1; c <= 2; c++) print "@SQ\tSN:chr" c "\tLN:" L > sam
    print "@RG\tID:rg1\tSM:demo" > sam
    for (i = 1; i <= 400; i++) {
        c = int(rand() * 2) + 1; p = int(rand() * (L - 100)) + 1
        printf "read%04d\t0\tchr%d\t%d\t60\t100M\t*\t0\t0\t%s\t%s\tRG:Z:rg1\n", i, c, p, substr(seq[c], p, 100), q > sam
    }
}'
S() { samtools "$@" 2>/dev/null; }

# Clean cases
S sort -o "$DIR/bams/sampleA.bam" "$T/reads.sam" && S index "$DIR/bams/sampleA.bam"   # coordinate-sorted, indexed
S sort -n -o "$DIR/bams/sampleB.bam" "$T/reads.sam"                                  # name-sorted
awk -F'\t' -v OFS='\t' '/^@SQ/{next} /^@HD/{print $0, "SO:queryname"; next} /^@/{print; next}
                        {$2=4; $3="*"; $4=0; $5=0; $6="*"; print}' "$T/reads.sam" \
    | S view -b -o "$DIR/bams/batch2/sampleC.bam" -                                  # unaligned reads, in a subfolder

# Problem cases
S view -h "$DIR/bams/sampleA.bam" | sed -e 's/SN:chr/SN:/' -e 's/\tchr\([12]\)\t/\t\1\t/' \
    | S view -b -o "$DIR/bams/sampleD.bam" -                                         # aligned to a "1, 2" reference
head -c $(( $(wc -c < "$DIR/bams/sampleA.bam") / 2 )) "$DIR/bams/sampleA.bam" > "$DIR/bams/sampleE.bam"   # truncated
S view -h "$DIR/bams/sampleB.bam" | sed -e 's/SO:queryname/SO:coordinate/' -e 's/\tSS:[^\t]*//' \
    | S view -b -o "$DIR/bams/sampleF.bam" -                                         # header says sorted, reads are not

rm -rf "$T"
