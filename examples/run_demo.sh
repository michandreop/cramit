#!/usr/bin/env bash
# Runs the whole cramit workflow on tiny demo data and prints every command
# with its output: build the data, scan, convert, check, fix a failed BAM and
# retry, delete the BAMs, read a CRAM back. examples/README.md walks through
# the same steps.
#
# Usage: bash examples/run_demo.sh [DIR]    (default ./cramit-demo; needs samtools)
# Takes a few seconds and needs no SLURM. Delete DIR when you are done.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="${1:-cramit-demo}"
if [ -n "${CRAMIT_SETUP:-}" ]; then eval "$CRAMIT_SETUP" >/dev/null 2>&1 || true; fi

cramit() { "$HERE/../cramit" "$@"; }
step()   { printf '\n========== %s\n' "$*"; }
show()   { printf '\n$ %s\n' "$*"; "$@"; }
tree_()  { printf '\n$ tree %s\n' "$*"; if command -v tree >/dev/null; then tree --noreport "$@"; else find "$@" | sort; fi; }

step "0. Build the demo data"
bash "$HERE/make_demo_data.sh" "$DIR"
cd "$DIR"
tree_ bams ref

step "1. Scan: check every BAM, set up the run folder"
show cramit scan -r ref/ref.fa -o run bams
tree_ run
show cat run/attention.txt

step "2. Preview a SLURM submit (dry run, works without SLURM)"
show cramit submit -n run

step "3. Convert and validate on this machine"
show cramit run run
tree_ bams
tree_ run/samples/sampleA__*
show cat run/samples/sampleA__*/convert.log

step "4. See what needs attention"
show cramit status -v run

step "5. Fix the failed BAM and retry it"
show samtools sort -o bams/sampleF.sorted.bam bams/sampleF.bam
show mv bams/sampleF.sorted.bam bams/sampleF.bam
show cramit run --retry-failed run

step "6. Delete the BAMs: dry run first, then for real"
show cramit delete-bams run
show cat run/delete-bams.*.log
show cramit delete-bams --force run
tree_ bams

step "7. Link the reference next to the CRAMs"
show cramit link-ref run
tree_ bams

step "8. Read a CRAM"
printf '\n$ samtools view -T bams/ref_for_bam_reconstruction.fa bams/sampleA.cram chr1:2000-2100 | head -3 | cut -f1-6\n'
samtools view -T bams/ref_for_bam_reconstruction.fa bams/sampleA.cram chr1:2000-2100 | head -3 | cut -f1-6

# The outcome the walkthrough describes; a failure here means the demo broke.
for f in bams/sampleA.cram bams/sampleA.cram.crai bams/sampleB.cram bams/batch2/sampleC.cram bams/sampleF.cram bams/sampleD.bam bams/sampleE.bam \
         bams/ref_for_bam_reconstruction.fa bams/batch2/ref_for_bam_reconstruction.fa.fai; do
    [ -s "$f" ] || { echo "demo FAILED: $f missing" >&2; exit 1; }
done
for f in bams/sampleA.bam bams/sampleB.bam bams/batch2/sampleC.bam bams/sampleF.bam; do
    [ ! -e "$f" ] || { echo "demo FAILED: $f was not deleted" >&2; exit 1; }
done
printf '\ndemo finished: 4 BAMs converted, validated and deleted; 2 set aside by scan. Files are in %s\n' "$PWD"
