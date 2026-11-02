#!/usr/bin/env bash
# Regression tests: builds small fake BAMs for every case in
# docs/SPECIAL_CASES.md, runs the real commands on them and checks the outcome.
# SLURM is faked (sbatch/squeue/scontrol stubs; array tasks run in-process).
#
# Usage: bash tests/run_tests.sh [-k]      -k  keep the test data afterwards
# Needs samtools and python3 (to build the data). Takes under a minute.
set -uo pipefail

C="$(cd "$(dirname "$0")/.." && pwd)"
KEEP=0; [ "${1:-}" = -k ] && KEEP=1
D=$(mktemp -d "${TMPDIR:-/tmp}/cramit_tests.XXXXXX")
if [ -n "${CRAMIT_SETUP:-}" ]; then eval "$CRAMIT_SETUP" >/dev/null 2>&1; fi
command -v samtools >/dev/null || { echo "samtools not available (put it on PATH or set CRAMIT_SETUP)"; exit 2; }

pass=0; fail=0
check() { local name="$1"; shift; if "$@" >/dev/null 2>&1; then echo "PASS  $name"; pass=$((pass+1)); else echo "FAIL  $name"; fail=$((fail+1)); fi; }
not()    { ! "$@"; }
cramit() { "$C/cramit" "$@"; }
key()    { echo "$(basename "$1" .bam)__$(printf '%s' "$1" | md5sum | cut -c1-8)"; }
pf()     { awk -F'\t' -v b="$1" '$8==b{print $1}' "$2/preflight.tsv"; }
state()  { awk -F'\t' -v b="$1" '$5==b{print $2 ($2=="failed" ? ":" $3 : "")}' "$2/status.tsv"; }
logged() { grep -q "^$1	$2	$3\$" "$4"; }
S()      { samtools "$@" 2>/dev/null; }

B=$D/bams; R1=$D/run1; mkdir -p $B/p $B/q $B/ro $D/ref
echo "== building test data in $D"
python3 -I - "$D" <<'EOF'
import random, sys
d = sys.argv[1]; random.seed(7); L = 50000
seqs = {c: "".join(random.choice("ACGT") for _ in range(L)) for c in ("chr1", "chr2")}
with open(f"{d}/ref/ref.fa", "w") as f:
    for c, s in seqs.items():
        f.write(f">{c}\n" + "\n".join(s[i:i+60] for i in range(0, L, 60)) + "\n")
def reads(n, tag):
    out = []
    for i in range(n):
        c = random.choice(list(seqs)); p = random.randint(1, L - 100)
        out.append((f"{tag}{i:06d}", c, p, seqs[c][p-1:p+99]))
    return out
def sam(path, rs, so=None, aligned=True, names=("chr1", "chr2"), m5=""):
    with open(path, "w") as f:
        f.write("@HD\tVN:1.6" + (f"\tSO:{so}" if so else "") + "\n")
        if aligned:
            for n in names: f.write(f"@SQ\tSN:{n}\tLN:{L}" + (f"\tM5:{m5}" if m5 else "") + "\n")
        f.write("@RG\tID:rg1\tSM:test\n")
        ren = dict(zip(seqs, names))
        for name, c, p, s in rs:
            if aligned: f.write(f"{name}\t0\t{ren[c]}\t{p}\t60\t100M\t*\t0\t0\t{s}\t{'I'*100}\tRG:Z:rg1\n")
            else:       f.write(f"{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{s}\t{'I'*100}\tRG:Z:rg1\n")
base = reads(40000, "r")
sam(f"{d}/coord.sam", base)
sam(f"{d}/unsorted.sam", base, so="unsorted")
sam(f"{d}/noso.sam", base)
sam(f"{d}/ubam.sam", base, so="queryname", aligned=False)
sam(f"{d}/wrongref.sam", base[:2000], names=("1", "2"))
sam(f"{d}/m5bad.sam", base[:2000], m5="0" * 32)
sam(f"{d}/other.sam", reads(5000, "s"))
EOF
samtools faidx $D/ref/ref.fa
REF=$D/ref/ref.fa
S sort -o $B/coord.bam $D/coord.sam && S index $B/coord.bam
S sort -n -o $B/qname.bam $D/coord.sam
S view -b -o $B/unsorted.bam $D/unsorted.sam
S view -b -o $B/noso.bam $D/noso.sam
S view -h $B/coord.bam | sed 's/\tSO:coordinate//' | S view -b -o $B/noso_sorted.bam -
S view -b -o $B/ubam.bam $D/ubam.sam
S sort -o $B/wrongref.bam $D/wrongref.sam
S sort -o $B/m5bad.bam $D/m5bad.sam
S sort -o $D/other.bam $D/other.sam
S view -h $D/other.bam | sed 's/LN:50000/LN:50001/' | S view -b -o $B/lenbad.bam -
S view -h $B/noso.bam | sed 's/^@HD\tVN:1.6$/@HD\tVN:1.6\tSO:coordinate/' | S view -b -o $B/liar.bam -
head -c $(( $(stat -c%s $B/coord.bam) - 4000 )) $B/coord.bam > $B/trunc.bam
cp $B/coord.bam $B/p/twin.bam; cp $B/coord.bam $B/q/twin.bam
cp $B/coord.bam $B/have_good.bam; S view -C -T $REF -o $B/have_good.cram $B/coord.bam
cp $B/coord.bam $B/have_bad.bam;  S view -C -T $REF -o $B/have_bad.cram $D/other.bam
cp $B/coord.bam $B/ro/locked.bam; chmod a-w $B/ro
ln -s $B/coord.bam $B/link.bam
rm -f $D/*.sam
HAVE_BAD_MD5=$(md5sum < $B/have_bad.cram)

echo "== scan (folder)"
cramit scan -r $REF -o $R1 -j 4 $B > $D/scan1.out 2>&1
for b in coord qname unsorted noso noso_sorted ubam liar p/twin q/twin have_good have_bad; do
    check "scan: $b is OK" [ "$(pf $B/$b.bam $R1)" = OK ]; done
check "scan: truncated BAM is CORRUPT"                 [ "$(pf $B/trunc.bam $R1)" = CORRUPT ]
check "scan: contigs not in reference → WRONG_REFERENCE" [ "$(pf $B/wrongref.bam $R1)" = WRONG_REFERENCE ]
check "scan: contig length differs → WRONG_REFERENCE"  [ "$(pf $B/lenbad.bam $R1)" = WRONG_REFERENCE ]
check "scan: M5 checksum differs → WRONG_REFERENCE"    [ "$(pf $B/m5bad.bam $R1)" = WRONG_REFERENCE ]
check "scan: read-only folder → NOT_WRITABLE"          [ "$(pf $B/ro/locked.bam $R1)" = NOT_WRITABLE ]
check "scan: symlinked BAM skipped in folder mode"     not grep -q link.bam $R1/preflight.tsv
check "scan: summary mentions the skipped symlink"     grep -q "1 symlinked" $D/scan1.out
check "scan: samples.tsv lists the 11 OK BAMs"         [ "$(( $(wc -l < $R1/samples.tsv) - 1 ))" -eq 11 ]
check "scan: attention.txt lists the 5 others"         [ "$(wc -l < $R1/attention.txt)" -eq 5 ]
check "scan: refuses a second reference for the run"   not cramit scan -r $D/other.bam -o $R1 $B

echo "== run (this machine)"
cramit run -j 4 -t 1 $R1 > $D/run1.out 2>&1
cramit status $R1 > /dev/null
for b in coord qname unsorted noso noso_sorted ubam p/twin q/twin have_good; do
    check "run: $b done" [ "$(state $B/$b.bam $R1)" = done ]; done
check "run: header claims coordinate, reads unsorted → not_sorted" [ "$(state $B/liar.bam $R1)" = failed:not_sorted ]
check "run:   ...and no CRAM published"                 not [ -e $B/liar.cram ]
check "run: different CRAM already there → existing_cram_differs" [ "$(state $B/have_bad.bam $R1)" = failed:existing_cram_differs ]
check "run:   ...and that CRAM is left untouched"       [ "$(md5sum < $B/have_bad.cram)" = "$HAVE_BAD_MD5" ]
check "run: coordinate-sorted CRAM is indexed"          [ -s $B/coord.cram.crai ]
check "run: no stated order but sorted → indexed"       [ -s $B/noso_sorted.cram.crai ]
for b in qname unsorted noso ubam; do check "run: no index for $b" not [ -e $B/$b.cram.crai ]; done
check "run: matching CRAM already there is kept and indexed" [ -s $B/have_good.cram.crai ]
check "run: CRAM gets M5 tags although the BAM has none" [ "$(samtools view -H $B/coord.cram | grep -c $'\tM5:')" -eq 2 ]
check "run: same-name BAMs get separate CRAMs and logs" bash -c "[ -e $B/p/twin.cram ] && [ -e $B/q/twin.cram ] && [ -d $R1/samples/$(key $B/p/twin.bam) ] && [ -d $R1/samples/$(key $B/q/twin.bam) ]"
check "run: status records the validated CRAM size"     [ "$(cut -f4 $R1/samples/$(key $B/coord.bam)/status)" = "$(stat -c%s $B/coord.cram)" ]
check "run: stats + flagstat kept for BAM and CRAM"     bash -c "ls $R1/samples/$(key $B/coord.bam)/{bam,cram}.{stats,flagstat}"
check "run: no temp files left anywhere"                [ -z "$(find $B -name '*cramit-tmp*')" ]
cramit run $R1 > $D/run1b.out 2>&1
check "run again: nothing to do (failed ones wait for a fix)" grep -q "nothing to do" $D/run1b.out
S sort -o $B/liar.sorted.bam $B/liar.bam && mv $B/liar.sorted.bam $B/liar.bam
cramit run --retry-failed -t 1 $R1 > $D/run1c.out 2>&1
cramit status $R1 > $D/status1.out
check "retry after fixing the BAM: done"                [ "$(state $B/liar.bam $R1)" = done ]
check "retry: unfixed BAM still fails"                  [ "$(state $B/have_bad.bam $R1)" = failed:existing_cram_differs ]
cramit status -v $R1 > $D/status1v.out
check "status: counts done BAMs"                        grep -qE "done +10" $D/status1.out
check "status -v: explains the failure and the fix"     grep -q "move or delete it" $D/status1v.out
check "status: points to delete-bams"                   grep -q "cramit delete-bams" $D/status1.out

echo "== submit (SLURM faked)"
B2=$D/b2; R2=$D/run2; FB=$D/fakebin; mkdir -p $B2 $FB
for n in a b c d e; do cp $B/coord.bam $B2/$n.bam; done; S index $B2/e.bam
ln -s $B/qname.bam $B2/link.bam
printf '%s\n' $B2/{a,b,c,d,e,link}.bam $B2/missing.bam > $D/list2.txt
cat > $FB/sbatch <<EOF
#!/bin/bash
echo "\$*" >> $D/sbatch.calls; echo "Submitted batch job 900"
EOF
printf '#!/bin/bash\ncat %s/squeue.out 2>/dev/null\n' $D > $FB/squeue
printf '#!/bin/bash\necho "MaxArraySize            = 4"\n' > $FB/scontrol
chmod +x $FB/*; : > $D/squeue.out
export PATH=$FB:$PATH
cramit scan -r $REF -o $R2 $D/list2.txt > /dev/null 2>&1
check "scan list: missing path → MISSING"               [ "$(pf $B2/missing.bam $R2)" = MISSING ]
check "scan list: symlinked BAM accepted"               [ "$(pf $B2/link.bam $R2)" = OK ]
CRAMIT_SBATCH_ARGS="--partition=x" cramit submit -n $R2 -- --account=lab > $D/sub_dry.out
check "submit -n: array capped by MaxArraySize, 50 at once" grep -q -- "--array=1-3%50" $D/sub_dry.out
check "submit -n: passes CRAMIT_SBATCH_ARGS and -- options" bash -c "grep -q -- --partition=x $D/sub_dry.out && grep -q -- --account=lab $D/sub_dry.out"
check "submit -n: says the rest needs another submit"   grep -q "3 more BAM" $D/sub_dry.out
check "submit -n: submitted nothing"                    not [ -e $D/sbatch.calls ]
cramit submit -p 7 -m 9G $R2 > $D/sub1.out
check "submit: calls sbatch with the resources"         grep -q -- "--array=1-3%7 --cpus-per-task=4 --mem=9G" $D/sbatch.calls
L=$(ls -t $R2/slurm/batch.*.tsv | head -1)
for i in 1 2 3; do SLURM_ARRAY_JOB_ID=900 SLURM_ARRAY_TASK_ID=$i SLURM_CPUS_PER_TASK=1 \
    "$C/cramit" _task $R2 $L > $R2/slurm/900_$i.out 2> $R2/slurm/900_$i.err; done
cramit status $R2 > /dev/null
for n in a b c; do check "array task converted $n" [ "$(state $B2/$n.bam $R2)" = done ]; done
check "array task records where it ran"                 grep -q "slurm:900_1" $R2/samples/$(key $B2/a.bam)/status
echo 901 > $D/squeue.out
check "submit refuses while the run's array is active"  not cramit submit $R2
: > $D/squeue.out
cramit submit $R2 > /dev/null
L=$(ls -t $R2/slurm/batch.*.tsv | head -1)
check "second submit takes only the BAMs not done"      [ "$(wc -l < $L)" -eq 3 ]
# simulate three tasks of job 901 that never finished
fake() { mkdir -p $R2/samples/$(key $1); printf 'running\t-\t-\t-\tslurm:%s\t-\n' $2 > $R2/samples/$(key $1)/status; }
fake $B2/d.bam 901_1; echo "slurmstepd: error: *** JOB 901 ON n1 CANCELLED AT 2026 DUE TO TIME LIMIT ***" > $R2/slurm/901_1.err
fake $B2/e.bam 901_2; echo "slurmstepd: error: Detected 1 oom_kill event in StepId=901.batch" > $R2/slurm/901_2.err
fake $B2/link.bam 901_3; echo "901_3" > $D/squeue.out
cramit status $R2 > $D/status2.out
check "status: SLURM time limit → timeout"              [ "$(state $B2/d.bam $R2)" = timeout ]
check "status: SLURM memory kill → out_of_memory"       [ "$(state $B2/e.bam $R2)" = out_of_memory ]
check "status: task still in squeue → running"          [ "$(state $B2/link.bam $R2)" = running ]
: > $D/squeue.out; cramit status $R2 > /dev/null
check "status: task gone without a result → incomplete" [ "$(state $B2/link.bam $R2)" = incomplete ]
check "status: suggests more time for timeouts"         grep -q "submit -T" $D/status2.out
cramit run -t 1 $R2 > /dev/null 2>&1; cramit status $R2 > /dev/null
check "killed BAMs are picked up again and finish"      bash -c "[ \"\$(awk -F'\t' 'NR>1 && \$2==\"done\"' $R2/status.tsv | wc -l)\" -eq 6 ]"

echo "== delete-bams"
rm -f $B2/a.cram.crai                                                  # coordinate-sorted, index gone
S view -C -T $REF -o $B2/b.cram $D/other.bam                           # CRAM replaced after validation
sz=$(stat -c%s $B2/c.cram); dd if=/dev/zero of=$B2/c.cram bs=1 seek=$((sz-38)) count=38 conv=notrunc 2>/dev/null  # EOF damaged, same size
samtools view -H $B2/d.cram | sed 's/LN:50000/LN:50009/' > $D/d.hdr; samtools reheader -i $D/d.hdr $B2/d.cram 2>/dev/null
cramit delete-bams $R2 > $D/del_dry.out
G=$(ls -t $R2/delete-bams.*.log | head -1)
check "dry run: would delete the one safe BAM"          logged WOULD_DELETE - $B2/e.bam $G
check "keeps: coordinate-sorted BAM whose CRAM lost its index" logged SKIP crai_missing $B2/a.bam $G
check "keeps: CRAM changed since validation"            logged SKIP cram_changed_since_validation $B2/b.bam $G
check "keeps: CRAM fails quickcheck"                    logged SKIP cram_quickcheck_failed $B2/c.bam $G
check "keeps: CRAM contigs differ from the BAM"         grep -qE "^SKIP	(contigs_differ|cram_changed_since_validation)	$B2/d.bam\$" $G
check "keeps: symlinked BAM"                            logged SKIP bam_is_symlink $B2/link.bam $G
check "dry run deleted nothing"                         bash -c "[ -e $B2/e.bam ] && grep -q 'DRY RUN' $D/del_dry.out"
cramit delete-bams --force $R2 > /dev/null
check "--force deletes the safe BAM and its .bai"       bash -c "[ ! -e $B2/e.bam ] && [ ! -e $B2/e.bam.bai ] && [ -e $B2/e.cram ]"
check "--force keeps every BAM that failed a check"     bash -c "[ -e $B2/a.bam ] && [ -e $B2/b.bam ] && [ -e $B2/c.bam ] && [ -e $B2/d.bam ]"
check "--force leaves the symlink target alone"         [ -e $B/qname.bam ]
cramit delete-bams --force $R1 > /dev/null
check "run1: validated BAMs deleted, CRAMs kept"        bash -c "[ ! -e $B/coord.bam ] && [ ! -e $B/coord.bam.bai ] && [ -e $B/coord.cram ] && [ ! -e $B/ubam.bam ]"
check "run1: unvalidated BAM kept"                      [ -e $B/have_bad.bam ]
check "run1: BAMs set aside by scan kept"               bash -c "[ -e $B/trunc.bam ] && [ -e $B/ro/locked.bam ]"
cramit run --retry-failed -t 1 $R1 > /dev/null 2>&1; cramit status $R1 > $D/status1d.out
check "after deletion, done BAMs stay done"             [ "$(state $B/coord.bam $R1)" = done ]
check "status counts deleted BAMs"                      grep -q "already deleted" $D/status1d.out

echo "== convert (one BAM)"
X=$D/single; mkdir -p $X/cwd; cp $D/other.bam $X/one.bam; cp $D/other.bam $X/clash.bam; cp $B/coord.cram $X/clash.cram
(cd $X/cwd && "$C/cramit" convert -r $REF -o $X/out -t 1 $X/one.bam > /dev/null 2>&1); rc=$?
check "convert: exit 0, CRAM in -o folder"              bash -c "[ $rc -eq 0 ] && [ -s $X/out/one.cram ] && [ -s $X/out/one.cram.crai ]"
check "convert: default log dir ./cramit_logs/<key>"    [ -f "$X/cwd/cramit_logs/$(key $X/one.bam)/convert.log" ]
cramit convert -r $REF -l $X/l2 -t 1 $X/one.bam > /dev/null 2>&1 && cramit convert -r $REF -l $X/l2 -t 1 $X/one.bam > $X/again.out 2>&1
check "convert: second call says already done"          grep -q "already done" $X/again.out
cramit convert -r $REF -l $X/l3 $B/wrongref.bam > /dev/null 2>&1; rc=$?
check "convert: wrong reference → exit 1, no CRAM"      bash -c "[ $rc -eq 1 ] && [ ! -e $B/wrongref.cram ]"
cramit convert -r $REF -l $X/l4 $B/trunc.bam > /dev/null 2>&1; rc=$?
check "convert: truncated BAM → exit 1, bam_corrupt"    bash -c "[ $rc -eq 1 ] && grep -q bam_corrupt $X/l4/status && [ ! -e $B/trunc.cram ]"
cramit convert -r $REF -l $X/l5 $X/clash.bam > /dev/null 2>&1; rc=$?
check "convert: mismatching existing CRAM → exit 2"     [ $rc -eq 2 ]

echo "== deep-check (aux tags)"
samtools calmd -b $D/other.bam $REF 2>/dev/null | samtools view -h - \
    | awk 'BEGIN{OFS="\t"} /^@/{print;next} {print $0, "RX:Z:ACGTACGT", "MI:Z:"NR}' | S view -b -o $X/tags.bam -
S view -C -T $REF -o $X/tags.cram $X/tags.bam
cramit deep-check -r $REF $X/tags.bam $X/tags.cram > /dev/null; rc=$?
check "deep-check: tags survive (reordered, MD/NM regenerated)" [ $rc -eq 0 ]
samtools view -h $X/tags.bam | awk 'BEGIN{OFS="\t"} /^@/{print;next} NR==100{sub(/RX:Z:ACGTACGT/,"RX:Z:TTTTTTTT")} {print}' \
    | S view -C -T $REF -o $X/tags_bad.cram -
cramit deep-check -r $REF $X/tags.bam $X/tags_bad.cram > /dev/null; rc=$?
check "deep-check: a changed tag is caught"             [ $rc -eq 1 ]
mkdir -p $X/t; cp $X/tags.bam $X/t/tags.bam; cp $X/tags_bad.cram $X/t/tags.cram
check "  ...which the standard validation does not see" cramit convert -r $REF -l $X/l6 $X/t/tags.bam

echo "== command line"
check "cramit help"                                     cramit help
check "unknown command fails"                           not cramit frobnicate
check "convert without a BAM fails"                     not cramit convert -r $REF

echo ""
echo "passed: $pass   failed: $fail"
chmod -R u+w $D 2>/dev/null
if [ "$fail" -eq 0 ] && [ "$KEEP" -eq 0 ]; then rm -rf "$D"; else echo "test data: $D"; fi
[ "$fail" -eq 0 ]
