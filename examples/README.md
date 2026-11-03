# Worked example

One cramit run on six tiny BAMs, from start to finish: what each command
prints, which files appear where, and what to do when a BAM has a problem.
It runs on any Linux machine in a few seconds, no SLURM needed.

```bash
bash cramit/examples/run_demo.sh ~/cramit-demo    # builds the demo data, runs every step below
rm -rf ~/cramit-demo                              # when you are done
```

The command output below is real output from that script, and the file trees
have notes added. Paths are shortened; your times, sizes and IDs will differ.

**Contents:**
[The demo data](#the-demo-data) ·
[1. Scan](#1-scan) ·
[2. Preview a SLURM submit](#2-preview-a-slurm-submit-optional) ·
[3. Convert and validate](#3-convert-and-validate) ·
[4. See what needs attention](#4-see-what-needs-attention) ·
[5. Fix and retry](#5-fix-and-retry) ·
[6. Delete the BAMs](#6-delete-the-bams) ·
[7. Read a CRAM](#7-read-a-cram)

## The demo data

A reference with two 5 kb contigs, and six BAMs of 400 reads each:

```
bams
├── batch2
│   └── sampleC.bam
├── sampleA.bam
├── sampleA.bam.bai
├── sampleB.bam
├── sampleD.bam
├── sampleE.bam
└── sampleF.bam
ref
└── ref.fa
```

| BAM | What it is | What happens |
|---|---|---|
| `sampleA` | coordinate-sorted, with a `.bai` | converted; the CRAM gets a `.crai` index |
| `sampleB` | sorted by read name | converted; no index (only position-sorted CRAMs can be indexed) |
| `batch2/sampleC` | unaligned reads, in a subfolder | converted; no index |
| `sampleD` | aligned to a reference that calls the contigs `1`, `2` instead of `chr1`, `chr2` | set aside by scan: `WRONG_REFERENCE` |
| `sampleE` | truncated: the file was cut in half | set aside by scan: `CORRUPT` |
| `sampleF` | header says coordinate-sorted, but the reads are not | fails during conversion: `not_sorted`; fixed and retried in step 5 |

## 1. Scan

```
$ cramit scan -r ref/ref.fa -o run bams
computing reference checksums (once per run)...
checking 6 BAMs in ~/cramit-demo/bams (8 at a time)...

cramit scan  2026-11-03 13:36:24
  BAMs found in ~/cramit-demo/bams: 6 (60.9 KiB)
  reference: ~/cramit-demo/ref/ref.fa

  ready to convert: 4    → ~/cramit-demo/run/samples.tsv
  set aside:        2    → ~/cramit-demo/run/attention.txt
      WRONG_REFERENCE  1
      CORRUPT          1
  note: 2 BAM(s) are queryname-sorted, unsorted or unaligned — their CRAMs will have no index

Next:  cramit run ~/cramit-demo/run        (this machine)
   or  cramit submit ~/cramit-demo/run     (SLURM)
```

Scan reads only each BAM's header and the end of the file, so it is quick
even on large BAMs. It searched `bams/` including subfolders, and created the
run folder `run/`:

```
run
├── attention.txt       BAMs set aside, and why
├── cramit.conf         the reference for this run
├── preflight.tsv       scan's verdict on every BAM found
├── reference.dict      contig checksums of the reference (made once)
├── samples             one folder per BAM, filled during conversion
├── samples.tsv         the BAMs to convert, and where each CRAM will go
├── scan_summary.txt    the summary printed above
└── slurm               job lists and job logs, if you use submit
```

The two BAMs set aside are left alone:

```
$ cat run/attention.txt
WRONG_REFERENCE	2 contig(s) not in the reference	~/cramit-demo/bams/sampleD.bam
CORRUPT	no end-of-file marker: truncated, or still being written	~/cramit-demo/bams/sampleE.bam
```

In a real run, you would scan `sampleD` into its own run with the reference
it was aligned to, and copy `sampleE` again from its source. See
[SPECIAL_CASES.md](../docs/SPECIAL_CASES.md).

## 2. Preview a SLURM submit (optional)

On a cluster, the next step would be `cramit submit run`. The dry run shows
what would be submitted, and works without SLURM:

```
$ cramit submit -n run
cramit submit — dry run, nothing is submitted
  BAMs            4 to convert, one array task each
  each BAM gets   4 CPUs, 8G memory, 24:00:00
  at once         at most 50 BAMs
  sbatch options  none, so SLURM uses your default account and partition
  samtools        samtools 1.21 from PATH (/usr/bin/samtools)
  job logs        ~/cramit-demo/run/slurm/<job>_<task>.out and .err

sbatch command:
 sbatch --job-name=cramit_run_364cb2 --array=1-4%50 --cpus-per-task=4 --mem=8G --time=24:00:00 ...

Change resources with -c -m -T -p, add sbatch options after -- (see: cramit submit -h).
```

Here we convert on this machine instead.

## 3. Convert and validate

```
$ cramit run run
converting 4 BAM(s), 1 at a time, 4 thread(s) each
13:36:25  done                     ~/cramit-demo/bams/batch2/sampleC.bam
13:36:25  done                     ~/cramit-demo/bams/sampleA.bam
13:36:26  done                     ~/cramit-demo/bams/sampleB.bam
13:36:26  failed: not_sorted       ~/cramit-demo/bams/sampleF.bam

cramit status: ~/cramit-demo/run
  BAMs in run           4
  done                  3
  failed                1   needs a fix — cramit status -v
      not_sorted             1
  set aside             2   by scan — see ~/cramit-demo/run/attention.txt

3 BAM(s) can be deleted now: cramit delete-bams ~/cramit-demo/run
```

Each CRAM is now next to its BAM. Nothing has been deleted yet:

```
bams
├── batch2
│   ├── sampleC.bam
│   └── sampleC.cram
├── sampleA.bam
├── sampleA.bam.bai
├── sampleA.cram
├── sampleA.cram.crai      only sampleA is coordinate-sorted, so only it has an index
├── sampleB.bam
├── sampleB.cram
├── sampleD.bam            set aside by scan
├── sampleE.bam            set aside by scan
└── sampleF.bam            failed, so no CRAM was kept
```

Each BAM also got a folder in `run/samples/`, named after the BAM plus a
short ID:

```
run/samples/sampleA__38b5a02a
├── status                       one line: done or failed, the reason, CRAM size, where and when it ran
├── convert.log                  every step, with times
├── bam.stats    cram.stats      samtools stats of the BAM and the CRAM: their checksums must match
├── bam.flagstat cram.flagstat   samtools flagstat of both: must be identical
└── *.err                        samtools messages for each step (empty when all went well)
```

```
$ cat run/samples/sampleA__38b5a02a/convert.log

=== cramit convert  2026-11-03 13:36:25
2026-11-03 13:36:25  BAM:  ~/cramit-demo/bams/sampleA.bam
2026-11-03 13:36:25  CRAM: ~/cramit-demo/bams/sampleA.cram
2026-11-03 13:36:25  REF:  ~/cramit-demo/ref/ref.fa
2026-11-03 13:36:25  samtools 1.21, 4 threads, local:mylaptop:2647291
2026-11-03 13:36:25  index: yes (coordinate-sorted)
2026-11-03 13:36:25  converting (BAM read once: CRAM + stats + flagstat) → sampleA.cramit-tmp.2647291.cram
2026-11-03 13:36:25  converted
2026-11-03 13:36:25  indexed
2026-11-03 13:36:25  validating (CRAM decoded once: stats + flagstat)
2026-11-03 13:36:25  OK: 10.1 KiB → 2.6 KiB (74% smaller) in 0h00m
```

The CRAM is written under a temporary name and renamed to `sampleA.cram`
only after it has been validated. A file called `NAME.cram` is therefore
always a finished, validated CRAM.

## 4. See what needs attention

```
$ cramit status -v run
cramit status: ~/cramit-demo/run
  BAMs in run           4
  done                  3
  failed                1   needs a fix — cramit status -v
      not_sorted             1
  set aside             2   by scan — see ~/cramit-demo/run/attention.txt

3 BAM(s) can be deleted now: cramit delete-bams ~/cramit-demo/run

not_sorted  sampleF__4dc14263
    the header says coordinate-sorted, but the reads are not in position order (samtools index failed)
    → sort it (samtools sort) or fix the header's SO tag, then retry with --retry-failed
```

`-v` lists every BAM that needs you, what went wrong and what to do.

## 5. Fix and retry

```
$ samtools sort -o bams/sampleF.sorted.bam bams/sampleF.bam
$ mv bams/sampleF.sorted.bam bams/sampleF.bam

$ cramit run --retry-failed run
converting 1 BAM(s), 1 at a time, 4 thread(s) each
13:36:27  done                     ~/cramit-demo/bams/sampleF.bam

cramit status: ~/cramit-demo/run
  BAMs in run           4
  done                  4
  set aside             2   by scan — see ~/cramit-demo/run/attention.txt

4 BAM(s) can be deleted now: cramit delete-bams ~/cramit-demo/run
```

Only `sampleF` was converted again; finished BAMs are always skipped. Without
`--retry-failed`, failed BAMs are left alone, so a BAM that cannot convert
is not retried on every run.

## 6. Delete the BAMs

The dry run first: it checks everything again on disk and writes a log, but
deletes nothing.

```
$ cramit delete-bams run

DRY RUN — nothing deleted. Run again with --force to delete.
would delete:    4 BAM(s), 43.4 KiB
kept (a check failed): 0
not validated yet:     0
deleted earlier:       0
log: ~/cramit-demo/run/delete-bams.20261103_133629.2647904.log

$ cat run/delete-bams.20261103_133629.2647904.log
# cramit delete-bams  2026-11-03 13:36:29  mode: DRY-RUN
# action <TAB> reason <TAB> bam
WOULD_DELETE	-	~/cramit-demo/bams/batch2/sampleC.bam
WOULD_DELETE	-	~/cramit-demo/bams/sampleA.bam
WOULD_DELETE	-	~/cramit-demo/bams/sampleB.bam
WOULD_DELETE	-	~/cramit-demo/bams/sampleF.bam
```

A BAM that failed a check would be listed as `SKIP` with the reason. Then for
real:

```
$ cramit delete-bams --force run

deleted:         4 BAM(s), 43.4 KiB freed
kept (a check failed): 0
not validated yet:     0
deleted earlier:       0
log: ~/cramit-demo/run/delete-bams.20261103_133630.2647979.log
```

What is left:

```
bams
├── batch2
│   └── sampleC.cram
├── sampleA.cram
├── sampleA.cram.crai
├── sampleB.cram
├── sampleD.bam            never in the run, so never touched
├── sampleE.bam            never in the run, so never touched
├── sampleF.cram
└── sampleF.cram.crai
```

The `.bai` index went with its BAM. `sampleD` and `sampleE` are still there
for you to deal with.

## 7. Read a CRAM

A CRAM is read like a BAM, but needs the reference it was made with:

```
$ samtools view -T ref/ref.fa bams/sampleA.cram chr1:2000-2100 | head -3 | cut -f1-6
read0294	0	chr1	1910	60	100M
read0069	0	chr1	1914	60	100M
read0222	0	chr1	1939	60	100M
```

To get a BAM back: `samtools view -b -T ref/ref.fa -o sampleA.bam bams/sampleA.cram`.

---

Back to the [main README](../README.md) for running cramit on your own data.
