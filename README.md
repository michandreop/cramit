<img src="docs/cramit-logo-beige.svg" alt="cramit" width="420">

[![tests](https://github.com/michandreop/cramit/actions/workflows/tests.yml/badge.svg)](https://github.com/michandreop/cramit/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![samtools](https://img.shields.io/badge/samtools-tested%201.21-brightgreen.svg)](https://www.htslib.org/)
[![bash](https://img.shields.io/badge/bash-%E2%89%A54.4-4EAA25.svg?logo=gnubash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Platform: Linux](https://img.shields.io/badge/platform-Linux-lightgrey.svg?logo=linux&logoColor=white)](#install)
[![SLURM: optional](https://img.shields.io/badge/SLURM-optional-orange.svg)](#convert-a-folder-of-bams)

Convert BAM files to CRAM, prove that every CRAM holds exactly the same reads
as its BAM, and only then delete the BAMs. Works on one BAM or on thousands,
on your own machine or on a SLURM cluster.

- **Safe:** a CRAM is kept only if it matches its BAM, and a BAM is deleted
  only after everything is checked again on disk.
- **Fast:** each BAM is read once. Conversion and validation share that one
  read.
- **Resumable:** stop at any point and run the same command again. Finished
  BAMs are skipped.

## Install

You need Linux, bash 4.4+ and [samtools](https://www.htslib.org/) (tested with
1.21) on your `PATH`.

```bash
git clone https://github.com/michandreop/cramit
cramit/tests/run_tests.sh          # optional self-test, under a minute
export PATH="$PWD/cramit:$PATH"    # or call it as /path/to/cramit/cramit
```

## Convert a folder of BAMs

You only need four commands. `RUN` is any folder you choose: cramit keeps its
lists, logs and progress there.

**1. Check the BAMs and set up the run**

```bash
cramit scan -r ref.fa.gz -o RUN /path/to/bams     # or a text file with one BAM path per line
```

Use the reference the BAMs were aligned to. Scan reads each BAM's header and
end of file, and sets aside any BAM that is truncated, aligned to a different
reference or in a folder you can't write to. The rest are listed as ready to
convert.

**2. Convert and validate**

```bash
cramit run RUN                      # on this machine (-j 4 converts 4 BAMs at a time)
cramit submit RUN                   # or as a SLURM job array
cramit submit RUN -- --account=mylab --partition=long    # extra sbatch options go after --
```

Each CRAM is written next to its BAM as `NAME.cram` (plus `NAME.cram.crai`
when the BAM is coordinate-sorted). It appears under that name only once it
has been validated.

**3. Check progress**

```bash
cramit status RUN        # counts: done, running, failed, ...
cramit status -v RUN     # every BAM that needs attention, why, and what to do
```

If jobs were killed or you stopped early, run step 2 again: it only picks up
the BAMs that are not done.

**4. Delete the BAMs**

```bash
cramit delete-bams RUN            # dry run: shows what would be deleted
cramit delete-bams --force RUN    # deletes
```

## Recommended: the safe way to run it

Deleting BAMs cannot be undone, so it pays to go slowly the first time.

1. **Start with a test run of 2–3 BAMs.** Put their paths in a text file,
   scan it into a separate run folder, and run it to the end.
   *Why:* a wrong reference, missing permissions or too little time and
   memory then show up on 3 BAMs, not on 3,000 queued jobs. A large one also
   tells you how much time and memory to request for the rest.
2. **Read the scan summary and `RUN/attention.txt` before converting.**
   *Why:* BAMs set aside there need a decision from you. If a folder mixes
   references, give each reference its own run.
3. **Delete only once `cramit status -v RUN` shows nothing still running,
   and you have looked at every failure.**
   *Why:* then you know exactly which BAMs are left, and why.
4. **Do the dry run of `delete-bams` first and read its log**, then run it
   with `--force`.
   *Why:* the log lists every BAM that would be deleted and the reason each
   kept BAM failed a check.
5. **Open one or two CRAMs with the tools you use downstream** (IGV, GATK,
   your pipeline) before deleting.
   *Why:* some tools need extra files next to the reference (GATK wants
   `.fai` and `.dict`). Better to find that out while the BAMs still exist.
6. **If aux tags matter** (UMIs, methylation), run `cramit deep-check` on a
   few samples before deleting.
   *Why:* the standard validation does not compare tags.
7. **Keep the reference safe:** at a stable path, backed up, never deleted.
   *Why:* a CRAM cannot be read without the exact reference it was made with.
   Losing the reference means losing the data.
8. **Back up the CRAMs if the BAMs were backed up.**
   *Why:* once the BAMs are gone, the CRAMs are the only copy.

## Convert one BAM

```bash
cramit convert -r ref.fa.gz sample.bam              # writes sample.cram next to it
cramit convert -r ref.fa.gz -o out/ sample.bam      # or into another folder
```

It exits with `0` if the CRAM is validated, `1` if it was not converted (an
input problem) and `2` if the CRAM does not match the BAM. This makes it easy
to call from Snakemake, Nextflow or your own scripts. It never deletes the
BAM.

## What is checked

**Before a CRAM is kept**, the BAM and the CRAM must give identical
`samtools flagstat` output and identical `samtools stats` checksums over every
read name, sequence and base quality. If they differ, the CRAM is discarded
and the BAM is left alone.

**Before a BAM is deleted**, all of these must hold:

- its CRAM was validated and has not changed since;
- the BAM is a real file, not a symlink;
- a coordinate-sorted BAM's CRAM has its index;
- the CRAM passes `samtools quickcheck`;
- the CRAM lists the same contigs as the BAM.

Aux tags (UMIs, methylation `MM`/`ML`, …) are kept in the CRAM but are not
part of the checksum. If they matter for your data, spot-check a few samples
with `cramit deep-check -r ref.fa.gz sample.bam sample.cram`.

## Settings

| What | How |
|---|---|
| Resources per BAM (SLURM) | `cramit submit -c 4 -m 8G -T 24:00:00 RUN` (these are the defaults, sized for whole-genome BAMs) |
| BAMs running at once (SLURM) | `cramit submit -p 50 RUN` (default 50) |
| Account, partition, other sbatch options | after `--`, or `export CRAMIT_SBATCH_ARGS="--account=mylab --partition=long"` |
| Parallel BAMs / threads on one machine | `cramit run -j 4 -t 4 RUN` |
| samtools comes from a module or conda | `export CRAMIT_SETUP="module load samtools/1.21"` (run before samtools is used) |
| Retry BAMs that failed, after fixing them | add `--retry-failed` to `run` or `submit` |

Put the `export` lines in a small file you `source`, or in your `~/.bashrc`.
SLURM jobs inherit them.

## When something goes wrong

`cramit status -v RUN` names the problem for each BAM and what to do.
[docs/SPECIAL_CASES.md](docs/SPECIAL_CASES.md) explains every case, including
BAMs whose CRAM gets no index, wrong references, truncated files and jobs
that hit their time limit.

## Reading the CRAMs later

A CRAM needs its reference to be decoded, so keep the reference at a stable
path. To get a BAM back:

```bash
samtools view -b -T ref.fa.gz -o sample.bam sample.cram
```

## Credits

The validation idea (comparing `samtools flagstat` and the `samtools stats`
CHK checksum of the BAM and the CRAM) comes from
[bam2cram-check](https://github.com/wtsi-hgi/bam2cram-check) by Irina Colgiu,
Human Genetics Informatics, Wellcome Sanger Institute. cramit is an
independent implementation, developed in the Park Lab at Harvard Medical
School. MIT license.
