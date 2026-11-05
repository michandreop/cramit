# Special cases

Most BAMs need nothing from you. This page covers the ones that do. You don't
have to read it up front: `cramit status -v RUN` tells you which case applies
to each BAM, and the table below tells you where to look.

Every case is decided per BAM. A problem with one BAM never stops the others.

## Quick lookup

| You see | Where | Meaning | Section |
|---|---|---|---|
| `CORRUPT` | scan | BAM is truncated, damaged or still being written | [Truncated BAM](#truncated-or-damaged-bam) |
| `WRONG_REFERENCE` | scan | BAM was aligned to a different reference | [Wrong reference](#wrong-reference) |
| `NOT_WRITABLE` | scan | You can't write the CRAM next to the BAM | [Folder not writable](#folder-not-writable) |
| `MISSING` | scan | Path in your list does not exist | — fix the list |
| `bam_corrupt` | status | BAM could not be read to the end | [Truncated BAM](#truncated-or-damaged-bam) |
| `wrong_reference` | status | Contigs don't match the reference | [Wrong reference](#wrong-reference) |
| `not_sorted` | status | Header says sorted, reads are not | [Not really sorted](#header-says-sorted-but-the-reads-are-not) |
| `existing_cram_differs` | status | A different CRAM is already next to the BAM | [CRAM already there](#a-cram-is-already-next-to-the-bam) |
| `cram_differs` | status | New CRAM didn't match its BAM (rare) | [Validation failed](#validation-failed) |
| `timeout`, `out_of_memory` | status | SLURM stopped the job | [Time and memory](#time-limit-or-out-of-memory) |
| `incomplete`, `cancelled` | status | Job ended without a result | [Stopped jobs](#stopped-jobs) |
| `convert_error`, `validate_error`, `internal_error` | status | Something else went wrong | read `RUN/samples/<key>/convert.log` |
| `SKIP ...` | delete-bams log | A safety check failed, so the BAM was kept | [Kept BAMs](#bams-that-delete-bams-keeps) |
| `expected differences` | deep-check | Changes CRAM makes on purpose | [What CRAMs change on purpose](#what-crams-change-on-purpose) |

---

## BAMs that convert fine, with a limitation

### No index: queryname-sorted, unsorted or unaligned BAMs

**What happens.** The CRAM is converted and validated, and the BAM can be
deleted as usual. But the CRAM has **no `.crai` index**, because only a CRAM
whose reads are in position order can be indexed. The read order is kept as
it was, so tools that need name-grouped reads keep working.

Typical sources: UMI/duplex consensus pipelines, `samtools sort -n`, STAR
`toTranscriptome`, unaligned BAMs (uBAMs, raw ONT/PacBio reads).

**Limitation.** No region queries (`samtools view x.cram chr1:1-1000`, IGV).
If you need them, sort the BAM first (`samtools sort`) and convert the sorted
file. Unaligned BAMs also shrink less, since there is no reference to
compress against.

If the header states no sort order at all, cramit tries to index and keeps
the index if the reads turn out to be sorted.

### Aux tags and other fields are not part of validation

**What happens.** CRAM keeps every field and aux tag (`RG`, `MD`/`NM`, UMI
tags like `RX`/`MI`, methylation `MM`/`ML`, ...). But validation compares
read names, sequences, qualities and flags only. A changed position, CIGAR
or tag would therefore not be noticed.

**What to do.** Before deleting, deep-check a few BAMs: at least one per
aligner or pipeline version in the folder (`samtools view -H x.bam | grep
'^@PG'` shows it). It matters most when tags carry data (UMI/duplex
consensus, methylation).

```bash
cramit deep-check RUN sample                               # a BAM of a run, every record
cramit deep-check -n 1000000 RUN sample                    # the first million records
cramit deep-check -r ref.fa.gz sample.bam sample.cram      # any BAM and CRAM
sbatch -c 3 --mem 4G -t 8:00:00 --wrap "cramit deep-check RUN sample"   # on SLURM
```

It reads both files to the end (1 to 1.5 hours per 100 million reads), counts
the differences by field, shows the first few, and exits with 0 (identical),
1 (different) or 2 (could not read a file). Tags are compared regardless of
their order. Two differences are expected and listed separately; see the
next section.

### What CRAMs change on purpose

deep-check lists these under "expected differences" and still reports
IDENTICAL. Neither one affects what downstream tools compute.

**Unmapped reads lose their CIGAR and MAPQ.** CRAM does not store the CIGAR
or MAPQ of a read flagged unmapped (`0x4`), so they come back as `*` and
`0`. Some aligners, notably older BWA, leave values there, e.g. on reads
hanging over the end of a contig that BWA then flags unmapped. The SAM
specification says that for unmapped reads no assumption can be made about
these fields; Picard's `ValidateSamFile` reports a MAPQ on an unmapped read
as an error, and `CleanSam` sets it to 0. Variant callers, coverage tools and realignment
(`samtools fastq`, `RevertSam`) do not use them. Expect anything from none
to a few percent of the unmapped reads in older data.

**`MD` and `NM` tags appear.** When samtools reads a CRAM, it computes `MD`
and `NM` from the reference and adds them to every mapped read that has
none. If the BAM had them, they come back unchanged.

### Symlinked BAMs

**What happens.** When scanning a folder, symlinks named `*.bam` are skipped
(scan says how many). If you list a symlink explicitly, it is converted and
the CRAM goes next to the link. But `delete-bams` never deletes through a
link: removing it would free no space and would leave the real BAM behind.

**What to do.** Scan the folder that holds the real files.

### BAMs without M5 tags (e.g. GTEx)

Nothing to do. samtools adds the reference checksums (`M5`) to the CRAM
header by itself, so the CRAM can find its reference by checksum later.

---

## BAMs cramit refuses until you fix them

A refused BAM is never deleted, and no CRAM is written for it.

### Truncated or damaged BAM

**How you see it.** Scan: `CORRUPT`. Status: `bam_corrupt`.

**What happens.** Scan checks only the header and the end of the file, so it
catches truncated files and files still being written or downloaded. Damage
in the middle of a file is found during conversion, which reads every record.
Either way, the BAM is set aside and no CRAM is kept.

**What to do.** Re-download or regenerate the BAM, then
`cramit run --retry-failed RUN` (or `submit --retry-failed`). If the file was
still being written, wait until it is finished and run `cramit scan` again.

### Wrong reference

**How you see it.** Scan: `WRONG_REFERENCE`, with the reason: contigs missing
from the reference, another contig length, or another sequence checksum.
Status: `wrong_reference`.

**Why cramit refuses.** With the wrong reference, samtools does not fail. It
quietly stores the reference sequence inside the CRAM, so the file comes out
larger than expected and is tied to the wrong reference. Typical causes:
GRCh37 vs GRCh38, `1` vs `chr1` naming, a folder that mixes references.

**Finding the right reference.** The reason in `attention.txt` names the
contigs that don't match, and the reference the BAM's header names (its
`AS`/`UR` tags), if any. Then:

- contigs missing: another naming (`1` vs `chr1`) or another build;
- another length: another build (GRCh37 vs GRCh38) or patch level;
- same names and lengths, another checksum (`M5`): the same build, prepared
  differently. Compare with `samtools dict ref.fa`. A common case is chrY:
  some GRCh37 references (hs37d5) mask its pseudo-autosomal regions with N,
  others don't; when only Y differs, try the other one.

**What to do.** Use one run per reference: put the set-aside BAMs in a list
and scan it into a new run folder with the right `-r`.

### Header says sorted, but the reads are not

**How you see it.** Status: `not_sorted`.

**What happens.** The header claims `SO:coordinate`, but the reads are out of
order, so the CRAM cannot be indexed. cramit refuses rather than leave a
"sorted" CRAM without an index.

**What to do.** Sort the BAM (`samtools sort`), or correct the header if the
tag is simply wrong. Then retry with `--retry-failed`.

### A CRAM is already next to the BAM

**What happens.** cramit never overwrites an existing `NAME.cram`. It checks
that file against the BAM instead:

- **it matches:** it is kept (and indexed if needed), and the BAM counts as
  done;
- **it differs or is incomplete:** status `existing_cram_differs`, and the
  file is left untouched.

**What to do (if it differs).** Find out where that CRAM came from. Move or
delete it, then retry with `--retry-failed`.

### Folder not writable

**How you see it.** Scan: `NOT_WRITABLE`.

**What to do.** CRAMs are written next to their BAMs, so you need write
permission there, and you need it to delete BAMs too. If you can't get it,
convert into a folder of your own, one BAM at a time:
`cramit convert -r REF -o OUTDIR BAM`.

---

## Jobs that stop

### Time limit or out of memory

**How you see it.** Status: `timeout` or `out_of_memory` (read from the
job's SLURM log in `RUN/slurm/`).

**What to do.** Submit again with more:

```bash
cramit submit -T 48:00:00 RUN    # more time
cramit submit -m 16G RUN         # more memory
```

Only the BAMs that are not done are submitted. A stopped conversion leaves no
CRAM behind, so the retry starts clean. `samtools stats` needs a few GB on
whole-genome data, so don't go below the default 8G for WGS.

### Stopped jobs

**How you see it.** Status: `incomplete` (the job ended without writing a
result: killed, node failure, Ctrl-C) or `cancelled` (`scancel`).

**What to do.** Run `cramit run RUN` or `cramit submit RUN` again. These BAMs
are picked up automatically.

### More BAMs than one SLURM array allows

Clusters cap the size of a job array (`MaxArraySize`). If a run has more BAMs
than that, `submit` sends as many as fit and tells you so. When that array
has finished, run the same `cramit submit` command again to send the rest.
While an array for the run is still queued or running, `submit` refuses, so
the same BAM never runs twice.

---

## BAMs that delete-bams keeps

`delete-bams` redoes every check on disk at the moment of deletion. If one
fails, the BAM is kept and the reason is written to
`RUN/delete-bams.<time>.<pid>.log`:

| Reason | Meaning |
|---|---|
| `bam_is_symlink` | See [symlinked BAMs](#symlinked-bams) |
| `cram_missing` | The CRAM was moved or deleted since validation |
| `cram_changed_since_validation` | The CRAM's size changed after it was validated: it is not the file that was checked |
| `crai_missing`, `crai_empty` | A coordinate-sorted BAM's CRAM has lost its index |
| `cram_quickcheck_failed` | The CRAM is damaged (e.g. its end-of-file block is missing) |
| `contigs_differ` | The CRAM header no longer lists the same contigs as the BAM |
| `bam_unreadable`, `cram_unreadable` | samtools cannot read the header |

To make a CRAM again (it changed, or it was made while the reference sat at
another path): `cramit redo RUN sample`, then `cramit run RUN` or
`cramit submit RUN`. redo deletes the CRAM only while its BAM is still there.

### Validation failed

**How you see it.** Status: `cram_differs`.

**What happens.** The new CRAM decoded to different reads than the BAM.
This should not happen with a correct reference; a transient storage error
is the most likely cause. The CRAM is discarded and the BAM is untouched.

**What to do.** Retry with `--retry-failed`. If it fails again on the same
BAM, please report it with the sample's `convert.log`.

---

## Good to know

**Re-running is always safe.** A CRAM gets its final name only after it has
been validated, so a `NAME.cram` written by cramit is never half-finished.
Finished BAMs are skipped, and stopped ones start over cleanly.

**Two BAMs with the same file name** in different folders don't clash: each
CRAM goes next to its own BAM, and the log folders are named
`<name>__<8 characters of a hash of the full path>`.

**What validation proves.** The BAM and the CRAM give identical
`samtools flagstat` output and identical `samtools stats` `CHK` checksums,
which cover every read name, sequence and base quality. So no read was lost,
added or changed in any of those. Not covered: aux tags (see above) and the
header text (samtools adds a `@PG` line).

**Reading the CRAMs later.** A CRAM needs its reference to be decoded. Keep
the reference at a stable path, or make it findable through htslib's
`REF_PATH`/`REF_CACHE`. To get a BAM back:
`samtools view -b -T ref.fa.gz -o sample.bam sample.cram`. cramit writes
CRAM 3.0, which samtools, GATK, Picard and IGV all read.

**Where the files are.** In `RUN/`: `samples.tsv` (BAMs to convert),
`attention.txt` (BAMs set aside by scan), `preflight.tsv` (scan details),
`status.tsv` (written by `cramit status`), `samples/<key>/` (each BAM's
`convert.log`, stats and status) and `slurm/` (job logs).
