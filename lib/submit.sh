#!/usr/bin/env bash
# cramit submit — convert and validate the run's BAMs as a SLURM job array.
#
# Usage: cramit submit [options] RUN [-- extra sbatch options]
#
# Each BAM that is not done yet becomes one array task. The task converts the
# BAM, validates the CRAM and puts NAME.cram next to the BAM. Job logs go to
# RUN/slurm/. Run the same command again to resume, or to resubmit BAMs whose
# job was killed (time limit, memory).
#
# Options (the defaults suit whole-genome BAMs):
#   -c CPUS         CPUs per BAM (default 4)
#   -m MEM          memory per BAM (default 8G)
#   -T TIME         time limit per BAM (default 24:00:00)
#   -p N            at most N BAMs running at once (default 50)
#   -n              dry run: show the plan and the sbatch command, submit nothing
#   --retry-failed  also resubmit BAMs that failed before (once you have fixed them)
#
# Fitting it to your cluster:
#   account, partition   sbatch options go after --, passed on as they are:
#                          cramit submit RUN -- --account=mylab --partition=long
#                        `sinfo -s` lists the partitions and their time limits.
#   samtools from a      export CRAMIT_SETUP="module load samtools/1.21"
#   module or conda      (every job runs this before it uses samtools)
#   same options always  export CRAMIT_SBATCH_ARGS="--account=mylab --partition=long"
#                        Put the exports in a file you source before cramit.
#   time and memory      submit 2-3 BAMs first, then see what they used:
#                          sacct -j JOBID -o JobID,Elapsed,MaxRSS
#
# Before the first real submit, run with -n to check the plan.
set -euo pipefail
source "$CRAMIT_HOME/lib/common.sh"

# ── Task mode: SLURM runs this once per array index ──────────
if [ "${1:-}" = --task ]; then
    load_run "$2"
    IFS=$'\t' read -r KEY BAM _ < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$3")
    [ -n "${BAM:-}" ] || die "no BAM at line $SLURM_ARRAY_TASK_ID of $3"
    export CRAMIT_WHERE="slurm:${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    exec bash "$CRAMIT_HOME/lib/convert.sh" -r "$REF" -l "$RUN/samples/$KEY" -t "${SLURM_CPUS_PER_TASK:-4}" "$BAM"
fi

# ── Submit mode ───────────────────────────────────────────────
CPUS=4; MEM=8G; TIME=24:00:00; CAP=50; DRY=0; RETRY=0; ARGS=(); EXTRA=()
read -r -a EXTRA <<< "${CRAMIT_SBATCH_ARGS:-}"
while [ $# -gt 0 ]; do
    case "$1" in
        -c) CPUS="$2"; shift 2 ;;  -m) MEM="$2"; shift 2 ;;  -T) TIME="$2"; shift 2 ;;
        -p) CAP="$2"; shift 2 ;;   -n) DRY=1; shift ;;       --retry-failed) RETRY=1; shift ;;
        --) shift; EXTRA+=("$@"); break ;;
        -h|--help) show_usage "$0" ;;
        -*) die "unknown option: $1 (sbatch options go after --)" ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ ${#ARGS[@]} -eq 1 ] || show_usage "$0" 1
[[ "$CPUS" =~ ^[1-9][0-9]*$ && "$CAP" =~ ^[1-9][0-9]*$ ]] || die "-c and -p must be positive integers"
load_run "${ARGS[0]}"
command -v sbatch >/dev/null 2>&1 || [ "$DRY" -eq 1 ] || die "sbatch not found — no SLURM here? Use: cramit run $RUN"

JOB="cramit_$(basename "$RUN")_$(printf '%s' "$RUN" | md5sum | cut -c1-6)"
if command -v squeue >/dev/null 2>&1; then
    ACTIVE=$(squeue -h -u "${USER:-$(id -un)}" -n "$JOB" -o '%A' 2>/dev/null | sort -u | paste -sd' ' - || true)
    [ -z "$ACTIVE" ] || die "an earlier array for this run is still queued or running (job $ACTIVE) — wait for it to finish, or scancel it"
fi

mkdir -p "$RUN/slurm"
LIST="$RUN/slurm/batch.$(date +%Y%m%d_%H%M%S).$$.tsv"
pending_samples "$RETRY" > "$LIST"
N=$(wc -l < "$LIST")
if [ "$N" -eq 0 ]; then
    rm -f "$LIST"
    echo "nothing to submit — every BAM is done, running, or failed (see: cramit status -v $RUN)"
    exit 0
fi

# One array can hold at most MaxArraySize-1 tasks; the rest wait for the next submit.
MAXA=$(scontrol show config 2>/dev/null | awk '$1=="MaxArraySize"{print $3}' || true)
MAXT=$(( ${MAXA:-1001} - 1 ))
LATER=0
if [ "$N" -gt "$MAXT" ]; then LATER=$((N - MAXT)); N=$MAXT; head -n "$N" "$LIST" > "$LIST.t" && mv "$LIST.t" "$LIST"; fi

SB=(sbatch --job-name="$JOB" --array="1-${N}%${CAP}" --cpus-per-task="$CPUS" --mem="$MEM" --time="$TIME"
    --output="$RUN/slurm/%A_%a.out" --error="$RUN/slurm/%A_%a.err" "${EXTRA[@]}"
    --wrap="bash $(printf '%q' "$CRAMIT_HOME/cramit") _task $(printf '%q' "$RUN") $(printf '%q' "$LIST")")

# The jobs inherit this environment, so samtools must be found here too.
if [ -n "${CRAMIT_SETUP:-}" ]; then eval "$CRAMIT_SETUP" >/dev/null 2>&1 || true; fi
if command -v samtools >/dev/null 2>&1; then
    ST="$(samtools --version 2>/dev/null | awk 'NR==1')"
    if [ -n "${CRAMIT_SETUP:-}" ]; then ST+=", loaded by CRAMIT_SETUP"; else ST+=" from PATH ($(command -v samtools))"; fi
else
    ST="NOT FOUND — every job would fail; set CRAMIT_SETUP or put samtools on PATH"
fi
if [ ${#EXTRA[@]} -gt 0 ]; then SO="${EXTRA[*]}"; else SO="none, so SLURM uses your default account and partition"; fi

if [ "$DRY" -eq 1 ]; then echo "cramit submit — dry run, nothing is submitted"; else echo "cramit submit"; fi
NB="$N to convert, one array task each"; [ "$LATER" -eq 0 ] || NB+=" ($LATER more after this array)"
printf '  %-15s %s\n' \
    "BAMs"            "$NB" \
    "each BAM gets"   "$CPUS CPUs, $MEM memory, $TIME" \
    "at once"         "at most $CAP BAMs" \
    "sbatch options"  "$SO" \
    "samtools"        "$ST" \
    "job logs"        "$RUN/slurm/<job>_<task>.out and .err"

if [ "$DRY" -eq 1 ]; then
    echo ""
    echo "sbatch command:"
    printf ' %q' "${SB[@]}"; echo
    echo ""
    echo "Change resources with -c -m -T -p, add sbatch options after -- (see: cramit submit -h)."
else
    command -v samtools >/dev/null 2>&1 || die "samtools not found — set CRAMIT_SETUP or put samtools on PATH (see: cramit submit -h)"
    "${SB[@]}"
    echo "progress:  cramit status $RUN"
fi
if [ "$LATER" -gt 0 ]; then
    echo ""
    echo "$LATER more BAM(s) did not fit in one array (cluster limit: $MAXT tasks)."
    echo "When this array has finished, run the same command again to submit them."
fi
