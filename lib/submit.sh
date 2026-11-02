#!/usr/bin/env bash
# cramit submit — convert and validate the run's BAMs as a SLURM job array.
#
# Usage: cramit submit [options] RUN [-- extra sbatch options]
#   -c CPUS         CPUs per BAM (default 4)
#   -m MEM          memory per BAM (default 8G)
#   -T TIME         time limit per BAM (default 24:00:00)
#   -p N            at most N BAMs running at once (default 50)
#   -n              dry run: print the sbatch command, submit nothing
#   --retry-failed  also resubmit BAMs that failed before (once you have fixed them)
#
# Account, partition and any other sbatch options go after `--`, or in the
# environment variable CRAMIT_SBATCH_ARGS, e.g.
#   cramit submit RUN -- --account=mylab --partition=long
#
# Submits only the BAMs that are not done yet, so the same command also
# resumes a run and retries BAMs whose job was killed (time limit, memory).
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

if [ "$DRY" -eq 1 ]; then
    echo "dry run — would submit $N BAM(s):"
    printf ' %q' "${SB[@]}"; echo
else
    OUT=$("${SB[@]}")
    echo "$OUT"
    echo "$N BAM(s), at most $CAP at a time ($CPUS CPUs, $MEM, $TIME each)"
    echo "progress:  cramit status $RUN"
fi
if [ "$LATER" -gt 0 ]; then
    echo ""
    echo "$LATER more BAM(s) did not fit in one array (cluster limit: $MAXT tasks)."
    echo "When this array has finished, run the same command again to submit them."
fi
