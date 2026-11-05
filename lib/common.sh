# Helpers shared by the cramit commands. Sourced, not run.
#
# A run directory (made by `cramit scan`) holds:
#   cramit.conf          the reference for this run
#   samples.tsv          BAMs to convert: key <TAB> bam <TAB> cram
#   preflight.tsv        the checks for every BAM found
#   attention.txt        BAMs set aside by scan, with the reason
#   samples/<key>/       per BAM: convert.log, stats, status
#   slurm/               job array lists and SLURM logs
# <key> is <name>__<8 hex of md5(BAM path)>, unique even when BAMs share a name.

die()  { echo "cramit: $*" >&2; exit 1; }
warn() { echo "cramit: $*" >&2; }

# Show the comment block at the top of a script (its usage text).
show_usage() { awk 'NR>1 && /^#/{sub(/^# ?/, ""); print; next} NR>1{exit}' "$1"; exit "${2:-0}"; }

# samtools: run $CRAMIT_SETUP first if set (e.g. "module load samtools/1.21"),
# then require samtools on PATH.
need_samtools() {
    if [ -n "${CRAMIT_SETUP:-}" ]; then
        eval "$CRAMIT_SETUP" >/dev/null 2>&1 || warn "CRAMIT_SETUP failed: $CRAMIT_SETUP"
    fi
    command -v samtools >/dev/null 2>&1 \
        || die "samtools not found — put it on PATH, or set CRAMIT_SETUP (see README)"
}

abspath() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s\n' "$PWD/${1#./}" ;; esac; }

sample_key() { printf '%s__%s\n' "$(basename "$1" .bam)" "$(printf '%s' "$1" | md5sum | cut -c1-8)"; }

# Non-blank, non-comment lines of a file.
list_lines() { { grep -v '^[[:space:]]*#' "$1" | grep -v '^[[:space:]]*$'; } || true; }

# Make sure REF has the .fai (and, if bgzipped, .gzi) index samtools needs.
index_reference() {
    [ -f "$1" ] || die "reference not found: $1"
    if [ ! -f "$1.fai" ] || { [[ "$1" == *.gz ]] && [ ! -f "$1.gzi" ]; }; then
        samtools faidx "$1" 2>/dev/null \
            || die "cannot index the reference (is it bgzip-compressed, and is its folder writable?) — run: samtools faidx $1"
    fi
}

load_run() {
    RUN=$(abspath "${1%/}")
    [ -f "$RUN/cramit.conf" ] || die "$RUN is not a cramit run directory — create it with: cramit scan"
    # shellcheck source=/dev/null
    source "$RUN/cramit.conf"
}

# samples/<key>/status holds one line:
#   state <TAB> reason <TAB> detail <TAB> cram_bytes <TAB> where <TAB> time
# state is running, ok or failed. `where` is slurm:<job>_<task> or local:<host>:<pid>.
write_status() {   # dir state [reason] [detail] [cram_bytes]
    local f="$1/status" detail
    detail=$(printf '%s' "${4:--}" | tr '\t\n' '  ' | cut -c1-300)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$2" "${3:--}" "$detail" "${5:--}" \
        "${CRAMIT_WHERE:-local:$(hostname):$$}" "$(date '+%F %T')" > "$f.tmp.$$"
    mv -f "$f.tmp.$$" "$f"
}

# One row per BAM in samples.tsv: key state reason detail bam cram, where state
# is one of: done failed running timeout out_of_memory cancelled incomplete not_started.
# A job killed by SLURM never writes its final status, so a "running" status
# whose job is gone is resolved from the job's SLURM log.
sample_states() {
    local active=""
    if command -v squeue >/dev/null 2>&1; then
        active=$(squeue -h -u "${USER:-$(id -un)}" -r -o '%i' 2>/dev/null | tr '\n' ' ' || true)
    fi
    awk -F'\t' -v OFS='\t' -v run="$RUN" -v active=" $active " -v host="$(hostname)" '
    function slurm_end(id,   f, l, why) {
        f = run "/slurm/" id ".err"; why = ""
        while ((getline l < f) > 0) {
            if (l ~ /DUE TO TIME LIMIT/) why = "timeout"
            else if (tolower(l) ~ /oom[-_ ]kill|out of memory|oom killed/) why = "out_of_memory"
            else if (l ~ /CANCELLED/ && why == "") why = "cancelled"
        }
        close(f); return why
    }
    $1 == "key" { next }
    {
        key = $1; bam = $2; cram = $3; f = run "/samples/" key "/status"; line = ""
        if ((getline line < f) <= 0) line = ""
        close(f)
        if (line == "") { print key, "not_started", "-", "-", bam, cram; next }
        split(line, s, "\t")
        if (s[1] == "ok")     { print key, "done", "-", "-", bam, cram; next }
        if (s[1] == "failed") { print key, "failed", s[2], s[3], bam, cram; next }
        where = s[5]
        if (where ~ /^slurm:/) {
            id = substr(where, 7)
            if (index(active, " " id " ")) { print key, "running", "-", "SLURM job " id, bam, cram; next }
            why = slurm_end(id)
            if (why != "") { print key, why, why, "SLURM job " id ", log: " run "/slurm/" id ".err", bam, cram; next }
            print key, "incomplete", "-", "SLURM job " id " ended without a result", bam, cram; next
        }
        split(where, w, ":")
        if (w[2] == host && system("kill -0 " w[3] " 2>/dev/null") == 0) { print key, "running", "-", "pid " w[3], bam, cram; next }
        print key, "incomplete", "-", "stopped without a result (on " w[2] ")", bam, cram
    }' "$RUN/samples.tsv"
}

# BAMs that `run` / `submit` should (re)start: everything not done and not
# running. Failed BAMs only with --retry-failed, since they need a fix first.
pending_samples() {   # retry_failed(0|1) → key <TAB> bam <TAB> cram
    sample_states | awk -F'\t' -v OFS='\t' -v retry="$1" '
        $2 == "done" || $2 == "running" { next }
        $2 == "failed" && !retry { next }
        { print $1, $5, $6 }'
}

# One BAM of the run, named by its key, BAM or CRAM path, or BAM file name
# (with or without .bam) → key <TAB> bam <TAB> cram. Fails if none or several match.
find_sample() {   # name
    local hits n
    hits=$(awk -F'\t' -v q="$1" -v qa="$(abspath "$1")" 'NR > 1 {
        name = $2; sub(/.*\//, "", name); base = name; sub(/\.bam$/, "", base)
        if ($1 == q || $2 == qa || $3 == qa || name == q || base == q) print }' "$RUN/samples.tsv")
    n=$(grep -c . <<< "$hits" || true)
    if [ "$n" -eq 0 ]; then warn "no BAM in $RUN matches '$1'"; return 1; fi
    if [ "$n" -gt 1 ]; then warn "'$1' matches $n BAMs — name one by its key: $(cut -f1 <<< "$hits" | tr '\n' ' ')"; return 1; fi
    printf '%s\n' "$hits"
}

human_bytes() { awk -v b="$1" 'BEGIN{ split("B KiB MiB GiB TiB PiB", u); i = 1
    while (b >= 1024 && i < 6) { b /= 1024; i++ } printf (i == 1 ? "%d %s" : "%.1f %s"), b, u[i] }'; }
