#!/usr/bin/env bash
#
# Run scripts/linelist/fit_joint.jl detached.
#
# The fit takes five to seven hours on two chains. That is long enough that the
# two things most likely to end it are not the model but the machine, a closed
# terminal and an idle laptop, so this wrapper deals with both and returns
# immediately.
#
#   caffeinate -is   holds off idle and system sleep for the fit's lifetime
#   nohup ... &      detaches from this terminal, so closing it sends no signal
#
# Usage, from anywhere:
#
#   LINELIST_INPUT_DIR=<dir> scripts/linelist/run_joint_fit.sh
#
# Then, at any point:
#
#   ps -o etime=,%cpu= -p "$(cat "$PIDFILE")"
#   tail -3 "$LOG"
#
# `nuts_sample` runs with `progress = false`, so this log says the fit started
# and, hours later, that it finished, and nothing in between. Per-iteration
# progress goes to logs/joint.log, written by the fit callback, which is the
# place to watch: it carries the iteration number, the log posterior and a
# running divergence count, with the two chains interleaved. Elapsed time and
# CPU from `ps` are the coarser indicator; about 200% CPU means both chains are
# sampling.
#
# There is no resume: the fit bypasses the fit cache, which stores completed
# fits only. An interrupted run starts again from zero.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${LINELIST_OUT_DIR:-$ROOT/ignore/linelist}"
mkdir -p "$OUT"
LOG="$OUT/fit_joint.log"
PIDFILE="$OUT/fit_joint.pid"

# Two chains under MCMCThreads, so two threads do the work and any more sit
# idle. The chain count comes from the fit registry and is not changed here.
THREADS=2

# The inputs live where the operator says, and the fit reads them from there.
# Checked up front: five hours is too long to find out at the end that one was
# missing. See scripts/linelist/README.md for both schemas.
if [ -z "${LINELIST_INPUT_DIR:-}" ]; then
    echo "set LINELIST_INPUT_DIR to the directory holding the inputs" >&2
    exit 1
fi
for f in linelist_streams.csv onset_curve_scanned.csv; do
    if [ ! -f "$LINELIST_INPUT_DIR/$f" ]; then
        echo "missing $LINELIST_INPUT_DIR/$f" >&2
        echo "see scripts/linelist/README.md for its schema" >&2
        exit 1
    fi
done

# A second fit would write the same outputs as the first.
if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "already running as PID $(cat "$PIDFILE")" >&2
    echo "stop it with: kill $(cat "$PIDFILE")" >&2
    exit 1
fi

cd "$ROOT"

{
    echo "== linelist joint fit started $(date '+%Y-%m-%d %H:%M:%S') =="
    echo "== threads $THREADS, expect five to seven hours =="
} > "$LOG"

nohup caffeinate -is julia -t "$THREADS" \
    --project=docs \
    scripts/linelist/fit_joint.jl >> "$LOG" 2>&1 &

echo $! > "$PIDFILE"

cat <<EOF
started PID $(cat "$PIDFILE"), five to seven hours

  progress   ps -o etime=,%cpu= -p \$(cat "$PIDFILE")
  log        tail -3 "$LOG"
  stop       kill \$(cat "$PIDFILE")

writes, on success:
  $OUT/linelist_forecast.csv
  $OUT/linelist_stream_estimates.csv
  $OUT/linelist_posterior_summary.csv
EOF
