#!/usr/bin/env bash
#
# Run scripts/linelist/fit_joint.jl detached, for docs/model-comparison.qmd.
#
# The fit took 6 hours 50 minutes when it was run on 2026-08-10: 500 samples on
# two chains, at roughly 46 seconds an iteration rather than the 22 first
# estimated. That is long enough that the two things most likely to end it are
# not the model but the machine — a closed terminal and an idle laptop — so this
# wrapper deals with both and returns immediately.
#
#   caffeinate -is   holds off idle and system sleep for the fit's lifetime
#   nohup ... &      detaches from this terminal, so closing it sends no signal
#
# A previous attempt died on `signal 15: Terminated`, which is an external kill
# rather than anything in the model. Started this way it survives.
#
# Usage, from anywhere:
#
#   scripts/linelist/run_joint_fit.sh
#
# Then, at any point:
#
#   ps -o etime=,%cpu= -p "$(cat "$PIDFILE")"
#   tail -3 "$LOG"
#
# `nuts_sample` runs with `progress = false`, so this log says the fit started
# and, seven hours later, that it finished, and nothing in between. Per-iteration
# progress goes to logs/joint.log, written by the registry's own fit callback,
# which is the place to watch: it carries the iteration number, the log posterior
# and a running divergence count, with the two chains interleaved. Elapsed time
# and CPU from `ps` are the coarser indicator; about 200% CPU means both chains
# are sampling.
#
# There is no resume: the fit bypasses the BVDOutbreakSize fit cache, which
# stores completed fits only. An interrupted run starts again from zero.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${LINELIST_OUT_DIR:-$ROOT/ignore/linelist}"
mkdir -p "$OUT"
LOG="$OUT/fit_joint.log"
PIDFILE="$OUT/fit_joint.pid"

# Two chains under MCMCThreads, so two threads do the work and any more sit
# idle. The chain count comes from the released report's own fit registry and
# is deliberately not changed here.
THREADS=2

# Inputs, with the script that writes each. Checked up front: four hours is too
# long to find out at the end that one was missing.
if [ ! -f "$OUT/bos_linelist_streams.csv" ]; then
    echo "missing $OUT/bos_linelist_streams.csv" >&2
    echo "run: Rscript bvd-analysis's stream builder 20260803" >&2
    exit 1
fi
if [ ! -f "$OUT/bos_observations.toml" ]; then
    echo "missing $OUT/bos_observations.toml" >&2
    echo "run: bvd-analysis's fetch_bos.sh" >&2
    exit 1
fi

# A second fit would write the same three outputs as the first.
if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "already running as PID $(cat "$PIDFILE")" >&2
    echo "stop it with: kill $(cat "$PIDFILE")" >&2
    exit 1
fi

cd "$ROOT"
mkdir -p "$OUT"

{
    echo "== linelist joint fit started $(date '+%Y-%m-%d %H:%M:%S') =="
    echo "== threads $THREADS, expect 7 hours or more =="
} > "$LOG"

nohup caffeinate -is julia -t "$THREADS" \
    --project=docs \
    scripts/linelist/fit_joint.jl >> "$LOG" 2>&1 &

echo $! > "$PIDFILE"

cat <<EOF
started PID $(cat "$PIDFILE"), 7 hours or more

  progress   ps -o etime=,%cpu= -p \$(cat "$PIDFILE")
  log        tail -3 "$LOG"
  stop       kill \$(cat "$PIDFILE")

writes, on success:
  $OUT/bos_linelist_forecast.csv
  $OUT/bos_linelist_stream_estimates.csv
  $OUT/bos_linelist_posterior_summary.csv
EOF
