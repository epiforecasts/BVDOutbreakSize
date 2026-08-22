#!/usr/bin/env bash
#
# Run a fixed grid of scripts/linelist/fit_single.jl calls, one after another,
# detached, each logged on its own.
#
# The grids exist to answer two separate questions with one axis moved at a
# time: stage 1 holds the delay priors at the package defaults and moves the
# data source (the situation reports against the line list), stage 2 holds the
# data source and moves the delay priors (the package defaults against cmmid's
# fresh estimates). Mixing the two into one pass would make a difference in
# R_t ambiguous between "the data changed" and "the delays changed", so the two
# stages are run, and read, separately. See scripts/linelist/README.md and
# scripts/linelist/delays.jl for what each axis means.
#
# Both fits condition on laboratory-confirmed cases alone, which is the only
# case definition the two sources share: the line list's reported stream counts
# every alert raised and the situation reports' counts triaged suspects, so a
# difference between those would not be attributable to the source.
#
# Runs are SERIAL rather than parallel. Each one already claims two threads for
# `nuts_sample`'s two chains (see fit_single.jl and the `-t 2` below); running
# several at once would have them contend for the same cores rather than
# finish sooner, and it would interleave five or six logs into one unreadable
# stream. A failing run is logged and skipped rather than aborting the rest,
# because the value of a comparison grid is in how many of its cells fill in,
# and one bad cell should not cost the others.
#
# An `onsets` run takes forty minutes to an hour and a quarter on this machine.
# A `confirmed` run has not been timed at full sample size; a ten-sample pilot
# takes about five minutes including compilation, which does not extrapolate.
# Either way a stage runs for hours, and a completed run is skipped on a rerun,
# so the cost of guessing the total wrong is only patience. That is long enough
# that, as with run_joint_fit.sh, the likeliest way a run ends early is the
# machine, not the model:
#
#   caffeinate -is   holds off idle and system sleep for the grid's lifetime
#   nohup ... &      detaches from this terminal, so closing it sends no signal
#
# What gets detached is not a single external command but this script's own
# grid loop, and a bash array (the grid itself) does not cross into a nohup'd
# child the way a plain argument does. So the detached process is this same
# script, re-invoked on itself with an internal `--worker` flag that only this
# script ever passes; the grid is defined once, in one place, and both the
# dry run and the real worker read it from there.
#
# Usage, from anywhere:
#
#   LINELIST_INPUT_DIR=<dir> LINELIST_AS_OF=2026-08-10 BVD_DELAY_DIR=<dir> \
#     scripts/linelist/run_grid.sh --stage=1
#
#   scripts/linelist/run_grid.sh --stage=2 --dry-run
#
# Then, at any point:
#
#   ps -o etime=,%cpu= -p "$(cat "$OUT/run_grid.pid")"
#   tail -f "$OUT/run_grid.log"
#
# See scripts/linelist/README.md for the inputs and scripts/linelist/delays.jl
# for BVD_DELAY_DIR and the delay configuration names.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

OUT="${LINELIST_OUT_DIR:-$ROOT/ignore/linelist}"
mkdir -p "$OUT"

GRIDLOG="$OUT/run_grid.log"
PIDFILE="$OUT/run_grid.pid"

# The path this script relaunches itself at, once validated, to get onto the
# far side of nohup. Built from ROOT rather than reused from $0, since $0 may
# be relative to a working directory this process has already left behind by
# the time it is needed.
SELF="$ROOT/scripts/linelist/run_grid.sh"

## -- arguments -------------------------------------------------------------

STAGE=""
DRY_RUN=false
WORKER=false
for a in "$@"; do
    case "$a" in
        --stage=1) STAGE=1 ;;
        --stage=2) STAGE=2 ;;
        --dry-run) DRY_RUN=true ;;
        --worker)
            # Internal only: set when this script relaunches itself detached.
            # Not part of the documented interface.
            WORKER=true
            ;;
        *)
            echo "unknown argument: $a" >&2
            echo "usage: $0 --stage=1|2 [--dry-run]" >&2
            exit 1
            ;;
    esac
done

if [ -z "$STAGE" ]; then
    echo "pass --stage=1 or --stage=2" >&2
    exit 1
fi

## -- the grid ----------------------------------------------------------
##
## fit:data:delays triples, run in this order. This is the one place the grid
## is written down; the file-requirement checks below and both the dry run and
## the worker read it from here rather than from a second copy.

if [ "$STAGE" = 1 ]; then
    # Data source moves, delays held at the package defaults. This is the
    # comparison: `onsets` estimates its own reporting delay from the triangle,
    # so at `repo` it takes no cmmid estimate at all and the only thing that
    # differs between the two runs is the data.
    #
    # `confirmed` was in this grid and has been removed. On linelist_known it
    # does not sample: both chains freeze (max R-hat 2.63, bulk ESS 2.4 of 1000
    # draws, zero divergences), one pinned at the `inv_sqrt_k = 0` truncation
    # boundary where `k` is `1/eps`, the other stuck at 1.327 with a
    # within-chain range of 0.02. The line-list confirmed history is 35 sparse
    # vintages whose increments are multi-day sums, so it carries less
    # dispersion than the negative-binomial observation model can express and
    # the dispersion parameter is driven to its boundary. That is a mismatch
    # between the known-by construction and the observation model, not a
    # sampler setting: four chains split two-and-two the same way, and neither
    # chain count nor initialisation moves a chain that is not exploring.
    # `fit_single.jl` still runs it; it is out of the comparison.
    GRID=(
        "onsets:sitrep:repo"
        "onsets:linelist_known:repo"
    )
else
    # Delay priors move, on both data sources. Each source gets the
    # generation-interval-only step as well as the full swap, because the two
    # do not act alike: the cmmid onset-to-diagnosis delay is about nine days
    # longer than the package prior, so it shifts the trajectory in calendar
    # time as well as changing its level, and a configuration that moves both
    # at once cannot say which did the work. The last confirmed entry brackets
    # the transmission-pair definition the generation interval is fitted from.
    #
    # cmmid_rep runs first and on both sources. It swaps the report delay and
    # leaves the generation interval at the package default, so it is the one
    # linelist-against-sitrep contrast in this stage that does not rest on a
    # generation interval fitted from line-list transmission pairs. Every other
    # cell here carries that generation interval.
    #
    # The onsets fit takes no report-delay override: it estimates that delay
    # from the reporting triangle rather than taking a prior, so only the
    # generation-interval configurations apply to it. delay_config refuses the
    # others there rather than ignoring them.
    # Only the generation-interval configs appear: `onsets` estimates its
    # reporting delay from the triangle, so `delay_config` refuses a
    # report-delay override there rather than ignoring it. All three
    # transmission-pair definitions are run on both sources, which brackets
    # what the cmmid generation interval does to the answer -- the estimate
    # stage 1 deliberately does not rest on.
    GRID=(
        "onsets:sitrep:cmmid_gi_any"
        "onsets:linelist_known:cmmid_gi_any"
        "onsets:sitrep:cmmid_gi_case"
        "onsets:linelist_known:cmmid_gi_case"
        "onsets:sitrep:cmmid_gi_diag"
        "onsets:linelist_known:cmmid_gi_diag"
    )
fi

## -- up-front checks -----------------------------------------------------
##
## Every check below runs before anything is launched, worker included: the
## worker is this same script re-invoked, and re-checking costs nothing next
## to the hours a run takes, whereas finding a missing input an hour in is the
## thing these checks exist to prevent.

if [ -z "${LINELIST_INPUT_DIR:-}" ]; then
    echo "set LINELIST_INPUT_DIR to the directory holding the inputs" >&2
    exit 1
fi

# Every run in both stages fits at least one stream that pulls the onset
# triangle (fit_single.jl calls place_onset_curve unconditionally), so this
# file is never optional.
if [ ! -f "$LINELIST_INPUT_DIR/onset_curve_scanned.csv" ]; then
    echo "missing $LINELIST_INPUT_DIR/onset_curve_scanned.csv" >&2
    echo "see scripts/linelist/README.md for its schema" >&2
    exit 1
fi

# The remaining requirements are only checked if a run in the chosen stage's
# grid actually needs them, worked out from the grid above rather than
# hard-coded per stage, so this cannot drift from the table it is checking.
NEEDS_LINELIST_KNOWN=false
NEEDS_AS_OF=false
NEEDS_DELAY_DIR=false
for entry in "${GRID[@]}"; do
    IFS=: read -r _fit _data _delays <<< "$entry"
    case "$_data" in
        linelist_known) NEEDS_LINELIST_KNOWN=true ;;
        sitrep) NEEDS_AS_OF=true ;;
    esac
    # delay_dir() in delays.jl is only called for a non-"repo" config; "repo"
    # is `nothing`, meaning no overrides and no read of BVD_DELAY_DIR.
    [ "$_delays" != "repo" ] && NEEDS_DELAY_DIR=true
done

if [ "$NEEDS_LINELIST_KNOWN" = true ] && [ ! -f "$LINELIST_INPUT_DIR/linelist_streams_known.csv" ]; then
    echo "missing $LINELIST_INPUT_DIR/linelist_streams_known.csv" >&2
    echo "see scripts/linelist/README.md for its schema" >&2
    exit 1
fi

if [ "$NEEDS_AS_OF" = true ] && [ -z "${LINELIST_AS_OF:-}" ]; then
    echo "set LINELIST_AS_OF: a run in stage $STAGE uses --data=sitrep, which" >&2
    echo "has no replacement streams of its own to take a cut-off from" >&2
    exit 1
fi

if [ "$NEEDS_DELAY_DIR" = true ] && [ -z "${BVD_DELAY_DIR:-}" ]; then
    echo "set BVD_DELAY_DIR: a run in stage $STAGE uses a cmmid delay" >&2
    echo "configuration. See scripts/linelist/delays.jl (delay_dir)." >&2
    exit 1
fi

## -- chain path, for the skip check ---------------------------------------
##
## Matches fit_single.jl's TAG and chain_path exactly: TAG is
## "<fit>_<data>_<delays>" with no "_pilot" suffix, since this grid never
## passes --pilot, and the cache key's date component is LINELIST_AS_OF
## because that check above guarantees it is set for every run that reaches
## here (sitrep needs it directly; the two line-list modes are pinned to the
## same cut-off for the comparison to sit on one grid, per fit_single.jl).
chain_path() {
    echo "$OUT/chains/$1_$2_$3_${LINELIST_AS_OF}.jls"
}

## -- dry run -----------------------------------------------------------
##
## Prints exactly what the worker would do, including which runs it would
## skip, and launches nothing.

if [ "$DRY_RUN" = true ]; then
    for entry in "${GRID[@]}"; do
        IFS=: read -r fit data delays <<< "$entry"
        cp="$(chain_path "$fit" "$data" "$delays")"
        if [ -f "$cp" ]; then
            echo "skip $fit/$data/$delays, chain exists: $cp"
        else
            echo "julia -t 2 --project=docs scripts/linelist/fit_single.jl $fit --data=$data --delays=$delays"
        fi
    done
    exit 0
fi

## -- worker ------------------------------------------------------------
##
## This branch is what actually runs the grid. It only ever runs inside the
## nohup'd process the real-launch branch below starts; nothing here checks
## for an already-running grid, because the one recorded in $PIDFILE is this
## very process.

if [ "$WORKER" = true ]; then
    GRID_STARTED=$(date +%s)
    FAILED=""
    for entry in "${GRID[@]}"; do
        IFS=: read -r fit data delays <<< "$entry"
        cp="$(chain_path "$fit" "$data" "$delays")"
        runlog="$OUT/run_${fit}_${data}_${delays}.log"

        if [ -f "$cp" ]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') skip $fit/$data/$delays, chain exists: $cp" >> "$GRIDLOG"
            continue
        fi

        echo "$(date '+%Y-%m-%d %H:%M:%S') start $fit/$data/$delays -> $runlog" >> "$GRIDLOG"
        run_started=$(date +%s)

        # A failing run must not take the rest of the grid down with it, so
        # its status is captured rather than left to set -e.
        if julia -t 2 --project=docs scripts/linelist/fit_single.jl "$fit" \
            --data="$data" --delays="$delays" > "$runlog" 2>&1; then
            status=ok
        else
            status=FAILED
            FAILED="${FAILED}${FAILED:+ }$fit/$data/$delays"
        fi

        run_elapsed=$(( $(date +%s) - run_started ))
        echo "$(date '+%Y-%m-%d %H:%M:%S') end $fit/$data/$delays status=$status elapsed=${run_elapsed}s" >> "$GRIDLOG"
    done

    total_elapsed=$(( $(date +%s) - GRID_STARTED ))
    {
        echo "== linelist grid (stage $STAGE) finished $(date '+%Y-%m-%d %H:%M:%S'), elapsed ${total_elapsed}s =="
        if [ -n "$FAILED" ]; then
            echo "== failed: $FAILED =="
        else
            echo "== all runs succeeded =="
        fi
    } >> "$GRIDLOG"

    # Empty FAILED -> success; non-empty -> the worker, and so the nohup'd
    # process ps reports on, exits non-zero.
    [ -z "$FAILED" ]
    exit 0
fi

## -- real launch ---------------------------------------------------------

# A second grid would write the same run logs and the same chains as the
# first, racing it rather than adding to it.
if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "a grid is already running as PID $(cat "$PIDFILE")" >&2
    echo "stop it with: kill $(cat "$PIDFILE")" >&2
    exit 1
fi

{
    echo "== linelist grid (stage $STAGE) started $(date '+%Y-%m-%d %H:%M:%S') =="
} >> "$GRIDLOG"

nohup caffeinate -is "$SELF" --stage="$STAGE" --worker >> "$GRIDLOG" 2>&1 &

echo $! > "$PIDFILE"

cat <<EOF
started grid PID $(cat "$PIDFILE") (stage $STAGE), running serially, detached

  progress   ps -o etime=,%cpu= -p \$(cat "$PIDFILE")
  grid log   tail -f "$GRIDLOG"
  run log    tail -f "$OUT/run_<fit>_<data>_<delays>.log"
  stop       kill \$(cat "$PIDFILE")
EOF
