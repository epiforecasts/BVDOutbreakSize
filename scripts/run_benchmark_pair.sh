#!/usr/bin/env bash
# Reproduce the CI benchmark comparison locally. Called by
# `task benchmark-pair` and re-usable from the command line:
#
#   ./scripts/run_benchmark_pair.sh [base-rev] [head-rev]
#
# Defaults to `main` against the current `HEAD`.
#
# It does what .github/workflows/benchmark.yaml does: materialise both
# revisions as git worktrees, then hand the two trees to
# benchmark/ci/run_pair.jl, which benchmarks them one after the other in one
# process. The worktrees are what let AirspeedVelocity run at all here; see
# the header of benchmark/ci/run_pair.jl.
#
# The comparison is only worth reading on a quiet machine. Two arms of the
# suite is roughly twice the cost of a single `task benchmark`, most of it
# compiling gradients.
set -euo pipefail

cd "$(dirname "$0")/.."

BASE_REV="${1:-main}"
HEAD_REV="${2:-HEAD}"

WORKTREE_DIR=".benchmark-worktrees"
RESULTS_DIR="benchmark-results"

# Worktrees are rebuilt from scratch each time: a stale tree left at an old
# revision would be benchmarked silently and reported as that revision.
for arm in main pr; do
    if [ -d "$WORKTREE_DIR/$arm" ]; then
        git worktree remove --force "$WORKTREE_DIR/$arm"
    fi
done
git worktree add --detach "$WORKTREE_DIR/main" "$BASE_REV"
git worktree add --detach "$WORKTREE_DIR/pr" "$HEAD_REV"

julia --project=benchmark/ci -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark/ci benchmark/ci/run_pair.jl \
    "$WORKTREE_DIR/main" "$WORKTREE_DIR/pr" "$RESULTS_DIR"
julia --project=benchmark/ci benchmark/ci/comment.jl \
    "$RESULTS_DIR" BVDOutbreakSize main pr "$RESULTS_DIR/comment.md"

echo
echo "Results in $RESULTS_DIR, comment in $RESULTS_DIR/comment.md."
echo "Worktrees left in $WORKTREE_DIR; remove with:"
echo "  git worktree remove --force $WORKTREE_DIR/main $WORKTREE_DIR/pr"
