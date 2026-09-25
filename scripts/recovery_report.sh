#!/usr/bin/env bash
# Summarise the parameter-recovery runs (scripts/recovery.jl) across seeds
# (scripts/recovery_report.jl) into the job summary. With `post`, on main,
# keep a tracking issue in step with the verdict: comment on it (opening it
# if needed) when recovery fails or a fit does not converge, and close it
# once recovery passes again. With `pr`, keep one comment on the pull
# request `PR_NUMBER` up to date with the summary.
#
# Usage: scripts/recovery_report.sh <dir with recovery_*.csv> [post|pr]
set -euo pipefail
dir=${1:?directory of recovery results}
post=${2:-}
title="Parameter recovery check failing"

shopt -s nullglob
files=("$dir"/recovery_[0-9]*.csv)
if [ ${#files[@]} -eq 0 ]; then
  echo "No recovery results found; the recovery runs did not finish." >&2
  exit 1
fi

julia --project=docs scripts/recovery_report.jl "$dir"
status=$(cat "$dir/status")
# `warn` passes with a note; `fail` and `unconverged` flag the run.
case "$status" in
  pass | warn) failed=0 ;;
  *) failed=1 ;;
esac

cat "$dir/report.md" >> "${GITHUB_STEP_SUMMARY:-/dev/stdout}"
footer=$'\nThis was opened by a bot. Please ping @seabbs for any questions.'

if [ "$post" = "pr" ]; then
  marker="<!-- bvd-parameter-recovery -->"
  body="${marker}"$'\n'"$(cat "$dir/report.md")"$'\n'"${footer}"
  existing=$(gh api "repos/${GITHUB_REPOSITORY}/issues/${PR_NUMBER}/comments" \
    --paginate --jq ".[] | select(.body | startswith(\"$marker\")) | .id" |
    head -n 1)
  if [ -n "$existing" ]; then
    gh api -X PATCH "repos/${GITHUB_REPOSITORY}/issues/comments/${existing}" \
      -f body="$body" > /dev/null
  else
    gh pr comment "$PR_NUMBER" --body "$body"
  fi
  exit 0
fi
[ "$post" = "post" ] || exit 0

run="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${GITHUB_RUN_ID}"
issue=$(gh issue list --state open --search "\"$title\" in:title" \
  --json number --jq '.[0].number // empty')
if [ "$failed" = "1" ]; then
  body="Recovery is ${status} on ${GITHUB_SHA:0:8} (${run})."$'\n\n'"$(cat "$dir/report.md")"$'\n'"${footer}"
  if [ -n "$issue" ]; then
    gh issue comment "$issue" --body "$body"
  else
    gh issue create --title "$title" --body "$body"
  fi
elif [ -n "$issue" ]; then
  gh issue comment "$issue" --body "Recovery passes again on ${GITHUB_SHA:0:8} (${run}). Closing.
${footer}"
  gh issue close "$issue"
fi
