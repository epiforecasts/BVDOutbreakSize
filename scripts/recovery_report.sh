#!/usr/bin/env bash
# Summarise the parameter-recovery runs (scripts/recovery.jl) into the job
# summary and, on main, keep a tracking issue in step with the verdict:
# comment on it (opening it if needed) when any seed fails, and close it
# once every seed passes again.
#
# Usage: scripts/recovery_report.sh <dir with verdict_*.txt> [post]
set -euo pipefail
dir=${1:?directory of verdict files}
post=${2:-}
title="Parameter recovery check failing"

shopt -s nullglob
files=("$dir"/verdict_*.txt)
if [ ${#files[@]} -eq 0 ]; then
  echo "No recovery verdicts found; the recovery runs did not finish." >&2
  exit 1
fi

failed=0
lines=""
for f in "${files[@]}"; do
  IFS=$'\t' read -r status line < "$f"
  [ "$status" = "pass" ] || failed=1
  lines+="- ${status^^}: ${line}"$'\n'
done

{
  echo "### Parameter recovery"
  echo
  printf '%s' "$lines"
} >> "${GITHUB_STEP_SUMMARY:-/dev/stdout}"

[ "$post" = "post" ] || exit 0

run="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${GITHUB_RUN_ID}"
issue=$(gh issue list --state open --search "\"$title\" in:title" \
  --json number --jq '.[0].number // empty')
if [ "$failed" = "1" ]; then
  body=$(printf 'Recovery fails badly on %s (%s):\n\n%s\nA seed fails when a checked quantity'"'"'s true value lies outside its 99%% posterior interval, or fewer than half are inside their 90%% interval.\n\nThis was opened by a bot. Please ping @seabbs for any questions.' \
    "${GITHUB_SHA:0:8}" "$run" "$lines")
  if [ -n "$issue" ]; then
    gh issue comment "$issue" --body "$body"
  else
    gh issue create --title "$title" --body "$body"
  fi
elif [ -n "$issue" ]; then
  gh issue comment "$issue" --body "Recovery passes again on ${GITHUB_SHA:0:8} (${run}). Closing.

This was opened by a bot. Please ping @seabbs for any questions."
  gh issue close "$issue"
fi
