#!/usr/bin/env bash
# Summarise the parameter-recovery runs (scripts/recovery.jl) into the job
# summary. With `post`, on main, keep a tracking issue in step with the
# verdict: comment on it (opening it if needed) when any seed fails, and
# close it once every seed passes again. With `pr`, keep one comment on the
# pull request `PR_NUMBER` up to date with the summary.
#
# Usage: scripts/recovery_report.sh <dir with verdict_*.txt> [post|pr]
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
  # `warn` passes with a note; `fail` and `unconverged` flag the run.
  case "$status" in
    pass | warn) ;;
    *) failed=1 ;;
  esac
  lines+="- ${status^^}: ${line}"$'\n'
done

{
  echo "### Parameter recovery"
  echo
  printf '%s' "$lines"
} >> "${GITHUB_STEP_SUMMARY:-/dev/stdout}"

if [ "$post" = "pr" ]; then
  marker="<!-- bvd-parameter-recovery -->"
  verdict=$([ "$failed" = "1" ] && echo "fails" || echo "passes")
  body=$(printf '%s\n### Parameter recovery %s\n\n%s\nThis was opened by a bot. Please ping @seabbs for any questions.' \
    "$marker" "$verdict" "$lines")
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
  body=$(printf 'Recovery fails badly on %s (%s):\n\n%s\nA seed fails when a checked quantity'"'"'s true value lies outside its 99%% posterior interval or fewer than 60%% are inside their 90%% interval, and is unconverged when its fit'"'"'s R-hat exceeds 1.05 or its bulk ESS is below 100.\n\nThis was opened by a bot. Please ping @seabbs for any questions.' \
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
