#!/usr/bin/env bash
# Download the model fits produced by a docs CI run into the local fit cache,
# so a local render loads them instead of refitting. This mirrors what the
# render jobs do in CI, where the per-fit matrix fits and the render only
# loads.
#
#   ./scripts/fetch_fits.sh              # latest successful docs run on main
#   ./scripts/fetch_fits.sh 12345678     # a specific run id
#   BVD_FIT_REF=my-branch ./scripts/fetch_fits.sh
#
# Fit cache keys are a content hash of the model source, the data and the
# sampler settings, so a downloaded fit is loaded only when the working tree
# matches the commit that produced it. Change either and the render reports a
# miss; fit those locally with `task fit-all`.
#
# The fit artifacts are kept for three days (docs.yml), so only recent runs
# can be fetched. Needs the `gh` CLI authenticated against the repo.
set -euo pipefail

REPO=${BVD_FIT_REPO:-epiforecasts/BVDOutbreakSize}
REF=${BVD_FIT_REF:-main}
CACHE=${BVD_FIT_CACHE:-logs/fit_cache}
RUN=${1:-}

if [ -z "$RUN" ]; then
  RUN=$(gh run list --repo "$REPO" --workflow docs.yml --branch "$REF" \
    --status success --limit 1 --json databaseId --jq '.[0].databaseId')
  if [ -z "$RUN" ] || [ "$RUN" = "null" ]; then
    echo "no successful docs.yml run found on $REF in $REPO" >&2
    exit 1
  fi
fi

mkdir -p "$CACHE"
echo "Fetching fits from $REPO run $RUN into $CACHE"
gh run download "$RUN" --repo "$REPO" --pattern 'fit-*' --dir "$CACHE"

# Each `fit-<id>` artifact unpacks into its own directory; the cache is flat
# and content-addressed, so the .jls files are moved up and the directories
# dropped.
find "$CACHE" -mindepth 2 -name '*.jls' -exec mv -f {} "$CACHE"/ \;
find "$CACHE" -mindepth 1 -maxdepth 1 -type d -name 'fit-*' -exec rmdir {} +

echo "Cached $(find "$CACHE" -maxdepth 1 -name '*.jls' | wc -l) fits in $CACHE"

# The fits were serialised under the manifest that run resolved. A local
# environment that resolved differently can fail to deserialise the chains
# (the chains library in particular). Pin the same versions when that happens.
cat <<'NOTE'

If the render fails to deserialise a chain, pin the run's package versions:
  gh run download <run> --repo <repo> --name manifest --dir .
  julia --project=docs -e 'using Pkg; Pkg.instantiate()'
NOTE
