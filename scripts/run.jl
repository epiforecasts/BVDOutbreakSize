# Entry point for regenerating the published results.
#
# Runs both report pages (analysis and sensitivity), which fit the models
# (loaded from the content-addressed cache when present) and write the summary
# tables, thinned posterior draws, per-stream and frozen-fit comparisons and a
# copy of the input data into `output/` at the repo root. The Release workflow
# bundles that directory into a GitHub Release on each push to `main`. Both
# pages share `docs/examples/_setup.jl`, which is loaded once per session, so
# the second include reuses the first's fitted chains rather than refitting.

using BVDOutbreakSize

const REPO_ROOT = pkgdir(BVDOutbreakSize)

include(joinpath(REPO_ROOT, "docs", "examples", "analysis.jl"))
include(joinpath(REPO_ROOT, "docs", "examples", "sensitivity.jl"))
