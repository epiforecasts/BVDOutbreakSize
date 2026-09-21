# Entry point for regenerating the published results.
#
# Runs every report page, which fits the models (loaded from the
# content-addressed cache when present) and writes the summary tables, thinned
# posterior draws, forecasts, per-stream and frozen-fit comparisons and a copy
# of the input data into `output/` at the repo root. The Release workflow
# bundles that directory into a GitHub Release on each push to `main`. The
# pages share `docs/pages/_setup.jl`, which is loaded once per session, so
# each include after the first reuses the fitted chains rather than refitting.

using BVDOutbreakSize

const REPO_ROOT = pkgdir(BVDOutbreakSize)

for page in
    ("analysis", "province", "forecast", "evaluation", "sensitivity")
    include(joinpath(REPO_ROOT, "docs", "pages", "$page.jl"))
end
