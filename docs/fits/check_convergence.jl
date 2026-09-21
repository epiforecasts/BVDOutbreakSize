# Gate the build on the fits having converged, from the same cached chains the
# render reads. Run after the fit matrix:
#
#   BVD_FIT_CACHE=logs/fit_cache julia --project=docs docs/fits/check_convergence.jl
#
# Select the fits with `BVD_CONVERGENCE_IDS` (comma separated, default
# `joint`), the cache directory with `BVD_FIT_CACHE`, and the thresholds with
# the `BVD_CONVERGENCE_FAIL_*` / `BVD_CONVERGENCE_WARN_*` variables
# `convergence.jl` documents.
#
# The report is written to stdout, to the GitHub Actions job summary and, when
# `BVD_CONVERGENCE_REPORT` names a path, to that file for the workflow's PR
# comment. The process exits non-zero when any gated fit failed, so the build
# goes red; nothing here touches the docs build, which runs in its own jobs and
# publishes its preview either way.
using Pkg: Pkg
Pkg.instantiate()

using BVDOutbreakSize
include(joinpath(@__DIR__, "registry.jl"))
include(joinpath(@__DIR__, "convergence.jl"))

const IDS = (
    ids = strip(get(ENV, "BVD_CONVERGENCE_IDS", ""));
    isempty(ids) ? CONVERGENCE_IDS :
        [strip(s) for s in split(ids, ",") if !isempty(strip(s))]
)
const CACHE = get(
    ENV, "BVD_FIT_CACHE",
    joinpath(pkgdir(BVDOutbreakSize), "logs", "fit_cache")
)
const REPORT = strip(get(ENV, "BVD_CONVERGENCE_REPORT", ""))
## Lets the workflow find the comment it posted last time and edit that one
## rather than leaving a new comment per build.
const MARKER = "<!-- bvd-fit-convergence -->"

## The frozen fits hand back `(; cutoff, o, chn)` rather than the chain.
_chain(x) = x isa NamedTuple && haskey(x, :chn) ? x.chn : x

## Where this verdict came from, so a comment read weeks later still says
## which commit and which run produced it.
function _context()
    sha = get(ENV, "GITHUB_SHA", "")
    isempty(sha) && return ""
    server = get(ENV, "GITHUB_SERVER_URL", "https://github.com")
    repo = get(ENV, "GITHUB_REPOSITORY", "")
    run = get(ENV, "GITHUB_RUN_ID", "")
    line = "Commit `$(first(sha, 7))`"
    isempty(run) || (line *= " — [build log]($server/$repo/actions/runs/$run)")
    return line * "."
end

checks = map(IDS) do id
    ## `strict` because the fit matrix has already produced every fit and the
    ## workflow downloads them here: a miss is a wiring bug, and refitting the
    ## joint model inside the gate would take hours.
    result = fit_or_load(
        fit_key(id), () -> nothing;
        cache_dir = CACHE, strict = true
    )
    fit_convergence(id, _chain(result))
end

md = convergence_markdown(checks; marker = MARKER, context = _context())
print(stdout, md)
for var in ("GITHUB_STEP_SUMMARY",)
    path = get(ENV, var, "")
    isempty(path) || open(io -> print(io, md), path, "a")
end
isempty(REPORT) || (mkpath(dirname(REPORT)); write(REPORT, md))

status = any(c -> c.status === :fail, checks) ? "fail" :
    (any(c -> c.status === :warn, checks) ? "warn" : "pass")
out = get(ENV, "GITHUB_OUTPUT", "")
isempty(out) || open(io -> println(io, "status=", status), out, "a")

if status == "fail"
    @error "gated fits have not converged" ids = IDS
    exit(1)
end
