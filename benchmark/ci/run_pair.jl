# Benchmark two revisions of this repository, one after the other, in one
# process.
#
# Usage, from the repository root:
#   julia --project=benchmark/ci benchmark/ci/run_pair.jl \
#       <main-worktree> <pr-worktree> <results-dir>
#
# Both arguments are working trees, not revisions. AirspeedVelocity's other
# entry points obtain a revision through `Pkg.add` on a git tree, which
# cannot check out a tree that declares a submodule, and this repository
# declares one. `rev = "dirty"` is its local mode and calls `Pkg.develop` on
# a path, so nothing is cloned. `benchmark/README.md` carries the reasoning.
#
# One process times both arms on one machine, so a ratio compares two
# revisions rather than two hosted runners, which differ by about a factor
# of two.
using AirspeedVelocity: benchmark
using Pkg: PackageSpec

const PACKAGE = "BVDOutbreakSize"

## AirspeedVelocity names the results file after the revision, and both arms
## are `dirty`, so each one is renamed as soon as it is written.
const DIRTY_RESULTS = "results_$(PACKAGE)@dirty.json"

## Package-load samples per arm. AirspeedVelocity always injects one
## in-process load benchmark; anything beyond the first relaunches Julia to
## take another. `comment.jl` drops the row either way, because the warmup
## BenchmarkTools runs before it samples performs the `using` and leaves the
## in-process sample timing a warm re-import. So ask for one and do not pay
## for the relaunches.
const LOAD_SAMPLES = 1

"""
    run_arm(label, worktree, results_dir) -> String

Benchmark `worktree` and write `results_<PACKAGE>@<label>.json` into
`results_dir`, returning that path.

The suite comes from the tree being timed, not from one revision applied to
both. So a pull request's own changes to `benchmark/` and to
`test/ad_fixtures.jl` are exercised by that pull request, and a component
present on one side only is reported as added or removed rather than
failing the other arm. `benchmark/Project.toml` is taken from the same tree,
which is what keeps the two arms resolving the versions each revision pins.
"""
function run_arm(
        label::AbstractString, worktree::AbstractString,
        results_dir::AbstractString
    )
    worktree = abspath(worktree)
    results_dir = abspath(results_dir)
    script = joinpath(worktree, "benchmark", "benchmarks.jl")
    project = joinpath(worktree, "benchmark", "Project.toml")
    out = joinpath(results_dir, "results_$(PACKAGE)@$(label).json")
    if !isfile(script)
        ## A revision from before the suite existed. An empty group reports
        ## every benchmark as new rather than failing the run.
        @info "no benchmark suite at $label, writing an empty result"
        write(out, "{\"data\":{}}")
        return out
    end
    benchmark(
        PackageSpec(; name = PACKAGE, path = worktree, rev = "dirty");
        output_dir = results_dir,
        script = script,
        project_toml = project,
        nsamples_load_time = LOAD_SAMPLES,
    )
    mv(joinpath(results_dir, DIRTY_RESULTS), out; force = true)
    return out
end

function main(args)
    length(args) == 3 ||
        error("usage: run_pair.jl <main-worktree> <pr-worktree> <results-dir>")
    main_worktree, pr_worktree, results_dir = args
    mkpath(results_dir)
    ## Fixed order, `main` then the pull request, so two runs of the same
    ## pair are comparable. Compilation is outside the timed region either
    ## way: AirspeedVelocity precompiles the environment and BenchmarkTools
    ## warms up before it samples.
    for (label, worktree) in (("main", main_worktree), ("pr", pr_worktree))
        @info "benchmarking $label from $worktree"
        println(run_arm(label, worktree, results_dir))
    end
    return nothing
end

main(ARGS)
