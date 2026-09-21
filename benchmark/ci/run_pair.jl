# Benchmark two revisions of this repository, one after the other, in one
# process on one machine.
#
# Usage, from the repository root:
#   julia --project=benchmark/ci benchmark/ci/run_pair.jl \
#       <main-worktree> <pr-worktree> <results-dir>
#
# Both arguments are working trees, not revisions, and that is the reason
# this file exists rather than a plain call to `benchpkg`. AirspeedVelocity's
# other entry points obtain a revision by handing it to `Pkg.add` as a git
# tree, and Pkg cannot check out a tree that declares a submodule:
#
#   GitError(Code:ERROR, Class:Submodule,
#            cannot get submodules without a working tree)
#
# This repository declares `external/bdbv-linelist-analysis`, so every route
# through `Pkg.add` fails, the `action-v1` action included, whether the
# package is named by url or by path. `rev = "dirty"` is AirspeedVelocity's
# supported local mode: it calls `Pkg.develop` on a path, which never clones
# and so never reaches the submodule. The workflow materialises the two
# revisions with `git worktree add`, where the submodule is left
# uninitialised because nothing under `src/` reads it, and this script
# benchmarks each tree in turn.
#
# One process on one machine is the point of the change. Each revision used
# to be timed in its own CI job, so every reported ratio divided one hosted
# runner's speed by another's. That pool is heterogeneous by about a factor
# of two and the difference landed whole on the ratio.
using AirspeedVelocity: benchmark
using Pkg: PackageSpec

const PACKAGE = "BVDOutbreakSize"

## AirspeedVelocity names the results file after the revision, and both arms
## are `dirty`, so each one is renamed as soon as it is written.
const DIRTY_RESULTS = "results_$(PACKAGE)@dirty.json"

## Package-load samples per arm. Each one relaunches Julia and loads the
## package, so the default of five is minutes of wall clock for a number
## this comment reads as context rather than as a result.
const LOAD_SAMPLES = 3

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
    ## Fixed order, `main` then the pull request. Compilation is outside
    ## the timed region either way, since AirspeedVelocity precompiles the
    ## environment and BenchmarkTools warms up before it samples, but the
    ## order is pinned so two runs of the same pair are comparable.
    for (label, worktree) in (("main", main_worktree), ("pr", pr_worktree))
        @info "benchmarking $label from $worktree"
        println(run_arm(label, worktree, results_dir))
    end
    return nothing
end

main(ARGS)
