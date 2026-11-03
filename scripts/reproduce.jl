#!/usr/bin/env julia
#
# Standalone reproduction runner. Fetches the package and re-fits the
# models without a manual `git clone`. Run it straight from the web:
#
#   curl -fsSL https://raw.githubusercontent.com/epiforecasts/BVDOutbreakSize/main/scripts/reproduce.jl | julia
#
# (Reading a script before piping it to a shell is good practice.)
#
# Outputs are written to `./bvd-output` by default; set
# `BVD_OUTPUT_DIR` to write them elsewhere. Set `BVD_REF` to a release
# tag or branch to reproduce a specific version (defaults to `main`).
#
# The analysis imports several of the package's dependencies directly,
# so it runs against the package's own project. We clone into a
# temporary directory and activate that: it is writable (the analysis
# writes its outputs into the package directory), always fresh, and
# leaves the caller's own Julia environments untouched. Requires `git`
# on the PATH.

using Pkg

const REPO_URL = "https://github.com/epiforecasts/BVDOutbreakSize"

ref = get(ENV, "BVD_REF", "main")
output_dir = abspath(
    get(
        ENV, "BVD_OUTPUT_DIR",
        joinpath(pwd(), "bvd-output")
    )
)

src = mktempdir()
run(`git clone --depth 1 --branch $ref $REPO_URL $src`)

Pkg.activate(src)
Pkg.instantiate()

withenv("BVD_OUTPUT_DIR" => output_dir) do
    include(joinpath(src, "scripts", "run.jl"))
end

@info "Reproduction complete" output_dir
