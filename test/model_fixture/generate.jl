# Writes `joint.toml` beside this script: the headline joint's log density,
# gradient and recorded `:=` quantities at a few fixed parameter points.
# `test/test_model_fixture.jl` checks the package against it, so a change
# meant to leave the model as it is (a performance refactor) cannot change
# it unnoticed.
#
# Run it against `main`'s code, never the branch under test. From this
# checkout, with a worktree of `main` at <main>:
#
#   julia --project=<main>/test test/model_fixture/generate.jl
#
# `<main>/test` resolves `BVDOutbreakSize` to `<main>`, so the values come
# from `main`'s model and data; the builder is this checkout's `model.jl`.
# Regenerate after a deliberate model change has merged, or if the data at
# `FIXTURE_CUTOFF` is corrected. One gradient needs a cold Mooncake compile
# of the joint, about 20 minutes.

using BVDOutbreakSize: BVDOutbreakSize, freeze_observations, production_joint,
    default_breakpoint, default_adtype
using Turing: DynamicPPL
using .DynamicPPL: VarInfo, InitFromParams, InitFromPrior, ParamsWithStats,
    LogDensityFunction, link, getlogjoint_internal,
    get_all_ranges_and_transforms
using LogDensityProblems: logdensity_and_gradient
using Random: seed!
using Dates: today
using TOML: TOML
import Mooncake

include(joinpath(@__DIR__, "model.jl"))

## The highest-density `NPOINTS` of `NDRAWS` prior draws, so no point sits
## where a likelihood is clamped at its floor.
const NDRAWS = 40
const NPOINTS = 3

model = fixture_joint()
draws = map(1:NDRAWS) do s
    seed!(s)
    pws = ParamsWithStats(InitFromPrior(), model; include_colon_eq = false)
    (; seed = s, params = pws.params, lj = pws.stats.logjoint)
end
points = sort(filter(d -> isfinite(d.lj), draws); by = d -> -d.lj)
points = points[1:NPOINTS]

vi0 = link(VarInfo(model, InitFromParams(points[1].params, nothing)), model)
ldf = fixture_ldf(model, vi0)
ranges = varname_ranges(ldf)
sampled = Set(keys(ranges))

out = map(points) do p
    vi = link(VarInfo(model, InitFromParams(p.params, nothing)), model)
    x = collect(vi[:])
    lp, g = logdensity_and_gradient(ldf, x)
    reported, stats = reported_values(ldf, x, sampled)
    Dict{String, Any}(
        "seed" => p.seed,
        "logdensity" => lp,
        "logjoint" => stats.logjoint,
        "x" => Dict(k => x[r] for (k, r) in ranges),
        "gradient" => Dict(k => g[r] for (k, r) in ranges),
        "reported" => reported,
    )
end

src = pkgdir(BVDOutbreakSize)
meta = Dict{String, Any}(
    "commit" => readchomp(`git -C $src rev-parse --short HEAD`),
    "cutoff" => FIXTURE_CUTOFF,
    "generated" => string(today()),
    "julia" => string(VERSION),
    "dynamicppl" => string(pkgversion(DynamicPPL)),
    "mooncake" => string(pkgversion(Mooncake)),
)
fixture = Dict{String, Any}("meta" => meta, "points" => out)
open(FIXTURE_PATH, "w") do io
    TOML.print(io, fixture; sorted = true)
end
@assert isequal(TOML.parsefile(FIXTURE_PATH)["points"], out) "round trip"
@info "wrote" FIXTURE_PATH meta length(ranges) length(first(out)["reported"])
