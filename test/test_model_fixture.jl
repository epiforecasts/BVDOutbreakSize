## The headline joint against `test/model_fixture/joint.toml`, which
## `test/model_fixture/generate.jl` writes from `main`'s code. A change that
## alters the model's log density, its gradient or a recorded `:=` quantity
## fails here; floating-point reassociation stays well inside the tolerance.
## Tagged `:ad`, since the gradient needs a cold Mooncake compile of the joint.

@testitem "model fixture: the headline joint matches its fixture" tags = [
    :ad,
] begin
    using BVDOutbreakSize: freeze_observations, production_joint,
        default_breakpoint, default_adtype
    using Turing: DynamicPPL
    using .DynamicPPL: VarInfo, ParamsWithStats, LogDensityFunction, link,
        getlogjoint_internal, get_all_ranges_and_transforms
    using LogDensityProblems: logdensity_and_gradient
    using TOML: TOML
    include(joinpath(@__DIR__, "model_fixture", "model.jl"))

    fixture = TOML.parsefile(FIXTURE_PATH)
    model = fixture_joint()
    vi = link(VarInfo(model), model)
    ldf = fixture_ldf(model, vi)
    ranges = varname_ranges(ldf)
    sampled = Set(keys(ranges))

    rtol = 1.0e-8
    for p in fixture["points"]
        ## The same sampled varnames, each the same length.
        @test sampled == Set(keys(p["x"]))
        x = zeros(length(vi[:]))
        for (k, r) in ranges
            x[r] = p["x"][k]
        end
        lp, g = logdensity_and_gradient(ldf, x)
        @test isapprox(lp, p["logdensity"]; rtol)
        ## Each gradient component to `rtol`, with a floor scaled to the
        ## largest component for entries near zero.
        gscale = maximum(v -> maximum(abs, v), values(p["gradient"]))
        grad_off = filter(collect(keys(ranges))) do k
            gref = p["gradient"][k]
            !all(
                @. abs(g[ranges[k]] - gref) <=
                    rtol * abs(gref) + 1.0e-10 * gscale
            )
        end
        @test isempty(grad_off)
        reported, stats = reported_values(ldf, x, sampled)
        @test isapprox(stats.logjoint, p["logjoint"]; rtol)
        ## Every quantity main records, to `rtol`. A vector is stored as its
        ## sum, which is near zero for a centred one, hence the floor.
        reported_off = filter(collect(keys(p["reported"]))) do k
            !isapprox(
                get(reported, k, NaN), p["reported"][k];
                rtol, atol = 1.0e-10, nans = true
            )
        end
        @test isempty(reported_off)
    end
end
