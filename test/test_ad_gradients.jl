## AD-gradient smoke check: the package default backend (Mooncake) must be
## able to differentiate the models — the one property the NUTS sampler
## actually relies on. A single unconstrained-space log-density gradient per
## component is the minimal, fast way to assert that, avoiding a full NUTS
## fit as a differentiability check: that takes ~25 min and can flake during
## sampler adaptation rather than in the gradient itself.
##
## The components come from `test/ad_fixtures.jl`, which
## `benchmark/benchmarks.jl` includes too. One source means a component
## cannot be benchmarked without also being asserted differentiable, or
## asserted without being timed. Adding a model there adds it to both.
##
## The full `bvd_joint` is not a component here: one gradient sits behind a
## cold Mooncake compile of roughly 18 minutes. Its gradient in the suite is
## the production joint with and without the rules in
## `test/test_mooncake_rules.jl`.
## The components below are the submodels, the latent process and every
## single-stream composer, so the joint is a composition of surfaces each
## asserted differentiable here. Its own gradient is timed
## by the benchmark suite and exercised end to end by the NUTS fits the docs
## build runs.
##
## Tagged `:ad`, the tag the downgrade-compat run skips (AD gradients drift
## below the package's pinned dependency versions), matching the other
## AD-sensitive items. The gradient pattern mirrors the Enzyme-extension
## check in `test/enzyme/runtests.jl`, which validates the same models.

@testitem "AD gradient: every component differentiates (Mooncake)" tags = [
    :ad,
] begin
    using LogDensityProblems: logdensity_and_gradient
    using BVDOutbreakSize: default_adtype
    include(joinpath(@__DIR__, "ad_fixtures.jl"))

    scenarios = ADFixtures.scenarios()
    @test length(scenarios) >= ADFixtures.MIN_SCENARIOS

    for scen in scenarios
        @testset "$(scen.group): $(scen.name)" begin
            ldf, x = ADFixtures.log_density_function(scen, default_adtype())
            logp, grad = logdensity_and_gradient(ldf, x)
            @test isfinite(logp)
            @test length(grad) == length(x)
            @test all(isfinite, grad)
            ## A non-trivial gradient (not identically zero) confirms AD
            ## actually ran.
            @test any(!iszero, grad)
        end
    end
end

@testitem "AD gradient: province_composition_model differentiates (Mooncake)" tags = [:ad] begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: province_composition_model, default_adtype

    ## The stick-breaking BetaBinomial composition is the likelihood the
    ## spatial data enter through, so its gradient has to exist.
    ##
    ## The reshaping that turns the per-province histories into this matrix
    ## (`province_increment_matrix`) looks provinces up by name in a
    ## `Dict{String}`. It must stay OUTSIDE the model body: a string compare
    ## on the tape is a `memcmp` foreigncall that Mooncake has no rule for,
    ## and it aborts the gradient of the whole joint.
    seed!(20260518)
    obs = [853 21 42; 77 2 5; 3 0 0]
    modelled = [800.0 20.0 40.0; 70.0 2.5 4.0; 2.0 0.1 0.2]
    model = province_composition_model(obs, modelled)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype()
    )
    logp, grad = logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
end
