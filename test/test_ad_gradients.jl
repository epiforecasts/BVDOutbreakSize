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
## The full `bvd_joint` is not a component here: differentiating it under
## Mooncake takes ~10 min and is unstable on the Julia LTS runner (the same
## path that makes the `bvd_joint` NUTS fits flaky there), so a dedicated
## check would reintroduce exactly the slowness and flakiness this file
## removes. The joint's gradient is still exercised end-to-end whenever the
## per-vintage predict/fit tests sample it (test_vintage_predict,
## test_lab_pipeline).
##
## Tagged `:ad`, the tag the downgrade-compat run skips (AD gradients drift
## below the package's pinned dependency versions), matching the other
## AD-sensitive items. The gradient pattern mirrors the Enzyme-extension
## check in `test/enzyme/runtests.jl`, which validates the same models.

@testitem "AD gradient: every component differentiates (Mooncake)" tags=[
    :ad
] begin
    using LogDensityProblems: logdensity_and_gradient
    using BVDOutbreakSize: default_adtype
    include(joinpath(@__DIR__, "ad_fixtures.jl"))

    scenarios = ADFixtures.scenarios()
    ## The fixtures are the benchmark suite's component list too, so a
    ## truncated list would pass this testset while quietly shrinking the
    ## benchmarks.
    @test length(scenarios) >= 16

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
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype())
    logp, grad = logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
end

@testitem "AD gradient: the testing covariate differentiates (Mooncake)" tags = [:ad] begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: province_composition_model, default_adtype

    ## The covariate samples `β_asc` and adds a fused broadcast over the
    ## patches, and that branch is taken only when the covariate is not all
    ## zeros. The test above passes the zero default, so it differentiates
    ## the model without the coefficient and leaves this path unproven.
    seed!(20260518)
    obs = [853 21 42; 77 2 5; 3 0 0]
    modelled = [800.0 20.0 40.0; 70.0 2.5 4.0; 2.0 0.1 0.2]
    covariate = [0.64, -0.23, -0.41]
    model = province_composition_model(obs, modelled;
        testing_covariate = covariate)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    ## The coefficient is in the parameter vector only on this path, so the
    ## gradient below covers it.
    @test any(k -> occursin("β_asc", string(k)), keys(vi))
    x0 = collect(vi[:])
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype())
    logp, grad = logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
end
