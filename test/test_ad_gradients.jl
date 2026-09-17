## AD-gradient smoke check: the package default backend (Mooncake) must be
## able to differentiate the models — the one property the NUTS sampler
## actually relies on. A single unconstrained-space log-density gradient per
## component is the minimal, fast way to assert that, avoiding a full NUTS
## fit as a differentiability check: that takes ~25 min and can flake during
## sampler adaptation rather than in the gradient itself.
##
## The components come from the `ADFixtures` path package at
## `test/ADFixtures`, which also drives `benchmark/src/ad_gradients.jl`. One
## source means a component cannot be benchmarked without also being
## asserted differentiable, or asserted without being timed. Adding a model
## there adds it to both.
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
    using ADFixtures
    using LogDensityProblems: logdensity_and_gradient
    using BVDOutbreakSize: default_adtype

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
