## AD-gradient smoke check: the package default backend (Mooncake) must be
## able to differentiate the models — the one property the NUTS sampler
## actually relies on. A single unconstrained-space log-density gradient is
## the minimal, fast way to assert that, avoiding a full NUTS fit as a
## differentiability check: that takes ~25 min and can flake during sampler
## adaptation rather than in the gradient itself.
##
## A composer (exports_only_model) exercises the renewal + onset + likelihood
## AD path quickly. The full bvd_joint's gradient is not checked here on its
## own: differentiating it under Mooncake takes ~10 min and is unstable on the
## Julia LTS runner (the same path that makes the bvd_joint NUTS fits flaky
## there), so a dedicated check would reintroduce exactly the slowness and
## flakiness this file removes. The joint's gradient is still exercised
## end-to-end whenever the per-vintage predict/fit tests sample it
## (test_vintage_predict, test_lab_pipeline).
##
## Tagged `:ad`, the tag the downgrade-compat run skips (AD gradients drift
## below the package's pinned dependency versions), matching the other
## AD-sensitive items. The gradient pattern mirrors the Enzyme-extension
## check in `test/enzyme/runtests.jl`, which validates the same model.

@testitem "AD gradient: exports_only_model differentiates (Mooncake)" tags=[
    :ad
] begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: exports_only_model, default_adtype

    ## Link to the unconstrained space and evaluate the log-density gradient
    ## at a fixed prior draw, exactly as the Enzyme check does.
    seed!(20260518)
    model = exports_only_model(40, 2)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype())
    logp, grad = logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    ## A non-trivial gradient (not identically zero) confirms AD actually ran.
    @test any(!iszero, grad)
end

@testitem "AD gradient: patch_infection_model differentiates (Mooncake)" tags = [:ad] begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: patch_infection_model, default_adtype

    ## The multi-patch renewal, the per-patch Rt walk and the implied-national
    ## Rt inversion are all new AD surface. Check the uncoupled default and
    ## the coupled (importation kernel) path, which brings epsilon and the
    ## between-patch term onto the tape.
    for kernel in (zeros(3, 3), [0.0 0.0 0.0; 1e-4 0.0 0.0; 1e-5 0.0 0.0])
        seed!(20260518)
        model = patch_infection_model(60, 3; importation_kernel = kernel)
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
