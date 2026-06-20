## AD-gradient smoke checks: the package default backend (Mooncake) must be
## able to differentiate the models — the one property the NUTS sampler
## actually relies on. A single unconstrained-space log-density gradient is
## the minimal, fast way to assert that, and it replaces the slow full NUTS
## fits that previously stood in as the "the models differentiate" coverage
## (those took ~25 min and flaked during sampler adaptation rather than in
## the gradient itself).
##
## Tagged `:ad`, the tag the downgrade-compat run skips (AD gradients drift
## below the package's pinned dependency versions), matching the other
## AD-sensitive items. The gradient pattern mirrors the Enzyme-extension
## check in `test/enzyme/runtests.jl`, which validates the same models.

@testitem "AD gradient: exports_only_model differentiates (Mooncake)" tags = [:ad] begin
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

@testitem "AD gradient: bvd_joint differentiates (Mooncake)" tags = [:ad] begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: bvd_joint, default_adtype

    ## Full joint with every stream conditioned; `confirmed_deaths` is passed
    ## as a count so the model carries no sampled discrete parameter and the
    ## whole posterior is differentiable.
    seed!(20260518)
    model = bvd_joint(40, 2, 18, 905, 0, 27; confirmed_deaths = 5)
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
