## AD-gradient smoke check: the package default backend (Mooncake) must be
## able to differentiate the models — the one property the NUTS sampler
## actually relies on. A single unconstrained-space log-density gradient is
## the minimal, fast way to assert that, and it replaces the slow full NUTS
## fits that previously stood in as the "the models differentiate" coverage
## (those took ~25 min and flaked during sampler adaptation rather than in
## the gradient itself).
##
## A composer (exports_only_model) exercises the renewal + onset + likelihood
## AD path quickly. The full bvd_joint's gradient is NOT checked here on its
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

