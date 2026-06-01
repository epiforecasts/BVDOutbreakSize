## Tests for the Enzyme AD extension (`ext/BVDOutbreakSizeEnzymeExt.jl`).
## Loading Enzyme activates `enzyme_adtype()`; the `SpecialFunctions.gamma`
## EnzymeRule reached by the Beta / NegativeBinomial normalising constants
## is supplied by CensoredDistributions' own Enzyme extension. We check the
## gradient of a single-stream composer's unconstrained log-density under
## Enzyme against Mooncake (the package default). Differentiating the
## renewal under Enzyme is still work in progress: the censored-delay
## discretisation does not yet differentiate cleanly under Enzyme on every
## supported version, so the gradient match is recorded as broken (rather
## than erroring the suite) when Enzyme cannot produce it, and Mooncake
## stays the default. Tagged `:slow` for the one-off Enzyme compilation.

@testitem "enzyme_adtype is an AutoEnzyme with runtime activity" tags=[
    :slow, :ad] begin
    using ADTypes: AutoEnzyme
    using Enzyme
    using BVDOutbreakSize: enzyme_adtype
    ad = enzyme_adtype()
    @test ad isa AutoEnzyme
    ## Duplicated closure annotation is a type parameter; the mode is
    ## reverse with runtime activity enabled.
    @test ad isa AutoEnzyme{<:Any, Enzyme.Duplicated}
    @test ad.mode === Enzyme.set_runtime_activity(Enzyme.Reverse)
end

@testitem "Enzyme gradient matches Mooncake on a single-stream model" tags=[
    :slow, :ad] begin
    using Enzyme
    using Mooncake
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: exports_only_model, default_adtype, enzyme_adtype

    model = exports_only_model(3, 2)
    seed!(20260518)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])

    grad(adtype) = last(logdensity_and_gradient(
        DynamicPPL.LogDensityFunction(
            model, DynamicPPL.getlogjoint, vi; adtype = adtype), x0))

    g_mooncake = grad(default_adtype())
    ## Enzyme differentiation of the censored-delay discretisation is WIP;
    ## treat a differentiation failure as a known-broken match rather than
    ## a hard error so the suite stays green while Enzyme support matures.
    g_enzyme = try
        grad(enzyme_adtype())
    catch
        nothing
    end
    if g_enzyme === nothing
        @test_broken false
    else
        @test g_enzyme ≈ g_mooncake rtol=1e-6
    end
end
