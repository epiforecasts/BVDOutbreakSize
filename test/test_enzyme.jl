## Tests for the Enzyme AD extension (`ext/BVDOutbreakSizeEnzymeExt.jl`).
## Loading Enzyme activates `enzyme_adtype()` and the EnzymeRules for the
## gamma CDF. The gradient of a composer's unconstrained log-density under
## Enzyme must match Mooncake (the package default) and, through it,
## finite differences. Tagged `:slow` for the one-off Enzyme compilation.

@testitem "enzyme_adtype is an AutoEnzyme with runtime activity" tags=[:slow, :ad] begin
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

@testitem "Enzyme gradient matches Mooncake on a single-stream model" tags=[:slow, :ad] begin
    using Enzyme
    using Mooncake
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: exports_only_model, default_adtype, enzyme_adtype

    ## With growth sampled as the rate `r`, the exports likelihood
    ## differentiates through `expm1(r·Δ)/r` and an `r`-capturing closure.
    ## Enzyme mis-handles that path on Julia LTS (1.10) and pre-release,
    ## disagreeing with Mooncake (the default backend, which matches
    ## finite differences and is correct on every version). Restrict this
    ## opt-in cross-AD check to stable Julia where Enzyme is reliable for
    ## this model; see issue #153.
    if VERSION < v"1.11" || !isempty(VERSION.prerelease)
        @test_skip "Enzyme gradient unreliable on Julia LTS / pre (#153)"
    else
        model = exports_only_model(3)
        seed!(20260518)
        vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
        x0 = collect(vi[:])

        grad(adtype) = last(logdensity_and_gradient(
            DynamicPPL.LogDensityFunction(
                model, DynamicPPL.getlogjoint, vi; adtype = adtype), x0))

        g_mooncake = grad(default_adtype())
        g_enzyme = grad(enzyme_adtype())

        @test g_enzyme ≈ g_mooncake rtol=1e-6
    end
end
