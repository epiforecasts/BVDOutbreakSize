## Tests for the Enzyme AD extension (`ext/BVDOutbreakSizeEnzymeExt.jl`).
## Loading Enzyme activates `enzyme_adtype()`; the `SpecialFunctions.gamma`
## EnzymeRule reached by the Beta / NegativeBinomial normalising constants
## is supplied by CensoredDistributions' own Enzyme extension. The renewal
## model differentiates cleanly under Enzyme, so the Enzyme gradient of a
## composer's unconstrained log-density must match Mooncake (the package
## default) for both a single-stream composer and the full joint. Tagged
## `:slow` for the one-off Enzyme compilation. Mooncake remains the default
## backend; Enzyme is the validated opt-in.
##
## The two gradient items run off Windows only: Enzyme reverse-mode
## differentiation crashes the Julia process on Windows with a native
## `EXCEPTION_ACCESS_VIOLATION` in the stack unwinder (an Enzyme/Windows
## issue, not a model issue — the gradients match Mooncake on Linux and
## macOS), and a process-level segfault cannot be caught in-test, so the
## items are skipped there. The adtype-construction item below runs
## everywhere.

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

@testsnippet EnzymeGrad begin
    using Enzyme
    using Mooncake
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: default_adtype, enzyme_adtype

    ## Unconstrained-space gradient of a composer's log-density under a
    ## given AD backend, at a fixed prior draw.
    function adgrad(model, adtype)
        seed!(20260518)
        vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
        x0 = collect(vi[:])
        return last(logdensity_and_gradient(
            DynamicPPL.LogDensityFunction(
                model, DynamicPPL.getlogjoint, vi; adtype = adtype), x0))
    end
end

@testitem "Enzyme gradient matches Mooncake on a single-stream model" tags=[
    :slow, :ad] setup=[EnzymeGrad] begin
    using BVDOutbreakSize: exports_only_model
    if Sys.iswindows()
        @test_skip "Enzyme reverse-mode segfaults on Windows"
    else
        model=exports_only_model(3, 2)
        @test adgrad(model, enzyme_adtype()) ≈
              adgrad(model, default_adtype()) rtol=1e-6
    end
end

@testitem "Enzyme gradient matches Mooncake on the joint" tags=[
    :slow, :ad] setup=[EnzymeGrad] begin
    using BVDOutbreakSize: bvd_joint
    if Sys.iswindows()
        @test_skip "Enzyme reverse-mode segfaults on Windows"
    else
        ## All streams plus the lab pipeline and the intervention breakpoint.
        model=bvd_joint(20, 2, 3, 5, 1, 4, 10; breakpoint = 14)
        @test adgrad(model, enzyme_adtype()) ≈
              adgrad(model, default_adtype()) rtol=1e-6
    end
end
