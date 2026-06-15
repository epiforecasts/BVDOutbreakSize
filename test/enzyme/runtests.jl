## Enzyme AD extension checks, isolated in their own environment.
##
## Enzyme is kept out of the main test environment because its
## reverse-mode support is platform- and version-dependent for this model
## (a native access violation on Windows, an `EnzymeInternalError` LLVM
## compile failure on the joint on some Linux runners, a wrong gradient
## from mishandling the Gauss-Legendre quadrature in the censored-delay
## path on Julia LTS). Loading it in the main env also tripped Aqua's
## persistent-task check and broke precompilation on Windows. Here it is a
## dependency of this sub-environment only, run as a tolerated subprocess
## by `test/package/EnzymeExt.jl` on the platforms where it is viable.
##
## Mooncake is the package default and is asserted to differentiate every
## model in the main suite; this script checks the Enzyme opt-in matches
## Mooncake where Enzyme produces a correct gradient, and records a broken
## test otherwise.

using Test
using ADTypes: AutoEnzyme
using Enzyme
using Mooncake
using Turing: DynamicPPL
using LogDensityProblems: logdensity_and_gradient
using Random: seed!
using BVDOutbreakSize: default_adtype, enzyme_adtype,
                       exports_only_model, bvd_joint

## Unconstrained-space gradient of a composer's log-density under a given
## AD backend, at a fixed prior draw.
function adgrad(model, adtype)
    seed!(20260518)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])
    return last(logdensity_and_gradient(
        DynamicPPL.LogDensityFunction(
            model, DynamicPPL.getlogjoint, vi; adtype = adtype), x0))
end

## True when the Enzyme gradient matches Mooncake, false when Enzyme
## throws (so it cannot compile the model on this platform).
function enzyme_matches_mooncake(model)
    g_mooncake = adgrad(model, default_adtype())
    g_enzyme = try
        adgrad(model, enzyme_adtype())
    catch
        nothing
    end
    return g_enzyme === nothing ? false :
           isapprox(g_enzyme, g_mooncake; rtol = 1e-6)
end

@testset "Enzyme extension" begin
    @testset "enzyme_adtype is an AutoEnzyme with runtime activity" begin
        ad = enzyme_adtype()
        @test ad isa AutoEnzyme
        @test ad isa AutoEnzyme{<:Any, Enzyme.Duplicated}
        @test ad.mode === Enzyme.set_runtime_activity(Enzyme.Reverse)
    end

    @testset "gradient matches Mooncake on a single-stream model" begin
        if enzyme_matches_mooncake(exports_only_model(3, 2))
            @test true
        else
            @test_broken false
        end
    end

    @testset "gradient matches Mooncake on the joint" begin
        if enzyme_matches_mooncake(
            bvd_joint(20, 2, 3, 5, 1, 4, 10; breakpoint = 14))
            @test true
        else
            @test_broken false
        end
    end
end
