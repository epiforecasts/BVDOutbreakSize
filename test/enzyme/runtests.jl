## Enzyme AD extension checks, isolated in their own environment because
## Enzyme is platform- and version-dependent for this model and must not
## enter the main test environment. Run as a tolerated subprocess by
## `test/package/EnzymeExt.jl`.
##
## Mooncake is the default and is asserted to differentiate every model in
## the main suite. This checks the Enzyme opt-in against it over the shared
## AD scenarios from `test/ad_fixtures.jl`, which also declare what is
## expected to fail and what is too slow to run.
##
## Exit code 2 means a scenario did not behave as declared; any other
## non-zero code means the script could not run. The wrapper fails on the
## first and tolerates the second.

using Test
using ADTypes: AutoEnzyme
using Enzyme
using Mooncake
using Turing: DynamicPPL
using LogDensityProblems: logdensity_and_gradient
using Random: seed!
using BVDOutbreakSize: default_adtype, enzyme_adtype,
    exports_only_model, bvd_joint

include(joinpath(@__DIR__, "..", "ad_fixtures.jl"))

## Unconstrained-space gradient of a model's log-density under a given
## AD backend, at a fixed prior draw.
function adgrad(model, adtype)
    seed!(20260518)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])
    return last(
        logdensity_and_gradient(
            DynamicPPL.LogDensityFunction(
                model, DynamicPPL.getlogjoint, vi; adtype = adtype
            ), x0
        )
    )
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
        isapprox(g_enzyme, g_mooncake; rtol = 1.0e-6)
end

## Assert one scenario. `@test_broken` is reached only on a real failure,
## so over-listing a name is safe and under-listing fails.
function check_scenario(name, matched, broken)
    return if matched
        @test matched
    elseif name in broken
        @test_broken matched
    else
        @test matched
    end
end

## Exit code 2 means a scenario did not behave as declared; any other
## non-zero code means the script could not run. `test/package/EnzymeExt.jl`
## fails on the first and tolerates the second.
failed = try
    @testset "Enzyme extension" begin
        @testset "enzyme_adtype is an AutoEnzyme with runtime activity" begin
            ad = enzyme_adtype()
            @test ad isa AutoEnzyme
            @test ad isa AutoEnzyme{<:Any, Enzyme.Duplicated}
            @test ad.mode === Enzyme.set_runtime_activity(Enzyme.Reverse)
        end

        @testset "gradient matches Mooncake on a single-stream model" begin
            @test enzyme_matches_mooncake(exports_only_model(3, 2))
        end

        @testset "every AD component matches Mooncake" begin
            broken = ADFixtures.enzyme_broken_scenarios()
            skipped = ADFixtures.enzyme_skip_scenarios()
            scenarios = ADFixtures.scenarios()
            @test length(scenarios) >= ADFixtures.MIN_SCENARIOS
            for scen in scenarios
                @testset "$(scen.group): $(scen.name)" begin
                    if scen.name in skipped
                        @test_skip "declared too slow to run under Enzyme"
                    else
                        check_scenario(
                            scen.name, enzyme_matches_mooncake(scen.model), broken
                        )
                    end
                end
            end
        end

        @testset "gradient matches Mooncake on the joint" begin
            check_scenario(
                "bvd_joint",
                enzyme_matches_mooncake(
                    bvd_joint(20, 2, 3, 5, 1, 4, 10; breakpoint = 14)
                ),
                ADFixtures.enzyme_broken_scenarios()
            )
        end
    end
    false
catch err
    ## A failing `@testset` throws this at the end of the block. Anything
    ## else is the script failing to run at all, which is the case the
    ## wrapper tolerates, so it propagates and leaves the exit code to
    ## Julia.
    err isa Test.TestSetException || rethrow()
    true
end

exit(failed ? 2 : 0)
