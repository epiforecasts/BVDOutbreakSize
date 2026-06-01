## Tests for the Enzyme AD backend (`src/enzyme.jl`). Enzyme reverse-mode
## with runtime activity is now the package default, so this is the
## correctness gate: the gradient of a composer's unconstrained
## log-density under Enzyme must match Mooncake and, through it, finite
## differences. Uses a single-stream composer (fast); the full joint is
## exercised by the docs analysis build, not the unit tests. Enzyme
## gradient correctness is validated on stable Julia (the supported
## versions); LTS/nightly are out of scope (#153). Tagged `:slow` for the
## one-off Enzyme compilation.

@testitem "enzyme_adtype is reverse with runtime activity" tags=[:slow, :ad] begin
    using ADTypes: AutoEnzyme
    using Enzyme
    using BVDOutbreakSize: enzyme_adtype
    ad = enzyme_adtype()
    @test ad isa AutoEnzyme
    ## Duplicated closure annotation is a type parameter; the mode is
    ## reverse with runtime activity (needed for the joint's incubation
    ## convolution, which defeats static activity analysis).
    @test ad isa AutoEnzyme{<:Any, Enzyme.Duplicated}
    @test ad.mode === Enzyme.set_runtime_activity(Enzyme.Reverse)
end

@testitem "Enzyme gradient matches Mooncake on a single-stream model" tags=[:slow, :ad] begin
    using Enzyme
    using Mooncake
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: exports_only_model, mooncake_adtype, enzyme_adtype
    using Random: MersenneTwister

    ## Validate Enzyme on stable Julia only. Enzyme is wrong/erroring on
    ## LTS and nightly under the growth reparam (#153), which are not
    ## supported versions, so skip there rather than turn those CI legs red.
    if VERSION < v"1.11" || !isempty(VERSION.prerelease)
        @test_skip true
    else
        model = exports_only_model(3)
        seed!(20260518)
        vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
        x0 = collect(vi[:])

        f_moon = DynamicPPL.LogDensityFunction(
            model, DynamicPPL.getlogjoint, vi; adtype = mooncake_adtype())
        f_enz = DynamicPPL.LogDensityFunction(
            model, DynamicPPL.getlogjoint, vi; adtype = enzyme_adtype())
        grad(f, x) = last(logdensity_and_gradient(f, x))

        ## Check the prior-mode draw plus two perturbed points so the
        ## runtime-activity-off config is exercised away from the mode.
        pts = [x0]
        for k in 1:2
            rng = MersenneTwister(20260518 + k)
            push!(pts, x0 .+ 0.5 .* randn(rng, length(x0)))
        end
        for x in pts
            @test grad(f_enz, x) ≈ grad(f_moon, x) rtol=1e-6
        end
    end
end
