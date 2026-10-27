## Julia does not specialise a method on a function argument that its body
## only passes on without calling. A composer hands its `infection` and
## `onset_incidence` submodel constructors straight to `_latent`, so the
## evaluator it compiles can take them as plain `Function`. Every submodel
## call downstream is then a dynamic dispatch that boxes its return value.
## Results stay correct, so only the allocations and the timings see it.
##
## The bounds sit between the bytes one log-density evaluation allocates
## with the composer's pass-through keywords typed (`where {I, O}`) and
## without: exports 11.1 kB against 15.0 kB, deaths 13.1 kB against 17.1 kB,
## cases 10.4 kB against 15.0 kB and confirmed 18.1 kB against 26.8 kB, on
## Julia 1.13. The treatment and onsets composers allocate the same either
## way, so they are not checked. If a bound fails, first check that the
## named composer types `infection::I` and `onset_incidence::O` with
## `where {I, O}`. If it does, the model allocates more for another reason
## and the bound needs moving to the new midpoint.
@testitem "composers compile specialised on their submodel arguments" tags = [
    :ad,
] begin
    using BVDOutbreakSize
    using LogDensityProblems: logdensity

    include(joinpath(@__DIR__, "ad_fixtures.jl"))
    using .ADFixtures: ADFixtures

    bounds = Dict(
        "exports_only_model" => 13_000,
        "deaths_only_model" => 15_000,
        "cases_only_model" => 12_500,
        "confirmed_only_model" => 22_000,
    )
    composers = filter(
        s -> haskey(bounds, s.name), ADFixtures.scenarios()
    )
    @test length(composers) == length(bounds)

    ## Inside a function so the measured call sees concrete argument types.
    function log_density_bytes(ldf, x)
        logdensity(ldf, x)
        return minimum(@allocated(logdensity(ldf, x)) for _ in 1:5)
    end
    for scen in composers
        ldf, x = ADFixtures.log_density_function(
            scen, BVDOutbreakSize.default_adtype()
        )
        bytes = log_density_bytes(ldf, x)
        @test bytes <= bounds[scen.name]
        bytes <= bounds[scen.name] ||
            @info "composer log density allocates over its bound" scen.name bytes
    end
end
