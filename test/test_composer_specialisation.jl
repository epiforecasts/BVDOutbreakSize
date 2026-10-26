## Julia does not specialise a method on a function argument that its body
## only passes on without calling. A composer hands its `infection` and
## `onset_incidence` submodel constructors straight to `_latent`, so the
## evaluator it compiles can take them as plain `Function`. Every submodel
## call downstream is then a dynamic dispatch that boxes its return value.
## Results stay correct, so only the benchmarks see it: on
## `confirmed_only_model` it cost 1.4 times the log-density allocations and
## 2.2 times the Mooncake gradient time.
##
## Whether the widened evaluator is the one that runs is up to the compiler,
## so the test asserts none exists. If it fails, give the named composer's
## pass-through function keywords a type parameter
## (`infection::I = infection_model` with `where {I}`), as the others have.
@testitem "composers compile specialised on their submodel arguments" tags = [
    :ad,
] begin
    using BVDOutbreakSize
    using LogDensityProblems: logdensity

    include(joinpath(@__DIR__, "ad_fixtures.jl"))
    using .ADFixtures: ADFixtures

    M = BVDOutbreakSize
    composers = filter(s -> s.group == "Composer", ADFixtures.scenarios())
    @test !isempty(composers)

    widened = String[]
    for scen in composers
        ldf, x = ADFixtures.log_density_function(scen, M.default_adtype())
        logdensity(ldf, x)
        ## The `@model` body is the gensym-named keyword-body method
        ## `#<name>#<n>`, not the user-facing constructor.
        prefix = "#" * scen.name * "#"
        for nm in names(M; all = true)
            startswith(string(nm), prefix) || continue
            for m in methods(getfield(M, nm)), mi in Base.specializations(m)
                mi isa Core.MethodInstance || continue
                ## Only the evaluator the log density runs.
                occursin("InitFromVector", string(mi.specTypes)) || continue
                any(p -> p === Function, mi.specTypes.parameters) &&
                    push!(widened, scen.name)
            end
        end
    end

    @test isempty(widened)
    isempty(widened) || @info "composers compiled on a bare Function" widened
end
