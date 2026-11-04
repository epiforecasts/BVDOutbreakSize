## Tests for the bound on the cryptic phase's share of the population. The
## cryptic exponential phase has no depletion, so `infection_model` and
## `patch_infection_model` give zero density to a seed whose cryptic
## infections exceed `max_cryptic_share` of the population, and leave the
## density unchanged below it.

@testitem "cryptic_infections: geometric sum back from the seed" begin
    using BVDOutbreakSize: cryptic_infections

    C_T, r = 40.0, 0.06
    @test cryptic_infections(C_T, r) ≈
        sum(C_T * exp(-r * k) for k in 0:100_000)
    @test cryptic_infections(2C_T, r) ≈ 2cryptic_infections(C_T, r)
end

@testitem "infection_model: zero density past the cryptic share" begin
    using BVDOutbreakSize: infection_model
    using Turing: fix, logjoint, @varname
    using Random: Xoshiro

    n = 60
    ## r and m at the prior centre, and at the stuck onsets chain (a seed
    ## of about 1e7 a day against 28 million people).
    ordinary = Dict(
        @varname(growth_state.r) => 0.06, @varname(growth_state.m) => 2.75
    )
    depleting = Dict(
        @varname(growth_state.r) => 0.16, @varname(growth_state.m) => 7.4
    )
    bounded = infection_model(n; rt_start = 10)
    unbounded = infection_model(n; rt_start = 10, max_cryptic_share = Inf)

    θ = rand(Xoshiro(1), fix(bounded, ordinary))
    lj = logjoint(fix(bounded, ordinary), θ)
    @test isfinite(lj)
    @test lj == logjoint(fix(unbounded, ordinary), θ)

    θ = rand(Xoshiro(1), fix(bounded, depleting))
    @test logjoint(fix(bounded, depleting), θ) == -Inf
    @test isfinite(logjoint(fix(unbounded, depleting), θ))
end

@testitem "patch_infection_model: zero density past the cryptic share" begin
    using BVDOutbreakSize: patch_infection_model
    using Turing: fix, logjoint, @varname
    using Random: Xoshiro

    n = 60
    ordinary = Dict(
        @varname(growth_state.r) => 0.06, @varname(growth_state.m) => 2.75
    )
    depleting = Dict(
        @varname(growth_state.r) => 0.16, @varname(growth_state.m) => 7.4
    )
    bounded = patch_infection_model(n, 2; rt_start = 10)
    unbounded = patch_infection_model(
        n, 2; rt_start = 10, max_cryptic_share = Inf
    )

    θ = rand(Xoshiro(1), fix(bounded, ordinary))
    lj = logjoint(fix(bounded, ordinary), θ)
    @test isfinite(lj)
    @test lj == logjoint(fix(unbounded, ordinary), θ)

    θ = rand(Xoshiro(1), fix(bounded, depleting))
    @test logjoint(fix(bounded, depleting), θ) == -Inf
    @test isfinite(logjoint(fix(unbounded, depleting), θ))
end
