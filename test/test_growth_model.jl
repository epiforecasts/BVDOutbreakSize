## Tests for the exponential growth submodel from `src/models/priors.jl`.
## The prior is placed on the growth rate `r`; the doubling time `τ`, the
## doubling count `m`, the elapsed time `T` and the cumulative size `C_T`
## are exposed as deterministics. The `r` prior is the exact reciprocal
## pushforward of the previous `τ ~ LogNormal(log 14, 0.4)`.

@testitem "exponential_growth_model exposes r, τ, m, T, C_T" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using BVDOutbreakSize: exponential_growth_model

    chn = sample(MersenneTwister(20260518), exponential_growth_model(),
        Prior(), 200; progress = false)
    r = vec(Array(chn[:r]))
    τ = vec(Array(chn[:τ]))
    m = vec(Array(chn[:m]))
    T = vec(Array(chn[:T]))
    C_T = vec(Array(chn[:C_T]))

    @test all(isfinite, r) && all(>(0), r)
    ## τ is the deterministic reciprocal log(2)/r, not a sampled variable.
    @test all(@. isapprox(τ, log(2) / r; rtol = 1e-8))
    ## T = m·τ and C_T = 2^m as documented.
    @test all(@. isapprox(T, m * τ; rtol = 1e-8))
    @test all(@. isapprox(C_T, 2.0^m; rtol = 1e-8))
end

@testitem "growth r-prior is the exact pushforward of τ" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    import Statistics
    using BVDOutbreakSize: exponential_growth_model

    ## Seed a version-stable MersenneTwister (the default RNG differs
    ## across Julia versions); assert on fast-converging moments rather
    ## than tail quantiles so the test is robust to any residual sampling
    ## variation. log(τ) must be Normal(log 14, 0.4) because r = log(2)/τ
    ## is a reciprocal, which preserves the log-scale SD.
    chn = sample(MersenneTwister(20260518), exponential_growth_model(),
        Prior(), 40_000; progress = false)
    logτ = log.(vec(Array(chn[:τ])))

    @test isapprox(Statistics.mean(logτ), log(14); atol = 0.02)
    @test isapprox(Statistics.std(logτ), 0.4; atol = 0.03)
end

@testitem "m_prior_centre advances with the cut-off date" begin
    using BVDOutbreakSize: m_prior_centre
    ## Base assumption: m = 9 at the 18 May 2026 report date.
    @test m_prior_centre("2026-05-18") ≈ 9.0
    ## Advances by one doubling per 14 days of elapsed time.
    @test m_prior_centre("2026-05-20") ≈ 9.0 + 2 / 14
    @test m_prior_centre("2026-06-01") ≈ 9.0 + 14 / 14
    ## Base value is configurable.
    @test m_prior_centre("2026-05-18"; m_base = 8.0) ≈ 8.0
end

@testitem "growth m-prior is recentred on McCabe central" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean
    using BVDOutbreakSize: exponential_growth_model

    ## Seeded, version-stable RNG; generous tolerance so it is robust to
    ## sampling variation. Centred on m = 9 (C_T = 2^9 = 512); truncation
    ## at 0 nudges the sample mean up only slightly from the location.
    chn = sample(MersenneTwister(20260518), exponential_growth_model(),
        Prior(), 40_000; progress = false)
    m = vec(Array(chn[:m]))
    @test isapprox(mean(m), 9.0; atol = 0.3)
end
