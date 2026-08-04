## Tests for the partially-pooled per-stream surveillance dispersion
## (`pooled_dispersion_model`) and its wiring into `bvd_joint`, where the
## suspected-case, suspected-death, confirmed-case and confirmed-death streams
## each draw their own negative-binomial dispersion from a shared population
## rather than sharing one global `k`.

@testitem "pooled_dispersion_model: per-stream draws from a shared population" begin
    using Turing: returned
    using Random: default_rng, seed!
    using BVDOutbreakSize: pooled_dispersion_model

    n = 6
    seed!(3)
    m = pooled_dispersion_model(n)
    rng = default_rng()
    for _ in 1:50
        s = returned(m, rand(rng, m))
        @test length(s.k) == n
        @test all(isfinite, s.k) && all(>(0), s.k)
        @test s.k_pop > 0 && isfinite(s.k_pop)
        @test s.τ >= 0
        ## `k_pop` is the dispersion at the population mean, and each stream's
        ## `k` is the reciprocal square of its `1/sqrt(k)` deviation.
        @test isapprox(s.k_pop, 1 / exp(s.μ_log)^2; rtol = 1e-6)
        @test all(isapprox.(s.k, 1.0 ./ s.inv_sqrt_k .^ 2; rtol = 1e-6))
    end
end

@testitem "pooled_dispersion_model: zero pooling SD collapses to shared k" begin
    using Turing: returned
    using Random: default_rng, seed!
    using BVDOutbreakSize: pooled_dispersion_model
    using Distributions: Normal, truncated

    ## With the pooling SD pinned at ~0 every stream takes the population
    ## value, recovering the shared-`k` model.
    seed!(5)
    m = pooled_dispersion_model(6;
        sd_prior = truncated(Normal(0, 1e-9); lower = 0))
    rng = default_rng()
    for _ in 1:20
        s = returned(m, rand(rng, m))
        @test all(isapprox.(s.k, s.k_pop; rtol = 1e-4))
    end
end

@testitem "pooled_dispersion_model: centred and non-centred both valid" begin
    using Turing: returned
    using Random: default_rng, seed!
    using BVDOutbreakSize: pooled_dispersion_model

    ## The default is the centred parameterisation; `centred = false` gives the
    ## non-centred form. Both yield a valid per-stream dispersion vector with
    ## the same `k = 1/inv_sqrt_k^2` structure.
    n = 6
    for centred in (true, false)
        seed!(7)
        m = pooled_dispersion_model(n; centred = centred)
        rng = default_rng()
        for _ in 1:30
            s = returned(m, rand(rng, m))
            @test length(s.k) == n
            @test all(isfinite, s.k) && all(>(0), s.k)
            @test all(isapprox.(s.k, 1.0 ./ s.inv_sqrt_k .^ 2; rtol = 1e-6))
        end
    end
end

@testitem "bvd_joint: exposes partially-pooled per-stream dispersions" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: load_observations, bvd_joint, genetic_seeding_model

    obs = load_observations()
    breakpoint = obs.n - obs.who_first_sitrep_days
    m = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        suspected_daily_history = obs.suspected_daily_history,
        isolation_history = obs.isolation_history,
        bed_capacity_history = obs.bed_capacity_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = breakpoint,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 30;
        chain_type = FlexiChains.VNChain, progress = false)

    ## Every per-stream dispersion (including the isolation and recovered
    ## streams, now pooled rather than independent in the joint), the
    ## population value and the pooling SD are exposed and positive.
    for key in (:k, :k_cases, :k_deaths, :k_confirmed, :k_confirmed_deaths,
        :isolation_dispersion, :recovered_dispersion)
        v = vec(Array(chn[key]))
        @test length(v) == 30
        @test all(isfinite, v)
        @test all(v .> 0)
    end
    sd = vec(Array(chn[:dispersion_sd]))
    @test all(isfinite, sd)
    @test all(sd .>= 0)
end
