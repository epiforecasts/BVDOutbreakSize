## Tests for the contact-tracing background covariate: the observed contact
## follow-up rate expanded onto the grid (`expand_covariate`), the coefficient
## submodel (`contact_background_model`) and its entry into the suspected-case
## background through `bvd_joint`. The covariate is a fixed offset on observed
## data, so it adds an outbreak-scaled non-BVD background component without the
## ascertainment/outbreak-size degeneracy a latent-scaled background carries.

@testitem "expand_covariate interpolates within the span and zeros outside" begin
    using BVDOutbreakSize: expand_covariate

    ## Linear interpolation between the reported knot days, zero before the
    ## first day and after the last (the covariate contributes only where it is
    ## observed).
    @test expand_covariate([3, 6], [1.0, 4.0], 8) ==
          [0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 0.0]
    ## A single knot holds flat on its own day only.
    @test expand_covariate([2], [0.5], 4) == [0.0, 0.5, 0.0, 0.0]
    ## An empty covariate is a zero vector (a no-op offset).
    @test expand_covariate(Int[], Float64[], 4) == zeros(4)
    ## The element type follows the values and the length follows the grid.
    out = expand_covariate([2, 5], [0.6, 0.9], 6)
    @test out isa Vector{Float64}
    @test length(out) == 6
    @test all(iszero, out[[1, 6]])
end

@testitem "contact_followup_history loads as fractions on the grid" begin
    using BVDOutbreakSize: load_observations, expand_covariate

    obs = load_observations()
    hist = obs.contact_followup_history
    ## The manifest carries the 7 June-1 July follow-up rate as fractions.
    @test !isempty(hist.values)
    @test length(hist.days) == length(hist.values)
    @test all(0 .< hist.values .< 1)
    @test issorted(hist.days)
    ## Expanded onto the grid it is a length-n series, nonzero only across the
    ## reported span.
    cov = expand_covariate(hist.days, hist.values, obs.n)
    @test length(cov) == obs.n
    @test count(!iszero, cov) > 0
    @test all(iszero, cov[1:(first(hist.days) - 1)])
end

@testitem "contact_background_model draws a non-negative coefficient" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using BVDOutbreakSize: contact_background_model

    ## The coefficient is a weakly-informative half-normal on the natural
    ## scale, so every draw is non-negative (more case-finding surfaces more,
    ## not fewer, non-BVD suspects).
    chn = sample(MersenneTwister(20260704), contact_background_model(),
        Prior(), 2_000; progress = false)
    β = vec(Array(chn[Symbol("β_contact")]))
    @test all(>=(0), β)
    @test any(>(0), β)
end

@testitem "contact follow-up enters the joint and adds sampled parameters" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: load_observations, bvd_joint, genetic_seeding_model

    obs = load_observations()
    breakpoint = obs.n - obs.who_first_sitrep_days
    empty_hist = (; days = Int[], values = Float64[])

    make(hist) = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
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
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = breakpoint,
        background_re = true,
        contact_followup_history = hist,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days)

    ## With the follow-up rate the joint samples the latent contact process and
    ## the background coefficient and stays finite and positive.
    chn = sample(make(obs.contact_followup_history), Prior(), 30;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("contact_state", k), ks)
    @test any(k -> occursin("contact_bg_state", k) || occursin("β_contact", k),
        ks)
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 30
    @test all(isfinite, C_T)
    @test all(C_T .> 0)

    ## Without the follow-up history no contact process or coefficient is
    ## sampled, so the fit is unchanged from the covariate-free model.
    chn0 = sample(make(empty_hist), Prior(), 30;
        chain_type = FlexiChains.VNChain, progress = false)
    ks0 = string.(collect(keys(chn0)))
    @test !any(k -> occursin("contact_state", k) || occursin("β_contact", k),
        ks0)
end
