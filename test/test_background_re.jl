## Tests for the non-BVD background walk (`background_walk_model`), its
## shared pooling SD (`background_pooling_model`) and the daily window
## expansion (`expand_vintage_rate`) from `src/models`. The walk is the
## time-varying generalisation of the scalar `λ_bg` background: it is
## tightly pooled toward an informative scalar baseline so it cannot
## out-explain the real suspected-case signal.

@testitem "expand_vintage_rate maps windows and carries the tail" begin
    using BVDOutbreakSize: expand_vintage_rate
    r = [1.0, 2.0, 3.0]
    days = [3, 6, 9]
    @test expand_vintage_rate(r, days, 10) ==
        [1, 1, 1, 2, 2, 2, 3, 3, 3, 3]
    ## Empty days recovers a flat scalar over the grid.
    @test expand_vintage_rate([0.5], Int[], 4) == fill(0.5, 4)
    @test all(iszero, expand_vintage_rate(Float64[], Int[], 3))
    ## A grid shorter than the last vintage edge clamps the windows.
    @test expand_vintage_rate([1.0, 2.0], [2, 8], 5) ==
        [1.0, 1.0, 2.0, 2.0, 2.0]
end

@testitem "background_pooling_model default σ_bg regularises the walk" tags = [
    :slow,
] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: median, quantile
    using BVDOutbreakSize: background_pooling_model

    ## The shared random-walk innovation SD is a half-normal of scale 0.3, so
    ## it is non-negative and its median is about 0.2. The daily background
    ## walk stays a gentle drift rather than per-day noise, while leaving room
    ## for the 0.17 to 0.22 the resumed daily new-suspect series pulls it to.
    chn = sample(
        MersenneTwister(20260604), background_pooling_model(),
        Prior(), 4_000; progress = false
    )
    σ = vec(Array(chn[:σ_bg]))
    @test all(>=(0), σ)
    @test median(σ) < 0.3
    ## Still regularised: the walk is not free to absorb whole windows.
    @test quantile(σ, 0.95) < 0.7
end

@testitem "gate_before zeroes a series before the onset day" begin
    using BVDOutbreakSize: gate_before
    v = collect(1.0:6.0)
    ## start ≤ 1 returns the series unchanged (no gating).
    @test gate_before(v, 1) == v
    @test gate_before(v, 0) == v
    ## start > 1 zeroes the entries before `start` and keeps the rest.
    @test gate_before(v, 4) == [0.0, 0.0, 0.0, 4.0, 5.0, 6.0]
    ## element type is preserved; the tail is untouched.
    g = gate_before(v, 3)
    @test eltype(g) == eltype(v)
    @test g[3:end] == v[3:end] && all(g[1:2] .== 0)
end

@testitem "background_walk_model default baseline is a half-normal SD 20" tags = [
    :slow,
] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean
    using BVDOutbreakSize: background_walk_model

    ## Default `truncated(Normal(0, 20); lower = 0)` on the window anchor: fold
    ## the half-normal back to its untruncated SD via E|X| = σ√(2/π). Wide
    ## enough that the joint posterior for `λ_mu` sits inside it rather than
    ## against its upper tail.
    chn = sample(
        MersenneTwister(20260604), background_walk_model(40, 0.03),
        Prior(), 40_000; progress = false
    )
    λ_mu = vec(Array(chn[:λ_mu]))
    @test isapprox(mean(λ_mu) * sqrt(pi / 2), 20.0; atol = 1.0)
    @test all(>=(0), λ_mu)
end

@testitem "background_walk_model edge cases (ungated, single day)" begin
    using BVDOutbreakSize: background_walk_model
    using Turing: returned
    using Random: MersenneTwister
    ## onset = 1 runs the walk over the whole grid (no leading zeros).
    s = returned(
        background_walk_model(10, 0.03; onset = 1),
        rand(MersenneTwister(2), background_walk_model(10, 0.03; onset = 1))
    )
    @test length(s.λ) == 10
    @test all(s.λ .> 0)
    ## a single-day window is just the baseline.
    s1 = returned(
        background_walk_model(5, 0.03; onset = 5),
        rand(MersenneTwister(2), background_walk_model(5, 0.03; onset = 5))
    )
    @test all(s1.λ[1:4] .== 0)
    @test s1.λ[5] ≈ s1.λ_mu
end

@testitem "background_walk_model ramps in across the onset boundary" begin
    using BVDOutbreakSize: background_walk_model
    using Turing: returned
    using Random: MersenneTwister
    ## With the default onset ramp the gated background grows in from zero
    ## rather than stepping straight to the baseline at the onset, so the first
    ## non-zero day is a small fraction of the baseline and the day-to-day rise
    ## across the boundary is gradual (no one-day jump to the full level).
    n, onset, σ_rw = 40, 18, 0.0
    st = returned(
        background_walk_model(n, σ_rw; onset = onset, onset_ramp = 7),
        rand(
            MersenneTwister(5),
            background_walk_model(n, σ_rw; onset = onset, onset_ramp = 7)
        )
    )
    @test all(st.λ[1:(onset - 1)] .== 0)
    @test st.λ[onset] ≈ st.λ_mu / 7  # first window day is 1/ramp of level
    @test st.λ[onset + 6] ≈ st.λ_mu  # reaches the level after the ramp
    @test st.λ[onset] < st.λ[onset + 1] < st.λ[onset + 6]  # monotone ramp-in
    ## onset_ramp = 1 recovers the old hard onset (full level on day one).
    hard = returned(
        background_walk_model(n, σ_rw; onset = onset, onset_ramp = 1),
        rand(
            MersenneTwister(5),
            background_walk_model(n, σ_rw; onset = onset, onset_ramp = 1)
        )
    )
    @test hard.λ[onset] ≈ hard.λ_mu
end

@testitem "background_walk_model is smooth, gated and bounded" begin
    using Turing: returned
    using Random: MersenneTwister
    using Statistics: mean
    using BVDOutbreakSize: background_walk_model

    ## Daily lognormal random walk over the surveillance window: zero before
    ## the onset, strictly positive after it, and (with a tight σ_rw) a smooth
    ## gentle drift around the half-normal baseline rather than per-vintage
    ## steps.
    n, onset, σ_rw = 30, 8, 0.04
    st = returned(
        background_walk_model(n, σ_rw; onset = onset),
        rand(
            MersenneTwister(11),
            background_walk_model(n, σ_rw; onset = onset)
        )
    )
    @test length(st.λ) == n
    @test all(st.λ[1:(onset - 1)] .== 0)            # gated before the onset
    @test all(st.λ[onset:end] .> 0)                 # positive after it
    @test st.σ_bg == σ_rw
    ## After the onset ramp (default 7 days) the walk is a tight gentle drift,
    ## so the log series has no large jumps over that stretch.
    logλ = log.(st.λ[(onset + 7):end])
    @test maximum(abs.(diff(logλ))) < 0.5
    ## σ_rw = 0 recovers a flat background at the baseline over the window once
    ## the onset ramp has completed.
    flat = returned(
        background_walk_model(n, 0.0; onset = onset),
        rand(MersenneTwister(11), background_walk_model(n, 0.0; onset = onset))
    )
    @test all(flat.λ[(onset + 7):end] .≈ flat.λ_mu)
end

@testitem "bvd_joint runs the pooled background branch" tags = [:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: load_observations, bvd_joint, genetic_seeding_model

    ## `background_re = true` is what every registry fit uses and the only
    ## path that samples the pooling SD, but nothing else in the suite sets it.
    obs = load_observations()
    breakpoint = obs.n - obs.who_first_sitrep_days
    m = bvd_joint(
        obs.n, obs.exported_cases, obs.total_deaths,
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
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days
    )
    chn = sample(
        m, Prior(), 20;
        chain_type = FlexiChains.VNChain, progress = false
    )

    ## The pooling SD reaches the chain, so the gated tilde ran.
    ks = collect(keys(chn))
    σ_key = only(filter(k -> occursin("σ_bg", string(k)), ks))
    σ = vec(Array(chn[σ_key]))
    @test length(σ) == 20
    @test all(isfinite, σ)
    @test all(>=(0), σ)

    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end
