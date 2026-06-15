## Tests for the per-vintage background random-effect submodel
## (`background_re_model`) and its daily expansion (`expand_vintage_rate`)
## from `src/models`. The random effect is the time-varying
## generalisation of the scalar `λ_bg` / `λ_bg_death` background: it is
## tightly pooled toward an informative scalar baseline so it cannot
## out-explain the real suspected-case signal, and `σ_bg → 0` recovers
## the scalar exactly.

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

@testitem "background_re_model σ_bg=0 recovers the scalar baseline" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: std
    using BVDOutbreakSize: background_re_model

    ## With the shared pooling SD passed as zero every per-vintage rate
    ## equals the sampled baseline, so the across-window spread is zero.
    nv = 5
    chn = sample(MersenneTwister(20260604),
        background_re_model(nv, 0.0), Prior(), 2_000; progress = false)
    ## `chn[:λ]` is a matrix (iter × chain) of length-`nv` vectors; each
    ## per-draw vector should be flat (zero across-window spread).
    draws = vec(Array(chn[:λ]))
    spreads = [std(d) for d in draws]
    @test maximum(spreads) < 1e-8
end

@testitem "background_pooling_model default σ_bg is a tight half-normal" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: median
    using BVDOutbreakSize: background_pooling_model

    ## The shared random-walk innovation SD is a half-normal SD 0.1, so its
    ## median is small (well under 0.2) and it is non-negative — the daily
    ## background walk is a gentle drift rather than per-day noise, though the
    ## data can pull it up to a modest rise over the surveillance window.
    chn = sample(MersenneTwister(20260604), background_pooling_model(),
        Prior(), 4_000; progress = false)
    σ = vec(Array(chn[:σ_bg]))
    @test all(>=(0), σ)
    @test median(σ) < 0.2
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

@testitem "background_walk_model edge cases (ungated, single day)" begin
    using BVDOutbreakSize: background_walk_model
    using Turing: returned
    using Random: MersenneTwister
    ## onset = 1 runs the walk over the whole grid (no leading zeros).
    s = returned(background_walk_model(10, 0.03; onset = 1),
        rand(MersenneTwister(2), background_walk_model(10, 0.03; onset = 1)))
    @test length(s.λ) == 10
    @test all(s.λ .> 0)
    ## a single-day window is just the baseline.
    s1 = returned(background_walk_model(5, 0.03; onset = 5),
        rand(MersenneTwister(2), background_walk_model(5, 0.03; onset = 5)))
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
    st = returned(background_walk_model(n, σ_rw; onset = onset, onset_ramp = 7),
        rand(MersenneTwister(5),
            background_walk_model(n, σ_rw; onset = onset, onset_ramp = 7)))
    @test all(st.λ[1:(onset - 1)] .== 0)
    @test st.λ[onset] ≈ st.λ_mu / 7         # first window day is 1/ramp of level
    @test st.λ[onset + 6] ≈ st.λ_mu          # reaches the level after the ramp
    @test st.λ[onset] < st.λ[onset + 1] < st.λ[onset + 6]   # monotone ramp-in
    ## onset_ramp = 1 recovers the old hard onset (full level on day one).
    hard = returned(background_walk_model(n, σ_rw; onset = onset, onset_ramp = 1),
        rand(MersenneTwister(5),
            background_walk_model(n, σ_rw; onset = onset, onset_ramp = 1)))
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
    st = returned(background_walk_model(n, σ_rw; onset = onset),
        rand(MersenneTwister(11), background_walk_model(n, σ_rw; onset = onset)))
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
    flat = returned(background_walk_model(n, 0.0; onset = onset),
        rand(MersenneTwister(11), background_walk_model(n, 0.0; onset = onset)))
    @test all(flat.λ[(onset + 7):end] .≈ flat.λ_mu)
end

@testitem "background_re_model is a positive perturbation of baseline" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using BVDOutbreakSize: background_re_model

    ## At a fixed shared σ_bg the per-vintage rates are a multiplicative
    ## log-normal deviation from the baseline: strictly positive.
    nv = 8
    chn = sample(MersenneTwister(20260604), background_re_model(nv, 0.3),
        Prior(), 4_000; progress = false)
    λ = reduce(vcat, vec(Array(chn[:λ])))
    @test all(>(0), λ)
end

@testitem "background_re_model baseline prior is overridable" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean
    using Distributions: truncated, Normal
    using BVDOutbreakSize: background_re_model

    ## A tight low baseline pulls every per-vintage rate down, confirming
    ## the baseline keyword is a real override point (used for the deaths
    ## background, which is far smaller than the cases background).
    chn = sample(MersenneTwister(20260604),
        background_re_model(4, 0.1;
            baseline_prior = truncated(Normal(0.0, 0.05); lower = 0)),
        Prior(), 4_000; progress = false)
    λ = reduce(vcat, vec(Array(chn[:λ])))
    @test mean(λ) < 0.15
end

@testitem "death_background_model default is a tight half-normal" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, std
    using BVDOutbreakSize: death_background_model

    ## Default `truncated(Normal(0, 0.25); lower = 0)`: fold the half-normal
    ## back to its untruncated SD via E|X| = σ√(2/π).
    chn = sample(MersenneTwister(20260604), death_background_model(),
        Prior(), 40_000; progress = false)
    λ = vec(Array(chn[:λ_bg_death]))
    @test isapprox(mean(λ) * sqrt(pi / 2), 0.25; atol = 0.02)
    @test all(>=(0), λ)
end
