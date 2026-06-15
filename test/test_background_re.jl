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

    ## The shared random-walk innovation SD is a tight half-normal SD 0.05, so
    ## its median is very small (well under 0.1) and it is non-negative — the
    ## daily background walk is a gentle drift, not per-day noise.
    chn = sample(MersenneTwister(20260604), background_pooling_model(),
        Prior(), 4_000; progress = false)
    σ = vec(Array(chn[:σ_bg]))
    @test all(>=(0), σ)
    @test median(σ) < 0.1
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
    ## The day-to-day multiplicative change is small (tight drift), so the log
    ## series has no large jumps.
    logλ = log.(st.λ[onset:end])
    @test maximum(abs.(diff(logλ))) < 0.5
    ## σ_rw = 0 recovers a flat background at the baseline over the window.
    flat = returned(background_walk_model(n, 0.0; onset = onset),
        rand(MersenneTwister(11), background_walk_model(n, 0.0; onset = onset)))
    @test all(flat.λ[onset:end] .≈ flat.λ_mu)
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
