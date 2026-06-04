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

@testitem "background_re_model pooling SD→0 recovers the scalar" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: std, mean
    using Distributions: truncated, Normal
    using BVDOutbreakSize: background_re_model

    ## With the pooling SD pinned at zero every per-vintage rate equals
    ## the sampled baseline, so the across-window spread is exactly zero.
    nv = 5
    chn = sample(MersenneTwister(20260604),
        background_re_model(nv;
            pooling_prior = truncated(Normal(0, 1e-8); lower = 0)),
        Prior(), 2_000; progress = false)
    λ = Array(chn[:λ])  # draws × nv
    spreads = [std(λ[i, :]) for i in axes(λ, 1)]
    @test maximum(spreads) < 1e-5
end

@testitem "background_re_model is a positive perturbation of baseline" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: median, mean
    using BVDOutbreakSize: background_re_model

    ## Default tight pooling: per-vintage rates stay close to the baseline
    ## (log-normal multiplicative deviation), are strictly positive, and
    ## their geometric centre tracks the half-normal baseline median.
    nv = 8
    chn = sample(MersenneTwister(20260604), background_re_model(nv),
        Prior(), 4_000; progress = false)
    λ = vec(Array(chn[:λ]))
    @test all(>(0), λ)
    σ = vec(Array(chn[:σ_bg]))
    ## The default pooling prior is half-normal SD 0.3, so the typical
    ## per-window log-deviation is modest (median σ_bg well under 0.3).
    @test median(σ) < 0.3
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
        background_re_model(4;
            baseline_prior = truncated(Normal(0.0, 0.05); lower = 0)),
        Prior(), 4_000; progress = false)
    λ = vec(Array(chn[:λ]))
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
