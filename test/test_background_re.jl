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

@testitem "gate_background zeros the pre-surveillance span" begin
    using BVDOutbreakSize: gate_background
    bg = fill(2.0, 10)
    ## A start of 4 zeros days 1-3 and keeps days 4-10.
    @test gate_background(bg, 4) == [0, 0, 0, 2, 2, 2, 2, 2, 2, 2]
    ## start ≤ 1 leaves the series unchanged (legacy ungated behaviour).
    @test gate_background(bg, 1) == bg
    @test gate_background(bg, 0) == bg
    ## A start past the grid zeros everything.
    @test all(iszero, gate_background(bg, 20))
    ## The element type follows the input and the BVD signal is untouched
    ## because gating only the background leaves the first bin's signal.
    @test eltype(gate_background(fill(1.5f0, 4), 2)) == Float32
    ## Cumulative over a gated window is strictly less than the ungated one
    ## whenever the start day is interior, the first-bin reduction.
    @test sum(gate_background(bg, 5)) < sum(bg)
end

@testitem "background_window gates the cumulative background lower" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: median
    import FlexiChains
    using BVDOutbreakSize: cases_only_model

    ## The first vintage sits well into the grid (day 30 of 40), so the
    ## ungated background accrues over 30 pre-surveillance days the gated
    ## one drops. Over matched prior draws the windowed background total is
    ## a strict reduction, the first-bin fix.
    history = (; days = [30, 35, 40], counts = [516, 905, 1077])
    common = (; reported_history = history, breakpoint = 10)

    gated = sample(
        cases_only_model(40, missing; common..., background_window = true),
        Prior(), 300; chain_type = FlexiChains.VNChain, progress = false)
    ungated = sample(
        cases_only_model(40, missing; common..., background_window = false),
        Prior(), 300; chain_type = FlexiChains.VNChain, progress = false)

    bg_g = vec(Array(gated[:background_total]))
    bg_u = vec(Array(ungated[:background_total]))
    @test all(isfinite, bg_g) && all(bg_g .>= 0)
    @test median(bg_g) < median(bg_u)
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

    ## The shared pooling SD is a half-normal SD 0.3, so its median is
    ## modest (well under 0.3) and it is non-negative.
    chn = sample(MersenneTwister(20260604), background_pooling_model(),
        Prior(), 4_000; progress = false)
    σ = vec(Array(chn[:σ_bg]))
    @test all(>=(0), σ)
    @test median(σ) < 0.3
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
