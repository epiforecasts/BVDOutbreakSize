## Tests for the dated per-day Uganda export likelihood: the
## `dated_event_bins` / `dated_poisson_model` helpers and the dated paths of
## `exports_model` / `exports_deaths_model` exercised through their joint
## composers.

@testitem "dated_event_bins groups events into unique days and counts" begin
    using BVDOutbreakSize: dated_event_bins

    ## Distinct days keep one term each.
    days, counts = dated_event_bins([3, 7, 12], 20)
    @test days == [3, 7, 12]
    @test counts == [1, 1, 1]

    ## Simultaneous events on one day share an edge with occupancy > 1.
    days2, counts2 = dated_event_bins([5, 5, 9], 20)
    @test days2 == [5, 9]
    @test counts2 == [2, 1]

    ## Out-of-grid days clamp into [1, n].
    days3, counts3 = dated_event_bins([0, 25], 10)
    @test days3 == [1, 10]
    @test counts3 == [1, 1]
end

@testitem "dated_poisson_model scores counts and generates when missing" begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: dated_poisson_model

    means = [1.0, 2.0, 0.5]

    ## Observed path: the supplied counts are conditioned, nothing sampled.
    chn = sample(dated_poisson_model(means, [1, 2, 0]), Prior(), 10;
        chain_type = FlexiChains.VNChain, progress = false)
    @test chn isa FlexiChains.VNChain

    ## Generator path: missing counts are sampled as one per-day vector,
    ## every entry a non-negative integer.
    gen = sample(dated_poisson_model(means, missing), Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)
    k = first(k for k in keys(gen) if occursin("counts", string(k)))
    draws = vec(Array(gen[k]))
    @test length(draws) == 50
    @test all(v -> length(v) == length(means), draws)
    @test all(v -> all(v .>= 0), draws)
end

@testitem "exports_only dated series prior draws are finite" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: exports_only_model

    ## Three dated imports on grid days 28, 33, 40 of a 40-day grid (the
    ## last import at the cut-off), fitted with the per-day Poisson.
    chn = sample(
        exports_only_model(40, 3; export_case_days = [28, 33, 40]),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 100
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "exports_joint_only fits cases and deaths together" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: exports_joint_only_model

    ## Conditions on both the dated import series and the dated import-death
    ## series over the one travel-gated prevalence; check_model is off (the
    ## predictive deaths/exports leave redundant discrete draws).
    chn = sample(
        exports_joint_only_model(40, 3, 1; export_case_days = [28, 33, 40],
            export_death_days = [30]),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 100
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "exports_only last_offset stops the clock before the cut-off" tags=[
    :slow
] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: exports_only_model

    ## The last import is on day 30 of a 40-day grid (last_offset = 10), so
    ## the export clock stops at day 30: the model still draws finite,
    ## positive sizes and does not accrue prevalence past the last import.
    chn = sample(
        exports_only_model(40, 3; export_case_days = [22, 26, 30]),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "dated export increments partition the cumulative expectation" begin
    using BVDOutbreakSize: bin_increments

    ## With the dated likelihood the pre-detection survival weight plus the
    ## per-day increments reconstruct the cumulative export intensity at the
    ## last import day, the same total a single cumulative Poisson scores.
    prevalence = collect(1.0:1.0:10.0)        # Λ(t) = cumsum
    days = [3, 6, 9]                          # detection days, t_last = 9
    d₁ = days[1]
    pre = sum(prevalence[1:(d₁ - 1)])         # Λ(d₁ - 1)
    μ = bin_increments(prevalence, days)
    μ[1] = μ[1] - pre
    Λ_tlast = sum(prevalence[1:days[end]])    # cumulative to t_last
    @test pre + sum(μ) ≈ Λ_tlast
    @test all(μ .> 0)
end
