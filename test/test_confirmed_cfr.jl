## Delay-corrected confirmed CFR: pure-logic tests for the residual-delay
## correction, a synthetic-chain contract test for the posterior wrapper, and
## a smoke test for the summary table.

@testitem "delay_corrected_cfr matches the naive ratio with no delay" begin
    using BVDOutbreakSize: delay_corrected_cfr
    ## Confirmation and death both at lag 0: every case has resolved by the
    ## cut-off, so no denominator shrinkage and the corrected ratio is the
    ## naive cumulative-deaths-over-cases.
    c_daily = [1.0, 1.0, 1.0]
    Kc = [1.0]
    Kd = [1.0]
    corrected = delay_corrected_cfr(c_daily, Kc, Kd, 0.9)
    @test corrected ≈ 0.9 / 3 atol = 1e-10
end

@testitem "delay_corrected_cfr lifts the naive ratio under a death delay" begin
    using BVDOutbreakSize: delay_corrected_cfr
    ## Deaths land five days after onset; confirmation is immediate. Only the
    ## earliest day's cases (five days before the cut-off) have had time to
    ## die, so the denominator collapses to that one day and the corrected
    ## ratio jumps well above the naive ratio.
    c_daily = ones(6)
    Kc = [1.0]
    Kd = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    corrected = delay_corrected_cfr(c_daily, Kc, Kd, 0.9)
    naive = 0.9 / sum(c_daily)
    @test corrected ≈ 0.9 atol = 1e-10
    @test corrected > naive
end

@testitem "delay_corrected_cfr grows with a longer death delay" begin
    using BVDOutbreakSize: delay_corrected_cfr
    c_daily = ones(8)
    Kc = [1.0]
    short = delay_corrected_cfr(c_daily, Kc, [0.0, 0.0, 1.0], 0.5)
    long = delay_corrected_cfr(c_daily, Kc, [0.0, 0.0, 0.0, 0.0, 0.0, 1.0], 0.5)
    @test long > short
end

@testitem "delay_corrected_cfr guards empty inputs" begin
    using BVDOutbreakSize: delay_corrected_cfr
    @test isnan(delay_corrected_cfr(Float64[], [1.0], [1.0], 1.0))
    @test isnan(delay_corrected_cfr([1.0], Float64[], [1.0], 1.0))
    @test isnan(delay_corrected_cfr([0.0, 0.0], [1.0], [1.0], 1.0))
end

@testsnippet ConfirmedCfrChain begin
    using Turing: @model, sample, Prior
    using Distributions: Beta, truncated, Normal
    import FlexiChains
    using BVDOutbreakSize: delay_corrected_confirmed_cfr

    ## Synthetic prior carrying the deterministics the posterior wrapper
    ## reads: the cumulative confirmed-case trajectory, the onset-to-
    ## confirmation and onset-to-death delay PMFs, the cut-off confirmed and
    ## confirmed-death totals, and the structural CFR. A late-weighted
    ## confirmed ramp guarantees recent unresolved cases, so the correction
    ## bites.
    @model function _ccfr_test()
        cfr ~ Beta(6.6, 13.4)
        CFR := cfr
        expected_confirmed_T ~ truncated(Normal(120.0, 20.0); lower = 1.0)
        expected_confirmed_deaths_T ~ truncated(Normal(8.0, 2.0); lower = 0.5)
        ## Rising daily confirmed-case incidence over a 30-day grid.
        daily = collect(range(1.0, 6.0; length = 30))
        cumulative_confirmed := cumsum(daily)
        ## Onset-to-confirmation concentrated early, onset-to-death-
        ## confirmation spread later, both normalised PMFs from lag 0.
        kc = [0.2, 0.5, 0.3]
        onset_to_confirmation_pmf := kc ./ sum(kc)
        kd = [0.0, 0.05, 0.1, 0.2, 0.3, 0.2, 0.15]
        onset_to_death_confirmation_pmf := kd ./ sum(kd)
        return nothing
    end

    _ccfr_chain(n) = sample(
        _ccfr_test(), Prior(), n;
        chain_type = FlexiChains.VNChain, progress = false
    )
end

@testitem "delay_corrected_confirmed_cfr returns aligned draw vectors" setup=[
    ConfirmedCfrChain
] begin
    chn = _ccfr_chain(300)
    res = delay_corrected_confirmed_cfr(chn;
        obs_confirmed = 210, obs_confirmed_deaths = 17)

    @test length(res.corrected) == 300
    @test length(res.modelled_naive) == 300
    @test length(res.structural) == 300
    @test res.naive_observed ≈ 17 / 210
    ## The correction shrinks the denominator below the cumulative confirmed
    ## cases, so the corrected ratio is at least the uncorrected modelled
    ## ratio on every finite draw.
    ok = isfinite.(res.corrected) .& isfinite.(res.modelled_naive)
    @test all(res.corrected[ok] .>= res.modelled_naive[ok] .- 1e-8)
    ## With a real death delay and a rising confirmed series the correction
    ## is strict in the median.
    using Statistics: median
    @test median(res.corrected[ok]) > median(res.modelled_naive[ok])
end

@testitem "confirmed_cfr_table summarises the four quantities" setup=[
    ConfirmedCfrChain
] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: confirmed_cfr_table
    chn = _ccfr_chain(200)
    res = delay_corrected_confirmed_cfr(chn;
        obs_confirmed = 210, obs_confirmed_deaths = 17)
    tbl = confirmed_cfr_table(res)
    @test tbl isa DataFrame
    @test nrow(tbl) == 4
end

@testitem "plot_confirmed_cfr returns a Makie figure" setup=[
    ConfirmedCfrChain
] begin
    using CairoMakie
    CairoMakie.activate!(type = "png")
    using BVDOutbreakSize: plot_confirmed_cfr
    chn = _ccfr_chain(200)
    res = delay_corrected_confirmed_cfr(chn;
        obs_confirmed = 210, obs_confirmed_deaths = 17)
    fig = plot_confirmed_cfr(res)
    @test fig isa CairoMakie.Makie.Figure
end
