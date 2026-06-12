## Tests for the confirmed-case window alignment, the per-window
## positivity random effect and the confirmed-death thinning stream.

@testitem "confirmed_positivity_windows splits early and observed windows" begin
    using BVDOutbreakSize: confirmed_positivity_windows

    ## The 28 May data: analysed cumulative stalls 24-25 May (flat at 295).
    ## Confirmed vintages up to the first lab date (18-23 May) become early
    ## windows scored against the modelled volume; the observed windows
    ## (24-28 May) merge the zero-denominator stall to positives
    ## [4, 16, 4, 85] of analysed [84, 108, 245, 107].
    confirmed = (; days = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        counts = [33, 51, 57, 79, 83, 101, 105, 106, 121, 125, 210])
    lab = (; days = [6, 7, 8, 9, 10, 11],
        counts = [211, 295, 295, 403, 648, 755])
    w = confirmed_positivity_windows(confirmed, lab)
    @test w.early_days == [1, 2, 3, 4, 5, 6]
    @test w.early_increments == [33, 18, 6, 22, 4, 18]
    @test w.obs_positives == [4, 16, 4, 85]
    @test w.obs_analysed == [84, 108, 245, 107]
    @test all(0 .<= w.obs_positives .<= w.obs_analysed)
    ## Early + observed partition the confirmed total, no double-count.
    @test sum(w.early_increments) + sum(w.obs_positives) ==
          confirmed.counts[end]
end

@testitem "confirmed_positivity_windows handles missing histories" begin
    using BVDOutbreakSize: confirmed_positivity_windows
    empty = (; days = Int[], counts = Int[])
    lab = (; days = [3, 5], counts = [10, 20])
    ## No confirmed history: everything empty.
    w0 = confirmed_positivity_windows(empty, lab)
    @test isempty(w0.obs_analysed) && isempty(w0.early_increments)
    ## No lab history: every confirmed vintage is an early window.
    conf = (; days = [3, 5], counts = [2, 6])
    w1 = confirmed_positivity_windows(conf, empty)
    @test isempty(w1.obs_analysed)
    @test w1.early_days == [3, 5]
    @test w1.early_increments == [2, 4]
end

@testitem "confirmed_positivity_model draws per-window probabilities" begin
    using BVDOutbreakSize: confirmed_positivity_model
    using Random: seed!
    seed!(1)
    s = confirmed_positivity_model(5)()
    @test length(s.p_pos) == 5
    @test all(0 .<= s.p_pos .<= 1)
    @test s.σ_q >= 0
end

@testitem "confirmed_deaths_model gives a bounded confirmation probability" begin
    using BVDOutbreakSize: confirmed_deaths_model
    using Turing: returned
    using Random: MersenneTwister

    ## 246 suspected deaths, 17 confirmed. The death-pool composition q_death
    ## is built from the death series' own BVD and background components, and
    ## the confirmation positivity is the assay transform of that composition.
    ## Modelled BVD and background suspected-death daily series; the suspected-
    ## death total is their sum. k is the dispersion.
    bvd_deaths = fill(5.5, 40)
    bg_death = fill(0.5, 40)
    deaths_daily = bvd_deaths .+ bg_death
    m = confirmed_deaths_model(17, 246, deaths_daily, bvd_deaths, bg_death, 5.0)
    st = returned(m, rand(MersenneTwister(2), m))
    @test 0 < st.q_death < 1
    @test 0 < st.p_death_conf < 1
    @test 0 < st.τ_death < 1
    @test st.expected_confirmed_deaths >= 0
end

@testitem "confirmed_deaths_model lags the series by the receipt delay" begin
    using BVDOutbreakSize: confirmed_deaths_model
    using Turing: returned
    using Random: MersenneTwister

    bvd_deaths = fill(5.5, 40)
    bg_death = fill(0.5, 40)
    deaths_daily = bvd_deaths .+ bg_death
    seed = 7
    ## Same args, one with no delay (identity PMF) and one with all mass at a
    ## five-day lag. The confirmation positivity is unchanged (the composition
    ## is flat, so a uniform delay leaves the BVD share fixed), but the delay
    ## pushes five of the forty flat days off the end of the grid, so the
    ## windowed expected confirmed deaths fall to 35/40 of the undelayed total.
    base = confirmed_deaths_model(17, 246, deaths_daily, bvd_deaths, bg_death,
        5.0)
    delayed = confirmed_deaths_model(17, 246, deaths_daily, bvd_deaths,
        bg_death, 5.0; receipt_pmf = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
    sb = returned(base, rand(MersenneTwister(seed), base))
    sd = returned(delayed, rand(MersenneTwister(seed), delayed))
    @test sd.p_death_conf ≈ sb.p_death_conf
    @test sd.expected_confirmed_deaths < sb.expected_confirmed_deaths
    @test sd.expected_confirmed_deaths ≈ (35 / 40) * sb.expected_confirmed_deaths rtol=1e-6
end

@testitem "confirmed_deaths_only_model conditions and stays finite" begin
    using BVDOutbreakSize: confirmed_deaths_only_model
    using Turing.DynamicPPL: logjoint
    using Random: MersenneTwister
    m = confirmed_deaths_only_model(40, 17, 246;
        deaths_history = (; days = [20, 40], counts = [120, 246]))
    draw = rand(MersenneTwister(1), m)
    @test isfinite(logjoint(m, draw))
end
