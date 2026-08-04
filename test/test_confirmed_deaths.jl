## Tests for the confirmed-case window alignment, the per-window
## positivity random effect and the confirmed-death thinning stream.

@testitem "confirmed_positivity_windows splits early and observed windows" begin
    using BVDOutbreakSize: confirmed_positivity_windows

    ## The 28 May data: analysed cumulative stalls 24-25 May (flat at 295).
    ## The first confirmed vintage (day 1, count 33) is the testing-onset
    ## BASELINE and is NOT scored, so the early windows are the vintages AFTER
    ## it up to the first lab date; `early_start` carries that baseline day so
    ## the model pins the early-window volume there. The observed windows
    ## (24-28 May) merge the zero-denominator stall to positives
    ## [4, 16, 4, 85] of analysed [84, 108, 245, 107].
    confirmed = (; days = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        counts = [33, 51, 57, 79, 83, 101, 105, 106, 121, 125, 210])
    lab = (; days = [6, 7, 8, 9, 10, 11],
        counts = [211, 295, 295, 403, 648, 755])
    w = confirmed_positivity_windows(confirmed, lab)
    @test w.early_start == 1
    @test w.early_days == [2, 3, 4, 5, 6]
    @test w.early_increments == [18, 6, 22, 4, 18]
    @test w.obs_positives == [4, 16, 4, 85]
    @test w.obs_analysed == [84, 108, 245, 107]
    @test all(0 .<= w.obs_positives .<= w.obs_analysed)
    ## The baseline (33) plus the fitted early increments and observed
    ## positives partition the confirmed total, no double-count.
    @test confirmed.counts[1] + sum(w.early_increments) +
          sum(w.obs_positives) == confirmed.counts[end]
end

@testitem "confirmed_positivity_windows handles missing histories" begin
    using BVDOutbreakSize: confirmed_positivity_windows
    empty = (; days = Int[], counts = Int[])
    lab = (; days = [3, 5], counts = [10, 20])
    ## No confirmed history: everything empty.
    w0 = confirmed_positivity_windows(empty, lab)
    @test isempty(w0.obs_analysed) && isempty(w0.early_increments)
    ## No lab history: the first confirmed vintage is the baseline and every
    ## vintage AFTER it is an early window.
    conf = (; days = [3, 5], counts = [2, 6])
    w1 = confirmed_positivity_windows(conf, empty)
    @test isempty(w1.obs_analysed)
    @test w1.early_start == 3
    @test w1.early_days == [5]
    @test w1.early_increments == [4]
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

@testitem "confirmed_deaths_model bounds the confirmation probability" begin
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
    @test sd.expected_confirmed_deaths ≈
          (35 / 40) * sb.expected_confirmed_deaths rtol=1e-6
end

@testitem "confirmed_deaths_model gates capacity before the testing onset" begin
    using BVDOutbreakSize: confirmed_deaths_model
    using Turing: returned
    using Random: MersenneTwister

    ## No deaths are confirmed before testing began: with a flat death series,
    ## gating the death "analysed" volume at day 21 of 40 zeroes the first 20
    ## days, so the expected confirmed deaths fall to 20/40 of the ungated
    ## total. The default `capacity_start = 0` leaves it ungated.
    bvd_deaths = fill(5.5, 40)
    bg_death = fill(0.5, 40)
    deaths_daily = bvd_deaths .+ bg_death
    seed = 4
    ung = confirmed_deaths_model(
        17, 246, deaths_daily, bvd_deaths, bg_death, 5.0)
    gat = confirmed_deaths_model(17, 246, deaths_daily, bvd_deaths, bg_death,
        5.0; capacity_start = 21)
    su = returned(ung, rand(MersenneTwister(seed), ung))
    sg = returned(gat, rand(MersenneTwister(seed), gat))
    @test sg.expected_confirmed_deaths < su.expected_confirmed_deaths
    @test sg.expected_confirmed_deaths ≈
          (20 / 40) * su.expected_confirmed_deaths rtol=1e-6
end

@testitem "confirmed_deaths_model scales the case analysed volume" begin
    using BVDOutbreakSize: confirmed_deaths_model
    using Turing: returned
    using Random: MersenneTwister

    ## In the joint the death analysed volume is scaling times the modelled
    ## case analysed volume carried at the per-day suspected death-to-case
    ## ratio, so the expected confirmed deaths are linear in the case analysed
    ## volume: doubling it doubles them (the scaling and assay draws are fixed
    ## by the seed, the suspected series are unchanged).
    bvd_deaths = fill(5.5, 40)
    bg_death = fill(0.5, 40)
    deaths_daily = bvd_deaths .+ bg_death
    susp_case = fill(20.0, 40)
    analysed_case = fill(6.0, 40)
    seed = 3
    base = confirmed_deaths_model(17, 246, deaths_daily, bvd_deaths, bg_death,
        5.0; case_analysed_daily = analysed_case,
        case_suspected_daily = susp_case)
    sb = returned(base, rand(MersenneTwister(seed), base))
    @test sb.scaling > 0
    @test 0 < sb.τ_death < 1
    @test sb.expected_confirmed_deaths >= 0
    twice = confirmed_deaths_model(17, 246, deaths_daily, bvd_deaths, bg_death,
        5.0; case_analysed_daily = 2 .* analysed_case,
        case_suspected_daily = susp_case)
    s2 = returned(twice, rand(MersenneTwister(seed), twice))
    @test s2.expected_confirmed_deaths ≈
          2 * sb.expected_confirmed_deaths rtol=1e-6
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

@testitem "confirmed deaths break day is sampled at composer level" begin
    using BVDOutbreakSize: confirmed_deaths_only_model, bvd_joint
    using Turing: logjoint
    using Random: MersenneTwister

    ## The composers rename the submodel's `confirmed_break_gross` to a
    ## per-stream `confirmed_break_gross_deaths`, so exercising the deaths break
    ## only through the submodel leaves the composer keyword unverified. A
    ## late vintage whose increment is dominated by a harmonisation: with the
    ## day declared, a `cdeath_step` is sampled and the log-density is finite.
    hist = (; days = [20, 30, 35, 40], counts = [200, 240, 476, 490])
    m = confirmed_deaths_only_model(40, 490, 800;
        deaths_history = (; days = [20, 40], counts = [400, 800]),
        confirmed_deaths_history = hist,
        confirmed_break_days = [35],
        confirmed_break_gross_deaths = [62])
    θ = rand(MersenneTwister(1), m)
    @test any(k -> occursin("cdeath_step", string(k)), keys(θ))
    @test isfinite(logjoint(m, θ))

    ## Undeclared, no step is sampled: the block is opt-in.
    m0 = confirmed_deaths_only_model(40, 490, 800;
        deaths_history = (; days = [20, 40], counts = [400, 800]),
        confirmed_deaths_history = hist)
    θ0 = rand(MersenneTwister(1), m0)
    @test !any(k -> occursin("cdeath_step", string(k)), keys(θ0))
    @test isfinite(logjoint(m0, θ0))

    ## `bvd_joint` forwards each confirmed stream its OWN gross vector, so the
    ## cases and deaths steps are both sampled and neither keyword is crossed.
    conf = (; days = [20, 30, 35, 40], counts = [500, 600, 969, 990])
    j = bvd_joint(40, missing, 800, missing, missing, 990, missing;
        confirmed_deaths = 490,
        deaths_history = (; days = [20, 40], counts = [400, 800]),
        confirmed_history = conf,
        confirmed_deaths_history = hist,
        lab_history = (; days = [10, 20], counts = [100, 300]),
        lab_daily_history = (; days = [30, 35, 40], counts = [50, 414, 60]),
        confirmed_break_days = [35],
        confirmed_break_gross_cases = [97],
        confirmed_break_gross_deaths = [62])
    θj = rand(MersenneTwister(1), j)
    ks = string.(keys(θj))
    @test any(k -> occursin("confirmed_step", k), ks)
    @test any(k -> occursin("cdeath_step", k), ks)
    @test isfinite(logjoint(j, θj))
end

@testitem "a deaths predictive replicates the break day" begin
    using BVDOutbreakSize: confirmed_deaths_only_model
    using Random: MersenneTwister

    ## Nulling the cut-off scalar is the generator gate, so `predict` resamples
    ## the increments while the dated history still supplies the vintage grid
    ## AND the published discrepancy the step is centred on. Emptying the
    ## history instead loses the discrepancy, which is why the scalar is the
    ## gate on this stream as it already is on the cases stream.
    hist = (; days = [20, 30, 35, 40], counts = [200, 240, 476, 490])
    gen = confirmed_deaths_only_model(40, missing, 800;
        deaths_history = (; days = [20, 40], counts = [400, 800]),
        confirmed_deaths_history = hist,
        confirmed_break_days = [35],
        confirmed_break_gross_deaths = [62])
    θ = rand(MersenneTwister(1), gen)

    ## The step survives into the generator, so a predictive carries the same
    ## dimensions as the fitted chain rather than dropping the break.
    @test any(k -> occursin("cdeath_step", string(k)), keys(θ))
end

@testitem "deaths_model background is the lagged scaled case background" begin
    using BVDOutbreakSize: deaths_model, background_walk_model, convolve_delay
    using Turing: returned
    using Random: MersenneTwister

    ## The joint passes the suspected-case background `case_bg_daily` (the
    ## smooth, gated, ramped daily random walk) into `deaths_model`; the death
    ## background is that background scaled by `cfr_bg` and lagged by the
    ## onset-to-death delay, so a background death follows its background case.
    ## It inherits the case background's gating (zero before the onset, since
    ## the delay only shifts mass later) and its smoothness.
    n = 60
    onset = 20
    seed = 5
    bgm = background_walk_model(n, 0.1; onset = onset)
    case_bg = returned(bgm, rand(MersenneTwister(seed), bgm)).λ
    onsets = [40.0 * exp(-((t - 35) / 10)^2) for t in 1:n]
    m = deaths_model((; days = Int[], counts = Int[]), missing, onsets, 5.0;
        case_bg_daily = case_bg)
    st = returned(m, rand(MersenneTwister(seed), m))

    ## The formula above holds pointwise, not just on the totals.
    @test st.bg_death_daily ≈ st.cfr_bg .* convolve_delay(case_bg, st.od_pmf)
    ## Gated: zero before the case-background onset (the lag only shifts
    ## later).
    @test all(st.bg_death_daily[1:(onset - 1)] .== 0)
    ## Smooth: no day-to-day jump exceeds the tight random-walk innovation
    ## scale (a per-vintage step background would jump abruptly).
    Δ = diff(st.bg_death_daily[onset:end])
    level = maximum(st.bg_death_daily)
    @test maximum(abs.(Δ)) <= 0.3 * level
    ## A positive case background gives a positive death background (when cfr_bg
    ## is non-trivial); both totals are finite.
    @test st.bg_death_total >= 0
    @test isfinite(st.bg_death_total)
end
