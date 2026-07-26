@testitem "lab pipeline daily likelihood runs and is finite" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing: logjoint
    using Random: MersenneTwister

    ## The confirmed-cases likelihood is exercised end to end through
    ## confirmed_only_model, which draws the shared report kernel, fits the
    ## analysed-specimen volume and scores the confirmed positives as a
    ## Binomial of the observed analysed denominator.
    m = confirmed_only_model(40, 8;
        confirmed_history = (; days = [20, 40], counts = [3, 8]),
        lab_history = (; days = [20, 40], counts = [5, 9]))
    lp = logjoint(m, rand(MersenneTwister(1), m))
    @test isfinite(lp)
end

@testitem "confirmed_positivity_windows anchors late days on 24h analysed" begin
    using BVDOutbreakSize: confirmed_positivity_windows

    ## Last cumulative laboratory date is day 20, so days 21-23 are late
    ## late windows. A published 24h analysed count on day 22 flags that
    ## window for the observed-denominator Binomial; the others stay unanchored.
    confirmed = (; days = [5, 20, 21, 22, 23], counts = [10, 50, 60, 75, 90])
    lab = (; days = [10, 20], counts = [100, 300])
    daily = (; days = [22], counts = [40])
    w = confirmed_positivity_windows(confirmed, lab, daily)
    @test w.late_days == [21, 22, 23]
    @test w.late_increments == [10, 15, 15]
    @test w.late_analysed == [0, 40, 0]

    ## With no daily series every late day is unanchored (no 24h count), and the
    ## field is present and aligned with `late_days`.
    w0 = confirmed_positivity_windows(confirmed, lab)
    @test length(w0.late_analysed) == length(w0.late_days)
    @test all(==(0), w0.late_analysed)

    ## A daily count on or before the last laboratory date is ignored: the
    ## cumulative series already covers it.
    w1 = confirmed_positivity_windows(confirmed, lab,
        (; days = [15, 22], counts = [99, 40]))
    @test w1.late_analysed == [0, 40, 0]
end

@testitem "lab pipeline daily anchor scores the 24h day as Binomial" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing: logjoint
    using Random: MersenneTwister

    ## Post-cutoff late vintages (days 30, 35, 40 after the last lab date
    ## 20) where day 35 publishes a 24h analysed count: it is scored as a
    ## Binomial of that observed denominator, the others as modelled-volume
    ## unanchored windows.
    m = confirmed_only_model(40, 20;
        confirmed_history =
        (; days = [20, 30, 35, 40], counts = [5, 9, 14, 20]),
        lab_history = (; days = [10, 20], counts = [12, 28]),
        lab_daily_history = (; days = [35], counts = [30]))
    lp = logjoint(m, rand(MersenneTwister(1), m))
    @test isfinite(lp)

    ## When the confirmed increment exceeds the observed denominator it is
    ## clamped into the Binomial support, so the likelihood stays finite.
    m2 = confirmed_only_model(40, 20;
        confirmed_history =
        (; days = [20, 30, 35, 40], counts = [5, 9, 14, 20]),
        lab_history = (; days = [10, 20], counts = [12, 28]),
        lab_daily_history = (; days = [35], counts = [3]))
    @test isfinite(logjoint(m2, rand(MersenneTwister(2), m2)))
end

@testitem "confirmed break days de-anchor the positivity denominator" begin
    using BVDOutbreakSize: confirmed_positivity_windows

    ## A retrospective harmonisation day carries a published 24h analysed
    ## count, but most of its confirmed increment is a reattached provincial
    ## base rather than positives out of that day's specimens. Declaring it a
    ## break day drops the denominator so the window falls back to the
    ## modelled-volume branch instead of entering the BetaBinomial as a
    ## near-100% positivity observation. Mirrors 22 July 2026 (SitRep 069):
    ## +369 confirmed against 414 analysed, an 89% implied positivity.
    confirmed = (; days = [5, 20, 21, 22, 23], counts = [10, 50, 60, 429, 444])
    lab = (; days = [10, 20], counts = [100, 300])
    daily = (; days = [21, 22, 23], counts = [50, 414, 60])

    ## Without the break day, day 22 is anchored on its 414 denominator and
    ## the harmonised +369 becomes the Binomial numerator.
    w = confirmed_positivity_windows(confirmed, lab, daily)
    @test w.late_days == [21, 22, 23]
    @test w.late_increments == [10, 369, 15]
    @test w.late_analysed == [50, 414, 60]

    ## Declaring day 22 a break day zeroes ONLY that denominator; the
    ## neighbouring anchored days are untouched, and the increments (the
    ## observed data) are not rewritten — the step absorbs them in the model.
    wb = confirmed_positivity_windows(confirmed, lab, daily, [22])
    @test wb.late_analysed == [50, 0, 60]
    @test wb.late_increments == w.late_increments
    @test wb.late_days == w.late_days

    ## Empty break days are a no-op, and a break day that matches no late
    ## window is ignored rather than shifting another day's denominator.
    @test confirmed_positivity_windows(confirmed, lab, daily, Int[]).late_analysed ==
          w.late_analysed
    @test confirmed_positivity_windows(confirmed, lab, daily, [99]).late_analysed ==
          w.late_analysed
end

@testitem "confirmed_break_offset places one step on each break window" begin
    using BVDOutbreakSize: confirmed_break_offset

    ## The confirmed likelihood fits INCREMENTS, so a one-off base
    ## integration inflates exactly one window: the offset is per-window, not
    ## cumulative (unlike the occupancy reclassification offset).
    late_days = [21, 22, 23, 24]
    @test confirmed_break_offset(late_days, [22], [250.0]) == [0.0, 250.0, 0.0, 0.0]

    ## Two break days each get their own step, and later windows stay at zero.
    @test confirmed_break_offset(late_days, [21, 24], [10.0, -5.0]) ==
          [10.0, 0.0, 0.0, -5.0]

    ## No break days (or none landing on a late window) is a zero no-op of the
    ## right length.
    @test confirmed_break_offset(late_days, Int[], Float64[]) == zeros(4)
    @test confirmed_break_offset(late_days, [99], [7.0]) == zeros(4)
    @test length(confirmed_break_offset(Int[], [22], [1.0])) == 0
end

@testitem "confirmed break day is sampled and keeps the likelihood finite" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing: logjoint
    using Random: MersenneTwister

    ## End-to-end: a late vintage whose increment is dominated by a
    ## harmonisation. With the day declared a break, a `confirmed_step` is
    ## sampled and the log-density stays finite.
    hist = (; days = [20, 30, 35, 40], counts = [50, 60, 429, 444])
    lab = (; days = [10, 20], counts = [100, 300])
    daily = (; days = [30, 35, 40], counts = [50, 414, 60])
    m = confirmed_only_model(40, 444;
        confirmed_history = hist, lab_history = lab,
        lab_daily_history = daily, confirmed_break_days = [35])
    θ = rand(MersenneTwister(1), m)
    @test any(k -> occursin("confirmed_step", string(k)), keys(θ))
    @test isfinite(logjoint(m, θ))

    ## Without a break day no step is sampled — the block is opt-in and empty
    ## by default, so an unlisted vintage keeps the previous behaviour.
    m0 = confirmed_only_model(40, 444;
        confirmed_history = hist, lab_history = lab,
        lab_daily_history = daily)
    θ0 = rand(MersenneTwister(1), m0)
    @test !any(k -> occursin("confirmed_step", string(k)), keys(θ0))
    @test isfinite(logjoint(m0, θ0))
end
