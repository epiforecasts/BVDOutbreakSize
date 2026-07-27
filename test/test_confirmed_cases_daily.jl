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

@testitem "break_step_centres derives the step from the printed gross count" begin
    using BVDOutbreakSize: break_step_centres

    ## The centre is `observed increment - printed 24h count`: the part of the
    ## vintage step the report attributes to base integration. 22 July 2026:
    ## 369 - 97 = 272 cases, 236 - 62 = 174 deaths.
    days, centres = break_step_centres([21, 22, 23], [10, 369, 15], [22], [97])
    @test days == [22]
    @test centres == [272.0]

    ## Deaths use the same helper with their own gross count.
    _, dcentres = break_step_centres([21, 22, 23], [10, 236, 15], [22], [62])
    @test dcentres == [174.0]

    ## A break day matching no window is dropped, so no inert step is sampled.
    @test break_step_centres([21, 23], [10, 15], [22], [97]) == (Int[], Float64[])

    ## A missing/short gross is read as a gross of zero, so the centre is the
    ## WHOLE increment: all of it attributed to the artefact rather than
    ## erroring. Conservative, not neutral — a centre near zero over a
    ## de-anchored day is the configuration that wrecks the geometry.
    @test break_step_centres([21, 22], [10, 369], [22], Int[])[2] == [369.0]

    ## No break days at all is empty.
    @test break_step_centres([21, 22], [10, 369], Int[], Int[]) ==
          (Int[], Float64[])

    ## Generator mode: no counts to difference, so there is no published
    ## discrepancy. The day still gets a step (a predictive keeps the fitted
    ## chain's dimensions) but the centre falls back to zero.
    @test break_step_centres([21, 22, 23], missing, [22], [97]) ==
          ([22], [0.0])
    @test break_step_centres([21, 23], missing, [22], [97]) ==
          (Int[], Float64[])
end

@testitem "break step centre and offset address the same window" begin
    using BVDOutbreakSize: break_step_centres, confirmed_break_offset

    ## The centre is read out of `increments` BY POSITION while the offset is
    ## written back BY DAY, so handing the two helpers day vectors that do not
    ## correspond to the same increments shifts a step onto its neighbour.
    ## Per-day increments that all differ make the pairing observable: only the
    ## break day's own increment can produce its centre, and the offset must be
    ## non-zero on exactly the index that increment came from.
    days = [30, 35, 40]
    increments = [11, 369, 22]
    bdays, centres = break_step_centres(days, increments, [35], [97])
    @test bdays == [35]
    @test centres == [272.0]
    pos = findfirst(==(35), days)
    Δ = confirmed_break_offset(days, bdays, centres)
    @test Δ == [0.0, 272.0, 0.0]
    @test Δ[pos] == increments[pos] - 97

    ## Shifting the declared day by one moves BOTH the centre and the offset,
    ## so a mismatch cannot hide behind a coincidentally equal step: the
    ## neighbour's centre is its own increment, not 22 July's.
    nbdays, ncentres = break_step_centres(days, increments, [40], [97])
    @test ncentres == [Float64(increments[3] - 97)]
    @test confirmed_break_offset(days, nbdays, ncentres) ==
          [0.0, 0.0, Float64(increments[3] - 97)]
end

@testitem "the listed break day is a live late window, not a silent no-op" begin
    using BVDOutbreakSize: load_observations, confirmed_positivity_windows

    ## Both mechanisms (the de-anchor and the fitted step) search `late_days`
    ## only, so if a later cumulative laboratory vintage ever moved the listed
    ## day out of the late group they would both vanish with no error and no
    ## warning. Assert the wiring against the real data rather than assume it.
    obs = load_observations()
    @test !isempty(obs.confirmed_break_days)
    w = confirmed_positivity_windows(obs.confirmed_history, obs.lab_history,
        obs.lab_daily_history, obs.confirmed_break_days)
    for d in obs.confirmed_break_days
        pos = findfirst(==(d), w.late_days)
        @test pos !== nothing

        ## De-anchored: the day publishes a 24h analysed count, and the break
        ## drops it so the window cannot enter the BetaBinomial as same-day
        ## positives at an implausible positivity.
        @test d in obs.lab_daily_history.days
        @test w.late_analysed[pos] == 0

        ## The step is centred on the vintage increment the model actually
        ## scores, which must exceed the printed gross count for the day to be
        ## a harmonisation at all.
        @test w.late_increments[pos] > 0
    end

    ## Left undeclared the same day keeps its published denominator, so the
    ## assertions above are testing the break and not a vacuous zero.
    w0 = confirmed_positivity_windows(obs.confirmed_history, obs.lab_history,
        obs.lab_daily_history)
    for d in obs.confirmed_break_days
        @test w0.late_analysed[findfirst(==(d), w0.late_days)] > 0
    end
end

@testitem "a zero break sd is a deterministic correction" begin
    using BVDOutbreakSize: confirmed_only_model
    using Turing: logjoint
    using Random: MersenneTwister

    ## `confirmed_break_sd = 0` takes the printed 24h count as exact: the step
    ## is pinned at the published discrepancy and no parameter is sampled, so
    ## the correction cannot cost any sampling geometry.
    hist = (; days = [20, 30, 35, 40], counts = [50, 60, 429, 444])
    lab = (; days = [10, 20], counts = [100, 300])
    daily = (; days = [30, 35, 40], counts = [50, 414, 60])
    build(sd, gross) = confirmed_only_model(40, 444;
        confirmed_history = hist, lab_history = lab,
        lab_daily_history = daily, confirmed_break_days = [35],
        confirmed_break_gross_cases = gross, confirmed_break_sd = sd)

    m = build(0.0, [97])
    θ = rand(MersenneTwister(1), m)
    @test !any(k -> occursin("confirmed_step", string(k)), keys(θ))
    @test isfinite(logjoint(m, θ))

    ## The pinned step still reaches the likelihood: the same draw scores
    ## differently once the declared gross changes the centre. Day 35's
    ## increment is 429 - 60 = 369, so a gross of 369 pins the step at zero
    ## while 97 pins it at 272.
    @test logjoint(build(0.0, [369]), θ) != logjoint(m, θ)
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
        lab_daily_history = daily, confirmed_break_days = [35],
        confirmed_break_gross_cases = [97])
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
