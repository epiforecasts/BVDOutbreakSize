## Tests for load_observations: returns the new grid-based fields from
## the bundled data/observations.toml file.

@testitem "load_observations returns the documented fields" begin
    using BVDOutbreakSize: load_observations, confirmed_positivity_windows
    using Dates: Date
    obs = load_observations()
    @test obs isa NamedTuple

    ## Grid dimensions
    @test obs.n isa Integer
    @test obs.n > 0
    @test obs.cutoff isa Date
    @test obs.seeding isa Date
    @test obs.cutoff >= obs.seeding
    @test obs.who_first_sitrep_days isa Integer
    @test obs.who_first_sitrep_days >= 1

    ## Cumulative stream totals
    @test obs.exported_cases isa Integer
    @test obs.exports_deaths isa Integer
    @test obs.total_deaths isa Integer
    @test obs.reported_cases isa Integer
    @test obs.confirmed_cases isa Integer
    @test obs.exported_cases >= 0
    @test obs.exports_deaths >= 0
    @test obs.total_deaths >= 0
    @test obs.reported_cases >= 0
    @test obs.confirmed_cases >= 0
    ## Tests-analysed is an optional scalar laboratory stream.
    @test ismissing(obs.tests_analysed) ||
          (obs.tests_analysed isa Integer && obs.tests_analysed >= 0)

    ## Per-vintage histories: named tuples with `days` and `counts`
    for key in (:deaths_history, :reported_history, :confirmed_history,
        :lab_history, :lab_daily_history, :suspected_daily_history,
        :suspected_daily_deaths_history,
        :isolation_history, :bed_capacity_history, :recovered_history)
        h = getproperty(obs, key)
        @test h isa NamedTuple
        @test hasproperty(h, :days)
        @test hasproperty(h, :counts)
        @test h.days isa AbstractVector{<:Integer}
        @test h.counts isa AbstractVector{<:Integer}
        @test length(h.days) == length(h.counts)
    end

    ## The post-cutoff 24h analysed series carries the late-window
    ## denominators, sorted oldest-first by day index and within the grid.
    @test !isempty(obs.lab_daily_history.counts)
    @test issorted(obs.lab_daily_history.days)
    @test all(d -> 1 <= d <= obs.n, obs.lab_daily_history.days)
    @test all(c -> c > 0, obs.lab_daily_history.counts)

    ## Validity invariant for the fitted stream rather than a literal echo of
    ## the data file: each day that anchors a late confirmed-positivity window
    ## scores the day's confirmed-case increment as a Binomial of its 24h
    ## analysed denominator, so the increment must not exceed the denominator
    ## (you cannot confirm more specimens than you analysed). A day breaching
    ## this would be silently clamped and is the one documented failure mode,
    ## so guarding it here catches a bad data update where the brittle literal
    ## would only have needed re-typing.
    let w = confirmed_positivity_windows(obs.confirmed_history,
            obs.lab_history, obs.lab_daily_history)
        anchored = findall(>(0), w.late_analysed)
        @test !isempty(anchored)
        @test all(i -> w.late_increments[i] <= w.late_analysed[i], anchored)
    end

    ## The post-cutoff daily new-suspect inflow ("nouveaux cas suspects du
    ## jour" / "cas suspects du jour"): a per-day incidence (never
    ## cumulative), so every count is a positive daily inflow, sorted
    ## oldest-first and within the grid. An invariant rather than a literal
    ## echo of the data file, which would break on every SitRep advance.
    @test !isempty(obs.suspected_daily_history.counts)
    @test all(c -> c > 0, obs.suspected_daily_history.counts)
    @test issorted(obs.suspected_daily_history.days)
    @test all(d -> 1 <= d <= obs.n, obs.suspected_daily_history.days)
    ## Strictly after the last cumulative suspected vintage (26 May): the two
    ## suspected likelihoods cover disjoint days, so they do not double-count.
    @test minimum(obs.suspected_daily_history.days) >
          maximum(obs.reported_history.days)

    ## The post-cutoff daily new suspected deaths ("cas suspects du jour N (M
    ## deces)"): the deaths analogue of the daily new-suspect inflow, a per-day
    ## count (never cumulative), so every count is positive, sorted oldest-first
    ## and within the grid. An invariant rather than a literal echo of the data
    ## file, which would break on every SitRep advance.
    @test !isempty(obs.suspected_daily_deaths_history.counts)
    @test all(c -> c > 0, obs.suspected_daily_deaths_history.counts)
    @test issorted(obs.suspected_daily_deaths_history.days)
    @test all(d -> 1 <= d <= obs.n, obs.suspected_daily_deaths_history.days)
    ## Strictly after the last cumulative suspected-death vintage (26 May): the
    ## two suspected-death likelihoods cover disjoint days, so they do not
    ## double-count.
    @test minimum(obs.suspected_daily_deaths_history.days) >
          maximum(obs.deaths_history.days)

    ## The daily isolation/treatment-bed occupancy ("Patients en isolement"):
    ## a per-day STOCK (point prevalence, never cumulative), so every count is
    ## a positive occupancy, sorted oldest-first and within the grid. An
    ## invariant rather than a literal echo, which would break on every
    ## SitRep advance. The fitted series starts at the all-patients column
    ## (1 June, SitRep 018), strictly after the frozen cumulative suspected
    ## vintages, so it shares no day with the cumulative reported series.
    @test !isempty(obs.isolation_history.counts)
    @test all(c -> c > 0, obs.isolation_history.counts)
    @test issorted(obs.isolation_history.days)
    @test all(d -> 1 <= d <= obs.n, obs.isolation_history.days)
    @test minimum(obs.isolation_history.days) >
          maximum(obs.reported_history.days)

    ## The implied bed-capacity series (occupancy / reported occupancy rate):
    ## positive bed counts, sorted oldest-first and within the grid. An
    ## invariant rather than a literal echo. Capacity must exceed the
    ## same-day occupancy (the occupancy rate is below 100%).
    @test !isempty(obs.bed_capacity_history.counts)
    @test all(c -> c > 0, obs.bed_capacity_history.counts)
    @test issorted(obs.bed_capacity_history.days)
    @test all(d -> 1 <= d <= obs.n, obs.bed_capacity_history.days)
    let isod = Dict(zip(obs.isolation_history.days,
            obs.isolation_history.counts))
        for (d, cap) in zip(obs.bed_capacity_history.days,
            obs.bed_capacity_history.counts)
            haskey(isod, d) && (@test cap >= isod[d])
        end
    end

    ## The cumulative recovered-among-confirmed series ("cumul guéris"): a
    ## CUMULATIVE count, so non-decreasing and positive, sorted oldest-first
    ## and within the grid, ending at the `recovered_cases` cut-off scalar.
    ## An invariant rather than a literal echo.
    @test !isempty(obs.recovered_history.counts)
    @test all(c -> c > 0, obs.recovered_history.counts)
    @test issorted(obs.recovered_history.counts)
    @test issorted(obs.recovered_history.days)
    @test all(d -> 1 <= d <= obs.n, obs.recovered_history.days)
    @test obs.recovered_cases == obs.recovered_history.counts[end]

    ## Tableau 6 treatment-centre patient-movement flows (per-day counts).
    ## Inflow/outflow streams are positive; absconded may be zero on a day.
    ## All sorted oldest-first and within the grid.
    for h in (obs.treatment_admissions_history, obs.treatment_deaths_history,
        obs.treatment_ruleout_history)
        @test !isempty(h.counts)
        @test all(c -> c > 0, h.counts)
        @test issorted(h.days)
        @test all(d -> 1 <= d <= obs.n, h.days)
    end
    @test !isempty(obs.treatment_absconded_history.counts)
    @test all(c -> c >= 0, obs.treatment_absconded_history.counts)
    @test issorted(obs.treatment_absconded_history.days)
    @test all(d -> 1 <= d <= obs.n, obs.treatment_absconded_history.days)
    ## Flow days fall within the occupancy window (they refine it).
    @test maximum(obs.treatment_admissions_history.days) <=
          maximum(obs.isolation_history.days)

    ## History day indices are in range
    dh = obs.deaths_history
    if !isempty(dh.days)
        @test all(1 .<= dh.days .<= obs.n)
        @test issorted(dh.days)
    end

    ## Suspected-death history (cumulative national, frozen 18-26 May): a
    ## non-decreasing positive series whose final vintage is the frozen
    ## `total_deaths` scalar. The 23 May zone-row correction is pinned (so a
    ## re-scan regression is caught) without echoing the whole vector: the
    ## 23 May deaths use the SitRep 009 zone-row sum (220), not the erroneous
    ## 119 headline. The 26 May revision to 246 (SitRep 012_v2, correcting
    ## the original 238) is pinned by the endpoint tie to `total_deaths`.
    @test !isempty(dh.counts)
    @test all(c -> c > 0, dh.counts)
    @test issorted(dh.counts)
    @test dh.counts[end] == obs.total_deaths
    @test 220 in dh.counts
    @test !(119 in dh.counts)

    ## Confirmed-case history (cumulative national, runs to the latest
    ## recorded vintage): non-decreasing and positive, with its final
    ## vintage equal to the `confirmed_cases` cut-off total. Grows with each
    ## SitRep, so an invariant rather than a literal every update re-types.
    @test !isempty(obs.confirmed_history.counts)
    @test all(c -> c > 0, obs.confirmed_history.counts)
    @test issorted(obs.confirmed_history.counts)
    @test obs.confirmed_history.counts[end] == obs.confirmed_cases

    ## Confirmed deaths (cumulative national, recorded for completeness):
    ## non-decreasing and positive, final vintage equal to the
    ## `confirmed_deaths` cut-off scalar. Grows with each SitRep.
    @test obs.confirmed_deaths isa Integer
    @test !isempty(obs.confirmed_deaths_history.counts)
    @test all(c -> c > 0, obs.confirmed_deaths_history.counts)
    @test issorted(obs.confirmed_deaths_history.counts)
    @test obs.confirmed_deaths_history.counts[end] == obs.confirmed_deaths

    ## Laboratory throughput histories (cumulative national, frozen
    ## 23-28 May): both non-decreasing and positive; specimens analysed
    ## cannot exceed specimens received vintage-by-vintage, and the analysed
    ## series ends at the cut-off `tests_analysed` total, the per-test
    ## positivity denominator.
    @test !isempty(obs.lab_history.counts)
    @test all(c -> c > 0, obs.tests_received_history.counts)
    @test all(c -> c > 0, obs.lab_history.counts)
    @test issorted(obs.tests_received_history.counts)
    @test issorted(obs.lab_history.counts)
    @test length(obs.lab_history.counts) ==
          length(obs.tests_received_history.counts)
    @test all(obs.lab_history.counts .<= obs.tests_received_history.counts)
    @test obs.lab_history.counts[end] == obs.tests_analysed

    ## Dated Uganda export series: grid day-indices of the three imports
    ## (11, 16, 23 May) and the single export death (14 May), sorted
    ## ascending and within the grid. Their lengths match the scalar totals.
    @test obs.export_case_days isa AbstractVector{<:Integer}
    @test obs.export_death_days isa AbstractVector{<:Integer}
    @test issorted(obs.export_case_days)
    @test issorted(obs.export_death_days)
    @test all(1 .<= obs.export_case_days .<= obs.n)
    @test all(1 .<= obs.export_death_days .<= obs.n)
    @test length(obs.export_case_days) == 3
    @test length(obs.export_death_days) == 1
    ## Detection days are strictly increasing (distinct import dates,
    ## 11→16→23 May).
    @test all(>(0), diff(obs.export_case_days))
    ## The export death (14 May) falls between imports #1 (11 May) and #2.
    @test obs.export_case_days[1] < obs.export_death_days[1] <
          obs.export_case_days[2]

    ## Manual occupancy break days: an opt-in dated list within the grid,
    ## sorted ascending. The 19 June DHIS2 reclassification is listed, so it
    ## resolves to that grid day and falls inside the isolation window.
    @test obs.occupancy_break_days isa AbstractVector{<:Integer}
    @test issorted(obs.occupancy_break_days)
    @test all(1 .<= obs.occupancy_break_days .<= obs.n)
    @test obs.occupancy_break_days ==
          [obs.n - (Date(obs.cutoff) - Date("2026-06-19")).value]
    @test minimum(obs.isolation_history.days) <=
          obs.occupancy_break_days[1] <= maximum(obs.isolation_history.days)

    ## Genetic TMRCA bound
    @test !ismissing(obs.tmrca_days)
    @test obs.tmrca_days isa Real
    @test obs.tmrca_days > 0

    ## Breakpoint day is consistent: n - who_first_sitrep_days
    breakpoint = obs.n - obs.who_first_sitrep_days
    @test breakpoint >= 1
    @test breakpoint <= obs.n
end

@testitem "load_observations histories have consistent counts" begin
    using BVDOutbreakSize: load_observations
    obs = load_observations()

    ## Cumulative counts in histories should be non-decreasing and bounded
    ## by the cut-off total
    dh = obs.deaths_history
    if length(dh.counts) > 1
        @test issorted(dh.counts)
    end
    if !isempty(dh.counts)
        @test dh.counts[end] <= obs.total_deaths
    end

    rh = obs.reported_history
    if length(rh.counts) > 1
        @test issorted(rh.counts)
    end
    if !isempty(rh.counts)
        @test rh.counts[end] <= obs.reported_cases
    end
end

@testitem "freeze_observations truncates to a past cut-off" begin
    using BVDOutbreakSize: load_observations, freeze_observations
    using Dates: Date

    full = load_observations()
    frozen = freeze_observations("2026-05-23")

    ## The cut-off moves to the freeze date; the grid shrinks by the
    ## number of days dropped.
    @test frozen.cutoff == Date("2026-05-23")
    @test frozen.n < full.n
    @test frozen.n == full.n - (Date(full.cutoff) - Date("2026-05-23")).value

    ## Every retained vintage is dated on or before the freeze date, so
    ## no history extends past the cut-off grid.
    for key in (:reported_history, :deaths_history, :confirmed_history,
        :lab_history)
        h = getproperty(frozen, key)
        isempty(h.days) && continue
        @test all(1 .<= h.days .<= frozen.n)
        @test maximum(h.days) <= frozen.n
    end

    ## The 23 May suspected streams: six vintages (18-23 May) of the nine
    ## in the full manifest. Truncation is a clean non-decreasing prefix of
    ## the full series (asserted against the loaded full data, not a hand-
    ## typed literal) and pins the cut-off scalars to the 23 May values.
    @test length(frozen.reported_history.counts) == 6
    @test length(frozen.deaths_history.counts) == 6
    @test frozen.reported_history.counts == full.reported_history.counts[1:6]
    @test frozen.deaths_history.counts == full.deaths_history.counts[1:6]
    @test issorted(frozen.reported_history.counts)
    @test issorted(frozen.deaths_history.counts)
    @test frozen.reported_history.counts[end] == frozen.reported_cases
    @test frozen.deaths_history.counts[end] == frozen.total_deaths
    @test frozen.reported_cases == 904
    @test frozen.total_deaths == 220

    ## The cut-off scalars come from the truncated history, not the
    ## manifest's full-data totals.
    @test frozen.reported_cases < full.reported_cases
    @test frozen.total_deaths < full.total_deaths
    @test frozen.confirmed_cases == frozen.confirmed_history.counts[end]

    ## The dated export series are truncated to detections on or before the
    ## freeze date. 23 May keeps all three imports (last detection 23 May)
    ## and the single export death (14 May); the export scalars come from
    ## the truncated dated count.
    @test length(frozen.export_case_days) == 3
    @test length(frozen.export_death_days) == 1
    @test frozen.exported_cases == 3
    @test frozen.exports_deaths == 1
    @test all(1 .<= frozen.export_case_days .<= frozen.n)
    @test frozen.export_case_days[end] == frozen.n   # last import = cut-off

    ## Freezing before the last import drops the post-cut-off detection.
    ## 18 May keeps imports #1 (11 May) and #2 (16 May) plus the death
    ## (14 May), but not import #3 (23 May).
    early = freeze_observations("2026-05-18")
    @test length(early.export_case_days) == 2        # 11 and 16 May imports
    @test early.exported_cases == 2
    @test length(early.export_death_days) == 1        # 14 May ≤ 18 May
    @test early.exports_deaths == 1
    @test all(1 .<= early.export_death_days .<= early.n)

    ## A Date argument is equivalent to the ISO string.
    @test freeze_observations(Date("2026-05-23")).n == frozen.n
end
