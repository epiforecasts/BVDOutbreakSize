## Tests for load_observations: returns the new grid-based fields from
## the bundled data/observations.toml file.

@testitem "load_observations returns the documented fields" begin
    using BVDOutbreakSize: load_observations
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
        :lab_history, :lab_daily_history, :suspected_daily_history)
        h = getproperty(obs, key)
        @test h isa NamedTuple
        @test hasproperty(h, :days)
        @test hasproperty(h, :counts)
        @test h.days isa AbstractVector{<:Integer}
        @test h.counts isa AbstractVector{<:Integer}
        @test length(h.days) == length(h.counts)
    end

    ## The post-cutoff 24h analysed series carries the trusted dark-window
    ## denominators (1, 4, 5, 6, 7, 8, 9 June), sorted oldest-first by day
    ## index.
    @test obs.lab_daily_history.counts == [76, 256, 126, 106, 67, 121, 68]
    @test issorted(obs.lab_daily_history.days)
    @test all(d -> d <= obs.n, obs.lab_daily_history.days)

    ## The post-cutoff daily new-suspect inflow ("nouveaux cas suspects du
    ## jour", 4-7 and 9 June; 8 June prints no figure), sorted oldest-first by
    ## day index and within the grid. A per-day incidence (never cumulative),
    ## so it begins where the frozen cumulative suspected series stops.
    @test obs.suspected_daily_history.counts == [153, 119, 117, 94, 119]
    @test issorted(obs.suspected_daily_history.days)
    @test all(d -> 1 <= d <= obs.n, obs.suspected_daily_history.days)
    ## Strictly after the last cumulative suspected vintage (26 May): the two
    ## suspected likelihoods cover disjoint days.
    @test minimum(obs.suspected_daily_history.days) >
          maximum(obs.reported_history.days)

    ## History day indices are in range
    dh = obs.deaths_history
    if !isempty(dh.days)
        @test all(1 .<= dh.days .<= obs.n)
        @test issorted(dh.days)
    end

    ## Suspected-death history, 18-26 May (nine vintages). The 23 May
    ## deaths use the SitRep 009 zone-row sum (220), not the erroneous 119
    ## headline; the 26 May value uses the revised SitRep 012_v2 (246).
    @test dh.counts == [131, 148, 160, 175, 204, 220, 223, 238, 246]

    ## Confirmed-case history runs to the 9 June recorded point (post-28
    ## May vintages have no analysed denominator); its final vintage equals
    ## the `confirmed_cases` total.
    @test obs.confirmed_history.counts ==
          [33, 51, 57, 79, 83, 101, 105, 106, 121, 125, 210,
        263, 282, 321, 344, 363, 381, 452, 488, 515, 550, 598, 635]
    @test obs.confirmed_history.counts[end] == obs.confirmed_cases

    ## Confirmed deaths: recorded, growing 17 → 127 over 26 May-9 June.
    @test obs.confirmed_deaths isa Integer
    @test obs.confirmed_deaths_history.counts ==
          [17, 17, 17, 42, 42, 48, 60, 62, 64, 82, 86, 91, 101, 115, 127]
    @test obs.confirmed_deaths_history.counts[end] == obs.confirmed_deaths

    ## Laboratory throughput histories (cumulative national, 23-28 May);
    ## the analysed series (`lab_history`) ends at the cut-off
    ## `tests_analysed` total and is the per-test positivity denominator.
    @test obs.tests_received_history.counts ==
          [418, 431, 431, 662, 774, 883]
    @test obs.lab_history.counts == [211, 295, 295, 403, 648, 755]
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
    ## Detection days are spaced 5 (11→16 May) then 7 (16→23 May) apart.
    @test diff(obs.export_case_days) == [5, 7]
    ## The export death (14 May) falls between imports #1 (11 May) and #2.
    @test obs.export_case_days[1] < obs.export_death_days[1] <
          obs.export_case_days[2]

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

    ## The 23 May suspected streams: six vintages (18-23 May) of the
    ## nine in the full manifest, ending at the frozen totals.
    @test frozen.reported_history.counts == [516, 575, 672, 745, 872, 904]
    @test frozen.deaths_history.counts == [131, 148, 160, 175, 204, 220]
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
