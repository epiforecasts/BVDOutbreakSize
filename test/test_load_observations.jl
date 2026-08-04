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
    ## a per-day stock (point prevalence, never cumulative), so every count is
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
    ## cumulative count, so non-decreasing and positive, sorted oldest-first
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
    ## sorted ascending. The 9 and 19 June and 14 July DHIS2
    ## reclassifications are listed, so each resolves to its grid day and
    ## falls inside the isolation window.
    @test obs.occupancy_break_days isa AbstractVector{<:Integer}
    @test issorted(obs.occupancy_break_days)
    @test all(1 .<= obs.occupancy_break_days .<= obs.n)
    @test obs.occupancy_break_days ==
          [obs.n - (Date(obs.cutoff) - Date(d)).value
           for d in (Date("2026-06-09"), Date("2026-06-19"),
        Date("2026-07-14"))]
    @test all(minimum(obs.isolation_history.days) .<=
              obs.occupancy_break_days .<= maximum(obs.isolation_history.days))

    ## Manual confirmed harmonisation-break days: the same opt-in dated-list
    ## shape, resolved onto the grid. 22 July 2026 (SitRep 069) is listed —
    ## the report itself states the cumulative step is mostly a provincial
    ## base integration rather than 24h notifications. An invariant plus the
    ## one pinned date, so a re-scan regression is caught without echoing a
    ## list that grows on every harmonisation.
    @test obs.confirmed_break_days isa AbstractVector{<:Integer}
    @test issorted(obs.confirmed_break_days)
    @test all(1 .<= obs.confirmed_break_days .<= obs.n)
    @test obs.confirmed_break_days ==
          [obs.n - (Date(obs.cutoff) - Date(d)).value
           for d in (Date("2026-07-22"),)]
    ## Every break day must land on a confirmed vintage, otherwise it fits an
    ## inert step against no observation.
    @test all(in(obs.confirmed_history.days), obs.confirmed_break_days)

    ## The printed 24h gross counts that make each step data-derived: aligned
    ## with the dates and filtered with them, positive, and strictly below the
    ## vintage increment they are subtracted from (otherwise the "discrepancy"
    ## the step centres on would be negative and the day is not a
    ## harmonisation at all).
    @test length(obs.confirmed_break_gross_cases) ==
          length(obs.confirmed_break_days)
    @test length(obs.confirmed_break_gross_deaths) ==
          length(obs.confirmed_break_days)
    @test all(>(0), obs.confirmed_break_gross_cases)
    @test all(>(0), obs.confirmed_break_gross_deaths)
    let cmap = Dict(zip(obs.confirmed_history.days,
            diff(vcat(0, collect(obs.confirmed_history.counts)))))
        for (d, g) in zip(obs.confirmed_break_days,
            obs.confirmed_break_gross_cases)
            @test g < cmap[d]
        end
    end

    ## Genetic TMRCA bound
    @test !ismissing(obs.tmrca_days)
    @test obs.tmrca_days isa Real
    @test obs.tmrca_days > 0

    ## Breakpoint day is consistent: n - who_first_sitrep_days
    breakpoint = obs.n - obs.who_first_sitrep_days
    @test breakpoint >= 1
    @test breakpoint <= obs.n
end

@testitem "confirmed break gross counts follow the sorted dates" begin
    using BVDOutbreakSize
    using BVDOutbreakSize: load_observations

    using Dates: Date, Day
    using TOML

    ## The grid days come back sorted while the TOML arrays keep the order
    ## they are written in, so one permutation has to carry the dates and both
    ## gross vectors together. Written out of date order, a per-vector filter
    ## would pair the earlier day with the later day's printed counts.
    ##
    ## The fixture is built by round-tripping the manifest through the TOML
    ## parser and rewriting only the break block, not by editing its text: text
    ## surgery matches neither a CRLF checkout nor a rewrapped block, and a
    ## fixture that silently fails to apply turns this into a test of nothing.
    path = joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml")
    raw = TOML.parsefile(path)
    blk = raw["confirmed_break_dates"]

    ## Derive both dates from the manifest's own listed day so the fixture does
    ## not pin a date that a later cut-off or data update invalidates, and
    ## write them descending so sorting has something to do.
    later = Date(String(blk["value"][1]))
    earlier = later - Day(1)
    blk["value"] = [string(later), string(earlier)]
    blk["gross_cases"] = [97, 11]
    blk["gross_deaths"] = [62, 22]
    @test Date(blk["value"][1]) > Date(blk["value"][2])

    tmp = joinpath(mktempdir(), "observations.toml")
    open(io -> TOML.print(io, raw), tmp, "w")
    obs = load_observations(tmp)

    ## The earlier day sorts first and must keep its own counts, not the later
    ## day's, which is exactly what a per-vector filter would hand it.
    @test length(obs.confirmed_break_days) == 2
    @test issorted(obs.confirmed_break_days)
    @test obs.confirmed_break_gross_cases == [11, 97]
    @test obs.confirmed_break_gross_deaths == [22, 62]
end

@testitem "a break day whose gross covers its increment is rejected" begin
    using BVDOutbreakSize
    using BVDOutbreakSize: load_observations
    using Dates: Date
    using TOML

    ## A day whose printed 24h count already covers its whole vintage step is
    ## not a harmonisation, and listing it is harmful rather than inert: the
    ## day is de-anchored from the positivity denominator while the step meant
    ## to absorb the backlog is centred at or below zero. Measured at 94
    ## divergences and a min bulk ESS of 15 against 20 and 522 with no break
    ## day, so the loader refuses it instead of sampling it.
    path = joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml")
    raw = TOML.parsefile(path)
    blk = raw["confirmed_break_dates"]
    day = String(blk["value"][1])

    ## That vintage's own increment, which the gross must stay below.
    ch = raw["confirmed_case_history"]
    idx = findfirst(==(day), String.(ch["dates"]))
    inc = Int(ch["values"][idx]) - Int(ch["values"][idx - 1])
    @test inc > 0

    write_toml(r) = (p = joinpath(mktempdir(), "observations.toml");
        open(io -> TOML.print(io, r), p, "w"); p)

    ## Equal to the increment is already too high: the centre would be zero.
    blk["gross_cases"] = [inc]
    @test_throws ErrorException load_observations(write_toml(raw))

    ## And above it, which would centre the step negatively.
    blk["gross_cases"] = [inc + 10]
    @test_throws ErrorException load_observations(write_toml(raw))

    ## Just below is accepted, so the bound is exactly `gross < increment`.
    blk["gross_cases"] = [inc - 1]
    obs = load_observations(write_toml(raw))
    @test obs.confirmed_break_gross_cases == [inc - 1]

    ## A zero gross is legal but warns: the whole increment becomes artefact.
    blk["gross_cases"] = [0]
    @test_logs (:warn, r"no printed 24h cases count") match_mode=:any begin
        load_observations(write_toml(raw))
    end

    ## The deaths stream is checked independently against its own increment, so
    ## a day can pass on cases and fail on deaths. The message names the stream.
    raw2 = TOML.parsefile(path)
    dh = raw2["confirmed_death_history"]
    didx = findfirst(==(day), String.(dh["dates"]))
    dinc = Int(dh["values"][didx]) - Int(dh["values"][didx - 1])
    raw2["confirmed_break_dates"]["gross_deaths"] = [dinc]
    derr = try
        load_observations(write_toml(raw2))
        nothing
    catch err
        err
    end
    @test derr isa ErrorException
    @test occursin("deaths", derr.msg)
    @test occursin(day, derr.msg)
end

@testitem "a break date matching no vintage is rejected" begin
    using BVDOutbreakSize
    using BVDOutbreakSize: load_observations
    using Dates: Date, Day
    using TOML

    ## The quietest failure of the three: a date that matches no vintage
    ## does nothing at all — no step, no de-anchor — and the gross bar
    ## cannot fire because there is no increment to compare against. A
    ## transposed digit or the wrong month therefore presents as silence
    ## while the user believes a harmonisation is being absorbed, so the
    ## loader refuses it.
    path = joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml")
    raw = TOML.parsefile(path)
    blk = raw["confirmed_break_dates"]
    listed = Date(String(blk["value"][1]))

    ## Find a date inside the history's own span that is not a vintage. There is
    ## at least one: SitRep 063 (16 July 2026) was never published upstream, so
    ## the confirmed series has no entry for it.
    dates = Set(String.(raw["confirmed_case_history"]["dates"]))
    first_vintage = Date(minimum(dates))
    absent = first(d for d in first_vintage:Day(1):listed
    if !(string(d) in dates))
    @test !(string(absent) in dates)

    blk["value"] = [string(absent)]
    blk["gross_cases"] = [10]
    blk["gross_deaths"] = [5]
    tmp = joinpath(mktempdir(), "observations.toml")
    open(io -> TOML.print(io, raw), tmp, "w")

    err = try
        load_observations(tmp)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin(string(absent), err.msg)
    @test occursin("matches no", err.msg)
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
