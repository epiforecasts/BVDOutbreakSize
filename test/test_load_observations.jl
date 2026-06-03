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
        :lab_history)
        h = getproperty(obs, key)
        @test h isa NamedTuple
        @test hasproperty(h, :days)
        @test hasproperty(h, :counts)
        @test h.days isa AbstractVector{<:Integer}
        @test h.counts isa AbstractVector{<:Integer}
        @test length(h.days) == length(h.counts)
    end

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

    ## Confirmed-case history runs to the 28 May cut-off; its final
    ## vintage equals the cut-off `confirmed_cases` total.
    @test obs.confirmed_history.counts ==
          [33, 51, 57, 79, 83, 101, 105, 106, 121, 125, 210]
    @test obs.confirmed_history.counts[end] == obs.confirmed_cases

    ## Confirmed deaths: recorded, flat at 17 over the fitted window.
    @test obs.confirmed_deaths isa Integer
    @test obs.confirmed_deaths_history.counts == [17, 17, 17]
    @test obs.confirmed_deaths_history.counts[end] == obs.confirmed_deaths

    ## Laboratory throughput histories (cumulative national, 23-28 May);
    ## the analysed series (`lab_history`) ends at the cut-off
    ## `tests_analysed` total and is the per-test positivity denominator.
    @test obs.tests_received_history.counts ==
          [418, 431, 431, 662, 774, 883]
    @test obs.lab_history.counts == [211, 295, 295, 403, 648, 755]
    @test obs.lab_history.counts[end] == obs.tests_analysed

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
