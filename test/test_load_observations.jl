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
