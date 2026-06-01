## Tests for load_observations: returns the documented fields from
## the bundled data/observations.toml file.

@testitem "load_observations returns the documented fields" begin
    using BVDOutbreakSize: load_observations
    obs = load_observations()
    @test obs isa NamedTuple
    @test obs.exported_cases isa Integer
    @test obs.exports_deaths isa Integer
    @test obs.total_deaths isa Integer
    @test obs.reported_cases isa Integer
    @test obs.confirmed_cases isa Integer
    @test obs.cumulative_tests_analysed isa Integer
    @test obs.source_population isa Integer
    @test obs.daily_outbound_travellers isa Real
    @test obs.daily_outbound_travellers_sd isa Real
    @test obs.genetic_tmrca_days isa Real
    @test obs.genetic_tmrca_days_sd isa Real
    @test obs.genetic_tmrca_alt_days isa Real
    @test obs.genetic_tmrca_alt_days_sd isa Real

    @test obs.exported_cases >= 0
    @test obs.exports_deaths >= 0
    @test obs.total_deaths >= 0
    @test obs.reported_cases >= 0
    @test obs.confirmed_cases >= 0
    @test obs.confirmed_cases <= obs.reported_cases
    @test obs.cumulative_tests_analysed >= obs.confirmed_cases
    @test obs.cumulative_tests_analysed <= obs.reported_cases
    @test obs.daily_outbound_travellers > 0
    @test obs.daily_outbound_travellers_sd > 0
    @test obs.source_population > 0
    @test obs.genetic_tmrca_days > 0
    @test obs.genetic_tmrca_days_sd > 0
    @test obs.genetic_tmrca_alt_days > 0
    @test obs.genetic_tmrca_alt_days_sd > 0
    ## The alternative (faster) clock dates the TMRCA more recently, so
    ## fewer days before the cut-off than the baseline estimate.
    @test obs.genetic_tmrca_alt_days < obs.genetic_tmrca_days

    @test obs.sources isa NamedTuple
    @test obs.sources.exported_cases isa String
    @test obs.sources.exports_deaths isa String
    @test obs.sources.total_deaths isa String
    @test obs.sources.reported_cases isa String
    @test obs.sources.confirmed_cases isa String
    @test obs.sources.cumulative_tests_analysed isa String
    @test obs.sources.daily_outbound_travellers isa String
    @test obs.sources.daily_outbound_travellers_sd isa String
    @test obs.sources.source_population isa String
    @test obs.sources.genetic_tmrca isa String

    @test !isempty(obs.sources.exported_cases)
    @test !isempty(obs.sources.genetic_tmrca)

    ## death_history: per-vintage cumulative suspected deaths.
    dh = obs.death_history
    @test dh isa NamedTuple
    @test hasproperty(dh, :dates)
    @test hasproperty(dh, :offsets)
    @test hasproperty(dh, :values)
    @test dh.values isa AbstractVector{<:Integer}
    @test dh.offsets isa AbstractVector{<:Integer}
    ## 18-26 May vintages: nine entries. The 23 May deaths use the
    ## SitRep 009 zone-row sum (220), not the erroneous 119 headline.
    ## The 26 May value uses the revised re-issue SitRep 012_v2 (246),
    ## correcting the original SitRep 012 headline of 238.
    @test dh.values == [131, 148, 160, 175, 204, 220, 223, 238, 246]
    @test length(dh.offsets) == 9
    ## Offsets are days before cut-off, sorted ascending (oldest first,
    ## largest offset first), so edges = T - offset are ascending.
    @test issorted(dh.offsets; rev = true)
    @test obs.sources.death_history isa String
    @test !isempty(obs.sources.death_history)
end

@testitem "export_deaths_daily is a daily series to the cut-off" begin
    using BVDOutbreakSize: load_observations

    function _write_obs(io; as_of, death_dates = nothing)
        write(io, "as_of_date = \"$as_of\"\n")
        death_dates === nothing || begin
            quoted = join(("\"$d\"" for d in death_dates), ", ")
            write(io, "[export_death_dates]\nvalue = [$quoted]\n",
                "source = \"x\"\n")
        end
        for k in ("exported_cases", "exports_deaths", "total_deaths",
            "reported_cases", "daily_outbound_travellers",
            "daily_outbound_travellers_sd", "source_population")
            write(io, "[$k]\nvalue = 1\nsource = \"x\"\n")
        end
    end

    mktempdir() do dir
        path = joinpath(dir, "obs.toml")
        open(
            io -> _write_obs(io; as_of = "2026-05-18",
                death_dates = ["2026-05-04", "2026-05-14"]),
            path, "w")
        daily = load_observations(path).export_deaths_daily
        ## Offsets 14 (2026-05-04) and 4 (2026-05-14); earliest = 14, so
        ## the series spans offsets 14..0 (length 15), with one death at
        ## index 1 (offset 14) and one at index 11 (offset 4).
        @test length(daily) == 15
        @test sum(daily) == 2
        @test daily[1] == 1
        @test daily[11] == 1
    end
end

@testitem "export_deaths_daily is empty when no dates are present" begin
    using BVDOutbreakSize: load_observations

    function _write_obs(io; as_of, death_dates = nothing)
        write(io, "as_of_date = \"$as_of\"\n")
        death_dates === nothing || begin
            quoted = join(("\"$d\"" for d in death_dates), ", ")
            write(io, "[export_death_dates]\nvalue = [$quoted]\n",
                "source = \"x\"\n")
        end
        for k in ("exported_cases", "exports_deaths", "total_deaths",
            "reported_cases", "daily_outbound_travellers",
            "daily_outbound_travellers_sd", "source_population")
            write(io, "[$k]\nvalue = 1\nsource = \"x\"\n")
        end
    end

    mktempdir() do dir
        path = joinpath(dir, "obs.toml")
        open(io -> _write_obs(io; as_of = "2026-05-18"), path, "w")
        @test load_observations(path).export_deaths_daily == Int[]
    end
end
