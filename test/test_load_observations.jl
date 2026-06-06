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

    ## confirmed-case history runs to the 3 June cut-off; its final
    ## vintage equals the cut-off `confirmed_cases` total.
    @test obs.confirmed_case_history.values ==
          [33, 51, 57, 79, 83, 101, 105, 106, 121, 125, 210, 263,
        282, 321, 344, 363, 381]
    @test obs.confirmed_case_history.values[end] == obs.confirmed_cases

    ## confirmed deaths: flat at 17 through 28 May, then the late-
    ## confirmation catch-up to 64 on 3 June.
    @test obs.confirmed_deaths isa Integer
    @test obs.confirmed_death_history isa NamedTuple
    @test obs.confirmed_death_history.values ==
          [17, 17, 17, 42, 42, 48, 60, 62, 64]
    @test obs.confirmed_death_history.values[end] == obs.confirmed_deaths
    @test obs.sources.confirmed_death_history isa String

    ## laboratory throughput histories (cumulative national, 23-28 May);
    ## the analysed series ends at the cut-off `cumulative_tests_analysed`.
    @test obs.tests_received_history.values ==
          [418, 431, 431, 662, 774, 883]
    @test obs.tests_analysed_history.values ==
          [211, 295, 295, 403, 648, 755]
    @test obs.tests_analysed_history.values[end] ==
          obs.cumulative_tests_analysed
    @test obs.sources.tests_received_history isa String
    @test obs.sources.tests_analysed_history isa String
end

@testitem "load_observations returns the lab-throughput histories" begin
    using BVDOutbreakSize: load_observations
    obs = load_observations()

    ## tests received / analysed span the six 23-28 May vintages whose
    ## lab Tableau IV is present, aligned with the tail of the confirmed
    ## history so the binomial denominator and numerator share dates.
    for h in (obs.tests_received_history, obs.tests_analysed_history)
        @test h isa NamedTuple
        @test hasproperty(h, :dates)
        @test hasproperty(h, :offsets)
        @test hasproperty(h, :values)
        @test h.values isa AbstractVector{<:Integer}
        @test length(h.values) == 6
        ## Oldest first, so edges = T - offset are ascending.
        @test issorted(h.offsets; rev = true)
    end

    @test obs.tests_received_history.values ==
          [418, 431, 431, 662, 774, 883]
    @test obs.tests_analysed_history.values ==
          [211, 295, 295, 403, 648, 755]
    ## Analysed never exceeds received; both align with the lab cut-off.
    @test all(obs.tests_analysed_history.values .<=
              obs.tests_received_history.values)
    @test obs.tests_analysed_history.values[end] ==
          obs.cumulative_tests_analysed
    @test obs.tests_received_history.dates ==
          obs.tests_analysed_history.dates
    ## Per-vintage confirmed positives never exceed tests analysed (the
    ## binomial denominator), checked on the shared 23-28 May lab dates
    ## (the confirmed history now runs past them to 3 June).
    ch = obs.confirmed_case_history
    aidx = [findfirst(==(d), ch.dates)
            for d in obs.tests_analysed_history.dates]
    @test all(ch.values[aidx] .<= obs.tests_analysed_history.values)

    @test obs.sources.tests_received_history isa String
    @test obs.sources.tests_analysed_history isa String
    @test !isempty(obs.sources.tests_analysed_history)
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

@testitem "exported_cases_daily is a daily series to the cut-off" begin
    using BVDOutbreakSize: load_observations

    function _write_obs(io; as_of, case_dates = nothing)
        write(io, "as_of_date = \"$as_of\"\n")
        case_dates === nothing || begin
            quoted = join(("\"$d\"" for d in case_dates), ", ")
            write(io, "[export_case_dates]\nvalue = [$quoted]\n",
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
            io -> _write_obs(io; as_of = "2026-05-26",
                case_dates = ["2026-05-11", "2026-05-16", "2026-05-23"]),
            path, "w")
        daily = load_observations(path).exported_cases_daily
        ## Offsets 15 (11 May), 10 (16 May), 3 (23 May); earliest = 15, so
        ## the series spans offsets 15..0 (length 16), with detections at
        ## index 1 (offset 15), index 6 (offset 10) and index 13 (offset 3).
        @test length(daily) == 16
        @test sum(daily) == 3
        @test daily[1] == 1
        @test daily[6] == 1
        @test daily[13] == 1
    end
end

@testitem "exported_cases_daily is empty when no dates are present" begin
    using BVDOutbreakSize: load_observations

    function _write_obs(io; as_of)
        write(io, "as_of_date = \"$as_of\"\n")
        for k in ("exported_cases", "exports_deaths", "total_deaths",
            "reported_cases", "daily_outbound_travellers",
            "daily_outbound_travellers_sd", "source_population")
            write(io, "[$k]\nvalue = 1\nsource = \"x\"\n")
        end
    end

    mktempdir() do dir
        path = joinpath(dir, "obs.toml")
        open(io -> _write_obs(io; as_of = "2026-05-26"), path, "w")
        @test load_observations(path).exported_cases_daily == Int[]
    end
end

@testitem "as_of_override truncates to the earlier cut-off" begin
    using BVDOutbreakSize: load_observations
    using Dates: Date, value

    ## The committed data file's own cut-off (3 June) and a one-week
    ## earlier as-of (28 May). The 28 May totals are the known
    ## INSP/WHO values: confirmed cases 210, confirmed deaths 17, samples
    ## analysed 755; the frozen suspected totals stay at their 26 May
    ## values (reported cases 1077, suspected deaths 246) and the dated
    ## exports are 3 cases / 1 death detected by 28 May.
    full = load_observations()
    cut = load_observations(; as_of_override = "2026-05-28")

    @test cut.as_of_date == "2026-05-28"
    @test cut.confirmed_cases == 210
    @test cut.confirmed_deaths == 17
    @test cut.cumulative_tests_analysed == 755
    @test cut.reported_cases == 1077
    @test cut.total_deaths == 246
    @test cut.exported_cases == 3
    @test cut.exports_deaths == 1

    ## Histories are truncated to entries on or before the cut-off.
    @test cut.confirmed_case_history.dates[end] == "2026-05-28"
    @test all(Date.(cut.confirmed_case_history.dates) .<= Date("2026-05-28"))
    @test cut.confirmed_death_history.dates[end] == "2026-05-28"
    @test cut.confirmed_case_history.values[end] == cut.confirmed_cases

    ## Elapsed-time offsets are recomputed relative to the new cut-off, so
    ## the genetic floor is one week closer than under the 3 June file.
    @test cut.genetic_tmrca_days ==
          full.genetic_tmrca_days -
          value(Date("2026-06-03") - Date("2026-05-28"))

    ## A `Date` argument is equivalent to the ISO string.
    @test load_observations(; as_of_override = Date("2026-05-28")).as_of_date ==
          "2026-05-28"

    ## The default (no override) loads the file's own cut-off unchanged.
    @test full.as_of_date == "2026-06-03"
    @test full.confirmed_cases == 381
    @test full.confirmed_deaths == 64
end
