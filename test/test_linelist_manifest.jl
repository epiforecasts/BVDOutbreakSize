@testitem "linelist manifest substitutes the case streams" begin
    using CSV: CSV
    using DataFrames: DataFrame
    using Dates: Date
    using TOML: TOML

    root = joinpath(@__DIR__, "..")
    include(joinpath(root, "scripts", "linelist", "manifest.jl"))

    fixture = joinpath(@__DIR__, "fixtures", "linelist")
    released = joinpath(root, "data", "observations.toml")
    out = mktempdir()

    manifest = write_manifest(
        released = released,
        streams = joinpath(fixture, "linelist_streams.csv"),
        out = joinpath(out, "linelist_observations.toml"))

    raw = TOML.parsefile(manifest)
    streams = CSV.read(joinpath(fixture, "linelist_streams.csv"), DataFrame)

    ## Every replaceable block now holds the fixture's own series, and nothing
    ## of the released one survives in it.
    for block in LINELIST_BLOCKS
        rows = sort(streams[streams.stream .== block, :], :date)
        @test raw[block]["values"] == rows.value
        @test raw[block]["dates"] == [string(d) for d in rows.date]
    end

    ## The cut-off scalar is rewritten from the replacement history rather than
    ## left at the released value, which the model would otherwise condition on
    ## while fitting a history that disagrees with it.
    reported = sort(streams[streams.stream .== "reported_case_history", :],
        :date)
    @test raw["reported_cases"]["value"] == last(reported.value)

    ## The situation reports' retrospective harmonisations are not events in the
    ## line list, and `load_observations` rejects them alongside a replaced
    ## confirmed history.
    @test !haskey(raw, "confirmed_break_dates")

    ## The cut-off moves to the last day the replacement streams cover.
    @test raw["as_of_date"] == string(maximum(streams.date))

    ## Streams with no line-list source are carried over untouched.
    original = TOML.parsefile(released)
    @test raw["recovered_history"] == original["recovered_history"]
    @test raw["occupancy_break_dates"] == original["occupancy_break_dates"]
end

@testitem "linelist manifest and triangle load together" begin
    using BVDOutbreakSize: load_observations
    using Dates: Date

    root = joinpath(@__DIR__, "..")
    include(joinpath(root, "scripts", "linelist", "manifest.jl"))

    fixture = joinpath(@__DIR__, "fixtures", "linelist")
    out = mktempdir()

    manifest = write_manifest(
        released = joinpath(root, "data", "observations.toml"),
        streams = joinpath(fixture, "linelist_streams.csv"),
        out = joinpath(out, "linelist_observations.toml"))

    ## `load_observations` reads the triangle from a fixed filename beside the
    ## manifest and degrades to an empty stream when it is absent, so placing it
    ## is what stops a fit silently dropping the onset likelihood.
    place_onset_curve(fixture, out)

    obs = load_observations(manifest)
    @test obs.cutoff == Date(2026, 7, 31)
    @test !isempty(obs.onset_curve_history.increments)
end

@testitem "linelist output directory refuses a tracked path" begin
    root = joinpath(@__DIR__, "..")
    include(joinpath(root, "scripts", "linelist", "manifest.jl"))

    ## This repository is public and these outputs derive from a line list, so
    ## the only writable location inside it is the ignored one.
    withenv("LINELIST_OUT_DIR" => joinpath(root, "data")) do
        @test_throws ErrorException linelist_output_dir(root)
    end
    withenv("LINELIST_OUT_DIR" => joinpath(root, "ignore", "elsewhere")) do
        @test isdir(linelist_output_dir(root))
    end
    withenv("LINELIST_OUT_DIR" => mktempdir()) do
        @test isdir(linelist_output_dir(root))
    end
end

@testitem "linelist input directory is required" begin
    root = joinpath(@__DIR__, "..")
    include(joinpath(root, "scripts", "linelist", "manifest.jl"))

    withenv("LINELIST_INPUT_DIR" => nothing) do
        @test_throws ErrorException linelist_input_dir()
    end
    withenv("LINELIST_INPUT_DIR" => joinpath(@__DIR__, "fixtures", "linelist")) do
        @test isdir(linelist_input_dir())
    end
end
