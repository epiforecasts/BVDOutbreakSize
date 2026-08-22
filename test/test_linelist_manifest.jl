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

@testitem "linelist baseline manifest only moves the cut-off" begin
    using TOML: TOML

    root = joinpath(@__DIR__, "..")
    include(joinpath(root, "scripts", "linelist", "manifest.jl"))

    released = joinpath(root, "data", "observations.toml")
    out = mktempdir()

    manifest = write_baseline_manifest(released = released,
        as_of = "2026-07-31", out = joinpath(out, "baseline.toml"))

    raw = TOML.parsefile(manifest)
    original = TOML.parsefile(released)

    ## The comparator differs from the released manifest in its cut-off and
    ## nothing else, so a fit of it against a line-list fit at the same cut-off
    ## differs in the data alone.
    @test raw["as_of_date"] == "2026-07-31"
    for k in keys(original)
        k == "as_of_date" && continue
        @test raw[k] == original[k]
    end

    ## Unlike a substituted manifest, this one is still the situation-report
    ## series, so the harmonisation breaks that series really has are kept.
    @test haskey(raw, "confirmed_break_dates")

    ## A cut-off past where the released data stop would leave the model
    ## fitting a tail of nothing.
    @test_throws ErrorException write_baseline_manifest(released = released,
        as_of = "2099-01-01", out = joinpath(out, "late.toml"))
end

@testitem "linelist onset triangle source is named, not inferred" begin
    using BVDOutbreakSize: load_observations

    root = joinpath(@__DIR__, "..")
    include(joinpath(root, "scripts", "linelist", "manifest.jl"))

    fixture = joinpath(@__DIR__, "fixtures", "linelist")

    ## Both triangles are called `onset_curve_scanned.csv` and they are
    ## different constructions of the same quantity: one digitised from the
    ## situation-report figures, one built from line-list records. A fit handed
    ## the wrong one is a mixture that nothing downstream reveals, so the source
    ## is named rather than taken from whichever file is to hand.
    ll = place_onset_curve(fixture, mktempdir(); source = :linelist,
        root = root)
    @test ll.source === :linelist
    @test ll.src == joinpath(fixture, "onset_curve_scanned.csv")
    @test ll.vintages > 0

    sr = place_onset_curve(fixture, mktempdir(); source = :sitrep, root = root)
    @test sr.source === :sitrep
    @test sr.src == joinpath(root, "data", "onset_curve_scanned.csv")
    @test sr.vintages > 0

    ## The two are not the same file, which is the whole point of naming them.
    @test ll.vintages != sr.vintages || ll.cells != sr.cells

    @test_throws ErrorException place_onset_curve(fixture, mktempdir();
        source = :nonsense, root = root)

    ## A triangle placed where the manifest can see it is the difference
    ## between fitting the onset stream and silently dropping it.
    out = mktempdir()
    manifest = write_manifest(
        released = joinpath(root, "data", "observations.toml"),
        streams = joinpath(fixture, "linelist_streams.csv"),
        out = joinpath(out, "linelist_observations.toml"))
    place_onset_curve(fixture, out; source = :linelist, root = root)
    @test !isempty(load_observations(manifest).onset_curve_history.increments)
end

@testitem "linelist manifest pins one cut-off across constructions" begin
    using CSV: CSV
    using DataFrames: DataFrame
    using Dates: Date
    using TOML: TOML

    root = joinpath(@__DIR__, "..")
    include(joinpath(root, "scripts", "linelist", "manifest.jl"))

    fixture = joinpath(@__DIR__, "fixtures", "linelist")
    out = mktempdir()

    ## The two stream constructions end on different days, so left to their own
    ## last day they sit on different grids and their trajectories are not on
    ## one axis.
    manifest = write_manifest(
        released = joinpath(root, "data", "observations.toml"),
        streams = joinpath(fixture, "linelist_streams.csv"),
        as_of = "2026-07-20",
        out = joinpath(out, "pinned.toml"))
    raw = TOML.parsefile(manifest)
    @test raw["as_of_date"] == "2026-07-20"

    ## The cut-off scalar is read at the cut-off rather than at the end of the
    ## replacement series. The known-by construction runs past a pinned
    ## cut-off, so taking its final value would tell the model a total that
    ## includes cases the history it fits has been truncated before, and it can
    ## only reconcile the two by distorting ascertainment.
    streams = CSV.read(joinpath(fixture, "linelist_streams.csv"), DataFrame)
    reported = sort(
        streams[(streams.stream .== "reported_case_history") .& (streams.date .<= Date(2026, 7, 20)), :], :date)
    @test raw["reported_cases"]["value"] == last(reported.value)
    @test raw["reported_cases"]["value"] != maximum(reported.value) ||
          last(reported.date) == Date(2026, 7, 20)
    ## and it is not the series end, which is what the pin is guarding against
    all_reported = sort(
        streams[streams.stream .== "reported_case_history", :], :date)
    @test last(all_reported.date) > Date(2026, 7, 20)
    @test raw["reported_cases"]["value"] != last(all_reported.value)

    @test_throws ErrorException write_manifest(
        released = joinpath(root, "data", "observations.toml"),
        streams = joinpath(fixture, "linelist_streams.csv"),
        as_of = "2099-01-01",
        out = joinpath(out, "late.toml"))
end

@testitem "linelist delay configs carry the estimates into priors" begin
    using Distributions: Normal, Truncated

    root = joinpath(@__DIR__, "..")
    include(joinpath(root, "scripts", "linelist", "delays.jl"))

    fixture = joinpath(@__DIR__, "fixtures", "linelist", "delays")

    withenv("BVD_DELAY_DIR" => fixture, "BVD_DELAY_PRIOR_INFLATE" => nothing) do
        ## `repo` is the absence of an override rather than a restatement of
        ## the package defaults, so it cannot drift from them.
        @test isnothing(delay_config("repo"))
        @test_throws ErrorException delay_config("not_a_config")

        d = delay_config("cmmid_any")

        ## The generation interval's shape and scale are taken as published,
        ## and their widths are the serial interval's own reported uncertainty
        ## on the mean carried across by the delta method. The fixture is built
        ## so those come out exact: mean 20 and SD 10 give shape 4 and scale 5,
        ## and a meanlog standard error of 0.05 gives 0.4 and 0.25.
        @test d.gi_alpha_prior.untruncated.μ ≈ 4
        @test d.gi_theta_prior.untruncated.μ ≈ 5
        @test d.gi_alpha_prior.untruncated.σ ≈ 0.4
        @test d.gi_theta_prior.untruncated.σ ≈ 0.25
        @test d.gi_nmax ≥ 5

        ## `censored_delay_model` builds `lognormal_meansd(mean, sd)`, so the
        ## priors are centred on the moment-matched values and the fitted
        ## lognormal is reproduced at the prior centre rather than approximated.
        μ, σ = 2.0, 1.0
        mean = exp(μ + σ^2 / 2)
        sd = mean * sqrt(exp(σ^2) - 1)
        @test d.report_mean_prior.untruncated.μ ≈ mean
        @test d.report_sd_prior.untruncated.μ ≈ sd
        ## With `sigma_Intercept` at zero the SD of the mean is the mean times
        ## the two standard errors in quadrature, which is an independent read
        ## of the delta method rather than a repeat of it.
        @test d.report_mean_prior.untruncated.σ / mean ≈ hypot(0.01, 0.02)
        @test d.report_nmax ≥ 5

        ## The generation-interval-only configs leave the report delay alone,
        ## which is what separates the two effects when both are run.
        gi_only = delay_config("cmmid_gi_any")
        @test isnothing(gi_only.report_mean_prior)
        @test gi_only.gi_alpha_prior.untruncated.μ ≈ 4

        ## Each config names a different pair definition.
        @test delay_config("cmmid_gi_diag").gi_theta_prior.untruncated.μ ≈ 6.25

        ## The onsets fit estimates its reporting delay from the triangle, so a
        ## report-delay override there would be silently dropped and the run
        ## would look like a result. It is refused instead.
        @test_throws ErrorException check_delay_config(d, "onsets")
        @test isnothing(check_delay_config(gi_only, "onsets"))
        @test isnothing(check_delay_config(d, "cases"))
        @test isnothing(check_delay_config(nothing, "onsets"))

        ## `cmmid_rep` is the mirror of `cmmid_gi_*`: the report delay on its
        ## own, no generation interval. It exists so a data-source comparison
        ## can take cmmid's report delay without also taking a generation
        ## interval fitted from line-list transmission pairs.
        rep = delay_config("cmmid_rep")
        @test isnothing(rep.gi)
        @test isnothing(rep.gi_alpha_prior)
        @test isnothing(rep.gi_theta_prior)
        @test isnothing(rep.gi_nmax)
        @test rep.has_report
        ## The same report delay `cmmid_any` carries, since there is one
        ## report prior rather than one per transmission-pair definition.
        @test rep.report_mean_prior.untruncated.μ ≈ d.report_mean_prior.untruncated.μ
        @test rep.report_sd_prior.untruncated.μ ≈ d.report_sd_prior.untruncated.μ

        ## It overrides the report delay, so the onsets fit refuses it for the
        ## same reason it refuses `cmmid_*`.
        @test_throws ErrorException check_delay_config(rep, "onsets")
        @test isnothing(check_delay_config(rep, "confirmed"))

        ## Provenance has to survive a config with no generation interval:
        ## it reads neither generation-interval table and writes no gi rows,
        ## rather than erroring on a `nothing`.
        provdir = mktempdir()
        provpath = write_delay_provenance(rep, provdir)
        @test isfile(provpath)
        prov = CSV.read(provpath, DataFrame)
        @test !any(startswith.(string.(prov.item), "gi_"))
        files = prov[prov.item .== "file", :value]
        @test "bayes_pooled_parameters.csv" in files
        @test !("gi_estimates.csv" in files)
        @test any(startswith.(string.(prov.item), "report_"))
    end

    ## The published standard errors come from large samples, so the priors are
    ## near-fixed; the inflation factor is how a run asks what the answer owes
    ## to the delay being exactly this.
    withenv("BVD_DELAY_DIR" => fixture, "BVD_DELAY_PRIOR_INFLATE" => "4") do
        @test delay_config("cmmid_any").gi_alpha_prior.untruncated.σ ≈ 1.6
    end

    withenv("BVD_DELAY_DIR" => nothing) do
        @test_throws ErrorException delay_config("cmmid_any")
    end
end
