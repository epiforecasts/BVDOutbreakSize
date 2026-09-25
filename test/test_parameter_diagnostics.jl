## Tests for the per-parameter fit diagnostics. The table and figure
## builders all accept a diagnostics frame in place of a chain, so most of
## them are exercised against a hand-built frame with no sampling at all.
## The few that need real sampler statistics fit a tiny model with NUTS.

@testsnippet DiagnosticsFrame begin
    using DataFrames: DataFrame

    ## Two scalar parameters and one four-element walk. `slow` and its
    ## deterministic copy `slow_alias` carry identical diagnostics, standing
    ## in for a model that names the same quantity twice.
    function _diag_frame()
        return DataFrame(
            parameter = [
                "fast", "slow", "slow_alias",
                "walk", "walk", "walk", "walk",
            ],
            index = [0, 0, 0, 1, 2, 3, 4],
            rhat = [1.001, 1.3, 1.3, 1.002, 1.05, 1.2, 1.4],
            ess_bulk = [900.0, 20.0, 20.0, 800.0, 300.0, 60.0, 15.0],
            ess_tail = [850.0, 25.0, 25.0, 700.0, 280.0, 70.0, 18.0]
        )
    end
end

@testitem "_split_index separates a parameter from its element index" begin
    using BVDOutbreakSize: _split_index
    @test _split_index("sigma") == ("sigma", 0)
    @test _split_index("rt_state.z[13]") == ("rt_state.z", 13)
    ## A multi-index name keeps the first index, the axis its elements
    ## run along.
    @test _split_index("beta[2, 3]") == ("beta", 2)
    ## A bracket that is not an index is left alone.
    @test _split_index("odd[a]") == ("odd[a]", 0)
end

@testitem "worst_parameters_table ranks by bulk effective sample size" setup = [
    DiagnosticsFrame,
] begin
    using DataFrames: nrow
    using BVDOutbreakSize: worst_parameters_table

    tbl = worst_parameters_table(_diag_frame(); n = 3)
    @test nrow(tbl) == 3
    ## Worst first, and the alias of the second-worst is collapsed away.
    @test tbl.parameter == ["walk[4]", "slow", "walk[3]"]
    @test tbl.ess_bulk == [15.0, 20.0, 60.0]

    kept = worst_parameters_table(
        _diag_frame(); n = 3,
        collapse_aliases = false
    )
    @test kept.parameter == ["walk[4]", "slow", "slow_alias"]

    labelled = worst_parameters_table(
        _diag_frame(); n = 1,
        labels = Dict(:walk => "the walk")
    )
    @test labelled.parameter == ["the walk[4]"]
end

@testitem "family_diagnostics_table groups a walk into one row" setup = [
    DiagnosticsFrame,
] begin
    using DataFrames: nrow
    using BVDOutbreakSize: family_diagnostics_table

    tbl = family_diagnostics_table(_diag_frame(); n = 10)
    ## The four walk elements collapse to one row, and the alias is gone.
    @test nrow(tbl) == 3
    @test tbl.parameter == ["walk", "slow", "fast"]
    @test tbl.elements == [4, 1, 1]
    @test tbl.max_rhat == [1.4, 1.3, 1.001]
    @test tbl.min_ess_bulk == [15.0, 20.0, 900.0]
    ## Two walk elements sit above an R-hat of 1.1, as does `slow`.
    @test tbl[!, "above_1.1"] == [2, 1, 0]
end

@testitem "diagnostic_spread_table counts bad parameters per fit" setup = [
    DiagnosticsFrame,
] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: diagnostic_spread_table

    good = DataFrame(
        parameter = ["a", "b"], index = [0, 0],
        rhat = [1.001, 1.002], ess_bulk = [900.0, 800.0],
        ess_tail = [850.0, 750.0]
    )
    tbl = diagnostic_spread_table("bad" => _diag_frame(), "good" => good)
    @test nrow(tbl) == 2
    @test tbl.parameters == [7, 2]
    @test tbl[!, "rhat_above_1.1"] == [4, 0]
    @test tbl[!, "ess_bulk_below_100"] == [4, 0]
    @test tbl.lowest_ess_parameter == ["walk", "b"]
end

@testitem "diagnostic_contrast pairs shared parameters" setup = [
    DiagnosticsFrame,
] begin
    using DataFrames: DataFrame, nrow
    using BVDOutbreakSize: diagnostic_contrast, diagnostic_contrast_table

    ## The comparison fit mixes the walk well and does not carry `fast`,
    ## so only the parameters both fits hold are paired.
    alone = DataFrame(
        parameter = ["slow", "walk", "walk", "walk", "walk"],
        index = [0, 1, 2, 3, 4],
        rhat = fill(1.001, 5),
        ess_bulk = [400.0, 900.0, 900.0, 900.0, 900.0],
        ess_tail = fill(900.0, 5)
    )
    df = diagnostic_contrast("joint" => _diag_frame(), "alone" => alone)
    @test nrow(df) == 5
    @test all(df.fit .== "alone")
    ## The reference fit's count over the comparison fit's.
    @test df.ess_ratio[df.parameter .== "slow"] == [20.0 / 400.0]

    tbl = diagnostic_contrast_table(df; n = 2)
    @test tbl.parameter == ["walk[4]", "slow"]
    @test tbl.ess_ratio == [round(15.0 / 900.0; digits = 2), 0.05]
end

@testitem "diagnostic figures build from a diagnostics frame" setup = [
    DiagnosticsFrame, HeadlessMakie,
] begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: plot_rhat_spread, diagnostic_contrast,
        plot_parameter_index_diagnostics,
        plot_diagnostic_contrast

    df = _diag_frame()
    @test plot_rhat_spread("joint" => df) isa CairoMakie.Makie.Figure
    ## The walk is the only vector-valued parameter and has four elements,
    ## so it is drawn only once the minimum element count allows it.
    fig = plot_parameter_index_diagnostics(df; min_elements = 4)
    @test fig isa CairoMakie.Makie.Figure
    ## With no group long enough the figure still builds, carrying the
    ## placeholder message instead of panels.
    @test plot_parameter_index_diagnostics(df; min_elements = 50) isa
        CairoMakie.Makie.Figure

    alone = DataFrame(
        parameter = ["slow"], index = [0], rhat = [1.001],
        ess_bulk = [400.0], ess_tail = [900.0]
    )
    contrast = diagnostic_contrast("joint" => df, "alone" => alone)
    @test plot_diagnostic_contrast(contrast) isa CairoMakie.Makie.Figure
    @test plot_diagnostic_contrast(contrast[1:0, :]) isa
        CairoMakie.Makie.Figure
end

@testitem "parameter diagnostics read a fitted chain" tags = [:slow] setup = [
    HeadlessMakie,
] begin
    using DataFrames: DataFrame, nrow
    using Distributions: Normal, product_distribution
    using Turing: @model
    using BVDOutbreakSize: nuts_sample, parameter_diagnostics,
        sampler_by_chain_table,
        divergence_location_table,
        plot_divergence_locations

    ## kept: a trivial Gaussian with one scalar and one three-element
    ## vector gives NUTS a fast target that still exercises the indexed
    ## names a vector-valued parameter produces.
    @model function _param_diag_synthetic()
        x ~ Normal(0, 1)
        z ~ product_distribution([Normal(0, 1) for _ in 1:3])
    end

    chn = nuts_sample(_param_diag_synthetic(); samples = 200, chains = 2)

    df = parameter_diagnostics(chn)
    @test df isa DataFrame
    @test nrow(df) == 4
    @test sort(unique(df.parameter)) == ["x", "z"]
    @test sort(df.index[df.parameter .== "z"]) == [1, 2, 3]
    @test df.index[df.parameter .== "x"] == [0]
    @test all(isfinite, df.rhat)
    @test all(df.ess_bulk .> 0)

    chains = sampler_by_chain_table(chn)
    @test nrow(chains) == 2
    @test chains.chain == [1, 2]
    @test all(chains.draws .== 200)
    @test all(chains.divergences .>= 0)

    ## A well-behaved target usually has no divergence at all, so the
    ## location table is only required to carry the right columns.
    loc = divergence_location_table(chn)
    @test sort(string.(propertynames(loc))) ==
        sort(["parameter", "all_draws", "divergent_draws", "separation"])
    @test plot_divergence_locations(chn, [:x]) isa CairoMakie.Makie.Figure
end
