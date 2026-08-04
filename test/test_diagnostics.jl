## Smoke tests for fit_diagnostics and diagnostics_table. A tiny model
## is fitted with NUTS (two chains) so the chain carries the R-hat,
## bulk-ESS and numerical-error (divergence) information the helpers
## read.

@testitem "fit_diagnostics summarises rhat, ess and divergences" tags=[
    :slow
] begin
    using Distributions: Normal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample, fit_diagnostics

    ## kept: a trivial 2-parameter Gaussian gives NUTS the fast, well-
    ## conditioned target needed to exercise rhat/ESS/divergence summaries.
    @model function _diag_synthetic()
        x ~ Normal(0, 1)
        y ~ Normal(0, 1)
    end

    chn = nuts_sample(_diag_synthetic(); samples = 200, chains = 2)
    d = fit_diagnostics(chn)
    @test isfinite(d.max_rhat)
    @test d.max_rhat > 0
    @test d.min_ess_bulk > 0
    @test d.n_divergent >= 0
end

@testitem "diagnostics_table has one row per fit" tags=[:slow] begin
    using DataFrames: DataFrame, nrow
    using Distributions: Normal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample, diagnostics_table

    ## kept: see the first @testitem; the diagnostics-table check just
    ## wants two cheap NUTS fits with no BVD-specific structure.
    @model function _diag_synthetic()
        x ~ Normal(0, 1)
        y ~ Normal(0, 1)
    end

    chn1 = nuts_sample(_diag_synthetic(); samples = 150, chains = 2)
    chn2 = nuts_sample(_diag_synthetic(); samples = 150, chains = 2)
    tbl = diagnostics_table("fit A" => chn1, "fit B" => chn2)
    @test tbl isa DataFrame
    @test nrow(tbl) == 2
    @test sort(string.(propertynames(tbl))) ==
          sort(["fit", "max_rhat", "min_ess_bulk", "divergences"])
    @test all(tbl.divergences .>= 0)
end
