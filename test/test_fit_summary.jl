## Tests for the per-fit job summary (`docs/fits/summary.jl`): the markdown
## carries the convergence diagnostics and the headline quantities the chain
## holds, and `write_fit_summary` appends it to the file GitHub names.

@testitem "fit_summary_markdown reports diagnostics and headlines" tags=[
    :slow
] begin
    using Distributions: Normal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample
    include(joinpath(@__DIR__, "..", "docs", "fits", "summary.jl"))

    ## kept: a trivial two-parameter Gaussian named after the reported
    ## quantities exercises the diagnostics row and both headline rows.
    @model function _summary_synthetic()
        C_T ~ Normal(100, 10)
        R_T ~ Normal(1, 0.1)
    end

    chn = nuts_sample(_summary_synthetic(); samples = 150, chains = 2)
    md = fit_summary_markdown("demo", chn)

    @test occursin("### Fit `demo`", md)
    @test occursin(
        "| max R-hat | min ESS bulk | min ESS tail | divergences |", md)
    @test occursin("| `C_T` |", md)
    @test occursin("| `R_T` |", md)
    ## Scannable in the run UI: a heading, the two tables and their rows.
    @test count("\n", md) <= 12

    ## The frozen fits hand back `(; cutoff, o, chn)` rather than a chain, and
    ## the summary is appended to the job-summary file when one is named.
    path = joinpath(mktempdir(), "step_summary.md")
    withenv("GITHUB_STEP_SUMMARY" => path) do
        write_fit_summary("frozen_demo", (; cutoff = "2026-05-20", o = 1, chn))
    end
    @test occursin("### Fit `frozen_demo`", read(path, String))
end
