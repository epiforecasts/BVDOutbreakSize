## Tests for the convergence gate (`docs/fits/convergence.jl`): which
## diagnostics fail a build, which only warn, that the thresholds are settable
## from the environment, and that the report names the fit and its reasons.

@testitem "convergence_verdict separates failure from poor mixing" begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "convergence.jl"))

    clean = (
        max_rhat = 1.01, min_ess_bulk = 800.0, min_ess_tail = 700.0,
        n_divergent = 0, n_draws = 3200,
    )
    @test convergence_verdict(clean).status === :pass
    @test isempty(convergence_verdict(clean).failures)

    ## Where the headline joint fit currently sits: converged, mixing badly.
    borderline = (
        max_rhat = 1.057, min_ess_bulk = 48.0, min_ess_tail = 120.0,
        n_divergent = 6, n_draws = 3200,
    )
    v = convergence_verdict(borderline)
    @test v.status === :warn
    @test isempty(v.failures)
    @test any(occursin("max R-hat", m) for m in v.warnings)
    @test any(occursin("min bulk ESS", m) for m in v.warnings)

    ## A chain stuck at its initialisation, as the v2.0.0 joint fit was.
    stuck = (
        max_rhat = 2.6, min_ess_bulk = 8.0, min_ess_tail = 9.0,
        n_divergent = 252, n_draws = 3200,
    )
    f = convergence_verdict(stuck)
    @test f.status === :fail
    @test any(occursin("max R-hat is 2.6", m) for m in f.failures)

    ## Divergences are judged as a share of the draws, not as a count: the
    ## same 252 out of far fewer draws is a worse fit, and out of far more is
    ## a better one.
    @test convergence_verdict(
        (;
            stuck..., max_rhat = 1.01, min_ess_bulk = 800.0,
            min_ess_tail = 700.0, n_draws = 200_000,
        )
    ).status === :pass

    ## A diagnostic that could not be computed is not evidence either way.
    @test convergence_verdict(
        (; clean..., max_rhat = NaN, min_ess_bulk = NaN, min_ess_tail = NaN)
    ).status === :pass
end

@testitem "convergence thresholds come from the environment" begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "convergence.jl"))

    borderline = (
        max_rhat = 1.057, min_ess_bulk = 48.0, min_ess_tail = 120.0,
        n_divergent = 6, n_draws = 3200,
    )
    @test convergence_verdict(borderline).status === :warn
    withenv("BVD_CONVERGENCE_FAIL_RHAT" => "1.01") do
        @test convergence_verdict(
            borderline; thresholds = convergence_thresholds()
        ).status === :fail
    end
    withenv("BVD_CONVERGENCE_WARN_ESS_BULK" => "10") do
        v = convergence_verdict(
            borderline; thresholds = convergence_thresholds()
        )
        @test !any(occursin("bulk ESS", m) for m in v.warnings)
    end
end

@testitem "convergence_markdown reports the verdict and its reasons" begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "convergence.jl"))

    stuck = (
        max_rhat = 2.6, min_ess_bulk = 8.0, min_ess_tail = 9.0,
        n_divergent = 252, n_draws = 3200,
    )
    check = (;
        id = "joint", diagnostics = stuck, per_parameter = nothing,
        chain = nothing, convergence_verdict(stuck)...,
    )
    md = convergence_markdown(
        [check]; marker = "<!-- mark -->", context = "Commit `abc1234`."
    )

    @test startswith(md, "<!-- mark -->")
    @test occursin("## Fit convergence:", md)
    @test occursin("should not be read", md)
    @test occursin("Commit `abc1234`.", md)
    @test occursin("| `joint` | 2.6 | 8 | 9 | 252 | 3200 |", md)
    @test occursin("**fail** — max R-hat is 2.6", md)
    @test occursin("Thresholds", md)

    ## A passing fit says so and carries no reasons.
    clean = (
        max_rhat = 1.01, min_ess_bulk = 800.0, min_ess_tail = 700.0,
        n_divergent = 0, n_draws = 3200,
    )
    passing = (;
        id = "joint", diagnostics = clean, per_parameter = nothing,
        chain = nothing, convergence_verdict(clean)...,
    )
    ok = convergence_markdown([passing])
    @test occursin("cleared the convergence thresholds", ok)
    @test !occursin("**fail**", ok)
end

@testitem "fit_convergence gates a real chain and names its parameters" tags = [
    :slow,
] begin
    using Distributions: Normal
    using Turing: @model
    using BVDOutbreakSize: nuts_sample
    include(joinpath(@__DIR__, "..", "docs", "fits", "convergence.jl"))

    ## kept: a trivial two-parameter Gaussian is the cheapest target that
    ## still carries R-hat, both effective sample sizes and a divergence
    ## flag, which is all the gate reads.
    @model function _gate_synthetic()
        x ~ Normal(0, 1)
        y ~ Normal(0, 1)
    end

    chn = nuts_sample(_gate_synthetic(); samples = 400, chains = 2)
    c = fit_convergence("demo", chn)
    @test c.id == "demo"
    @test c.diagnostics.n_draws == 800
    @test c.status in (:pass, :warn)

    ## Forcing a failure makes the per-parameter breakdown available, and the
    ## report then names the parameters rather than only the headline number.
    bad = withenv("BVD_CONVERGENCE_FAIL_ESS_BULK" => "1e9") do
        fit_convergence("demo", chn; thresholds = convergence_thresholds())
    end
    @test bad.status === :fail
    @test bad.per_parameter !== nothing
    md = convergence_markdown([bad])
    @test occursin("Worst-mixing parameters of `demo`", md)
    @test occursin("ess_bulk", md)
end

@testitem "the report uses the thresholds it was handed throughout" begin
    using DataFrames: DataFrame
    include(joinpath(@__DIR__, "..", "docs", "fits", "convergence.jl"))

    ## The grouped table counts elements past an R-hat threshold, and that
    ## threshold has to be the one the rest of the report is quoting. A
    ## report whose footer and whose tables disagree about what counts as a
    ## bad R-hat says two things at once.
    per_parameter = DataFrame(
        parameter = ["a", "b"], index = [0, 0],
        rhat = [1.3, 1.01], ess_bulk = [12.0, 900.0],
        ess_tail = [15.0, 800.0]
    )
    d = (
        max_rhat = 1.3, min_ess_bulk = 12.0, min_ess_tail = 15.0,
        n_divergent = 0, n_draws = 3200,
    )
    thresholds = (
        fail = (
            rhat = 1.25, ess_bulk = 25.0, ess_tail = 25.0,
            divergent_fraction = 0.05,
        ),
        warn = (
            rhat = 1.2, ess_bulk = 100.0, ess_tail = 100.0,
            divergent_fraction = 0.01,
        ),
    )
    check = (;
        id = "joint", diagnostics = d, per_parameter, chain = nothing,
        convergence_verdict(d; thresholds = thresholds)...,
    )
    md = convergence_markdown([check]; thresholds = thresholds)
    @test occursin("above_1.2", md)
    @test !occursin("above_1.05", md)
    @test occursin("past the failure threshold 1.25", md)
end

@testitem "convergence_summary gives the verdict in plain prose" begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "convergence.jl"))

    clean = (
        max_rhat = 1.01, min_ess_bulk = 800.0, min_ess_tail = 700.0,
        n_divergent = 0, n_draws = 3200,
    )
    ok = convergence_summary("joint", clean)
    @test startswith(ok, "The joint fit passes the convergence checks. ")
    @test occursin("worst R-hat is 1.01", ok)
    @test occursin("sample size is 800,", ok)
    @test occursin("0 divergent transitions in 3200 draws", ok)
    @test occursin("at most 1.05,", ok)
    @test occursin("at least 100 and at most 1% of its draws", ok)

    borderline = (; clean..., max_rhat = 1.057, min_ess_bulk = 48.0)
    @test occursin(
        "passes the convergence checks with warnings.",
        convergence_summary("joint", borderline)
    )

    stuck = (; clean..., max_rhat = 2.6, n_divergent = 252)
    @test occursin(
        "The joint fit fails the convergence checks.",
        convergence_summary("joint", stuck)
    )

    ## An undefined diagnostic is printed as such rather than as NaN.
    @test occursin(
        "worst R-hat is n/a",
        convergence_summary("joint", (; clean..., max_rhat = NaN))
    )

    ## The thresholds quoted are the ones handed in.
    thresholds = (
        fail = convergence_thresholds().fail,
        warn = (
            rhat = 1.2, ess_bulk = 40.0, ess_tail = 40.0,
            divergent_fraction = 0.025,
        ),
    )
    t = convergence_summary("joint", borderline; thresholds = thresholds)
    @test occursin("passes the convergence checks.", t)
    @test occursin("at most 1.2,", t)
    @test occursin("at least 40 and at most 2.5% of its draws", t)
end

@testitem "a malformed threshold variable names itself" begin
    include(joinpath(@__DIR__, "..", "docs", "fits", "convergence.jl"))

    ## The gate dies before writing its report, so the parse failure is all
    ## anyone sees. It has to say which variable was wrong.
    withenv("BVD_CONVERGENCE_FAIL_RHAT" => "1,1") do
        e = try
            convergence_thresholds()
            nothing
        catch err
            err
        end
        @test e isa ErrorException
        @test occursin("BVD_CONVERGENCE_FAIL_RHAT", e.msg)
    end

    ## `nan` and `inf` parse as numbers and then turn the threshold off in
    ## silence, since nothing compares greater than either. A gate a typo can
    ## disable without saying so is the failure this gate exists to prevent.
    for bad in ("nan", "inf", "-inf")
        withenv("BVD_CONVERGENCE_FAIL_ESS_BULK" => bad) do
            e = try
                convergence_thresholds()
                nothing
            catch err
                err
            end
            @test e isa ErrorException
            @test occursin("BVD_CONVERGENCE_FAIL_ESS_BULK", e.msg)
            @test occursin("finite", e.msg)
        end
    end
end
