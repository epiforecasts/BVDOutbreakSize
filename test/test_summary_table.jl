## Tests for summary_table: builds a DataFrame with one row per
## parameter and the documented quantile columns. We sample from
## Prior() on a trivial Turing model so the test does not depend on
## NUTS warm-up.

@testitem "summary_table returns expected columns and rows" tags=[:slow] begin
    using DataFrames: DataFrame, nrow
    using Distributions: Normal
    using Turing: @model, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: summary_table

    ## kept: summary_table only needs two named parameters with sensible
    ## quantiles; the real models drag in BVD-specific structure that the
    ## test does not need.
    @model function _summary_model()
        a ~ Normal(0.0, 1.0)
        b ~ Normal(2.0, 0.5)
    end

    chn = sample(_summary_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    params = [:a, :b]
    df = summary_table(chn, params)

    @test df isa DataFrame
    @test names(df) ==
          ["Quantity", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%"]
    @test nrow(df) == length(params)
    @test df[!, "Quantity"] == ["a", "b"]

    # Each row's quantile columns are monotone (an internal sanity
    # check that the rounded entries still respect the ordering).
    for r in eachrow(df)
        @test r["Lower 90%"] <= r["Lower 60%"] <= r["Lower 30%"]
        @test r["Upper 30%"] <= r["Upper 60%"] <= r["Upper 90%"]
        @test r["Lower 30%"] <= r["Upper 30%"]
    end
end

## `doubling_time = log(2) / r` is unbounded at zero growth, so quantiling
## its own draws returns endpoints that bound nothing once the posterior for
## `r` spans zero. The row is built from `r`'s interval instead.

@testitem "summary_table maps doubling_time through r's interval" tags=[
    :slow
] begin
    using DataFrames: DataFrame
    using Distributions: Normal
    using Statistics: quantile
    using Turing: @model, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: summary_table, posterior_summary

    ## An `r` centred on zero is the case the direct quantiles get wrong.
    @model function _growth_model()
        r ~ Normal(0.0, 0.02)
        doubling_time := log(2) / r
    end

    chn = sample(_growth_model(), Prior(), 2000;
        chain_type = FlexiChains.VNChain, progress = false)
    df = summary_table(chn, [:r, :doubling_time])
    row = df[2, :]
    @test row["Quantity"] == "doubling_time"

    ## Every endpoint is the image of the matching quantile of `r`.
    r_draws = vec(Array(chn[:r]))
    for (col, p) in [("Lower 90%", 0.05), ("Lower 60%", 0.20),
        ("Lower 30%", 0.35), ("Upper 30%", 0.65),
        ("Upper 60%", 0.80), ("Upper 90%", 0.95)]
        @test row[col] == round(log(2) / quantile(r_draws, p); digits = 2)
    end

    ## Ordering is induced by `r`, so the row runs from the fastest decline
    ## through the zero-growth pole to the fastest growth and is not sorted.
    @test row["Lower 90%"] < 0 < row["Upper 90%"]
    @test abs(row["Lower 30%"]) > abs(row["Lower 90%"])
    @test row["Upper 30%"] > row["Upper 90%"]

    ## Quantiling doubling_time's own draws puts a large negative number at
    ## the low end when r spans zero, so the mapped row must differ from it.
    direct = posterior_summary(vec(Array(chn[:doubling_time])))
    @test round(direct.lo90; digits = 2) != row["Lower 90%"]
end

@testitem "summary_table falls back when the source is absent" tags=[
    :slow
] begin
    using Distributions: Normal
    using Turing: @model, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: summary_table

    ## A chain carrying the derived quantity but not `r` still summarises.
    @model function _no_r_model()
        doubling_time ~ Normal(30.0, 5.0)
    end

    chn = sample(_no_r_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    df = summary_table(chn, [:doubling_time])
    @test df[1, "Lower 90%"] <= df[1, "Upper 90%"]
end

## The review flagged that only the zero-spanning case was covered. With a
## one-signed growth rate the transform is strictly decreasing, so the row is
## ordered by growth rate and runs from the longest doubling time to the
## shortest. That ordering is the point of the change, not an accident, so
## pin it.

@testitem "summary_table orders doubling_time by r when r is one-signed" tags=[
    :slow
] begin
    using Statistics: quantile
    using Distributions: Normal, truncated
    using Turing: @model, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: summary_table

    @model function _positive_growth()
        r ~ truncated(Normal(0.05, 0.01); lower = 0.005)
        doubling_time := log(2) / r
    end

    chn = sample(_positive_growth(), Prior(), 2000;
        chain_type = FlexiChains.VNChain, progress = false)
    df = summary_table(chn, [:r, :doubling_time])
    row = df[2, :]

    r_draws = vec(Array(chn[:r]))
    for (col, p) in [("Lower 90%", 0.05), ("Upper 90%", 0.95)]
        @test row[col] == round(log(2) / quantile(r_draws, p); digits = 2)
    end

    ## Every value is a positive doubling time, and the slowest growth sits
    ## in the Lower column because the row is ordered by growth rate.
    @test row["Lower 90%"] > row["Upper 90%"] > 0
    @test row["Lower 30%"] > row["Upper 30%"] > 0
end
