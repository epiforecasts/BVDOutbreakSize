## Tests for summary_table: builds a DataFrame with one row per
## parameter and the documented quantile columns. We sample from
## Prior() on a trivial Turing model so the test does not depend on
## NUTS warm-up.

@testitem "summary_table returns expected columns and rows" tags = [:slow] begin
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

    chn = sample(
        _summary_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false
    )
    params = [:a, :b]
    df = summary_table(chn, params)

    @test df isa DataFrame
    @test names(df) ==
        [
        "Quantity", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%",
    ]
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

@testitem "summary_table maps doubling_time through r's interval" tags = [
    :slow,
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

    chn = sample(
        _growth_model(), Prior(), 2000;
        chain_type = FlexiChains.VNChain, progress = false
    )
    df = summary_table(chn, [:r, :doubling_time])
    row = df[2, :]
    @test row["Quantity"] == "doubling_time"

    ## Every endpoint is the image of the matching quantile of `r`.
    r_draws = vec(Array(chn[:r]))
    for (col, p) in [
            ("Lower 90%", 0.05), ("Lower 60%", 0.2),
            ("Lower 30%", 0.35), ("Upper 30%", 0.65),
            ("Upper 60%", 0.8), ("Upper 90%", 0.95),
        ]
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

@testitem "summary_table falls back when the source is absent" tags = [
    :slow,
] begin
    using Distributions: Normal
    using Turing: @model, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: summary_table

    ## A chain carrying the derived quantity but not `r` still summarises.
    @model function _no_r_model()
        doubling_time ~ Normal(30.0, 5.0)
    end

    chn = sample(
        _no_r_model(), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false
    )
    df = summary_table(chn, [:doubling_time])
    @test df[1, "Lower 90%"] <= df[1, "Upper 90%"]
end

## The review flagged that only the zero-spanning case was covered. With a
## one-signed growth rate the transform is strictly decreasing, so the row is
## ordered by growth rate and runs from the longest doubling time to the
## shortest. That ordering is the point of the change, not an accident, so
## pin it.

@testitem "summary_table orders doubling_time by r when r is one-signed" tags = [
    :slow,
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

    chn = sample(
        _positive_growth(), Prior(), 2000;
        chain_type = FlexiChains.VNChain, progress = false
    )
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

@testsnippet PatchHeadlineDraws begin
    using Random: MersenneTwister
    using BVDOutbreakSize: PROVINCE_LABELS

    rng = MersenneTwister(7)
    nd = 400
    np = 3
    ndays = 10
    draws(centre, sd) = [centre .+ sd .* randn(rng, np) for _ in 1:nd]
    base = (;
        C_T_patch = draws([900.0, 300.0, 80.0], 20.0),
        R_T_patch = draws([1.2, 0.9, 0.7], 0.05),
    )
    full = (;
        base...,
        CFR_patch = draws([0.4, 0.3, 0.2], 0.01),
        province_ascertainment = draws([1.4, 0.8, 0.6], 0.05),
        ## One import a day into every patch, flattened patch-fastest.
        importation_patch = [ones(np * ndays) for _ in 1:nd],
    )
end

@testitem "patch_headline compares the provinces" setup = [
    PatchHeadlineDraws,
] begin
    using BVDOutbreakSize: patch_headline

    md = patch_headline(base, np)
    bullets = filter(startswith("- **"), split(md, "\n"))
    ## Shares of infections and the reproduction number only.
    @test length(bullets) == 2
    @test !occursin("median", md)
    ## Every province's share is given, and Ituri, the largest in every
    ## draw, is named as carrying the most.
    for p in 1:np
        @test occursin(Regex("$(PROVINCE_LABELS[p]) \\d+–\\d+%"), md)
    end
    @test occursin("Ituri has the most infections", md)
    @test occursin("over 99%", md)
    ## The reproduction number runs from the lowest to the highest province,
    ## with the probability each is above one read from the draws.
    @test occursin(
        r"from [\d.]+–[\d.]+ in Haut-Uele to [\d.]+–[\d.]+ in Ituri", md
    )
    @test occursin("under 1% in Haut-Uele", md)
    @test occursin("1 of 3 provinces", md)

    ## The optional comparisons appear only when the chain carries them.
    fmd = patch_headline(full, np)
    @test length(filter(startswith("- **"), split(fmd, "\n"))) == 5
    @test occursin(
        r"from 1\d\.\d–2\d\.\d% in Haut-Uele to 3\d\.\d–4\d\.\d% in Ituri",
        fmd
    )
    @test occursin(
        r"national average:\*\* from [\d.]+–[\d.]+ in Haut-Uele to " *
            r"[\d.]+–[\d.]+ in Ituri",
        fmd
    )
    ## Thirty imports against about 1 280 infections.
    @test occursin(r"another province make up 2\.\d–2\.\d%", fmd)

    @test_throws ErrorException patch_headline((; base.C_T_patch), np)
end

@testitem "patch_headline reads importation for the requested patches only" setup = [
    PatchHeadlineDraws,
] begin
    using BVDOutbreakSize: patch_headline

    ## One import a day into the first two patches and a hundred into the
    ## third, flattened patch-fastest.
    skewed = (;
        base...,
        importation_patch = [
            repeat([1.0, 1.0, 100.0], ndays) for _ in 1:nd
        ],
    )
    ## Twenty imports against about 1 200 infections in the first two
    ## patches. Summing the third patch's imports would give over 80%.
    md = patch_headline(skewed, 2)
    @test occursin(r"another province make up 1\.\d–1\.\d%", md)
    ## An importation vector that is not a whole number of days of the
    ## chain's patches is an error rather than a silent misread.
    bad = (; base..., importation_patch = [ones(3 * ndays + 1) for _ in 1:nd])
    @test_throws ErrorException patch_headline(bad, np)
end

@testitem "patch_detail_headline gives each province its intervals" setup = [
    PatchHeadlineDraws,
] begin
    using BVDOutbreakSize: patch_detail_headline

    md = patch_detail_headline(base, np)
    ## A bold lead per province, in patch order, each with its own bullets.
    leads = filter(startswith("**"), split(md, "\n"))
    @test leads == ["**$(l)**" for l in PROVINCE_LABELS[1:np]]
    @test count(startswith("- **"), split(md, "\n")) == 2 * np
    ## Equal-tailed intervals at every level, as on the National page.
    @test count("30% ", md) == 2 * np
    @test count("60% ", md) == 2 * np
    @test count("90% ", md) == 2 * np
    @test !occursin("median", md)
    @test !occursin(r"\d\.\d+ infections", md)
    @test !occursin("Case-fatality", md)

    fmd = patch_detail_headline(full, np)
    @test count(startswith("- **"), split(fmd, "\n")) == 4 * np
    @test count("30% ", fmd) == 4 * np
    ## The case-fatality ratio is written as a percentage.
    @test occursin(r"90% 3\d\.\d–4\d\.\d%", fmd)

    @test_throws ErrorException patch_detail_headline(
        (; base.C_T_patch), np
    )
end

@testitem "median_interval_text reads as plain words" begin
    using BVDOutbreakSize: median_interval_text

    draws = collect(0.0:100.0)
    @test median_interval_text(draws) ==
        "about 50 (90% credible interval 5 to 95)"
    @test median_interval_text(draws ./ 100; digits = 2) ==
        "about 0.50 (90% credible interval 0.05 to 0.95)"
    ## Decimals are padded, so a reproduction number reads 0.80, not 0.8.
    @test median_interval_text(fill(0.8, 10); digits = 2) ==
        "about 0.80 (90% credible interval 0.80 to 0.80)"
    ## A share reads as a percentage through the same phrase.
    @test median_interval_text(draws ./ 100; scale = 100, suffix = "%") ==
        "about 50% (90% credible interval 5% to 95%)"
end
