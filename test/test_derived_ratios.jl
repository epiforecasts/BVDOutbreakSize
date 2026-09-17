## Tests for the ratios a nested model inherits from the joint fit: the
## delay-weighted window totals they divide by, the per-province
## confirmed-case ascertainment built from the case composition's own
## weights, and the chain accessors that hand the draws on.

@testitem "window totals: a delay CDF weighting matches the convolution" begin
    using BVDOutbreakSize: _window_total_from_cdf, convolve_delay

    ## The total of a delayed series to the cut-off is the series weighted
    ## by the delay CDF at the days remaining, so the helper must agree
    ## with summing the convolution it avoids forming.
    x = [2.0, 5.0, 11.0, 7.0, 3.0, 1.0, 9.0, 4.0]
    pmf = [0.2, 0.5, 0.2, 0.1]
    @test _window_total_from_cdf(x, cumsum(pmf)) ≈ sum(convolve_delay(x, pmf))

    ## A point mass at lag zero loses nothing to censoring.
    @test _window_total_from_cdf(x, cumsum([1.0])) ≈ sum(x)

    ## A delay longer than the window keeps only what has cleared it.
    late = _window_total_from_cdf(x, cumsum([0.0, 0.0, 1.0]))
    @test late ≈ sum(x[1:(end - 2)])
    @test late < sum(x)
end

@testitem "province ascertainment: the parts reconstruct the national total" begin
    using BVDOutbreakSize: _province_confirmed_ascertainment, convolve_delay

    ## Three patches on different epidemic phases, so the delay censoring
    ## differs between them and the ratio cannot be a common rescaling of
    ## the cumulative counts.
    n = 30
    onsets = reduce(
        vcat,
        [
            [10.0 * exp(0.05 * t) for t in 1:n]',
            [40.0 + 0.0 * t for t in 1:n]',
            [5.0 * exp(0.12 * t) for t in 1:n]',
        ]
    )
    receipt = [0.3, 0.4, 0.2, 0.1]
    confirmation = [0.0, 0.1, 0.2, 0.3, 0.25, 0.15]
    asc = [1.4, 0.8, 1.0 / (1.4 * 0.8)]      ## geometric mean one
    total = 750.0

    alpha = _province_confirmed_ascertainment(
        onsets, receipt, confirmation,
        asc, total
    )
    @test length(alpha) == 3
    @test all(>(0), alpha)

    ## The infections each province has had time to have confirmed.
    exposed = [
        sum(convolve_delay(vec(onsets[p, :]), confirmation))
            for p in 1:3
    ]
    ## Every confirmed case is attributed to some province, so weighting the
    ## province ratios by their own denominators returns the national total.
    @test sum(alpha .* exposed) ≈ total

    ## A province that finds more of its infections has a higher ratio, and
    ## the split is a partition, so the others must fall.
    raised = _province_confirmed_ascertainment(
        onsets, receipt, confirmation,
        [2.0 * asc[1], asc[2], asc[3]], total
    )
    @test raised[1] > alpha[1]
    @test raised[2] < alpha[2]
    @test raised[3] < alpha[3]

    ## Scaling every province by the same factor is what the composition
    ## cannot see, so it must leave the ratios where they were.
    @test _province_confirmed_ascertainment(
        onsets, receipt, confirmation,
        3.0 .* asc, total
    ) ≈ alpha
end

@testitem "province ascertainment: one kernel makes it a rescaled contrast" begin
    using BVDOutbreakSize: _province_confirmed_ascertainment

    ## When the confirmation delay is the composition's own receipt delay,
    ## the province volumes cancel and the ratio is the sampled contrast on
    ## a common scale. That is the identity the report quotes, so it is
    ## worth pinning rather than leaving to the general case above.
    n = 25
    onsets = reduce(
        vcat,
        [
            [3.0 * exp(0.08 * t) for t in 1:n]',
            [20.0 + 0.0 * t for t in 1:n]',
            [1.0 * exp(0.15 * t) for t in 1:n]',
        ]
    )
    pmf = [0.5, 0.3, 0.2]
    asc = [1.3, 0.9, 1.0 / (1.3 * 0.9)]

    alpha = _province_confirmed_ascertainment(onsets, pmf, pmf, asc, 500.0)
    ratios = alpha ./ asc
    @test all(r -> isapprox(r, ratios[1]; rtol = 1.0e-10), ratios)
end

@testitem "bvd_joint: the inherited ratios reach the chain" tags = [:slow] begin
    using BVDOutbreakSize
    using BVDOutbreakSize: derived_ratio_draws, derived_ratio_table
    using Turing: sample, Prior
    import FlexiChains
    using DataFrames: DataFrame, nrow

    ## A nested health-zone model reads these off the parent chain, so they
    ## have to survive the trip through the sampler under these names.
    obs = load_observations()
    np = length(PROVINCE_NAMES)
    prov = province_increment_matrix(
        obs.province_confirmed_history, PROVINCE_NAMES, np
    )
    prov_deaths = province_increment_matrix(
        obs.province_death_history, PROVINCE_NAMES, np
    )
    m = bvd_joint(
        obs.n,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        n_patches = np,
        province_increments = prov.increments,
        province_days = prov.days,
        province_death_increments = prov_deaths.increments,
        province_death_days = prov_deaths.days,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days
    )

    chn = sample(
        m, Prior(), 50; chain_type = FlexiChains.VNChain, progress = false
    )

    ## The per-province fatality ratio is the per-province fatality
    ## parameter: the incubation map thins nothing, so the deaths a
    ## province's infections go on to cause are that fraction of them, and
    ## no ascertainment step relates the two. A change that broke the
    ## identity would be a change in what the ratio means.
    ifr_patch = vec(collect(chn[:IFR_patch]))
    cfr_patch = vec(collect(chn[:CFR_patch]))
    @test all(v -> length(v) == np, ifr_patch)
    @test ifr_patch == cfr_patch
    ## It pools through the death composition's sum-to-zero contrast, so it
    ## is the pooled ratio times a contrast with geometric mean one.
    @test vec(Array(chn[:IFR])) == vec(Array(chn[:CFR]))
    sev = vec(collect(chn[:province_cfr_relative]))
    @test all(v -> abs(sum(log.(v))) < 1.0e-8, sev)

    asc = vec(Array(chn[:confirmed_ascertainment]))
    @test length(asc) == 50
    @test all(isfinite, asc)
    @test all(>(0), asc)

    per_patch = vec(collect(chn[:province_confirmed_ascertainment]))
    @test all(v -> length(v) == np, per_patch)
    @test all(v -> all(x -> isfinite(x) && x > 0, v), per_patch)

    d = derived_ratio_draws(chn)
    @test length(d.IFR_patch) == np
    @test d.IFR_patch[2] == [v[2] for v in ifr_patch]
    @test length(d.ascertainment_patch) == np
    @test d.ascertainment_patch[2] == [v[2] for v in per_patch]
    @test d.confirmed_ascertainment == asc
    ## The pooling scales come back so a reader can tell a provincial spread
    ## the data found from one the prior put there.
    @test length(d.cfr_pooling_sd) == 50
    @test all(>=(0), d.cfr_pooling_sd)
    @test all(>=(0), d.ascertainment_pooling_sd)

    df = derived_ratio_table(chn)
    @test df isa DataFrame
    ## Two per-province blocks, the two pooled values, the two scales.
    @test nrow(df) == 2 * np + 4
    ## Intervals only: the table reports no central estimate.
    @test names(df) ==
        [
        "Quantity", "Lower 90%", "Lower 60%", "Lower 30%",
        "Upper 30%", "Upper 60%", "Upper 90%",
    ]
    for r in eachrow(df)
        @test r["Lower 90%"] <= r["Lower 60%"] <= r["Lower 30%"]
        @test r["Upper 30%"] <= r["Upper 60%"] <= r["Upper 90%"]
    end
end
