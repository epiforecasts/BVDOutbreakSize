## Province-level isolation occupancy and bed capacity. The national tile
## and the national implied capacity keep their likelihoods on every day;
## the province figures add a split conditional on the printed sum of the
## provinces present that day, over per-patch bed demand and per-patch
## capacity walks.

@testitem "province_split_logpdf: a full day matches the composition model" begin
    using BVDOutbreakSize: province_split_logpdf, province_composition_model
    using Turing: DynamicPPL
    using Random: Xoshiro

    obs = [80 40; 15 8; 5 2]
    level = [8.0 4.0; 1.5 0.8; 0.5 0.2]
    ρ = 0.05
    ## Long format: one row per (day, patch) printed.
    days = [1, 1, 1, 2, 2, 2]
    patches = [1, 2, 3, 1, 2, 3]
    counts = [80, 15, 5, 40, 8, 2]
    lp = province_split_logpdf(days, patches, counts, level, ρ)

    m = DynamicPPL.fix(
        province_composition_model(
            obs, level; ascertainment_sd_prior = nothing
        ); ρ = ρ
    )
    vi = DynamicPPL.VarInfo(Xoshiro(1), m)
    @test lp ≈ DynamicPPL.loglikelihood(m, vi) rtol = 1.0e-10

    ## Invariant to the scale of the modelled level: only the split counts.
    @test province_split_logpdf(days, patches, counts, 10 .* level, ρ) ≈ lp

    ## A day with a province absent conditions on the provinces present.
    part_days = [1, 1, 2, 2, 2]
    part_patches = [1, 3, 1, 2, 3]
    part_counts = [80, 5, 40, 8, 2]
    lp_part = province_split_logpdf(
        part_days, part_patches, part_counts, level, ρ
    )
    two = DynamicPPL.fix(
        province_composition_model(
            [80; 5;;], level[[1, 3], 1:1]; ascertainment_sd_prior = nothing
        ); ρ = ρ
    )
    one_day = DynamicPPL.fix(
        province_composition_model(
            obs[:, 2:2], level[:, 2:2]; ascertainment_sd_prior = nothing
        ); ρ = ρ
    )
    expected = DynamicPPL.loglikelihood(two, DynamicPPL.VarInfo(Xoshiro(1), two)) +
        DynamicPPL.loglikelihood(one_day, DynamicPPL.VarInfo(Xoshiro(1), one_day))
    @test lp_part ≈ expected rtol = 1.0e-10

    ## A day with one province present carries no split and scores zero.
    @test province_split_logpdf([3], [2], [7], level, ρ) == 0.0
    @test province_split_logpdf(Int[], Int[], Int[], level, ρ) == 0.0
end

@testitem "province_care_observations: pooled patches need every member" begin
    using BVDOutbreakSize: province_care_observations

    ## Two source provinces pooled into one patch, one on its own.
    names = ["a", "b", "other"]
    members = Dict("other" => ["c", "d"])
    hist = Dict(
        "a" => (; days = [1, 2, 3], counts = [10, 11, 12]),
        "b" => (; days = [1, 3], counts = [5, 6]),
        "c" => (; days = [2, 3], counts = [2, 3]),
        "d" => (; days = [3], counts = [1]),
    )
    obs = province_care_observations(hist, names; members)
    @test issorted(obs.days)
    rows = Set(zip(obs.days, obs.patches, obs.counts))
    ## Day 1: a and b print; the pooled patch has no member yet, so no row.
    @test (1, 1, 10) in rows
    @test (1, 2, 5) in rows
    @test !any(d == 1 && p == 3 for (d, p, _) in rows)
    ## Day 2: b is silent; c prints for the first time, d has never printed,
    ## so the pooled patch is c alone.
    @test (2, 1, 11) in rows
    @test !any(d == 2 && p == 2 for (d, p, _) in rows)
    @test (2, 3, 2) in rows
    ## Day 3: every member prints, the pooled patch is their sum.
    @test (3, 3, 4) in rows
    @test length(rows) == 7

    ## Day 4: d has printed before and is silent now, so the pooled patch is
    ## absent rather than a partial sum.
    hist["c"] = (; days = [2, 3, 4], counts = [2, 3, 5])
    obs4 = province_care_observations(hist, names; members)
    @test !any(d == 4 && p == 3 for (d, p) in zip(obs4.days, obs4.patches))

    ## Nothing in, nothing out.
    none = Dict{String, @NamedTuple{days::Vector{Int}, counts::Vector{Int}}}()
    @test isempty(province_care_observations(none, names; members).days)
end

@testitem "treatment_flow_model: patch demands sum to the national" begin
    using BVDOutbreakSize: treatment_flow_model, treatment_only_model
    using Random: seed!

    n = 60
    t = 1:n
    bvd = 2.0 .* exp.(0.05 .* t)
    bg = fill(6.0, n)
    reports = vcat(bvd' .* 0.5, bvd' .* 0.5)
    iso = (; days = collect(30:40), counts = fill(40, 11))

    seed!(11)
    res = treatment_flow_model(
        iso, bvd, bg, 0.5, 0.3;
        bvd_reports_matrix = reports,
        background_split = [0.5, 0.5]
    )()
    @test size(res.demand_patch) == (2, n)
    @test size(res.capacity_patch) == (2, n)
    ## An equal split is two half-sized copies of the national pipeline.
    @test vec(sum(res.demand_patch; dims = 1)) ≈ res.demand
    @test res.demand_patch[1, :] ≈ res.demand ./ 2
    @test vec(sum(res.capacity_patch; dims = 1)) ≈ res.capacity_series
    @test all(>(0), res.capacity_patch)

    ## One patch is the national model and exposes no split.
    single = treatment_flow_model(iso, bvd, bg, 0.5, 0.3)()
    @test size(single.demand_patch) == (1, n)
    @test single.demand_patch[1, :] == single.demand
end

@testitem "bvd_joint: province occupancy and beds enter as splits" begin
    using BVDOutbreakSize
    using Turing: DynamicPPL
    using Random: Xoshiro

    obs = load_observations()
    np = length(PROVINCE_NAMES)

    function build(iso_h, cap_h; n_patches = np)
        care = province_care_observations(iso_h, PROVINCE_NAMES)
        beds = province_care_observations(cap_h, PROVINCE_NAMES)
        return bvd_joint(
            obs.n,
            obs.exported_cases, obs.total_deaths, obs.reported_cases,
            obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
            confirmed_deaths = obs.confirmed_deaths,
            recovered_cases = obs.recovered_cases,
            deaths_history = obs.deaths_history,
            reported_history = obs.reported_history,
            confirmed_history = obs.confirmed_history,
            confirmed_deaths_history = obs.confirmed_deaths_history,
            lab_history = obs.lab_history,
            lab_daily_history = obs.lab_daily_history,
            suspected_daily_history = obs.suspected_daily_history,
            isolation_history = obs.isolation_history,
            bed_capacity_history = obs.bed_capacity_history,
            occupancy_break_days = obs.occupancy_break_days,
            export_case_days = obs.export_case_days,
            export_death_days = obs.export_death_days,
            breakpoint = obs.who_first_sitrep_days,
            n_patches,
            province_isolation = care,
            province_capacity = beds,
            tmrca_days = obs.tmrca_days
        )
    end

    ## Synthetic province rows on real isolation days, so the test does not
    ## wait on the transcription: Ituri takes most of the tile, Nord-Kivu the
    ## rest, Haut-Uélé a little from late July.
    iso = obs.isolation_history
    late = iso.days .>= iso.days[end] - 40
    iso_h = Dict(
        "ituri" => (; days = iso.days, counts = round.(Int, 0.7 .* iso.counts)),
        "nord_kivu" => (; days = iso.days, counts = round.(Int, 0.25 .* iso.counts)),
        "haut_uele" => (;
            days = iso.days[late],
            counts = round.(Int, 0.05 .* iso.counts[late]),
        ),
    )
    cap = obs.bed_capacity_history
    cap_h = Dict(
        "ituri" => (; days = cap.days, counts = round.(Int, 0.7 .* cap.counts)),
        "nord_kivu" => (; days = cap.days, counts = round.(Int, 0.2 .* cap.counts)),
    )

    m = build(iso_h, cap_h)
    vi = DynamicPPL.VarInfo(Xoshiro(7), m)
    base = DynamicPPL.logjoint(m, vi)
    @test isfinite(base)

    ## Moving occupied beds between provinces at a fixed printed sum moves
    ## the density.
    shifted = Dict(k => v for (k, v) in iso_h)
    mv = round.(Int, 0.2 .* iso.counts)
    shifted["ituri"] = (; days = iso.days, counts = iso_h["ituri"].counts .- mv)
    shifted["nord_kivu"] = (;
        days = iso.days, counts = iso_h["nord_kivu"].counts .+ mv,
    )
    @test !isapprox(
        DynamicPPL.logjoint(build(shifted, cap_h), vi), base; rtol = 1.0e-8
    )

    ## So does moving beds.
    shifted_cap = Dict(k => v for (k, v) in cap_h)
    mvc = round.(Int, 0.3 .* cap.counts)
    shifted_cap["ituri"] = (; days = cap.days, counts = cap_h["ituri"].counts .- mvc)
    shifted_cap["nord_kivu"] = (;
        days = cap.days, counts = cap_h["nord_kivu"].counts .+ mvc,
    )
    @test !isapprox(
        DynamicPPL.logjoint(build(iso_h, shifted_cap), vi), base; rtol = 1.0e-8
    )

    ## Without province rows the terms are absent and the density is finite.
    none = Dict{String, @NamedTuple{days::Vector{Int}, counts::Vector{Int}}}()
    @test isfinite(DynamicPPL.logjoint(build(none, none), vi))

    ## Province care data with one patch would be dropped silently.
    @test_throws ErrorException build(iso_h, cap_h; n_patches = 1)()
end

@testitem "province_bed_table: beds, demand and shortfall by province" tags = [
    :slow,
] begin
    using BVDOutbreakSize
    using Turing: sample, Prior
    import FlexiChains
    using DataFrames: nrow, names

    obs = load_observations()
    np = length(PROVINCE_NAMES)
    m = bvd_joint(
        obs.n,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        isolation_history = obs.isolation_history,
        bed_capacity_history = obs.bed_capacity_history,
        n_patches = np,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days
    )
    chn = sample(
        m, Prior(), 50; chain_type = FlexiChains.VNChain, progress = false
    )
    df = province_bed_table(chn, np)
    @test nrow(df) == np
    @test names(df) == [
        "Province", "Beds", "Bed demand", "Occupied beds",
        "Utilisation (%)", "Shortfall",
    ]
    @test df.Province == collect(PROVINCE_LABELS[1:np])
    ## Every cell is a median with an interval.
    @test all(contains("("), df[!, "Beds"])

    ## A single-population chain carries no per-province beds.
    single = sample(
        bvd_joint(
            obs.n,
            obs.exported_cases, obs.total_deaths, obs.reported_cases,
            obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
            reported_history = obs.reported_history,
            isolation_history = obs.isolation_history,
            breakpoint = obs.who_first_sitrep_days,
            tmrca_days = obs.tmrca_days
        ),
        Prior(), 5; chain_type = FlexiChains.VNChain, progress = false
    )
    @test_throws ErrorException province_bed_table(single, np)
end

@testitem "province care blocks: sparse, dated and inside the grid" begin
    using BVDOutbreakSize

    obs = load_observations()
    iso = obs.province_isolation_history
    beds = obs.province_bed_capacity_history
    @test length(iso) >= 5
    @test length(beds) >= 5
    ## The province tables start before the national tile does (the May
    ## facility tables), so the rows are bounded by the grid, not the tile.
    hi = maximum(obs.isolation_history.days)
    for (name, h) in iso
        @test issorted(h.days) && allunique(h.days)
        @test all(>=(0), h.counts)
        @test all(d -> 1 <= d <= hi, h.days)
    end
    for (name, h) in beds
        @test issorted(h.days) && allunique(h.days)
        @test all(>(0), h.counts)
    end
    ## The two large provinces print almost every day.
    @test length(iso["ituri"].days) > 60
    @test length(iso["nord_kivu"].days) > 60
    ## The long format the model scores has a split on most days.
    care = province_care_observations(iso, PROVINCE_NAMES)
    days = unique(care.days)
    @test count(d -> count(==(d), care.days) >= 2, days) > 50
end
