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

    ## Province rows switch the per-patch demand and capacity on.
    rows = (; days = [35, 35, 36, 36], patches = [1, 2, 1, 2], counts = [20, 20, 21, 19])
    seed!(11)
    res = treatment_flow_model(
        iso, bvd, bg, 0.5, 0.3;
        bvd_reports_matrix = reports,
        background_split = [0.5, 0.5],
        province_isolation = rows
    )()
    @test size(res.demand_patch) == (2, n)
    @test size(res.capacity_patch) == (2, n)
    ## An equal split is two half-sized copies of the national pipeline.
    @test vec(sum(res.demand_patch; dims = 1)) ≈ res.demand
    @test res.demand_patch[1, :] ≈ res.demand ./ 2
    @test vec(sum(res.capacity_patch; dims = 1)) ≈ res.capacity_series
    @test all(>(0), res.capacity_patch)
    ## The capacity share is daily, a simplex on each day.
    @test size(res.capacity_shares) == (2, n)
    @test vec(sum(res.capacity_shares; dims = 1)) ≈ ones(n)
    ## The pooled deviation scale is returned for diagnostics.
    @test res.capacity_pooling_sd >= 0

    ## One patch is the national model and exposes no split.
    single = treatment_flow_model(iso, bvd, bg, 0.5, 0.3)()
    @test size(single.demand_patch) == (1, n)
    @test single.demand_patch[1, :] == single.demand
    ## Patches without province rows stay on the national walk and demand.
    quiet = treatment_flow_model(
        iso, bvd, bg, 0.5, 0.3;
        bvd_reports_matrix = reports, background_split = [0.5, 0.5]
    )()
    @test size(quiet.demand_patch) == (1, n)
    @test size(quiet.capacity_patch) == (1, n)
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

    ## Without province rows the terms are absent, the stream is back on the
    ## national walk, and the density is finite at that model's own draw.
    none = Dict{String, @NamedTuple{days::Vector{Int}, counts::Vector{Int}}}()
    m_none = build(none, none)
    vi_none = DynamicPPL.VarInfo(Xoshiro(7), m_none)
    @test isfinite(DynamicPPL.logjoint(m_none, vi_none))
    @test !any(contains("occupancy_split_rho"), string.(keys(vi_none)))

    ## Province care data with one patch would be dropped silently.
    @test_throws ErrorException build(iso_h, cap_h; n_patches = 1)()
end

@testitem "bvd_joint: province admissions enter as a split" begin
    using BVDOutbreakSize
    using Turing: DynamicPPL
    using Random: Xoshiro

    obs = load_observations()
    adm = obs.treatment_admissions_history
    function build(adm_h)
        return bvd_joint(
            obs.n,
            obs.exported_cases, obs.total_deaths, obs.reported_cases,
            obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
            reported_history = obs.reported_history,
            confirmed_history = obs.confirmed_history,
            isolation_history = obs.isolation_history,
            bed_capacity_history = obs.bed_capacity_history,
            treatment_admissions_history = adm,
            breakpoint = obs.who_first_sitrep_days,
            n_patches = length(PROVINCE_NAMES),
            province_admissions = province_care_observations(
                adm_h, PROVINCE_NAMES
            ),
            tmrca_days = obs.tmrca_days
        )
    end
    ## Synthetic province rows on the national admission days.
    adm_h = Dict(
        "ituri" => (; days = adm.days, counts = round.(Int, 0.7 .* adm.counts)),
        "nord_kivu" => (;
            days = adm.days, counts = round.(Int, 0.3 .* adm.counts),
        ),
    )
    m = build(adm_h)
    vi = DynamicPPL.VarInfo(Xoshiro(7), m)
    base = DynamicPPL.logjoint(m, vi)
    @test isfinite(base)
    @test any(contains("admissions_split_rho"), string.(keys(vi)))
    ## Moving admissions between provinces at a fixed sum moves the density.
    mv = round.(Int, 0.2 .* adm.counts)
    shifted = Dict(
        "ituri" => (; days = adm.days, counts = adm_h["ituri"].counts .- mv),
        "nord_kivu" => (;
            days = adm.days, counts = adm_h["nord_kivu"].counts .+ mv,
        ),
    )
    @test !isapprox(
        DynamicPPL.logjoint(build(shifted), vi), base; rtol = 1.0e-8
    )
end

@testitem "_median_ci: drops non-finite draws" begin
    using BVDOutbreakSize: _median_ci

    ## Median 2, and the 5% and 95% quantiles of 1, 2, 3 are 1.1 and 2.9.
    @test _median_ci([1.0, NaN, 2.0, Inf, 3.0, -Inf]; digits = 1) ==
        "2.0 (1.1–2.9)"
    @test _median_ci([1.0, NaN, 2.0, 3.0]; digits = 0) == "2 (1–3)"
    @test _median_ci([NaN, Inf]) == "—"
end

@testitem "province_bed_table: beds, demand and occupancy by province" tags = [
    :slow,
] begin
    using BVDOutbreakSize
    using Turing: sample, Prior
    using Random: Xoshiro
    import FlexiChains
    using DataFrames: nrow, names

    ## The headline joint carries the province occupancy and bed splits.
    obs = load_observations()
    np = length(PROVINCE_NAMES)
    m = production_joint(obs; breakpoint = default_breakpoint(obs))
    chn = sample(
        Xoshiro(20261006), m, Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    df = province_bed_table(chn, np)
    @test nrow(df) == np
    @test names(df) == [
        "Province", "Beds", "Bed demand", "Occupied beds",
        "Utilisation (%)", "Patients above beds",
    ]
    @test df.Province == collect(PROVINCE_LABELS[1:np])
    ## Every cell is a median with an interval.
    @test all(contains("("), df[!, "Beds"])
    ## The national occupancy and beds are the provinces' sums, and a
    ## province's patients above beds are its occupancy over its beds.
    ## Nothing caps the occupancy at the beds. The national shortfall is
    ## the bed demand above the national beds.
    occ = BVDOutbreakSize._per_patch(chn, :province_expected_isolation, np)
    beds = BVDOutbreakSize._per_patch(chn, :province_bed_capacity, np)
    short = BVDOutbreakSize._per_patch(chn, :province_bed_shortfall, np)
    @test all(isapprox.(sum(occ), vec(Array(chn[:expected_isolation_T])); rtol = 1.0e-6, atol = 1.0e-8))
    @test all(isapprox.(sum(beds), vec(Array(chn[:bed_capacity])); rtol = 1.0e-6))
    demand = vec(Array(chn[:expected_bed_demand_T]))
    national_short = vec(Array(chn[:bed_shortfall_T]))
    @test all(
        isapprox.(
            national_short, max.(demand .- sum(beds), 0.0);
            rtol = 1.0e-6, atol = 1.0e-8
        )
    )
    for p in 1:np
        @test short[p] ≈ max.(occ[p] .- beds[p], 0.0)
    end

    ## A single-population chain carries no per-province beds.
    single = sample(
        Xoshiro(20261006),
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
    adm = obs.province_admissions_history
    @test length(adm) >= 5
    for (name, h) in adm
        @test issorted(h.days) && allunique(h.days)
        @test all(>=(0), h.counts)
    end
    ## On a day every province that has printed prints, the provinces'
    ## admissions add up to the national admissions.
    nat = Dict(
        zip(
            obs.treatment_admissions_history.days,
            obs.treatment_admissions_history.counts
        )
    )
    active(d) = [h for h in values(adm) if first(h.days) <= d]
    full = filter(
        d -> haskey(nat, d) && all(h -> d in h.days, active(d)),
        sort!(unique!(reduce(vcat, [h.days for h in values(adm)])))
    )
    @test length(full) >= 20
    for d in full
        @test sum(h.counts[findfirst(==(d), h.days)] for h in active(d)) ==
            nat[d]
    end
    ## The two large provinces print almost every day.
    @test length(iso["ituri"].days) > 60
    @test length(iso["nord_kivu"].days) > 60
    ## The long format the model scores has a split on most days.
    care = province_care_observations(iso, PROVINCE_NAMES)
    days = unique(care.days)
    @test count(d -> count(==(d), care.days) >= 2, days) > 50
end

@testitem "province_care_observations: weekly and change-day thinning" begin
    using BVDOutbreakSize: province_care_observations

    names = ["a", "b"]
    members = Dict{String, Vector{String}}()
    days = collect(1:21)
    hist = Dict(
        "a" => (; days = days, counts = fill(10, 21)),
        "b" => (; days = days, counts = vcat(fill(5, 10), fill(7, 11))),
    )
    weekly = province_care_observations(hist, names; members, every = 7)
    @test unique(weekly.days) == [1, 8, 15]
    @test length(weekly.days) == 6

    ## Change days: the first day, and the day b moves from 5 to 7.
    changes = province_care_observations(hist, names; members, changes_only = true)
    @test unique(changes.days) == [1, 11]
    @test length(changes.days) == 4

    ## Days with a single province carry no split and do not count in the
    ## weekly spacing.
    hist["b"] = (; days = [1, 2, 9, 10], counts = [5, 5, 5, 5])
    weekly2 = province_care_observations(hist, names; members, every = 7)
    @test unique(weekly2.days) == [1, 9]
    @test all(d -> count(==(d), weekly2.days) == 2, unique(weekly2.days))
end

@testitem "per-patch background stock matches the national running balance" begin
    using BVDOutbreakSize: accumulate_occupancy, abscond_thinned,
        convolve_delay, _background_stay_survival

    ## With no BVD admissions and no confirmation, the national balance holds
    ## only the background stock, so the per-cohort survival the province
    ## split convolves with must reproduce it day by day.
    n = 40
    A_bg = [5.0 + 3 * sin(t / 4) for t in 1:n]
    f = [0.1, 0.3, 0.25, 0.2, 0.1, 0.05]
    κ = 0.03
    z = zeros(n)
    ruleout = convolve_delay(A_bg, abscond_thinned(f, κ))
    national = accumulate_occupancy(z, A_bg, z, z, ruleout, κ, z).demand
    patch = convolve_delay(A_bg, _background_stay_survival(f, κ))
    @test patch ≈ national rtol = 1.0e-10
end

@testitem "background admission shares nest the constant share" begin
    using BVDOutbreakSize: background_admission_share_model, knot_days
    using Turing: returned, fix, @varname
    using Random: Xoshiro

    w = [0.5, 0.3, 0.15, 0.05]
    n = 80
    days = collect(20:60)
    m = background_admission_share_model(w, n, days)
    θ = rand(Xoshiro(3), m)
    r = returned(m, θ)
    ## Daily simplex columns that move only over the printed days.
    @test size(r.W) == (4, n)
    @test all(c -> isapprox(sum(c), 1; rtol = 1.0e-12), eachcol(r.W))
    @test all(>(0), r.W)
    @test all(t -> r.W[:, t] ≈ r.W[:, 20], 1:20)
    @test all(t -> r.W[:, t] ≈ r.W[:, 60], 60:n)
    @test !(r.W[:, 20] ≈ r.W[:, 60])
    @test r.level_sd > 0 && r.walk_sd > 0
    nk = length(knot_days(60; start = 20))
    z = θ[@varname(z_bga)]
    @test length(z) == 3 * (nk - 1)
    ## A zero level and zero walk give back the background suspect shares
    ## on every day.
    flat = fix(m, (; y_bga = zeros(3), σ_bga = 0.0, z_bga = z))
    @test all(c -> c ≈ w, eachcol(returned(flat, rand(Xoshiro(4), flat)).W))
    ## A level with no walk moves the shares by the same amount every day.
    lvl = fix(
        m, (; y_bga = [0.4, -0.2, 0.1], σ_bga = 0.0, z_bga = z)
    )
    W = returned(lvl, rand(Xoshiro(5), lvl)).W
    @test all(t -> W[:, t] ≈ W[:, 1], 1:n)
    @test !(W[:, 1] ≈ w)
    ## No printed days, or one patch, sample nothing.
    @test background_admission_share_model(w, n, Int[])().W == repeat(w, 1, n)
    @test background_admission_share_model([1.0], n, days)().W == ones(1, n)
end

@testitem "background admission shares hold provinces without admissions" begin
    using BVDOutbreakSize: background_admission_share_model, knot_days
    using Turing: returned, fix, @varname
    using Random: Xoshiro

    ## Only patches 1 and 2 print admissions, so the deviations are one
    ## contrast between them, and patches 3 and 4 keep the ratio of their
    ## background suspect shares on every day.
    w = [0.5, 0.3, 0.15, 0.05]
    n = 80
    days = collect(20:60)
    m = background_admission_share_model(w, n, days; patches = [1, 2, 1, 2])
    θ = rand(Xoshiro(13), m)
    nk = length(knot_days(60; start = 20))
    @test length(θ[@varname(y_bga)]) == 1
    @test length(θ[@varname(z_bga)]) == nk - 1
    W = returned(m, θ).W
    @test all(c -> isapprox(sum(c), 1; rtol = 1.0e-12), eachcol(W))
    @test all(t -> W[3, t] / W[4, t] ≈ w[3] / w[4], 1:n)
    ## The printing patches move against each other.
    @test !(W[1, 30] / W[2, 30] ≈ w[1] / w[2])
    ## A zero level and walk give back the background suspect shares.
    flat = fix(m, (; y_bga = [0.0], σ_bga = 0.0, z_bga = θ[@varname(z_bga)]))
    @test all(c -> c ≈ w, eachcol(returned(flat, rand(Xoshiro(14), flat)).W))
    ## One printing patch carries no contrast, so nothing is sampled.
    one_patch = background_admission_share_model(w, n, days; patches = [2])
    @test isempty(keys(rand(Xoshiro(15), one_patch)))
    @test one_patch().W == repeat(w, 1, n)
end

@testitem "background admission shares on one printed day have a level and no walk" begin
    using BVDOutbreakSize: background_admission_share_model
    using Turing: returned, fix, @varname
    using Random: Xoshiro

    ## A single printed day gives a single knot, so the level is drawn and
    ## the walk is not.
    w = [0.6, 0.3, 0.1]
    n = 40
    m = background_admission_share_model(w, n, [25])
    θ = rand(Xoshiro(7), m)
    @test haskey(θ, @varname(y_bga))
    @test !haskey(θ, @varname(σ_bga))
    @test !haskey(θ, @varname(z_bga))
    r = returned(m, θ)
    @test size(r.W) == (3, n)
    @test r.walk_sd == 0
    @test r.level_sd > 0
    @test all(c -> isapprox(sum(c), 1; rtol = 1.0e-12), eachcol(r.W))
    @test all(t -> r.W[:, t] ≈ r.W[:, 1], 1:n)
    ## A zero level gives back the background suspect shares, and a nonzero
    ## one moves them by the same amount on every day.
    flat = fix(m, (; y_bga = zeros(2)))
    @test all(c -> c ≈ w, eachcol(returned(flat, rand(Xoshiro(8), flat)).W))
    lvl = fix(m, (; y_bga = [0.5, -0.3]))
    W = returned(lvl, rand(Xoshiro(9), lvl)).W
    @test !(W[:, 1] ≈ w)
    @test all(t -> W[:, t] ≈ W[:, 1], 1:n)
end

@testitem "background admission shares carry into the patch admissions and demand" begin
    using BVDOutbreakSize: _patch_admissions, _patch_demand,
        _background_stay_survival, convolve_delay

    n = 30
    A_bvd = [fill(4.0, 1, n); fill(1.0, 1, n)]
    A_bg = fill(6.0, n)
    w = [0.25, 0.75]
    S = [1.0, 0.8, 0.5, 0.2]
    f = [0.4, 0.4, 0.2]
    κ = 0.02
    demand = fill(50.0, n)
    ## Constant daily shares match the vector form.
    Wc = repeat(w, 1, n)
    @test _patch_admissions(A_bvd, A_bg, Wc) ≈ _patch_admissions(A_bvd, A_bg, w)
    @test _patch_demand(A_bvd, A_bg, Wc, S, f, κ, demand) ≈
        _patch_demand(A_bvd, A_bg, w, S, f, κ, demand)
    ## Moving background to patch 1 from day 15 moves its admissions at once
    ## and its demand share after, through the rule-out stay.
    Wm = copy(Wc)
    Wm[:, 15:end] .= [0.75, 0.25]
    adm = _patch_admissions(A_bvd, A_bg, Wm)
    @test vec(sum(adm; dims = 1)) ≈ vec(sum(A_bvd; dims = 1)) .+ A_bg
    @test adm[1, 15] ≈ 4.0 + 4.5
    dm = _patch_demand(A_bvd, A_bg, Wm, S, f, κ, demand)
    d0 = _patch_demand(A_bvd, A_bg, w, S, f, κ, demand)
    @test vec(sum(dm; dims = 1)) ≈ demand
    @test dm[:, 1:14] ≈ d0[:, 1:14]
    @test dm[1, 20] > d0[1, 20]
end

@testitem "background admission shares: Mooncake gradient matches ForwardDiff" tags = [
    :ad,
] begin
    using BVDOutbreakSize: background_admission_share_model,
        province_split_model, _patch_admissions, default_adtype
    using Turing: @model, to_submodel, DynamicPPL
    using Distributions: Normal, truncated
    using LogDensityProblems: logdensity_and_gradient
    using ADTypes: AutoForwardDiff
    import ForwardDiff
    using Random: Xoshiro

    ## The walk feeding the admissions split, as in the treatment model.
    n = 30
    days = collect(5:28)
    A_bvd = [fill(4.0, 1, n); fill(1.0, 1, n); fill(0.5, 1, n)]
    A_bg = fill(6.0, n)
    rows = (;
        days = repeat(days; inner = 3), patches = repeat(1:3, length(days)),
        counts = repeat([7, 3, 2], length(days)),
    )
    @model function split_demo()
        st ~ to_submodel(background_admission_share_model([0.5, 0.3, 0.2], n, days))
        ρ ~ truncated(Normal(0, 0.1); lower = 0, upper = 1)
        sp ~ to_submodel(province_split_model(rows, _patch_admissions(A_bvd, A_bg, st.W), ρ))
    end
    model = split_demo()
    vi = DynamicPPL.link(DynamicPPL.VarInfo(Xoshiro(11), model), model)
    x = collect(vi[:])
    grad_with(adtype) = logdensity_and_gradient(
        DynamicPPL.LogDensityFunction(model, DynamicPPL.getlogjoint, vi; adtype), x
    )
    logp, grad = grad_with(default_adtype())
    @test isfinite(logp)
    @test all(isfinite, grad)
    _, grad_fd = grad_with(AutoForwardDiff())
    @test grad ≈ grad_fd rtol = 1.0e-6
end

@testitem "province effective beds: printed or rate-implied" begin
    using BVDOutbreakSize
    include(
        joinpath(
            pkgdir(BVDOutbreakSize), "scripts", "province_effective_beds.jl"
        )
    )

    ## Nord-Kivu, 15 July: 171 patients at 118.8% against 141 printed beds.
    ## The rate is computed on 144 beds, and the patients are not raised
    ## into beds, so the province is over its beds.
    @test effective_beds(141, 171, 118.8) == 144
    ## 21 September: 338 patients at 66.2% against a stale 308.
    @test effective_beds(308, 338, 66.2) == 511
    ## A province under its printed beds keeps them, whatever its rounded
    ## rate implies (Ituri, 15 July: 686 at 98.99% would give 693 + 2).
    @test effective_beds(25, 10, 40.0) == 25
    @test effective_beds(693, 686, 98.7) == 693
    ## Haut-Uélé, 6 September: a rate that disagrees with the printed beds.
    @test effective_beds(128, 68, 50.8) == 128
    ## No rate, or no patients, printed.
    @test effective_beds(20, 24, nothing) == 20
    @test effective_beds(20, nothing, nothing) == 20

    ## The recorded beds are not raised to the patients: Nord-Kivu prints
    ## more patients than beds on some recorded day.
    obs = load_observations()
    h = obs.province_bed_capacity_history["nord_kivu"]
    iso = obs.province_isolation_history["nord_kivu"]
    held = Dict(zip(iso.days, iso.counts))
    @test any(get(held, d, 0) > b for (d, b) in zip(h.days, h.counts))
end

@testitem "incomplete_capacity_days: a silent province's beds count" begin
    using BVDOutbreakSize: incomplete_capacity_days

    ## Province A prints every day. B prints on days 1, 2 and 5, and C
    ## prints 5 beds on day 1 only.
    provinces = Dict(
        "a" => (; days = [1, 2, 3, 4, 5], counts = [100, 100, 100, 100, 100]),
        "b" => (; days = [1, 2, 5], counts = [50, 60, 60]),
        "c" => (; days = [1], counts = [5]),
    )
    ## Day 0 predates every province print, so nothing is missing.
    national = (;
        days = [0, 1, 2, 3, 4, 5], counts = [90, 155, 165, 105, 160, 160],
    )
    ## Day 3 lacks B's 60 beds (36%). Day 4 carries B's beds although the
    ## province block does not, so only C's 5 beds (3%) are missing, under
    ## the tolerance.
    @test incomplete_capacity_days(national, provinces) == [3]
    @test incomplete_capacity_days(national, provinces; tolerance = 0.01) ==
        [3, 4, 5]
    @test isempty(incomplete_capacity_days(national, Dict()))

    ## Two small provinces silent together miss 8% of the beds, but neither
    ## holds 5% on its own, so the day is kept.
    small = Dict(
        "a" => (; days = [1, 2], counts = [92, 92]),
        "d" => (; days = [1], counts = [4]),
        "e" => (; days = [1], counts = [4]),
    )
    @test isempty(
        incomplete_capacity_days((; days = [1, 2], counts = [100, 92]), small)
    )

    ## A large province that stops printing counts as silent for 14 days
    ## after its last print, then as no longer reporting.
    stopped = Dict(
        "a" => (; days = collect(1:30), counts = fill(80, 30)),
        "b" => (; days = [1], counts = [20]),
    )
    gone = (; days = collect(2:30), counts = fill(80, 29))
    @test incomplete_capacity_days(gone, stopped) == collect(2:15)
    @test incomplete_capacity_days(gone, stopped; max_silence = 5) ==
        collect(2:6)
    ## A tolerance outside [0, 1] would silently disable the check.
    @test_throws ArgumentError incomplete_capacity_days(
        gone, stopped; tolerance = 1.5
    )
    @test_throws ArgumentError incomplete_capacity_days(
        gone, stopped; tolerance = -0.1
    )
end

@testitem "load_observations drops national capacity days a province misses" begin
    using BVDOutbreakSize: load_observations, grid_date
    using Dates: Date

    obs = load_observations()
    kept = Set(grid_date.(Ref(obs), obs.bed_capacity_history.days))
    ## Days on which Nord-Kivu is silent and the national figure is its
    ## provinces' partial sum are among those dropped.
    dips = vcat(
        Date(2026, 8, 4):Date(2026, 8, 12),
        Date(2026, 9, 10):Date(2026, 9, 15),
        [Date(2026, 9, 23), Date(2026, 9, 25)],
    )
    @test isempty(intersect(kept, dips))
    ## Days on which every large province prints are kept, as is 26 August,
    ## whose national figure carries Haut-Uélé's beds although the province
    ## block does not.
    complete = (
        Date(2026, 8, 26), Date(2026, 9, 9), Date(2026, 9, 16),
        Date(2026, 9, 26),
    )
    @test all(in(kept), complete)
end

@testitem "province care rules: occupancy after the last national tile" begin
    using BVDOutbreakSize, Dates
    include(joinpath(pkgdir(BVDOutbreakSize), "scripts", "province_care_rules.jl"))

    manifest = Dict(
        "isolation_history" => Dict("dates" => [Date(2026, 10, 23), Date(2026, 10, 24)])
    )
    @test last_isolation_date(manifest) == Date(2026, 10, 24)
    @test_throws ErrorException last_isolation_date(Dict{String, Any}())

    blocks = Dict(
        "nord_kivu" => [(Date(2026, 10, 24), 495), (Date(2026, 10, 25), 497)],
        "ituri" => [(Date(2026, 10, 24), 422)],
        "bas_uele" => [(Date(2026, 10, 25), 38)],
    )
    dropped = drop_after!(blocks, Date(2026, 10, 24))
    @test dropped == [
        ("bas_uele", Date(2026, 10, 25), 38), ("nord_kivu", Date(2026, 10, 25), 497),
    ]
    @test blocks["nord_kivu"] == [(Date(2026, 10, 24), 495)]
    @test blocks["ituri"] == [(Date(2026, 10, 24), 422)]
    @test !haskey(blocks, "bas_uele")
end

@testitem "province care data: SitRep 164 cells held back or left out" begin
    using BVDOutbreakSize, TOML
    manifest = TOML.parsefile(
        joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml")
    )
    adm = manifest["province_admissions_history"]
    @test !("2026-10-25" in string.(adm["nord_kivu"]["dates"]))
    @test "2026-10-25" in string.(adm["ituri"]["dates"])
    last_tile = maximum(string.(manifest["isolation_history"]["dates"]))
    for (prov, h) in manifest["province_isolation_history"]
        h isa Dict || continue
        @test all(<=(last_tile), string.(h["dates"]))
    end
end
