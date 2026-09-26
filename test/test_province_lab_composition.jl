@testitem "background_split_model: one patch is the whole background" begin
    using BVDOutbreakSize: background_split_model

    res = background_split_model(1)()
    @test res.w == [1.0]
end

@testitem "background_split_model: a simplex centred on population share" begin
    using BVDOutbreakSize: background_split_model
    using Turing: DynamicPPL

    pops = [4.0, 2.0, 1.0, 1.0]
    m = background_split_model(4; populations = pops)
    res = m()
    @test length(res.w) == 4
    @test sum(res.w) ≈ 1.0
    @test all(>(0), res.w)

    ## With every log-ratio at zero each patch sits on its population share,
    ## so the prior centre is the population split.
    flat = DynamicPPL.fix(m; bg_log_ratio = zeros(3))()
    @test flat.w ≈ pops ./ sum(pops)

    ## A log-ratio moves its patch against the reference by exactly that
    ## factor, and a wide one cannot overflow the normalisation.
    moved = DynamicPPL.fix(m; bg_log_ratio = [log(2.0), 0.0, 0.0])()
    @test moved.w[2] / moved.w[1] ≈ 2 * pops[2] / pops[1]
    wide = DynamicPPL.fix(m; bg_log_ratio = [800.0, 0.0, -800.0])()
    @test all(isfinite, wide.w)
    @test sum(wide.w) ≈ 1.0
end

@testitem "patch_capacity_share_model: fixed-scale log-ratios, no pooling scale" begin
    using BVDOutbreakSize: patch_capacity_share_model
    using Turing: DynamicPPL
    using Random: Xoshiro

    @test patch_capacity_share_model(1)().s == [1.0]

    pops = [4.0, 2.0, 1.0, 1.0]
    m = patch_capacity_share_model(4; populations = pops)
    names = string.(keys(DynamicPPL.VarInfo(Xoshiro(3), m)))
    @test names == ["cap_log_ratio"]
    flat = DynamicPPL.fix(m; cap_log_ratio = zeros(3))()
    @test flat.s ≈ pops ./ sum(pops)
    ## Log-ratios within two prior sds move the reference patch from half
    ## the population to most of the beds.
    s = DynamicPPL.fix(m; cap_log_ratio = [-2.0, -2.5, -4.0])().s
    @test sum(s) ≈ 1.0
    @test s[1] > 0.7
end

@testitem "province_lab_increment_matrix: pools members, matches the national" begin
    using BVDOutbreakSize

    obs = load_observations()
    np = length(PROVINCE_NAMES)
    lab = province_lab_increment_matrix(
        obs.province_lab_daily_history, PROVINCE_NAMES, np
    )
    @test size(lab.increments) == (np, length(lab.days))
    @test issorted(lab.days)

    ## The per-province analysed counts partition the national daily
    ## analysed series, so each column sums to the national count that day.
    nat = Dict(zip(obs.lab_daily_history.days, obs.lab_daily_history.counts))
    shared = [i for (i, d) in enumerate(lab.days) if haskey(nat, d)]
    @test length(shared) > 50
    for i in shared
        @test sum(@view lab.increments[:, i]) == nat[lab.days[i]]
    end

    ## A pooled patch is the sum of its members.
    members = PROVINCE_MEMBERS[PROVINCE_NAMES[end]]
    pooled = sum(
        obs.province_lab_daily_history["$(m)_analysed"].counts
            for m in members
    )
    @test lab.increments[end, :] == pooled

    ## Weekly bins sum the daily columns, calendar weeks from the first day.
    weekly = province_lab_increment_matrix(
        obs.province_lab_daily_history, PROVINCE_NAMES, np; every = 7
    )
    @test weekly.days == lab.days
    @test length(weekly.bins) == length(lab.days)
    @test size(weekly.increments, 2) == maximum(weekly.bins)
    @test size(weekly.increments, 2) < length(lab.days) ÷ 4
    @test sum(weekly.increments) == sum(lab.increments)
    for (i, b) in enumerate(weekly.bins)
        @test fld(lab.days[i] - lab.days[1], 7) ==
            fld(lab.days[findfirst(==(b), weekly.bins)] - lab.days[1], 7)
    end
    @test lab.bins == 1:length(lab.days)

    ## No data, no term.
    none = Dict{String, @NamedTuple{days::Vector{Int}, counts::Vector{Int}}}()
    @test isempty(province_lab_increment_matrix(none, PROVINCE_NAMES, np).days)
end

@testitem "province_composition_model: the ascertainment contrast can be off" begin
    using BVDOutbreakSize: province_composition_model
    using Turing: DynamicPPL
    using Random: Xoshiro

    obs = [80 40; 15 8; 5 2]
    modelled = [8.0 4.0; 1.5 0.8; 0.5 0.2]
    m = province_composition_model(
        obs, modelled; ascertainment_sd_prior = nothing
    )
    vi = DynamicPPL.VarInfo(Xoshiro(1), m)
    names = string.(keys(vi))
    @test !any(contains("τ_asc"), names)
    @test !any(contains("z_asc"), names)
    @test isfinite(DynamicPPL.loglikelihood(m, vi))

    res = m()
    @test res.province_ascertainment == ones(3)
    @test res.ascertainment_sd == 0.0

    ## Off, the shares are the modelled split and nothing else.
    expected = modelled ./ sum(modelled; dims = 1)
    @test res.shares ≈ expected
end

@testitem "bvd_joint: the province analysed volume is scored as a composition" begin
    using BVDOutbreakSize
    using Turing: DynamicPPL
    using Random: Xoshiro

    obs = load_observations()
    np = length(PROVINCE_NAMES)

    function build(labh; n_patches = np)
        lab = province_lab_increment_matrix(labh, PROVINCE_NAMES, n_patches)
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
            suspected_daily_deaths_history = obs.suspected_daily_deaths_history,
            isolation_history = obs.isolation_history,
            bed_capacity_history = obs.bed_capacity_history,
            recovered_history = obs.recovered_history,
            occupancy_break_days = obs.occupancy_break_days,
            export_case_days = obs.export_case_days,
            export_death_days = obs.export_death_days,
            breakpoint = obs.who_first_sitrep_days,
            n_patches,
            province_lab_increments = lab.increments,
            province_lab_days = lab.days,
            province_lab_bins = lab.bins,
            tmrca_days = obs.tmrca_days
        )
    end

    labh = obs.province_lab_daily_history
    m = build(labh)
    vi = DynamicPPL.VarInfo(Xoshiro(7), m)
    base = DynamicPPL.logjoint(m, vi)
    @test isfinite(base)

    ## The background split is sampled once the country has patches.
    names = string.(keys(vi))
    @test any(contains("bg_log_ratio"), names)
    @test any(contains("lab_composition_state"), names)

    ## Moving analysed specimens between provinces at a fixed daily total
    ## must move the density: the split is what this term adds.
    shifted = Dict(k => v for (k, v) in labh)
    it = labh["ituri_analysed"]
    nk = labh["nord_kivu_analysed"]
    mv = min.(it.counts, 20)
    shifted["ituri_analysed"] = (; days = it.days, counts = it.counts .- mv)
    shifted["nord_kivu_analysed"] = (; days = nk.days, counts = nk.counts .+ mv)
    @test !isapprox(DynamicPPL.logjoint(build(shifted), vi), base; rtol = 1.0e-8)

    ## Without the province series the term is absent and nothing else moves.
    none = Dict{String, @NamedTuple{days::Vector{Int}, counts::Vector{Int}}}()
    m_none = build(none)
    @test isfinite(DynamicPPL.logjoint(m_none, vi))
    @test !any(contains("lab_composition_state"), string.(keys(DynamicPPL.VarInfo(Xoshiro(7), m_none))))

    ## One patch has no split to sample.
    single = build(none; n_patches = 1)
    @test !any(contains("bg_log_ratio"), string.(keys(DynamicPPL.VarInfo(Xoshiro(7), single))))

    ## Province laboratory data with one patch would be dropped silently.
    @test_throws ErrorException build(labh; n_patches = 1)()
end

@testitem "_patch_analysed_increments: ascertainment scales the BVD part only" begin
    using BVDOutbreakSize: _patch_analysed_increments

    carried = [10.0 20.0 30.0; 5.0 5.0 5.0]
    bg = [100.0, 100.0, 100.0]
    w = [0.5, 0.5]
    days = [1, 2, 3]
    flat = _patch_analysed_increments(carried, 0.5, [1.0, 1.0], bg, w, days)
    ## Ascertainment re-splits the national BVD volume between the patches
    ## and leaves the column sums, BVD plus background, unchanged.
    up = _patch_analysed_increments(carried, 0.5, [2.0, 1.0], bg, w, days)
    @test vec(sum(up; dims = 1)) ≈ vec(sum(flat; dims = 1))
    @test all(up[1, :] .> flat[1, :])
    @test all(up[2, :] .< flat[2, :])
    ## Day 1: BVD 0.5 * 15 = 7.5 split 20:5 with ascertainment 2:1.
    @test up[1, 1] ≈ 7.5 * 20 / 25 + 50
    @test up[2, 1] ≈ 7.5 * 5 / 25 + 50
    ## Bins sum the printed days.
    binned = _patch_analysed_increments(carried, 0.5, [1.0, 1.0], bg, w, days, [1, 1, 2])
    @test binned[:, 1] ≈ flat[:, 1] .+ flat[:, 2]
    @test binned[:, 2] ≈ flat[:, 3]
end

@testitem "the ascertainment prior carries no testing coefficient" begin
    using BVDOutbreakSize
    using Turing: DynamicPPL
    using Random: Xoshiro

    ## The province laboratory volumes are scored by the lab composition, so
    ## no covariate built from the same series sits on the ascertainment
    ## prior, and the model has no coefficient for one.
    obs = load_observations()
    m = production_joint(obs; breakpoint = default_breakpoint(obs))
    vi = DynamicPPL.VarInfo(Xoshiro(1), m)
    ks = Set(string(k) for k in keys(vi))
    @test !any(k -> occursin("β_asc", k), ks)
    @test !any(k -> occursin("testing_coefficient", k), ks)
end

@testitem "province case-fatality contrasts pool with a scale of 0.1" begin
    using BVDOutbreakSize
    using Turing: Prior, sample
    using Random: Xoshiro
    using Statistics: median

    ## The scale prior is half-normal with sd 0.1, so a typical province sits
    ## within about ten percent of the national case-fatality ratio. Its
    ## median is 0.0674.
    obs = load_observations()
    m = production_joint(obs; breakpoint = default_breakpoint(obs))
    chn = sample(Xoshiro(3), m, Prior(), 400; progress = false)
    τ = vec(Array(chn[:province_cfr_sd]))
    @test 0.05 < median(τ) < 0.09
end
