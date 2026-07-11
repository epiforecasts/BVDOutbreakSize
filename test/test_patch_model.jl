## Unit tests for the patch (meta-population) model: the multi-patch renewal
## primitives in src/renewal.jl, the per-patch Rt and infection models, the
## per-province composition likelihood, and the bvd_patch_joint composer.

@testitem "patch_infections: uncoupled patches match single-patch renewal" begin
    using BVDOutbreakSize: patch_infections, renewal_infections

    ## With a zero importation kernel each patch is an independent renewal
    ## process, so every row must reproduce `renewal_infections` exactly.
    g = [0.2, 0.3, 0.3, 0.2]
    n = 40
    Rt = [fill(1.3, n)'; fill(0.8, n)'; fill(2.0, n)']
    seeds = [1.0 2.0; 0.5 0.6; 0.1 0.2]
    K = zeros(3, 3)

    I = patch_infections(Rt, g, seeds, K, 0.0)
    @test size(I) == (3, n)
    for p in 1:3
        single = renewal_infections(vec(Rt[p, :]), g, vec(seeds[p, :]))
        @test I[p, :] ≈ single
    end
end

@testitem "patch_infections: importation moves infections between patches" begin
    using BVDOutbreakSize: patch_infections

    ## Patch 2 has Rt = 0 and no seed, so it can only ever have infections
    ## that arrive by importation from patch 1.
    g = [0.5, 0.5]
    n = 12
    Rt = [fill(1.5, n)'; fill(0.0, n)']
    seeds = [1.0 1.0; 0.0 0.0]
    K = [0.0 0.0; 0.1 0.0]   ## K[2, 1]: flow from patch 1 into patch 2

    off = patch_infections(Rt, g, seeds, K, 0.0)
    @test all(iszero, off[2, :])           ## epsilon = 0: no importation

    on = patch_infections(Rt, g, seeds, K, 0.5)
    @test all(>(0), on[2, 3:n])            ## epsilon > 0: patch 2 seeded
    ## Patch 1 is upstream and must be untouched by the coupling.
    @test on[1, :] ≈ off[1, :]
    ## The imported amount is exactly epsilon * K[2,1] * I_1[t-1].
    for t in 3:n
        @test on[2, t] ≈ 0.5 * 0.1 * on[1, t - 1]
    end
end

@testitem "importation_from_kernel: matches the explicit sum" begin
    using BVDOutbreakSize: importation_from_kernel

    K = [0.0 0.2 0.3; 0.1 0.0 0.4; 0.5 0.6 0.0]
    I_prev = [10.0, 20.0, 30.0]
    ε = 0.25
    imp = importation_from_kernel(K, I_prev, ε)
    for p in 1:3
        @test imp[p] ≈ ε * sum(K[p, q] * I_prev[q] for q in 1:3)
    end
    ## A zero kernel imports nothing.
    @test all(iszero, importation_from_kernel(zeros(3, 3), I_prev, ε))
end

@testitem "implied_national_Rt: recovers the incidence-weighted patch Rt" begin
    using BVDOutbreakSize: patch_infections, implied_national_Rt

    ## Inverting the renewal equation on the SUMMED infections must return
    ## the force-of-infection-weighted mean of the per-patch Rts. Build two
    ## uncoupled patches with different constant Rts and check the identity
    ## Rt_implied(t) = sum_p Rt_p * force_p(t) / sum_p force_p(t).
    g = [0.3, 0.4, 0.3]
    n = 30
    R1, R2 = 1.8, 0.6
    Rt = [fill(R1, n)'; fill(R2, n)']
    seeds = [5.0 5.0; 1.0 1.0]
    I = patch_infections(Rt, g, seeds, zeros(2, 2), 0.0)
    total = vec(sum(I; dims = 1))
    implied = implied_national_Rt(total, g)

    force(p, t) = sum(I[p, t - s] * g[s] for s in 1:min(t - 1, length(g)))
    for t in 6:n
        want = (R1 * force(1, t) + R2 * force(2, t)) /
               (force(1, t) + force(2, t))
        @test implied[t] ≈ want rtol = 1e-10
    end
    ## The aggregate must sit between the two patch Rts, never outside.
    @test all(R2 - 1e-9 .<= implied[6:n] .<= R1 + 1e-9)
end

@testitem "implied_national_Rt: zero force gives zero, never NaN" begin
    using BVDOutbreakSize: implied_national_Rt

    @test implied_national_Rt(zeros(10), [0.5, 0.5]) == zeros(10)
    @test all(isfinite, implied_national_Rt([0.0, 0.0, 3.0, 4.0], [0.5, 0.5]))
end

@testitem "patch_rt_model: primary patch is the reference (delta_1 = 0)" begin
    using BVDOutbreakSize: patch_rt_model
    using Random: seed!

    n, np = 60, 3
    m = patch_rt_model(n, np, log(1.5); rt_start = 1)
    for s in 1:5
        seed!(s)
        res = m()
        ## The primary patch carries no modifier, so its Rt IS the national
        ## walk; the secondary patches are offset from it by exp(delta_p).
        @test res.δ_patch[1] == 0
        @test res.Rt_matrix[1, :] ≈ res.Rt_national
        for p in 2:np
            offset = exp(res.δ_patch[p])
            @test res.Rt_matrix[p, :] ≈ res.Rt_national .* offset
        end
        ## Rt is strictly positive everywhere: it is exp of the log scale.
        @test all(>(0), res.Rt_matrix)
    end
end

@testitem "province_increment_matrix: differences cumulative province counts" begin
    using BVDOutbreakSize: province_increment_matrix

    hist = Dict(
        "ituri" => (; days = [10, 20, 30], counts = [100, 150, 220]),
        "nord_kivu" => (; days = [10, 20, 30], counts = [10, 12, 20]),
        "sud_kivu" => (; days = [10, 20, 30], counts = [3, 3, 3]))
    names = ["ituri", "nord_kivu", "sud_kivu"]

    got = province_increment_matrix(hist, names, 3)
    @test got.days == [10, 20, 30]
    ## The first increment is the cumulative to the first vintage day, which
    ## is what `bin_increments` produces for the modelled series.
    @test got.increments == [100 50 70; 10 2 8; 3 0 0]

    ## No data -> empty, so the caller skips the composition term.
    @test isempty(province_increment_matrix(
        Dict{String, @NamedTuple{days::Vector{Int}, counts::Vector{Int}}}(),
        names, 3).days)

    ## Provinces on mismatched vintages are an error, not a silent reshape:
    ## the composition allocates each vintage's total across all provinces.
    bad = Dict(
        "ituri" => (; days = [10, 20], counts = [100, 150]),
        "nord_kivu" => (; days = [10, 30], counts = [10, 12]),
        "sud_kivu" => (; days = [10, 20], counts = [3, 3]))
    @test_throws ErrorException province_increment_matrix(bad, names, 3)
end

@testitem "province_composition_model: scores shares, not the total" begin
    using BVDOutbreakSize: province_composition_model
    using Turing: DynamicPPL
    using Random: Xoshiro

    ## The whole point of the composition likelihood is that it adds ONLY the
    ## spatial split and leaves the national total to the national confirmed
    ## stream. So its log-density must be invariant to rescaling the modelled
    ## confirmed level: doubling every patch's modelled cases leaves the
    ## shares unchanged and must leave the log-density unchanged.
    obs = [80 40; 15 8; 5 2]
    modelled = [8.0 4.0; 1.5 0.8; 0.5 0.2]

    logp(mod) = begin
        m = province_composition_model(obs, mod)
        vi = DynamicPPL.VarInfo(Xoshiro(1), m)
        DynamicPPL.loglikelihood(m, vi)
    end

    base = logp(modelled)
    @test logp(2.0 .* modelled) ≈ base rtol = 1e-10
    @test logp(100.0 .* modelled) ≈ base rtol = 1e-10

    ## But it MUST respond to a change in the split.
    skewed = [8.0 4.0; 5.0 3.0; 0.5 0.2]
    @test !isapprox(logp(skewed), base; rtol = 1e-6)

    ## The best-fitting shares are the observed ones: a modelled split equal
    ## to the observed split must beat a badly wrong one.
    matched = Float64.(obs)
    wrong = Float64.(reverse(obs; dims = 1))
    @test logp(matched) > logp(wrong)
end

@testitem "province_composition_model: predictive path fills every patch" begin
    using BVDOutbreakSize: province_composition_model
    using Random: seed!

    modelled = [8.0 4.0; 1.5 0.8; 0.5 0.2]
    seed!(3)
    res = province_composition_model(missing, modelled)()
    inc = res.obs_increments
    @test size(inc) == (3, 2)
    @test !any(ismissing, inc)
    @test all(>=(0), inc)
    ## Each vintage's generated counts must sum to the modelled total, since
    ## the composition is conditioned on that total.
    for i in 1:2
        @test sum(inc[:, i]) == round(Int, sum(modelled[:, i]))
    end
end

@testitem "province_composition_model: does not mutate observed data" begin
    using BVDOutbreakSize: province_composition_model
    using Turing: DynamicPPL
    using Random: Xoshiro

    ## The last patch is the stick-breaking remainder. On the fitting path it
    ## is an observation and must never be written over.
    obs = [80 40; 15 8; 5 2]
    before = copy(obs)
    m = province_composition_model(obs, [8.0 4.0; 1.5 0.8; 0.5 0.2])
    DynamicPPL.VarInfo(Xoshiro(1), m)
    @test obs == before
end

@testitem "patch_infection_model: epsilon only enters with a coupling kernel" begin
    using BVDOutbreakSize: patch_infection_model
    using Turing: DynamicPPL
    using Random: Xoshiro

    ## An all-zero kernel makes the importation term identically zero, so a
    ## sampled epsilon would be a dimension the likelihood never touches.
    ## It must not be in the parameter space at all.
    n, np = 80, 3
    has_eps(m) = any(
        k -> occursin("ε", string(k)) || occursin("epsilon", string(k)),
        keys(DynamicPPL.VarInfo(Xoshiro(1), m)))

    @test !has_eps(patch_infection_model(n, np))
    K = [0.0 0.0 0.0; 1e-4 0.0 0.0; 1e-5 0.0 0.0]
    @test has_eps(patch_infection_model(n, np; importation_kernel = K))
end

@testitem "patch_infection_model: national aggregates are the patch sums" begin
    using BVDOutbreakSize: patch_infection_model
    using Random: seed!

    n, np = 90, 3
    m = patch_infection_model(n, np; rt_start = 10)
    for i in 1:3
        seed!(i)
        s = m()
        @test size(s.infections_matrix) == (np, n)
        @test size(s.onsets_matrix) == (np, n)
        @test s.infections_total ≈ vec(sum(s.infections_matrix; dims = 1))
        @test s.C_T ≈ sum(s.C_T_patch)
        @test s.C_T ≈ s.cumulative_total[n]
        ## Headline quantities present and finite, matching infection_model.
        for q in (s.R0, s.r, s.r0, s.R_T, s.T, s.doubling_time, s.C_T)
            @test isfinite(q)
        end
        @test s.R_T > 0
        ## `r` is derived from `R_T` through Euler-Lotka, so their signs must
        ## agree by construction: r < 0 iff R_T < 1.
        @test (s.r < 0) == (s.R_T < 1)
    end
end

@testitem "bvd_patch_joint: fits real per-province data without double counting" begin
    using BVDOutbreakSize
    using Turing: DynamicPPL
    using Random: Xoshiro

    obs = load_observations()
    n, np = obs.n, 3

    function build(pch)
        prov = province_increment_matrix(pch, PROVINCE_NAMES, np)
        return bvd_patch_joint(n, np,
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
            treatment_admissions_history = obs.treatment_admissions_history,
            treatment_deaths_history = obs.treatment_deaths_history,
            treatment_ruleout_history = obs.treatment_ruleout_history,
            treatment_absconded_history = obs.treatment_absconded_history,
            occupancy_break_days = obs.occupancy_break_days,
            export_case_days = obs.export_case_days,
            export_death_days = obs.export_death_days,
            breakpoint = obs.who_first_sitrep_days,
            province_increments = prov.increments,
            province_days = prov.days,
            tmrca_days = obs.tmrca_days)
    end

    pch = obs.province_confirmed_history
    m = build(pch)
    vi = DynamicPPL.VarInfo(Xoshiro(7), m)
    base = DynamicPPL.logjoint(m, vi)
    @test isfinite(base)

    ## Reallocating cases BETWEEN provinces, at a fixed national total, must
    ## move the log-density: the spatial split is what this model adds.
    shifted = Dict(k => v for (k, v) in pch)
    shifted["ituri"] = (; days = pch["ituri"].days,
        counts = pch["ituri"].counts .- 40)
    shifted["nord_kivu"] = (; days = pch["nord_kivu"].days,
        counts = pch["nord_kivu"].counts .+ 40)
    @test !isapprox(DynamicPPL.logjoint(build(shifted), vi), base; rtol = 1e-8)

    ## Dropping the per-province data must leave the national streams alone:
    ## the composition term is the only thing it can remove.
    none = Dict{String, @NamedTuple{days::Vector{Int}, counts::Vector{Int}}}()
    @test isfinite(DynamicPPL.logjoint(build(none), vi))
end

@testitem "bvd_patch_joint: carries the bvd_joint headline quantities" tags=[:slow] begin
    using BVDOutbreakSize
    using Turing: sample, Prior
    import FlexiChains
    using DataFrames: nrow

    ## A patch chain must drop into the existing summary machinery, which
    ## keys off these names. Missing any of them silently breaks analysis.jl,
    ## so assert on a real chain rather than on the model's return value.
    obs = load_observations()
    prov = province_increment_matrix(obs.province_confirmed_history,
        PROVINCE_NAMES, 3)
    m = bvd_patch_joint(obs.n, 3,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        province_increments = prov.increments,
        province_days = prov.days,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)

    chn = sample(m, Prior(), 50; chain_type = FlexiChains.VNChain,
        progress = false)

    ## The headline quantities analysis.jl summarises, under the same names
    ## bvd_joint uses.
    for q in (:C_T, :R_T, :r, :r0, :T, :CFR, :R0, :doubling_time)
        draws = vec(Array(chn[q]))
        @test length(draws) == 50
        @test all(isfinite, draws)
    end
    ## summary_table must work on a patch chain unchanged.
    @test nrow(summary_table(chn, [:C_T, :R_T, :r, :T, :CFR])) == 5

    ## Per-patch vector deterministics, one entry per patch.
    for q in (:C_T_patch, :R_T_patch, :delta_patch, :infections_T_patch)
        @test all(v -> length(v) == 3, vec(collect(chn[q])))
    end
    ## The primary patch is the reference, so its modifier is exactly zero in
    ## every draw; the aggregate C_T is the sum over patches.
    @test all(v -> v[1] == 0, vec(collect(chn[:delta_patch])))
    C_T = vec(Array(chn[:C_T]))
    per_patch = vec(collect(chn[:C_T_patch]))
    @test all(i -> isapprox(C_T[i], sum(per_patch[i]); rtol = 1e-8),
        eachindex(C_T))
end

@testitem "patch_summary_table: one block per patch, ordered quantiles" tags=[:slow] begin
    using BVDOutbreakSize
    using Turing: sample, Prior, @model
    using Distributions: Normal
    import FlexiChains
    using DataFrames: DataFrame, nrow

    obs = load_observations()
    prov = province_increment_matrix(obs.province_confirmed_history,
        PROVINCE_NAMES, 3)
    m = bvd_patch_joint(obs.n, 3,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        province_increments = prov.increments,
        province_days = prov.days,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 100; chain_type = FlexiChains.VNChain,
        progress = false)

    df = patch_summary_table(chn, 3)
    @test df isa DataFrame
    ## Four quantities per patch, three patches.
    @test nrow(df) == 12
    @test unique(df[!, "Patch"]) == ["Ituri", "Nord-Kivu", "Sud-Kivu"]
    ## The reported quantiles must be ordered, which the previous
    ## implementation's invented "median" (the midpoint of the 20-80
    ## interval) did not guarantee.
    for r in eachrow(df)
        @test r["Lower 90%"] <= r["Lower 60%"] <= r["Lower 30%"]
        @test r["Upper 30%"] <= r["Upper 60%"] <= r["Upper 90%"]
        @test r["Lower 30%"] <= r["Upper 30%"]
    end

    ## A chain without the per-patch deterministics (e.g. a single-patch
    ## bvd_joint chain) must be rejected, not silently summarised.
    @model _no_patches() = x ~ Normal(0.0, 1.0)
    plain = sample(_no_patches(), Prior(), 5;
        chain_type = FlexiChains.VNChain, progress = false)
    @test_throws ErrorException patch_summary_table(plain, 3)
end
