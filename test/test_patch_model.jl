## Unit tests for the patch (meta-population) model: the multi-patch renewal
## primitives in src/renewal.jl, the per-patch Rt and infection models, the
## per-province composition likelihood, and the bvd_joint composer.

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

@testitem "patch_rt_model: deviations sum to zero, no privileged patch" begin
    using BVDOutbreakSize: patch_rt_model
    using Random: seed!

    n, np = 60, 3
    m = patch_rt_model(n, np, log(1.5); rt_start = 1)
    for s in 1:5
        seed!(s)
        res = m()
        ## Sum-to-zero at EVERY time: the national walk is the common trend
        ## and the deviations are contrasts around it.
        @test maximum(abs, sum(res.δ_patch; dims = 1)) < 1e-10
        ## Every patch, INCLUDING the primary, carries its own deviation.
        ## Reference coding (delta_1 = 0) also identifies the model, but it
        ## forces Ituri to BE the national trend while the other provinces
        ## carry noise -- an artefact of the identifiability fix, not
        ## epidemiology.
        @test any(!iszero, res.δ_patch[1, :])
        ## Rt is the common trend times the patch deviation.
        for p in 1:np
            @test res.Rt_matrix[p, :] ≈
                  res.Rt_national .* exp.(res.δ_patch[p, :])
        end
        @test all(>(0), res.Rt_matrix)
        @test size(res.δ_patch) == (np, n)
        ## Omega is a genuine correlation matrix: symmetric, unit diagonal,
        ## entries in [-1, 1]. It is reconstructed as L * L', so allow the
        ## round-off that puts the diagonal a few ulps either side of 1.
        for i in 1:np
            @test res.Ω[i, i] ≈ 1
            for j in 1:np
                @test res.Ω[i, j] ≈ res.Ω[j, i]
                @test -1 - 1e-12 <= res.Ω[i, j] <= 1 + 1e-12
            end
        end
    end
end

@testitem "patch_rt_model: Rt may vary across space AND over time" begin
    using BVDOutbreakSize: patch_rt_model
    using Distributions: Normal, truncated
    using Random: seed!

    ## The deviations follow a multivariate-normal random walk whose scale is
    ## itself sampled, so the model NESTS both hypotheses and lets the data
    ## choose. Pin both limits: the point of this parameterisation is that
    ## neither is hard-coded.
    n, np = 120, 3

    ## sigma_delta -> 0: every province collapses onto the common national
    ## trend, so they share one temporal Rt shape at a fixed ratio.
    flat = patch_rt_model(n, np, log(1.5); rt_start = 20, breakpoint = 60.0,
        region_sd_prior = truncated(Normal(0, 1e-12); lower = 0),
        region_drift_sd_prior = truncated(Normal(0, 1e-12); lower = 0))
    seed!(9)
    rf = flat()
    @test maximum(abs, rf.δ_patch) < 1e-8
    for p in 1:np
        @test rf.Rt_matrix[p, :] ≈ rf.Rt_national
    end

    ## A large drift scale lets the provincial trajectories genuinely
    ## separate over time, which a constant modifier could not represent.
    wide = patch_rt_model(n, np, log(1.5); rt_start = 20, breakpoint = 60.0,
        region_drift_sd_prior = truncated(Normal(0.5, 0.01); lower = 0))
    seed!(3)
    rw = wide()
    ## The Ituri / Nord-Kivu contrast must actually move over the window.
    contrast(t) = rw.δ_patch[2, t] - rw.δ_patch[1, t]
    @test abs(contrast(n) - contrast(20)) > 0.2
    ## Sum-to-zero survives however wide the walk.
    @test maximum(abs, sum(rw.δ_patch; dims = 1)) < 1e-10
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

@testitem "patch_infection_model: importation is on by default, off if uncoupled" begin
    using BVDOutbreakSize
    using BVDOutbreakSize: patch_infection_model
    using Turing: DynamicPPL
    using Random: Xoshiro

    n, np = 80, 3
    has_eps(m) = any(
        k -> occursin("ε", string(k)) || occursin("epsilon", string(k)),
        keys(DynamicPPL.VarInfo(Xoshiro(1), m)))

    ## Provinces are COUPLED by default: the outbreak spread from Ituri into
    ## Nord-Kivu and Sud-Kivu, so a meta-population model that cannot move
    ## infections between patches is not describing what happened.
    @test has_eps(patch_infection_model(n, np))

    ## The kernel is a gravity kernel weighted by DESTINATION population, with
    ## a zero diagonal (no self-importation).
    K = province_importation_kernel()
    @test size(K) == (3, 3)
    @test all(iszero, [K[p, p] for p in 1:3])
    @test all(>(0), [K[p, q] for p in 1:3, q in 1:3 if p != q])
    ## A bigger destination province absorbs proportionally more.
    @test K[2, 1] > K[1, 2]          ## Nord-Kivu (6.7M) > Ituri (4.4M)

    ## Passing an all-zero kernel switches the coupling off, and then epsilon
    ## must NOT be sampled: it would multiply zero and be a dimension the
    ## likelihood never touches, whose posterior is exactly its prior.
    @test !has_eps(patch_infection_model(n, np;
        importation_kernel = zeros(np, np)))
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

@testitem "bvd_joint: fits real per-province data without double counting" begin
    using BVDOutbreakSize
    using Turing: DynamicPPL
    using Random: Xoshiro

    obs = load_observations()
    n, np = obs.n, 3

    function build(pch)
        prov = province_increment_matrix(pch, PROVINCE_NAMES, np)
        return bvd_joint(n,
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
            n_patches = 3,
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

@testitem "bvd_joint: carries the bvd_joint headline quantities" tags=[:slow] begin
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
    m = bvd_joint(obs.n,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        n_patches = 3,
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
    for q in (:C_T_patch, :R_T_patch, :delta_patch, :infections_T_patch,
        :region_drift_sd, :log_rt_contrast)
        @test all(v -> length(v) == 3, vec(collect(chn[q])))
    end
    ## The deviations are contrasts around the common trend, so they sum to
    ## zero in every draw.
    @test all(v -> abs(sum(v)) < 1e-10, vec(collect(chn[:delta_patch])))
    ## The contrast is measured against the primary patch, so its own entry
    ## is identically zero and the others are log Rt relative to Ituri.
    @test all(v -> v[1] == 0, vec(collect(chn[:log_rt_contrast])))
    ## The spatial diagnostics the fit is read off.
    @test all(isfinite, vec(Array(chn[:region_sd])))
    @test all(c -> -1 <= c <= 1,
        vec(Array(chn[:region_corr_primary_secondary])))

    ## The aggregate C_T is the sum over patches.
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
    m = bvd_joint(obs.n,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        n_patches = 3,
        province_increments = prov.increments,
        province_days = prov.days,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 100; chain_type = FlexiChains.VNChain,
        progress = false)

    df = patch_summary_table(chn, 3)
    @test df isa DataFrame
    ## Seven quantities per patch, three patches: cumulative infections, Rt,
    ## daily infections, the log-Rt deviation from the national trend, the
    ## log-Rt contrast against the primary patch, and the relative case
    ## ascertainment. The last two must BOTH be present: the case composition
    ## identifies only their product, so reporting a provincial Rt without the
    ## ascertainment beside it invites a case-finding artefact to be read as
    ## epidemiology.
    @test nrow(df) == 21
    quantities = unique(df[!, "Quantity"])
    @test "Relative case ascertainment" in quantities
    @test "log-Rt vs primary patch" in quantities
    @test "Rt deviation drift" in quantities
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

@testitem "patch reporting: one table per province, and an overview" tags=[:slow] begin
    using BVDOutbreakSize
    using Turing: sample, Prior
    import FlexiChains
    using DataFrames: DataFrame, nrow, names

    obs = load_observations()
    prov = province_increment_matrix(obs.province_confirmed_history,
        PROVINCE_NAMES, 3)
    m = bvd_joint(obs.n,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        n_patches = 3,
        province_increments = prov.increments,
        province_days = prov.days,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 100; chain_type = FlexiChains.VNChain,
        progress = false)

    ## The cross-province overview is one ROW per province, not one row per
    ## (province, quantity). This is the whole point of it: the long-format
    ## table is unreadable as a comparison across provinces.
    ov = patch_overview_table(chn, 3)
    @test ov isa DataFrame
    @test nrow(ov) == 3
    @test ov[!, "Province"] == ["Ituri", "Nord-Kivu", "Sud-Kivu"]
    @test "Reproduction number" in names(ov)
    @test "Share of infections (%)" in names(ov)
    ## Shares are computed per draw and must therefore still sum to 100 in the
    ## median only approximately, but every entry must be a parseable
    ## `median (lower–upper)` cell rather than a raw number.
    for v in ov[!, "Reproduction number"]
        @test occursin("(", v) && occursin("–", v)
    end

    ## Selecting one province gives that province's rows only, and drops the
    ## Patch column, which would otherwise repeat one value down every row.
    full = patch_summary_table(chn, 3)
    one = patch_summary_table(chn, 3; patch = "Nord-Kivu")
    @test nrow(one) == nrow(full) / 3
    @test !("Patch" in names(one))
    @test "Quantity" in names(one)
    ## Selecting by index and by label must agree.
    @test patch_summary_table(chn, 3; patch = 2) == one
    ## The selected rows must be the SAME numbers the full table reports for
    ## that province, not a re-summary of a different patch.
    nk = full[full[!, "Patch"] .== "Nord-Kivu", :]
    @test one[!, "Lower 90%"] == nk[!, "Lower 90%"]
    @test one[!, "Upper 90%"] == nk[!, "Upper 90%"]

    @test_throws ErrorException patch_summary_table(chn, 3; patch = "Kinshasa")
    @test_throws ErrorException patch_summary_table(chn, 3; patch = 9)
end

@testitem "reconstruct_patch_rt: matches the chain's own per-patch Rt" tags=[:slow] begin
    using BVDOutbreakSize
    using Turing: sample, Prior
    import FlexiChains
    using Statistics: median

    obs = load_observations()
    prov = province_increment_matrix(obs.province_confirmed_history,
        PROVINCE_NAMES, 3)
    m = bvd_joint(obs.n,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        n_patches = 3,
        province_increments = prov.increments,
        province_days = prov.days,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 40; chain_type = FlexiChains.VNChain,
        progress = false)

    ## The same renewal start and walk start the model derives internally.
    rt_start = clamp(obs.n - round(Int, obs.tmrca_days) + RENEWAL_START_LEAD,
        1, obs.n)
    rt_walk_start = clamp(
        round(Int, obs.who_first_sitrep_days) - RT_WALK_LEAD, rt_start, obs.n)

    rt = reconstruct_patch_rt(chn; n = obs.n,
        breakpoint = obs.who_first_sitrep_days, n_patches = 3,
        rt_start = rt_start, rt_walk_start = rt_walk_start)
    @test length(rt) == 3
    @test all(size(r) == (40, obs.n) for r in rt)

    ## THE test that makes the figure trustworthy: rebuilding the provincial
    ## trajectory from the deviation knots must reproduce, at the cut-off, the
    ## `R_T_patch` the model itself computed. Without this the panels could
    ## drift from the tables and nothing would catch it.
    rtp = [collect(v) for v in vec(collect(chn[:R_T_patch]))]
    for p in 1:3, i in 1:length(rtp)

        @test rt[p][i, obs.n] ≈ rtp[i][p] rtol=1e-8
    end

    ## The deviations sum to zero, so the incidence-unweighted geometric mean
    ## of the provincial Rt is the national trajectory.
    nat = reconstruct_rt(chn; n = obs.n,
        breakpoint = obs.who_first_sitrep_days,
        rt_start = rt_start, rt_walk_start = rt_walk_start)
    for i in 1:5, d in (obs.n, obs.n - 7)

        gm = exp(sum(log(rt[p][i, d]) for p in 1:3) / 3)
        @test gm ≈ nat[i, d] rtol=1e-8
    end

    ## A chain with no patch structure carries no deviation knots, so the
    ## provincial trajectories cannot be rebuilt. That must be an error rather
    ## than three copies of the national trajectory.
    single = bvd_joint(obs.n,
        obs.exported_cases, obs.total_deaths, obs.reported_cases,
        obs.exports_deaths, obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)
    chn1 = sample(single, Prior(), 5; chain_type = FlexiChains.VNChain,
        progress = false)
    @test_throws ErrorException reconstruct_patch_rt(chn1; n = obs.n,
        breakpoint = obs.who_first_sitrep_days, n_patches = 3,
        rt_start = rt_start, rt_walk_start = rt_walk_start)
end

@testitem "province lab data: an exact partition of the national analysed" begin
    using BVDOutbreakSize

    ## The per-province laboratory throughput (sitrep section 4.3) sums, on
    ## every date, to the national `tests_analysed_daily_history` the model
    ## already fits. That is what makes it a COMPOSITION rather than a new
    ## count stream, exactly like the per-province confirmed cases. If this
    ## ever fails, either the scan drifted or the national series moved, and
    ## the two must be reconciled before the data is used.
    obs = load_observations()
    lab = obs.province_lab_daily_history
    @test !isempty(lab)
    for k in ("ituri_analysed", "ituri_positive", "nord_kivu_analysed",
        "nord_kivu_positive", "sud_kivu_analysed", "sud_kivu_positive")
        @test haskey(lab, k)
    end

    nat = Dict(zip(obs.lab_daily_history.days, obs.lab_daily_history.counts))
    an = [lab["$(p)_analysed"] for p in ("ituri", "nord_kivu", "sud_kivu")]
    days = an[1].days
    @test all(a -> a.days == days, an)
    for (i, d) in enumerate(days)
        haskey(nat, d) || continue
        @test sum(a.counts[i] for a in an) == nat[d]
    end

    ## Positives never exceed the samples analysed in the same province.
    for p in ("ituri", "nord_kivu", "sud_kivu")
        @test all(lab["$(p)_positive"].counts .<= lab["$(p)_analysed"].counts)
    end

    ## The finding that motivates using this data at all: the provinces test
    ## very differently-selected pools, so confirmed-case share is NOT
    ## infection share. Ituri runs ~32% positivity against Nord-Kivu's ~6%.
    it_pos = sum(lab["ituri_positive"].counts) /
             sum(lab["ituri_analysed"].counts)
    nk_pos = sum(lab["nord_kivu_positive"].counts) /
             sum(lab["nord_kivu_analysed"].counts)
    @test it_pos > 3 * nk_pos
end

@testitem "patch_rt_model: the drift prior permits real divergence" tags=[:slow] begin
    using BVDOutbreakSize
    using BVDOutbreakSize: patch_rt_model, knot_days
    using Distributions: Normal, truncated
    using Statistics: median, mean
    using Random: seed!

    ## The headline spatial claim is read off `region_drift_sd`: a posterior
    ## near zero is reported as "the provinces share one temporal Rt shape".
    ## That claim is only meaningful if the PRIOR would have allowed them to
    ## diverge. If the prior were tight, a near-zero posterior would be an
    ## artefact of the prior, not a finding about the outbreak.
    ##
    ## So pin the prior predictive: under the default drift prior, the change
    ## in the Ituri / Nord-Kivu Rt RATIO across the fitted window must have a
    ## substantial chance of exceeding 25%.
    obs = load_observations()
    n = obs.n
    bp = obs.who_first_sitrep_days
    rt_start = clamp(n - round(Int, obs.tmrca_days) + RENEWAL_START_LEAD, 1, n)
    rt_walk_start = clamp(round(Int, bp) - RT_WALK_LEAD, rt_start, n)

    m = patch_rt_model(n, 3, log(1.5); breakpoint = bp, rt_start, rt_walk_start)
    seed!(1)
    ratios = Float64[]
    for _ in 1:600
        r = m()
        c0 = r.δ_patch[2, rt_walk_start] - r.δ_patch[1, rt_walk_start]
        c1 = r.δ_patch[2, n] - r.δ_patch[1, n]
        push!(ratios, exp(abs(c1 - c0)))
    end

    ## Not a straitjacket: a real chance of >25% divergence in the Rt ratio.
    @test mean(ratios .> 1.25) > 0.15
    ## And not vacuous either: the prior is still centred near "no divergence",
    ## so it shrinks toward a shared shape rather than assuming divergence.
    @test median(ratios) < 1.4
end

@testitem "patch_infection_model: the seed prior can reach the observed split" begin
    using BVDOutbreakSize
    using BVDOutbreakSize: patch_infections, seed_infections,
                           seed_at_renewal_start, generation_interval_model,
                           r_to_R0, cdf_nmax
    using Distributions: Gamma, LogNormal, quantile
    using Random: seed!

    ## With importation off (the default), a secondary patch has only two
    ## routes to infections: its own seed and its own Rt. So the seed prior
    ## MUST be able to reach the observed provincial split on its own, at
    ## zero Rt difference. If it cannot, the log-Rt deviation is forced to
    ## absorb the level difference and the reported provincial Rt gap becomes
    ## an artefact of the seed prior -- which is the one thing the patch model
    ## exists to estimate.
    ##
    ## An earlier version used an ABSOLUTE seed prior N+(0.01, 0.01) against
    ## the primary patch's 2^m ~ 164: a seed ratio of ~13,700:1, reaching only
    ## 0.007% of infections where the data want ~9%. Reaching 9% then needed a
    ## 3.4-sigma draw on the deviation prior. The seed is now a FRACTION of the
    ## primary seed, so the level is explained by the seed and the deviation is
    ## identified by the time trend.
    obs = load_observations()
    n = obs.n
    rt_start = clamp(n - round(Int, obs.tmrca_days) + RENEWAL_START_LEAD, 1, n)

    ## The observed Nord-Kivu share of confirmed cases.
    prov = province_increment_matrix(obs.province_confirmed_history,
        PROVINCE_NAMES, 3)
    it = sum(prov.increments[1, :])
    nk = sum(prov.increments[2, :])
    obs_share = nk / (it + nk)
    @test 0.05 < obs_share < 0.15          ## ~9%, guards the fixture

    seed!(1)
    g = generation_interval_model(cdf_nmax(Gamma(2.71, 5.65)))().g
    r = 0.0593
    Rtv = fill(r_to_R0(r, g), n)
    seed0 = seed_at_renewal_start(2.0^m_prior_centre(obs.cutoff))

    ## Share of infections in patch 2 at ZERO Rt difference, for a given seed
    ## fraction of the primary patch's seed.
    function share(frac)
        Rtm = [Rtv'; Rtv'; Rtv']
        seeds = zeros(3, rt_start)
        seeds[1, :] = seed_infections(seed0, r, rt_start)
        seeds[2, :] = seed_infections(frac * seed0, r, rt_start)
        I = patch_infections(Rtm, g, seeds, zeros(3, 3), 0.0)
        a, b = sum(I[1, :]), sum(I[2, :])
        return b / (a + b)
    end

    ## The default prior on the seed fraction must BRACKET the observed share
    ## at zero Rt difference: the 5th percentile below it, the 95th above.
    prior = LogNormal(log(0.05), 1.0)
    lo, hi = quantile(prior, 0.05), quantile(prior, 0.95)
    @test share(lo) < obs_share
    @test share(hi) > obs_share

    ## And the prior median must be within an order of magnitude of what the
    ## data need, so the seed is not fighting the likelihood from the start.
    @test 0.1 * obs_share < share(quantile(prior, 0.5)) < 10 * obs_share
end

@testitem "province deaths: the split that identifies ascertainment" begin
    using BVDOutbreakSize

    ## The per-province CASE count is the product of ascertainment and
    ## incidence, and only the product is observed, so the case split alone
    ## can never separate "more infections" from "better case-finding".
    ##
    ## The DEATH split breaks that. Deaths are far harder to miss than cases,
    ## and the case-fatality ratio and death-confirmation probability belong to
    ## the virus and to a national laboratory, not to a province -- so they are
    ## common factors and cancel out of the normalised death shares, leaving
    ## the provincial INCIDENCE split free of case-ascertainment.
    obs = load_observations()
    cases = province_increment_matrix(obs.province_confirmed_history,
        PROVINCE_NAMES, 3)
    deaths = province_increment_matrix(obs.province_death_history,
        PROVINCE_NAMES, 3)

    ## Both are exact partitions of their national totals (the scanner gates
    ## on this), and they share a vintage grid so the two compositions line up.
    @test deaths.days == cases.days
    @test size(deaths.increments) == size(cases.increments)
    natc = Dict(zip(obs.confirmed_history.days, obs.confirmed_history.counts))
    natd = Dict(zip(obs.confirmed_deaths_history.days,
        obs.confirmed_deaths_history.counts))
    for (i, d) in enumerate(cases.days)
        haskey(natc, d) && @test sum(cases.increments[:, i]) ==
              (i == 1 ? natc[d] : natc[d] - natc[cases.days[i - 1]])
    end

    ## THE identifying signal: Nord-Kivu holds a far larger share of the
    ## confirmed DEATHS than of the confirmed CASES, at every vintage. If this
    ## gap ever vanishes the deaths stop identifying ascertainment, and the
    ## per-patch split silently falls back to being prior-driven -- so pin it.
    ct = sum(cases.increments; dims = 2)
    dt = sum(deaths.increments; dims = 2)
    nk_case_share = ct[2] / sum(ct)
    nk_death_share = dt[2] / sum(dt)
    @test 0.05 < nk_case_share < 0.12          ## ~9%
    @test 0.12 < nk_death_share < 0.22         ## ~14%
    @test nk_death_share > 1.4 * nk_case_share ## the gap that does the work
end

@testitem "bvd_joint is a drop-in for bvd_joint (headline model)" tags=[:slow] begin
    using BVDOutbreakSize
    using Turing: sample, Prior
    import FlexiChains

    ## The patch model is the HEADLINE joint, not a side analysis. So a patch
    ## chain must carry every quantity a single-patch chain does: analysis.jl,
    ## the forecast machinery and the plots all key off these names, and a
    ## missing one is a silent failure at report-render time rather than a
    ## test failure here.
    ##
    ## `forecast_reported` in particular reads a long list of `expected_*_T`
    ## deterministics off the chain; if any is absent the one-week-ahead
    ## forecast cannot be produced from a patch fit at all.
    obs = load_observations()
    prov = province_increment_matrix(obs.province_confirmed_history,
        PROVINCE_NAMES, 3)
    provd = province_increment_matrix(obs.province_death_history,
        PROVINCE_NAMES, 3)

    m = bvd_joint(obs.n,
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
        n_patches = 3,
        province_increments = prov.increments, province_days = prov.days,
        province_death_increments = provd.increments,
        province_death_days = provd.days,
        tmrca_days = obs.tmrca_days)

    chn = sample(m, Prior(), 40; chain_type = FlexiChains.VNChain,
        progress = false)

    ## THE CHAIN READ CONTRACT. Every key that src/ or docs/ pulls off a joint
    ## chain. Two CI render failures in a row came from keys the model tests
    ## could not see:
    ##
    ##   - `cumulative_infections` was surfaced by the OLD composer's `_latent`
    ##     submodel, not by its body, so a diff of the two function bodies
    ##     showed every other deterministic as missing but not that one.
    ##   - the `rt_state.*` keys existed but under BARE names, because the
    ##     patch model attached the Rt walk unprefixed. The same parameters are
    ##     sampled either way, so a parameter COUNT check passes and only a
    ##     render fails.
    ##
    ## Both are invisible to a test that asserts what the author remembers.
    ## Assert against what is actually READ.
    for q in (:bed_capacity, :CFR, :C_T, :cumulative_confirmed,
        :cumulative_infections, :cumulative_onsets,
        :cumulative_expected_deaths, :doubling_time, :expected_admissions_T,
        :expected_bed_demand_T, :expected_confirmed_deaths_T,
        :expected_confirmed_T, :expected_deaths_T, :expected_incare_deaths_T,
        :expected_infections_T, :expected_recovered_T, :expected_reports_T,
        :expected_ruleouts_T, :isolation_dispersion, :k,
        :onset_to_confirmation_pmf, :onset_to_death_confirmation_pmf,
        :p_drc, :p_uganda, :r, :r0, :R0, :recovered_dispersion, :R_T, :T,
        :lambda_bg, :lambda_bg_death, :tau_death, :tau_test)
        @test chn[q] !== nothing
    end

    ## Submodel-PREFIXED keys. These are the ones that bite, because a prefix
    ## change leaves the parameter set identical and only breaks the read.
    for q in ("rt_state.sigma_rw", "rt_state.log_R0", "rt_state.z",
        "rt_state.intervention_effect", "gi_state.α", "gi_state.θ",
        "inc_state.delay_mean", "inc_state.delay_sd",
        "cases_state.report_state.α", "cases_state.report_state.θ",
        "confirmed_state.receipt_state.d.delay_mean",
        "confirmed_state.receipt_state.d.delay_sd",
        "deaths_state.od_state.oa.α", "deaths_state.od_state.oa.θ",
        "deaths_state.od_state.ad.α", "deaths_state.od_state.ad.θ",
        "exports_state.detect_state.α", "exports_state.detect_state.θ",
        "exports_state.travel_state.daily_travellers")
        @test chn[Symbol(q)] !== nothing
    end

    ## The scalars must be finite; the trajectories must be non-empty.
    for q in (:C_T, :R_T, :r, :r0, :R0, :doubling_time, :CFR, :p_drc,
        :p_uganda, :k)
        @test all(isfinite, vec(Array(chn[q])))
    end
    for q in (:cumulative_infections, :cumulative_onsets,
        :cumulative_confirmed, :onset_to_confirmation_pmf)
        @test !isempty(vec(collect(chn[q]))[1])
    end
end

@testitem "bvd_joint: n_patches = 1 turns the spatial structure off cleanly" begin
    using BVDOutbreakSize
    using Turing: DynamicPPL
    using Random: Xoshiro

    ## There is ONE model, not two. The patch machinery is a strict
    ## generalisation: with a single patch the sum-to-zero deviations vanish,
    ## the importation kernel has nothing to couple, and there are no
    ## composition terms -- so the model collapses exactly onto the
    ## single-population one. That is what lets the old bvd_joint be retired
    ## rather than maintained alongside.
    ##
    ## The failure this guards against is a QUIET one: if the deviation
    ## machinery still sampled with one patch, it would add prior-only
    ## dimensions the likelihood never touches, and nothing would look wrong.
    obs = load_observations()
    m1 = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        breakpoint = obs.who_first_sitrep_days,
        tmrca_days = obs.tmrca_days)

    ks = Set(string(k) for k in keys(DynamicPPL.VarInfo(Xoshiro(1), m1)))
    has(s) = any(k -> occursin(s, k), ks)

    ## None of the patch machinery may be sampled with a single patch.
    @test !has("ε") && !has("epsilon")     ## nothing to import between
    @test !has("seed_fraction")            ## no secondary patch to seed
    @test !has("σ_δ")                      ## deviations are identically zero
    @test !has("σ_level")
    @test !has("Ω_L")                      ## no cross-patch correlation
    @test !has("z_drift")
    @test !has("z_level")
    @test !has("composition")              ## no per-province likelihood
end

@testitem "bvd_joint: province data with n_patches = 1 is an error" begin
    using BVDOutbreakSize
    using Turing: DynamicPPL
    using Random: Xoshiro

    ## The silent failure mode: per-province data supplied but n_patches left
    ## at its default of 1. The compositions would be scored against a single
    ## patch holding the entire national total, the spatial structure would
    ## quietly vanish, and the fit would look perfectly healthy. This exact
    ## mistake was made once already, by a rename that dropped the argument.
    obs = load_observations()
    prov = province_increment_matrix(obs.province_confirmed_history,
        PROVINCE_NAMES, 3)

    ## The guard lives in the model body, so it fires on EVALUATION, not on
    ## construction -- which is the right place: it is the fit that would be
    ## silently wrong.
    bad = bvd_joint(obs.n, obs.exported_cases,
        obs.total_deaths, obs.reported_cases, obs.exports_deaths,
        obs.confirmed_cases, obs.tests_analysed;
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        deaths_history = obs.deaths_history,
        breakpoint = obs.who_first_sitrep_days,
        province_increments = prov.increments,
        province_days = prov.days,
        tmrca_days = obs.tmrca_days)
    @test_throws ErrorException DynamicPPL.VarInfo(Xoshiro(1), bad)
end
