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
