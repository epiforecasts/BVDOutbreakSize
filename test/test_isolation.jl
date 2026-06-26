## Tests for the isolation/treatment-bed occupancy stream ("Patients en
## isolement"), a prevalence (length-of-stay) stream: the suspect inflow
## (BVD treatment stay plus non-BVD rule-out stay) carried through a
## length-of-stay survival into a daily stock, scored against the modelled
## occupancy on each report day. Exercised through the `convolve_survival`
## helper, `treatment_only_model` and `bvd_joint`.

@testitem "convolve_survival: same-day discharge returns the inflow" begin
    using BVDOutbreakSize: convolve_survival
    x = [1.0, 2.0, 3.0, 4.0]
    ## A length-of-stay point mass at 0 (`pmf = [1]`) gives survival
    ## `S(0) = 1`, `S(τ>0) = 0`, so the admission day is the only occupancy
    ## day and the occupancy equals the inflow.
    @test convolve_survival(x, [1.0]) == x
end

@testitem "convolve_survival: fixed stay accumulates the right occupancy" begin
    using BVDOutbreakSize: convolve_survival
    ## A length-of-stay fixed at 2 days (`pmf = [0, 0, 1]`) means a patient
    ## occupies a bed on the admission day and the next two days, so the
    ## survival weights are `S(0)=S(1)=S(2)=1`, `S(τ≥3)=0`. With one
    ## admission per day the occupancy ramps 1, 2, 3 then holds at 3.
    los = [0.0, 0.0, 1.0]
    x = ones(5)
    occ = convolve_survival(x, los)
    @test occ == [1.0, 2.0, 3.0, 3.0, 3.0]
    ## Total occupancy equals the total inflow times `E[LOS] + 1` (here 3).
    @test sum(convolve_survival([0.0, 0.0, 1.0, 0.0, 0.0], los)) ≈ 3.0
end

@testitem "convolve_survival: survival weights are non-increasing" begin
    using BVDOutbreakSize: convolve_survival, discretise_censored,
                           lognormal_meansd
    ## A single unit admission on day 1 traces the survival curve directly:
    ## occupancy[t] = S(t-1), which must be non-increasing and start at 1.
    los = discretise_censored(lognormal_meansd(6.0, 4.0), 30)
    x = zeros(40)
    x[1] = 1.0
    occ = convolve_survival(x, los)
    @test occ[1] ≈ 1.0
    @test all(diff(occ) .<= 1e-10)
    @test all(occ .>= -1e-12)
end

@testitem "censoring_cap: fixed bound from recorded capacity, never below obs" begin
    using BVDOutbreakSize: censoring_cap

    iso_days = [9, 10, 11]
    iso_obs = [260, 262, 315]
    capacity_history = (; days = [9, 10, 11], counts = [400, 410, 452])
    ## Each day takes its recorded capacity (all above the observed count), so
    ## the bound is the data, fixed, and never below the occupancy.
    cap = censoring_cap(iso_days, iso_obs, capacity_history)
    @test cap == [400.0, 410.0, 452.0]
    @test all(cap .>= iso_obs)

    ## A capacity below the count is floored at the count, so the censored NB
    ## never has zero probability on the observation (the -Inf wall is removed).
    low_cap = (; days = [9], counts = [200])
    @test censoring_cap([9], [260], low_cap) == [260.0]

    ## No recorded capacity gives a large finite no-op bound (not Inf, which
    ## safe_rate would map to eps), so the censoring is inactive.
    nocap = censoring_cap([9], [260], (; days = Int[], counts = Int[]))
    @test only(nocap) > 1.0e5

    ## The predictive generator passes the capacity history days-only (counts
    ## emptied) and the occupancy missing; the guard is on the counts, so this
    ## must not index the empty counts vector and returns the no-op cap.
    daysonly = (; days = [9, 10, 11], counts = Int[])
    @test only(censoring_cap([10], missing, daysonly)) > 1.0e5

    ## Nearest-day fill for a report day without its own capacity entry.
    @test censoring_cap([12], [300], capacity_history) == [452.0]
end

@testitem "bed_capacity_walk: positive capacity path over the grid" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bed_capacity_walk_model

    ## The walk returns a positive bed-capacity path; with a tight innovation
    ## SD it stays a gentle drift around the baseline rather than blowing up.
    chn = sample(bed_capacity_walk_model(30), Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("C0", k), ks)
    C0 = vec(Array(chn[:C0]))
    @test all(C0 .> 0)
end

@testitem "isolation occupancy: conditioned fit stays positive" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## A daily occupancy stock on later days, supplied as observed counts.
    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    chn = sample(
        treatment_only_model(33; isolation_history),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 100
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "isolation occupancy: predictive path samples the counts" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## Days but no counts: the occupancy is a predictive generator, so its
    ## per-day counts are sampled under the `treatment_state.isolation`
    ## submodel rather than conditioned.
    isolation_history = (; days = [28, 29, 30, 31, 32, 33], counts = Int[])
    chn = sample(
        treatment_only_model(33; isolation_history),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("isolation", k), ks)
end

@testitem "isolation occupancy: empty history is a no-op" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## With no isolation history the occupancy submodel scores nothing and
    ## adds no sampled occupancy keys.
    chn = sample(
        treatment_only_model(33),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ks = string.(collect(keys(chn)))
    @test !any(k -> occursin("isolation.increments", k), ks)
    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "accumulate_occupancy: balance closes and stays non-negative" begin
    using BVDOutbreakSize: accumulate_occupancy
    ## A constant admission inflow with simple discharge series builds a
    ## non-negative running balance; the confirmed sub-stock is a subset of the
    ## occupied BVD stock and the suspect sub-stock is the complement, so the two
    ## sum to the total demand on every day.
    n = 30
    A_bvd = fill(4.0, n)
    A_bg = fill(6.0, n)
    death_pmf = [0.0, 0.1, 0.2, 0.3, 0.2, 0.2]
    recover_pmf = vcat(zeros(6), [0.2, 0.3, 0.3, 0.2])
    ruleout_pmf = [0.0, 0.3, 0.4, 0.3]
    using BVDOutbreakSize: convolve_delay
    CFR_iso = 0.35
    deaths = convolve_delay(CFR_iso .* A_bvd, death_pmf)
    recover = convolve_delay((1 - CFR_iso) .* A_bvd, recover_pmf)
    ruleout = convolve_delay(A_bg, ruleout_pmf)
    conf_hazard = fill(0.55, n)
    acc = accumulate_occupancy(A_bvd, A_bg, deaths, recover, ruleout, 0.01,
        conf_hazard)
    @test all(acc.demand .>= -1e-9)
    @test all(acc.O_conf .>= -1e-9)
    @test all(acc.O_susp .>= -1e-9)
    ## Confirmed is a subset of the occupied BVD stock.
    @test all(acc.O_conf .<= acc.O_bvd .+ 1e-9)
    ## Suspect + confirmed = total demand.
    @test all(abs.(acc.O_conf .+ acc.O_susp .- acc.demand) .< 1e-6)
    ## A zero confirmation hazard leaves the confirmed sub-stock empty.
    acc0 = accumulate_occupancy(A_bvd, A_bg, deaths, recover, ruleout, 0.01,
        zeros(n))
    @test all(acc0.O_conf .== 0)
    @test all(acc0.O_susp .≈ acc0.demand)
end

@testitem "admission_headroom: fixed bound above obs, never on the boundary" begin
    using BVDOutbreakSize: admission_headroom
    ## Capacity 400 with previous-day occupancy 260 leaves 140 free beds; a
    ## slack admission count is censored at the headroom, a saturated count is
    ## kept strictly inside the support (obs + 0.5), so it never sits on the
    ## non-differentiable NegBin-CDF boundary.
    capacity_history = (; days = [9, 10], counts = [400, 410])
    isolation_history = (; days = [9, 10], counts = [260, 262])
    adm_days = [10, 11]
    adm_obs = [50, 300]
    head = admission_headroom(adm_days, adm_obs, capacity_history,
        isolation_history)
    @test head[1] ≈ 410 - 260            # day 10: capacity 410 less prev occ 260
    @test head[2] > adm_obs[2]           # day 11: saturated, strictly above obs
    @test all(head .> adm_obs .- 1e-9)
    ## No capacity record gives a large no-op headroom.
    nocap = admission_headroom([10], [50], (; days = Int[], counts = Int[]),
        isolation_history)
    @test only(nocap) > 1.0e5
end

@testitem "occupancy split: sub-stock parameters sampled, fit stays positive" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## Occupancy on days 28-33 with a published split on the last three days.
    ## On split days the two sub-stock censuses are scored instead of the total
    ## (a per-day total-OR-split switch); the abscond fraction and the two
    ## reporting breaks (suspect↔confirmed reclassification and the overnight
    ## total offset) are sampled and stay in range.
    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    confirmed_incare = (; days = [31, 32, 33], counts = [120, 130, 140])
    suspect_incare = (; days = [31, 32, 33], counts = [147, 153, 169])
    chn = sample(
        treatment_only_model(33; isolation_history,
            treatment_confirmed_incare_history = confirmed_incare,
            treatment_suspect_incare_history = suspect_incare,
            treatment_reclass_break_days = [33]),
        Prior(), 60;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("reclass_break", k), ks)
    @test any(k -> occursin("total_break", k), ks)
    @test any(k -> occursin("abscond_frac", k), ks)
    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
    ## The abscond fraction stays a small non-negative loss-to-follow-up rate
    ## under the prior. Index by the chain key object itself (FlexiChains
    ## resolves the submodel-prefixed varname), not a reconstructed Symbol.
    ab_key = first(k for k in keys(chn) if occursin("abscond_frac", string(k)))
    ab = vec(Array(chn[ab_key]))
    @test all(ab .>= 0)
end

@testitem "occupancy split: predictive path samples the sub-stock censuses" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## Days but no counts on the split histories: the two sub-stock censuses are
    ## predictive generators, so their per-day increments are sampled under the
    ## `confirmed_incare_obs` / `suspect_incare_obs` submodels.
    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    chn = sample(
        treatment_only_model(33; isolation_history,
            treatment_confirmed_incare_history = (; days = [31, 32, 33],
                counts = Int[]),
            treatment_suspect_incare_history = (; days = [31, 32, 33],
                counts = Int[])),
        Prior(), 40;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    @test any(k -> occursin("confirmed_incare_obs.increments", k), ks)
    @test any(k -> occursin("suspect_incare_obs.increments", k), ks)
end

@testitem "occupancy split: empty split history is a no-op" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: treatment_only_model

    ## With no split history the two sub-stock census likelihoods score nothing
    ## (no days → no scored or sampled increments), so no split-observation keys
    ## appear, while the total-occupancy backbone still runs.
    isolation_history = (; days = [28, 29, 30, 31, 32, 33],
        counts = [206, 233, 258, 267, 283, 309])
    chn = sample(
        treatment_only_model(33; isolation_history),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)
    ks = string.(collect(keys(chn)))
    @test !any(k -> occursin("confirmed_incare_obs.increments", k), ks)
    @test !any(k -> occursin("suspect_incare_obs.increments", k), ks)
    C_T = vec(Array(chn[:C_T]))
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "isolation occupancy: joint prior runs with the live data" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: load_observations, bvd_joint, genetic_seeding_model

    obs = load_observations()
    @test !isempty(obs.isolation_history.counts)
    breakpoint = obs.n - obs.who_first_sitrep_days
    m = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
        obs.tests_analysed;
        confirmed_deaths = obs.confirmed_deaths,
        deaths_history = obs.deaths_history,
        reported_history = obs.reported_history,
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        suspected_daily_history = obs.suspected_daily_history,
        isolation_history = obs.isolation_history,
        bed_capacity_history = obs.bed_capacity_history,
        treatment_confirmed_incare_history =
        obs.treatment_confirmed_incare_history,
        treatment_suspect_incare_history =
        obs.treatment_suspect_incare_history,
        treatment_reclass_break_days = obs.treatment_reclass_break_days,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = breakpoint,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days)
    chn = sample(m, Prior(), 30;
        chain_type = FlexiChains.VNChain, progress = false)
    C_T = vec(Array(chn[:C_T]))
    iso = vec(Array(chn[:expected_isolation_T]))
    dem = vec(Array(chn[:expected_bed_demand_T]))
    cap = vec(Array(chn[:bed_capacity]))
    @test length(C_T) == 30
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
    @test all(isfinite, iso)
    @test all(iso .> 0)
    ## Occupancy is `min(demand, C)`, so it never exceeds the latent demand,
    ## and the supply-limited occupancy never exceeds the bed capacity.
    @test all(iso .<= dem .+ 1e-6)
    @test all(iso .<= cap .+ 1e-6)
    ## The severity skew is non-negative and admits BVD suspects at least as
    ## readily as the base (non-BVD rule-out) rate.
    skew = vec(Array(chn[:isolation_severity]))
    p_base = vec(Array(chn[:isolation_admission]))
    p_bvd = vec(Array(chn[:isolation_bvd_admission]))
    @test all(skew .>= -1e-9)
    @test all(p_bvd .>= p_base .- 1e-9)
    @test all(0 .<= p_bvd .<= 1)
end
