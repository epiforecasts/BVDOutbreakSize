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
    ## Occupancy never exceeds the latent demand (the soft cap), and the
    ## supply-limited occupancy never exceeds the bed capacity.
    @test all(iso .<= dem .+ 1e-6)
    @test all(iso .<= cap .+ 1e-6)
end
