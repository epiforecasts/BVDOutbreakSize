## A forecast is the fitted model run past its cut-off. `with_horizon`
## lengthens the latent grid, draws the future walk innovations and the
## future observations as new variables, and leaves everything up to the
## cut-off as it was. These items check that the fitted density does not
## move, that `predict` keeps each draw's parameters, and that the future
## counts come out of the model's own delays and ascertainment.

@testsnippet HorizonFixtures begin
    using BVDOutbreakSize
    using BVDOutbreakSize: ForecastHorizon, with_horizon, horizon_days
    using Turing: DynamicPPL, logjoint, fix, sample, Prior, predict,
        returned, @varname
    using Random: Xoshiro
    import FlexiChains

    ## Cut-off day 40 with the walk starting at day 2 (breakpoint 30 less
    ## the 28-day lead), so the last fitted knot segment is three days
    ## long. A horizon of 10 days is not a whole number of weeks either.
    ## Both put a misplaced knot on a fitted day if the extension got the
    ## grid wrong.
    const N = 40
    const H = 10
    const BP = 30

    const DH = (; days = [13, 18, 40], counts = [10, 14, 18])
    const RH = (; days = [13, 18, 40], counts = [340, 516, 905])
    const CH = (; days = [13, 18, 28, 35, 40], counts = [9, 17, 27, 33, 40])
    const LH = (; days = [18, 28], counts = [30, 50])
    const LDH = (; days = [35], counts = [12])
    const CDH = (; days = [18, 28, 40], counts = [1, 3, 5])
    const SDH = (; days = [35, 36, 37], counts = [12, 9, 7])
    const SDDH = (; days = [35, 36], counts = [1, 2])
    const ISO = (; days = [30, 32, 34, 36, 38, 40], counts = [20, 22, 25, 24, 26, 28])
    const CAP = (; days = [30, 40], counts = [60, 60])
    const REC = (; days = [30, 40], counts = [2, 6])
    const ADM = (; days = [36, 38, 40], counts = [4, 5, 3])
    const TDTH = (; days = [36, 38, 40], counts = [1, 0, 1])
    const RO = (; days = [36, 38, 40], counts = [2, 3, 2])
    const OC = (;
        onset_days = [20, 21, 22, 23, 20, 21, 22, 23, 24],
        report_days = [25, 25, 25, 25, 32, 32, 32, 32, 32],
        prev_report_days = [0, 0, 0, 0, 25, 25, 25, 25, 0],
        increments = [2, 3, 1, 0, 1, 2, 3, 4, 5],
    )
    ## Four patches, the shape of the headline fit. The province tables are
    ## partitions of the national confirmed and confirmed-death increments.
    const NP = length(PROVINCE_NAMES)
    const PROV = [4 5 3; 2 1 3; 1 0 1; 1 0 0]
    const PROV_DAYS = [28, 35, 40]
    const PROVD = [1 1; 1 1; 0 0; 0 1]
    const PROVD_DAYS = [28, 40]

    joint_args() = (;
        confirmed_deaths = 5, recovered_cases = 6,
        deaths_history = DH, reported_history = RH, confirmed_history = CH,
        confirmed_deaths_history = CDH, lab_history = LH,
        lab_daily_history = LDH, suspected_daily_history = SDH,
        suspected_daily_deaths_history = SDDH, isolation_history = ISO,
        bed_capacity_history = CAP, recovered_history = REC,
        treatment_admissions_history = ADM,
        treatment_deaths_history = TDTH, treatment_ruleout_history = RO,
        export_case_days = [20, 25], export_death_days = [30],
        onset_curve_history = OC, breakpoint = BP,
        background_pooling = background_pooling_model,
        genetic = genetic_seeding_model, tmrca_days = 60.0,
    )
    patch_args() = (;
        n_patches = NP, province_increments = PROV,
        province_days = PROV_DAYS, province_death_increments = PROVD,
        province_death_days = PROVD_DAYS,
    )
    joint(; kw...) = bvd_joint(N, 2, 18, 905, 1, 40, 50; joint_args()..., kw...)
    patch_joint(; kw...) = joint(; patch_args()..., kw...)

    ## Every composer a forecast is drawn from, at the fixture cut-off.
    composers() = [
        "joint" => joint(),
        "patch joint" => patch_joint(),
        "cases" => cases_only_model(
            N, 905; reported_history = RH,
            suspected_daily_history = SDH, breakpoint = BP
        ),
        "deaths" => deaths_only_model(
            N, 18; deaths_history = DH,
            suspected_daily_deaths_history = SDDH, breakpoint = BP
        ),
        "confirmed" => confirmed_only_model(
            N, 40; confirmed_history = CH, lab_history = LH,
            lab_daily_history = LDH, breakpoint = BP
        ),
        "confirmed deaths" => confirmed_deaths_only_model(
            N, 5, 18; deaths_history = DH,
            confirmed_deaths_history = CDH, breakpoint = BP
        ),
        "treatment" => treatment_only_model(
            N; isolation_history = ISO, bed_capacity_history = CAP,
            treatment_admissions_history = ADM,
            treatment_deaths_history = TDTH,
            treatment_ruleout_history = RO, breakpoint = BP
        ),
        "onsets" => onsets_only_model(
            N; onset_curve_history = OC, breakpoint = BP
        ),
        "exports" => exports_joint_only_model(
            N, 2, 1; export_case_days = [20, 25],
            export_death_days = [30], breakpoint = BP
        ),
    ]

    ## The variables the extended model adds: those a draw from it carries
    ## and a draw from the fitted model does not.
    function future_draws(m0, mh, rng)
        θ0 = rand(rng, m0)
        θh = rand(rng, mh)
        k0 = collect(keys(θ0))
        kh = collect(keys(θh))
        fut = filter(k -> !(k in k0), kh)
        return (; θ0, θh, k0, kh, fut)
    end
    is_future(vn) = occursin(r"future|forecast", string(vn))

    ## A prior draw's parameters with the future variables fixed, so the
    ## extended model's density can be evaluated at the fitted draw.
    fix_future(mh, θh, fut) = fix(mh, Dict(k => θh[k] for k in fut))
end

@testitem "with_horizon rebuilds the model with a forecast horizon" setup = [
    HorizonFixtures,
] begin
    m0 = joint()
    mh = with_horizon(m0, H)
    @test horizon_days(nothing) == 0
    @test horizon_days(ForecastHorizon(H)) == H
    @test mh.defaults.forecast == ForecastHorizon(H)
    @test values(mh.args) == values(m0.args)
    ## Every other keyword is carried over unchanged.
    for k in keys(m0.defaults)
        k === :forecast && continue
        @test getproperty(mh.defaults, k) === getproperty(m0.defaults, k)
    end
    @test with_horizon(mh, 0).defaults.forecast === nothing
end

@testitem "a forecast horizon leaves the fitted density unchanged" setup = [
    HorizonFixtures,
] begin
    for (name, m0) in composers()
        mh = with_horizon(m0, H)
        for seed in 1:3
            d = future_draws(m0, mh, Xoshiro(seed))
            ## The extension adds variables and drops none, and every one it
            ## adds is a future quantity.
            @test all(k -> k in d.kh, d.k0)
            @test !isempty(d.fut)
            @test all(is_future, d.fut)
            ## Fixed variables add nothing to the density, so with the
            ## future ones fixed the extended model's density at a fitted
            ## draw is the fitted model's, to the bit.
            @test logjoint(fix_future(mh, d.θh, d.fut), d.θ0) ==
                logjoint(m0, d.θ0)
        end
    end
end

@testitem "the fitted days of the latent series do not move" setup = [
    HorizonFixtures,
] begin
    m0 = patch_joint()
    mh = with_horizon(m0, H)
    d = future_draws(m0, mh, Xoshiro(4))
    r0 = returned(m0, d.θ0)
    rh = returned(fix_future(mh, d.θh, d.fut), d.θ0)
    ## Every convolution and the renewal are causal, so the first `N` days of
    ## each daily series are the fitted model's.
    @test size(rh.patch_state.infections_matrix, 2) == N + H
    @test rh.patch_state.infections_matrix[:, 1:N] ==
        r0.patch_state.infections_matrix
    @test rh.cases_state.reports_daily[1:N] == r0.cases_state.reports_daily
    @test rh.confirmed_state.confirmed_daily[1:N] ==
        r0.confirmed_state.confirmed_daily
    @test rh.deaths_state.deaths_daily[1:N] == r0.deaths_state.deaths_daily
    @test rh.treatment_state.demand[1:N] == r0.treatment_state.demand
end

@testitem "predict keeps each draw and draws only the future" setup = [
    HorizonFixtures,
] begin
    m0 = patch_joint()
    chn = sample(
        Xoshiro(5), m0, Prior(), 6;
        chain_type = FlexiChains.VNChain, progress = false
    )
    pp = predict(Xoshiro(6), with_horizon(m0, H), chn)
    for key in (:C_T, :R_T, Symbol("rt_state.sigma_rw"))
        @test vec(Array(pp[key])) == vec(Array(chn[key]))
    end
    ## One future count per day for the daily streams.
    for key in (
            "forecast_reports.increments", "forecast_deaths.increments",
            "forecast_confirmed.increments",
            "forecast_confirmed_deaths.increments",
            "forecast_recovered.increments", "forecast_isolation.obs",
            "forecast_admissions.obs", "forecast_incare_deaths.increments",
            "forecast_ruleouts.increments", "forecast_infections",
            "forecast_onsets", "forecast_rt",
        )
        draws = [collect(v) for v in vec(collect(pp[Symbol(key)]))]
        @test length(draws) == 6
        @test all(v -> length(v) == H, draws)
        @test all(v -> all(isfinite, v) && all(>=(0), v), draws)
    end
end

@testitem "forecast counts follow the model's own daily series" setup = [
    HorizonFixtures,
] begin
    m0 = patch_joint()
    mh = with_horizon(m0, H)
    d = future_draws(m0, mh, Xoshiro(7))
    rh = returned(fix_future(mh, d.θh, d.fut), d.θ0)
    fd = (N + 1):(N + H)
    fm = rh.forecast_means
    @test fm.reports == rh.cases_state.reports_daily[fd]
    @test fm.deaths == rh.deaths_state.deaths_daily[fd]
    @test fm.confirmed == rh.confirmed_state.confirmed_daily[fd]
    @test fm.confirmed_deaths ==
        rh.confirmed_deaths_state.confirmed_death_daily[fd]
    @test fm.recovered == rh.recovered_state.recovered_daily[fd]
    @test fm.isolation == rh.treatment_state.occupancy_mean[fd]
end

@testitem "forecast counts carry the fitted delays past the cut-off" setup = [
    HorizonFixtures,
] begin
    m0 = patch_joint()
    mh = with_horizon(m0, H)
    d = future_draws(m0, mh, Xoshiro(8))
    ## Drive the reproduction number to zero from the cut-off: no infection
    ## happens after day `N`. Counts still arrive, from infections before the
    ## cut-off working through the reporting, laboratory and death delays.
    nz = length(d.θh[@varname(rt_state.z_future)])
    fut = Dict(k => d.θh[k] for k in d.fut)
    fut[@varname(rt_state.z_future)] = fill(-5.0e3, nz)
    rh = returned(fix(mh, fut), d.θ0)
    @test all(<(1.0e-6), rh.patch_state.infections_total[(N + 2):(N + H)])
    fm = rh.forecast_means
    @test fm.confirmed[1] > 0
    @test fm.confirmed_deaths[1] > 0
    @test fm.deaths[1] > 0
    ## The first future day's counts come almost wholly from infections
    ## before the cut-off, so stopping transmission barely moves them. Later
    ## days lose the infections that no longer happen.
    base = returned(fix_future(mh, d.θh, d.fut), d.θ0).forecast_means
    @test fm.confirmed[1] ≈ base.confirmed[1] rtol = 0.05
    @test fm.confirmed[end] < base.confirmed[end]
end

@testitem "province forecasts add up to the national forecast" setup = [
    HorizonFixtures,
] begin
    m0 = patch_joint()
    chn = sample(
        Xoshiro(9), m0, Prior(), 8;
        chain_type = FlexiChains.VNChain, progress = false
    )
    pp = predict(Xoshiro(10), with_horizon(m0, 14), chn)
    draws(key) = [collect(v) for v in vec(collect(pp[Symbol(key)]))]
    nat = draws("forecast_confirmed.increments")
    natd = draws("forecast_confirmed_deaths.increments")
    prov = draws("forecast_province_confirmed")
    provd = draws("forecast_province_deaths")
    for i in eachindex(nat)
        ## One column per future week, one row per patch, flattened
        ## column-major.
        p = reshape(prov[i], NP, :)
        pd = reshape(provd[i], NP, :)
        @test size(p, 2) == 2
        @test vec(sum(p; dims = 1)) == [sum(nat[i][1:7]), sum(nat[i][8:14])]
        @test vec(sum(pd; dims = 1)) ==
            [sum(natd[i][1:7]), sum(natd[i][8:14])]
        @test all(>=(0), p)
    end
end

@testitem "the onset forecast is the fitted reporting hazard run forward" setup = [
    HorizonFixtures,
] begin
    using BVDOutbreakSize: onset_report_expected_total
    m0 = joint()
    mh = with_horizon(m0, 14)
    d = future_draws(m0, mh, Xoshiro(11))
    rh = returned(fix_future(mh, d.θh, d.fut), d.θ0)
    st = rh.onset_report_state
    total(as_of) = onset_report_expected_total(
        rh.onsets, st.logit_h0, st.γ, st.grid_start, st.alpha, as_of
    )
    ## One future vintage a week, each the reported total it should print
    ## less the total at the cut-off.
    @test rh.forecast_means.onset_reports ≈
        [total(N + 7) - total(N), total(N + 14) - total(N)]
end

@testitem "production joint keeps its density with a horizon" tags = [
    :slow,
] begin
    using BVDOutbreakSize
    using BVDOutbreakSize: with_horizon
    using Turing: logjoint, fix, sample, Prior
    using Random: Xoshiro
    import FlexiChains

    ## The headline model on the live data, at prior draws: the fixtures
    ## above cannot rule out a data shape only the real observations have.
    obs = load_observations()
    m0 = production_joint(obs; breakpoint = default_breakpoint(obs))
    mh = with_horizon(m0, 28)
    for seed in 1:3
        θ0 = rand(Xoshiro(seed), m0)
        θh = rand(Xoshiro(seed + 10), mh)
        k0 = collect(keys(θ0))
        fut = filter(k -> !(k in k0), collect(keys(θh)))
        @test all(k -> occursin(r"future|forecast", string(k)), fut)
        @test logjoint(fix(mh, Dict(k => θh[k] for k in fut)), θ0) ==
            logjoint(m0, θ0)
    end
end
