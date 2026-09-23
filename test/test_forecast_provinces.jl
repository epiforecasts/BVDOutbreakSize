## Tests for the per-province forward projection, `forecast_provinces`.

@testsnippet ProvinceProjection begin
    using BVDOutbreakSize: patch_infections, province_importation_kernel,
        PROVINCE_POPULATIONS, cdf_nmax, convolve_delay, discretise_censored,
        lognormal_meansd
    using Distributions: Gamma
    using DataFrames: DataFrame
    using LinearAlgebra: cholesky, I as Id

    ## A three-patch outbreak run for `n + h` days at a constant reproduction
    ## number per patch. The chain carries the first `n` days, so the
    ## projection of days `n + 1 … n + h` has a known answer.
    const NP = 3
    const N = 80
    const H = 7
    const ND = 4
    const ALPHA = 2.71
    const THETA = 5.65
    const G = BVDOutbreakSize._gi_pmf(
        ALPHA, THETA; nmax = cdf_nmax(Gamma(2.71, 5.65))
    )
    const RT = [1.1, 0.9, 1.0]
    const EPS = [0.01, 0.02, 0.005]
    const KERNEL = province_importation_kernel(PROVINCE_POPULATIONS[1:NP])
    ## The model's incubation period and the two onset-to-confirmation
    ## kernels, as the chain carries them.
    const INC = (6.3, 3.5)
    const INC_PMF = discretise_censored(
        lognormal_meansd(INC...), cdf_nmax(lognormal_meansd(6.3, 3.5))
    )
    const KC = discretise_censored(lognormal_meansd(9.0, 4.0), 40)
    const KD = discretise_censored(lognormal_meansd(16.0, 6.0), 50)
    const ASC = [1.4, 0.8, 0.9]
    const ASC_D = [1.1, 0.9, 1.0]
    const SEV = [1.2, 0.7, 1.2]

    function full_run(eps; rt = RT)
        seeds = fill(5.0, NP, 10)
        Rt = repeat(rt, 1, N + H)
        ε = repeat(eps, 1, N + H)
        return patch_infections(Rt, G, seeds, KERNEL, ε).infections
    end

    ## Expected share of each patch in the forecast window, built the way the
    ## fitted composition builds its shares: onsets through the incubation
    ## period, then the onset-to-confirmation kernel, summed over the window
    ## and weighted.
    function window_shares(I, kernel, weight)
        m = [
            sum(
                convolve_delay(
                    convolve_delay(vec(I[p, :]), INC_PMF), kernel
                )[(N + 1):(N + H)]
            ) for p in 1:NP
        ]
        π = weight .* m
        return π ./ sum(π)
    end

    function projection_chain(
            I; delta = zeros(NP), drift_sd = zeros(NP),
            halflife = 42.0, sigma_rw = 0.0, eps = EPS,
            omega = Matrix{Float64}(Id, NP, NP), nd = ND,
            cases = false, deaths = false
        )
        base = (;
            R_T_patch = [copy(RT) for _ in 1:nd],
            delta_patch = [copy(delta) for _ in 1:nd],
            region_drift_sd = [copy(drift_sd) for _ in 1:nd],
            region_halflife = fill(halflife, nd),
            Ω_L = [cholesky(omega) for _ in 1:nd],
            var"rt_state.sigma_rw" = fill(sigma_rw, nd),
            var"rt_state.intervention_effect" = fill(0.0, nd),
            var"gi_state.α" = fill(ALPHA, nd),
            var"gi_state.θ" = fill(THETA, nd),
            var"inc_state.delay_mean" = fill(INC[1], nd),
            var"inc_state.delay_sd" = fill(INC[2], nd),
            infections_patch = [vec(I[:, 1:N]) for _ in 1:nd],
        )
        eps === nothing || (
            base = merge(
                base, (;
                    importation_epsilon_patch = [copy(eps) for _ in 1:nd],
                    importation_epsilon_effect = fill(0.0, nd),
                )
            )
        )
        cases && (
            base = merge(
                base, (;
                    onset_to_confirmation_pmf = [copy(KC) for _ in 1:nd],
                    province_ascertainment = [copy(ASC) for _ in 1:nd],
                    province_composition_rho = fill(1.0e-6, nd),
                )
            )
        )
        deaths && (
            base = merge(
                base, (;
                    onset_to_death_confirmation_pmf = [copy(KD) for _ in 1:nd],
                    province_death_ascertainment = [copy(ASC_D) for _ in 1:nd],
                    province_cfr_relative = [copy(SEV) for _ in 1:nd],
                    province_death_composition_rho = fill(1.0e-6, nd),
                )
            )
        )
        return base
    end

    ## A chain carrying both confirmed compositions, and the national
    ## forecast whose totals the provinces split.
    full_chain(; nd = ND) = projection_chain(
        full_run(EPS); cases = true, deaths = true, nd
    )
    national(; nd = ND) = DataFrame(
        confirmed_new = fill(1_000_000, nd),
        confirmed_deaths_new = fill(200_000, nd)
    )
end

@testitem "forecast_provinces continues the patch renewal" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces, PROVINCE_LABELS

    I = full_run(EPS)
    fc = forecast_provinces(
        projection_chain(I), national();
        horizon = H, n_patches = NP
    )
    @test sort(unique(fc.patch)) == 1:NP
    @test fc.province == PROVINCE_LABELS[fc.patch]
    @test count(==(1), fc.patch) == ND
    for p in 1:NP
        rows = fc[fc.patch .== p, :]
        ## With the walks switched off the projection is the renewal the
        ## model itself runs, so it reproduces the held-back days.
        @test all(rows.infections_new .≈ sum(I[p, (N + 1):(N + H)]))
        @test all(rows.rt_forecast .≈ RT[p])
    end
    ## A national forecast with a draw count other than the chain's is an
    ## error, since the provinces split its draws one for one.
    @test_throws ArgumentError forecast_provinces(
        projection_chain(I), national(; nd = ND + 1); n_patches = NP
    )
end

@testitem "forecast_provinces runs uncoupled patches" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces

    ## A chain sampled without the importation kernel carries no intensity,
    ## and the patches then renew on their own.
    I = full_run(zeros(NP))
    fc = forecast_provinces(
        projection_chain(I; eps = nothing), national();
        horizon = H, n_patches = NP
    )
    for p in 1:NP
        @test all(
            fc[fc.patch .== p, :infections_new] .≈
                sum(I[p, (N + 1):(N + H)])
        )
    end
end

@testitem "forecast_provinces reverts provincial deviations" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces

    I = full_run(EPS)
    delta = [0.2, -0.1, -0.1]
    ## A one-week half-life halves the deviation over the week, so each
    ## province moves back toward the national trend by half its gap.
    fc = forecast_provinces(
        projection_chain(I; delta, halflife = 7.0), national();
        horizon = H, n_patches = NP
    )
    for p in 1:NP
        @test all(
            fc[fc.patch .== p, :rt_forecast] .≈ RT[p] * exp(-0.5 * delta[p])
        )
    end

    ## With fresh deviation innovations the deviations stay centred, so with
    ## the national walk off the provinces' log changes sum to zero.
    fc = forecast_provinces(
        projection_chain(I; delta, drift_sd = [0.3, 0.2, 0.1]), national();
        horizon = H, n_patches = NP
    )
    for d in 1:ND
        change = [
            log(only(fc[(fc.patch .== p) .& (fc.draw .== d), :rt_forecast]) / RT[p])
                for p in 1:NP
        ]
        @test abs(sum(change)) < 1.0e-10
    end
end

@testitem "forecast_provinces draws deviations with the fitted correlation" setup = [
    ProvinceProjection,
] begin
    using Statistics: cov
    using LinearAlgebra: Diagonal, I as Id
    using BVDOutbreakSize: forecast_provinces

    ## Each knot's innovations are the fitted scales times the fitted
    ## correlation's Cholesky factor times standard normals, centred across
    ## the provinces. Centring changes the covariance, so the target is
    ## P D Ω D P rather than Ω itself.
    nd = 4000
    sd = [0.3, 0.2, 0.1]
    omega = [1.0 0.8 0.3; 0.8 1.0 0.5; 0.3 0.5 1.0]
    fc = forecast_provinces(
        projection_chain(
            full_run(EPS); drift_sd = sd, omega, nd, halflife = 1.0e9
        ),
        national(; nd); horizon = H, n_patches = NP
    )
    ## With the national walk off and the deviations starting at zero, the
    ## week's log change in each province is its innovation.
    innov = [
        log(fc.rt_forecast[(fc.patch .== p) .& (fc.draw .== d)][1] / RT[p])
            for d in 1:nd, p in 1:NP
    ]
    P = Matrix{Float64}(Id, NP, NP) .- 1 / NP
    D = Diagonal(sd)
    target = P * D * omega * D * P
    independent = P * D * D * P
    got = cov(innov)
    @test all(abs.(got .- target) .< 0.006)
    ## The same check separates the fitted correlation from independent
    ## draws, so it would catch the correlation being dropped.
    @test abs(got[1, 2] - independent[1, 2]) > 0.03
end

@testitem "forecast_provinces splits the national total by the fitted composition" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces

    ## Diverging reproduction numbers, so the provinces' mix of recent
    ## infections moves across the delay and a share that ignored it would
    ## differ from the model's.
    rt = [1.25, 0.8, 1.0]
    I = full_run(EPS; rt)
    chn = merge(
        projection_chain(I; cases = true, deaths = true),
        (; R_T_patch = [copy(rt) for _ in 1:ND])
    )
    fc = forecast_provinces(chn, national(); horizon = H, n_patches = NP)
    for d in 1:ND
        rows = fc.draw .== d
        ## The provinces partition the national forecast, draw by draw, as
        ## the composition partitions the national count.
        @test sum(fc.confirmed_new[rows]) == 1_000_000
        @test sum(fc.confirmed_deaths_new[rows]) == 200_000
    end
    ## Expected shares from the model's own pipeline: incubation, then the
    ## onset-to-confirmation kernel, weighted by relative ascertainment, and
    ## for deaths the onset-to-death-confirmation kernel weighted by death
    ## ascertainment times relative severity.
    cases_share = window_shares(I, KC, ASC)
    deaths_share = window_shares(I, KD, ASC_D .* SEV)
    undelayed = let m = [sum(I[p, (N + 1):(N + H)]) for p in 1:NP]
        (ASC .* m) ./ sum(ASC .* m)
    end
    ## The delay matters in this fixture, so matching the delayed shares is a
    ## real check.
    @test maximum(abs.(cases_share .- undelayed)) > 0.02
    for p in 1:NP
        got_c = fc.confirmed_new[fc.patch .== p] ./ 1_000_000
        got_d = fc.confirmed_deaths_new[fc.patch .== p] ./ 200_000
        @test all(abs.(got_c .- cases_share[p]) .< 0.005)
        @test all(abs.(got_d .- deaths_share[p]) .< 0.005)
    end
end

@testitem "forecast_provinces writes no deaths without a death composition" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces

    ## The model generates province deaths only through the death
    ## composition, so a chain without one projects cases alone.
    fc = forecast_provinces(
        projection_chain(full_run(EPS); cases = true), national();
        n_patches = NP
    )
    @test :confirmed_new in propertynames(fc)
    @test !(:confirmed_deaths_new in propertynames(fc))
end

@testitem "province forecast summaries read a projection frame" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: province_forecast_table, PROVINCE_LABELS

    ## A projection frame carries each province's own draws, so the table
    ## reads them directly and needs no shares on the chain.
    fc = DataFrame(
        patch = repeat(1:2; inner = 100),
        province = repeat(PROVINCE_LABELS[1:2]; inner = 100),
        draw = repeat(1:100; outer = 2),
        confirmed_new = vcat(fill(10.0, 100), fill(50.0, 100))
    )
    tbl = province_forecast_table((;), fc; n_patches = 2)
    @test tbl.Province == PROVINCE_LABELS[1:2]
    @test tbl[!, "Lower 90%"] == [10.0, 50.0]
    @test tbl[!, "Upper 90%"] == [10.0, 50.0]

    ## The latent quantities follow each province's observed streams, the
    ## reproduction number to two decimals.
    fc.infections_new = vcat(fill(30.0, 100), fill(90.0, 100))
    fc.rt_forecast = vcat(fill(1.234, 100), fill(0.876, 100))
    tbl = province_forecast_table((;), fc; n_patches = 2)
    @test tbl.Province == repeat(PROVINCE_LABELS[1:2]; inner = 3)
    @test tbl.Quantity == repeat(
        [
            "New confirmed cases by T+7", "New infections by T+7",
            "Reproduction number at T+7",
        ], 2
    )
    @test tbl[!, "Upper 90%"] == [10.0, 30.0, 1.23, 50.0, 90.0, 0.88]
end

@testitem "plot_province_forecast captions a projection" setup = [
    HeadlessMakie,
] begin
    using DataFrames: DataFrame
    using CairoMakie: Makie as Mk
    using BVDOutbreakSize: plot_province_forecast, PROVINCE_LABELS

    fc = DataFrame(
        patch = repeat(1:2; inner = 50),
        province = repeat(PROVINCE_LABELS[1:2]; inner = 50),
        draw = repeat(1:50; outer = 2),
        confirmed_new = collect(1.0:100.0)
    )
    fig = plot_province_forecast((;), fc; n_patches = 2)
    captions = [x.text[] for x in fig.content if x isa Mk.Label]
    @test any(c -> occursin("own renewal equation", c), captions)
    @test !any(c -> occursin("modelled share", c), captions)
end

@testitem "province summaries project a national forecast" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces, province_forecast_table

    ## A national forecast is replaced by the projection from the chain, so
    ## the table reads the same draws whichever it is given.
    chn = full_chain()
    proj = forecast_provinces(chn, national(); n_patches = NP)
    @test province_forecast_table(chn, national(); n_patches = NP) ==
        province_forecast_table(chn, proj; n_patches = NP)

    ## A chain with no per-patch state cannot be projected by province.
    @test_throws ErrorException province_forecast_table(
        (; a = 1), national(); n_patches = NP
    )
end

@testitem "province_forecast_archive archives the projection" setup = [
    ProvinceProjection,
] begin
    using Dates: Date, Day
    using DataFrames: nrow
    using BVDOutbreakSize: forecast_provinces, province_forecast_archive,
        PROVINCE_NAMES

    chn = full_chain()
    made = Date("2026-09-06")
    arch = province_forecast_archive(
        chn, [(7, national()), (14, national())];
        made_date = made, n_patches = NP
    )
    @test names(arch) == [
        "made_date", "horizon", "target_date", "province", "stream",
        "draw", "value", "method",
    ]
    ## The method is recorded so scoring can leave out older archives.
    @test all(==("projection-v2"), arch.method)
    @test nrow(arch) == NP * 2 * 2 * ND
    @test Set(arch.province) == Set(PROVINCE_NAMES[1:NP])
    @test all(arch.target_date .== made .+ Day.(arch.horizon))
    ## Each horizon's values are that horizon's projection, not the national
    ## forecast passed in.
    for h in (7, 14)
        proj = forecast_provinces(
            chn, national(); horizon = h, n_patches = NP,
            patch_labels = PROVINCE_NAMES
        )
        got = arch[
            (arch.horizon .== h) .& (arch.province .== PROVINCE_NAMES[1]) .&
                (arch.stream .== "confirmed cases"), :value,
        ]
        @test got == float.(proj[proj.patch .== 1, :confirmed_new])
    end

    thinned = province_forecast_archive(
        chn, [(7, national())];
        made_date = made, n_patches = NP, thin = 2
    )
    @test nrow(thinned) == NP * 2 * cld(ND, 2)
end

@testitem "province_forecast_vs_truth scores the projection" setup = [
    ProvinceProjection,
] begin
    using Statistics: quantile
    using BVDOutbreakSize: forecast_provinces, province_forecast_vs_truth

    chn = full_chain()
    kw = (;
        observed = [900, 100, 30], baseline = [800, 60, 20],
        death_observed = [40, 20, 5], death_baseline = [20, 10, 5],
        n_patches = NP,
    )
    df = province_forecast_vs_truth(chn, national(); kw...)
    @test size(df, 1) == 2 * NP
    cases = df[df[!, "Stream"] .== "Confirmed cases", :]
    @test cases[!, "Observed"] == [100, 40, 10]
    deaths = df[df[!, "Stream"] .== "Confirmed deaths", :]
    @test deaths[!, "Observed"] == [20, 10, 0]
    ## The interval is the projection's, so the national forecast passed in
    ## does not reach it.
    proj = forecast_provinces(chn, national(); n_patches = NP)
    v = float.(proj[proj.patch .== 1, :confirmed_new])
    @test cases[1, "Upper 90%"] == round(quantile(v, 0.95))
    @test province_forecast_vs_truth(chn, proj; kw...) == df
    @test eltype(df[!, "Within 90% PI"]) == Bool
    @test !("Central estimate" in names(df))

    ## A chain with no per-patch state cannot be scored by province.
    @test_throws ErrorException province_forecast_vs_truth(
        (; a = 1), national(); kw...
    )
end

@testitem "province_share_draws splits each draw's total" begin
    using DataFrames: DataFrame
    using BVDOutbreakSize: province_share_draws

    fc = DataFrame(
        patch = [1, 2, 1, 2], draw = [1, 1, 2, 2],
        confirmed_new = [30.0, 10.0, 0.0, 0.0]
    )
    sh = province_share_draws(fc, :confirmed_new; n_patches = 2)
    ## Shares are per draw, and a draw with nothing projected anywhere has
    ## no share to give, so it is left out rather than divided by zero.
    @test sh == [[0.75], [0.25]]
end

@testitem "province_forecast_table projects a national forecast at its horizon" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces, province_forecast_table

    chn = full_chain()
    proj = forecast_provinces(chn, national(); horizon = 14, n_patches = NP)
    tbl = province_forecast_table(chn, national(); horizon = 14, n_patches = NP)
    @test tbl == province_forecast_table(chn, proj; horizon = 14, n_patches = NP)
    @test "New confirmed cases by T+14" in tbl.Quantity
    @test "Reproduction number at T+14" in tbl.Quantity
end

@testitem "forecast_provinces shares the national walk across provinces" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces

    ## With the deviations held at zero, each province's reproduction number
    ## moves only with the national walk, so every province moves by the same
    ## factor within a draw, and the walk does move it.
    fc = forecast_provinces(
        projection_chain(full_run(EPS); sigma_rw = 0.3), national();
        horizon = H, n_patches = NP
    )
    moves = [
        [
            log(only(fc[(fc.patch .== p) .& (fc.draw .== d), :rt_forecast]) / RT[p])
                for p in 1:NP
        ] for d in 1:ND
    ]
    @test all(m -> all(x -> isapprox(x, m[1]; atol = 1.0e-12), m), moves)
    @test any(m -> abs(m[1]) > 1.0e-6, moves)
end

@testitem "the national walk keeps its draws" begin
    using Random: MersenneTwister
    using Distributions: Gamma
    using BVDOutbreakSize: cdf_nmax, euler_lotka_r

    ## The walk continuation as `_evolving_rates` wrote it before it was
    ## shared with `forecast_provinces`, kept here so the national forecast
    ## is shown to draw exactly the same values with the same seed.
    function reference(chn, horizon; rng, week = 7)
        sigma = chn[Symbol("rt_state.sigma_rw")]
        R_T = chn[:R_T]
        α = chn[Symbol("gi_state.α")]
        θ = chn[Symbol("gi_state.θ")]
        nknots = cld(horizon, week)
        paths = Vector{Vector{Float64}}(undef, length(R_T))
        rt_term = Vector{Float64}(undef, length(R_T))
        for i in eachindex(R_T)
            innov = sigma[i] .* randn(rng, nknots)
            cum_innov = cumsum(innov)
            g = BVDOutbreakSize._gi_pmf(α[i], θ[i])
            rs = Vector{Float64}(undef, horizon)
            log_R = log(max(R_T[i], 1.0e-6))
            log_rt = log_R
            for d in 1:horizon
                weeks = d / week
                j = floor(Int, weeks)
                whole = j == 0 ? 0.0 : cum_innov[j]
                part = j < nknots ? (weeks - j) * innov[j + 1] : 0.0
                log_rt = log_R + whole + part
                rs[d] = euler_lotka_r(
                    max(exp(log_rt), BVDOutbreakSize._RT_EULER_FLOOR), g
                )
            end
            paths[i] = rs
            rt_term[i] = exp(log_rt)
        end
        return (; paths, rt_term)
    end

    nd = 25
    chn = (;
        R_T = collect(range(0.6, 1.6; length = nd)),
        var"rt_state.z" = zeros(nd),
        var"rt_state.sigma_rw" = fill(0.2, nd),
        var"gi_state.α" = fill(2.71, nd),
        var"gi_state.θ" = fill(5.65, nd),
    )
    for horizon in (7, 10, 14, 28)
        got = BVDOutbreakSize._evolving_rates(
            chn, horizon; rng = MersenneTwister(11)
        )
        want = reference(chn, horizon; rng = MersenneTwister(11))
        @test got.paths == want.paths
        @test got.rt_term == want.rt_term
    end
end
