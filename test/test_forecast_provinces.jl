## Tests for the per-province forward projection, `forecast_provinces`.

@testsnippet ProvinceProjection begin
    using BVDOutbreakSize: patch_infections, province_importation_kernel,
        PROVINCE_POPULATIONS, cdf_nmax
    using Distributions: Gamma
    using DataFrames: DataFrame

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

    function full_run(eps)
        seeds = fill(5.0, NP, 10)
        Rt = repeat(RT, 1, N + H)
        ε = repeat(eps, 1, N + H)
        return patch_infections(Rt, G, seeds, KERNEL, ε).infections
    end

    function projection_chain(
            I; delta = zeros(NP), drift_sd = zeros(NP),
            halflife = 42.0, sigma_rw = 0.0, eps = EPS,
            drift_factor = nothing,
            conf_daily = nothing, shares = nothing,
            deaths_daily = nothing, death_shares = nothing
        )
        base = (;
            R_T_patch = [copy(RT) for _ in 1:ND],
            delta_patch = [copy(delta) for _ in 1:ND],
            region_drift_sd = [copy(drift_sd) for _ in 1:ND],
            region_halflife = fill(halflife, ND),
            var"rt_state.sigma_rw" = fill(sigma_rw, ND),
            var"gi_state.α" = fill(ALPHA, ND),
            var"gi_state.θ" = fill(THETA, ND),
            infections_patch = [vec(I[:, 1:N]) for _ in 1:ND],
        )
        eps === nothing ||
            (base = merge(base, (; importation_epsilon_patch = [copy(eps) for _ in 1:ND])))
        drift_factor === nothing || (
            base = merge(
                base,
                (; region_drift_factor = [vec(drift_factor) for _ in 1:ND])
            )
        )
        if conf_daily !== nothing
            ## A cumulative whose last increment is the cut-off daily rate.
            cum = collect(range(0.0, conf_daily * N; length = N))
            cum[end] = cum[end - 1] + conf_daily
            base = merge(
                base, (;
                    cumulative_confirmed = [copy(cum) for _ in 1:ND],
                    k_confirmed = fill(1.0e8, ND),
                    province_shares = [copy(shares) for _ in 1:ND],
                )
            )
        end
        if deaths_daily !== nothing
            cum = collect(range(0.0, deaths_daily * N; length = N))
            cum[end] = cum[end - 1] + deaths_daily
            base = merge(
                base, (;
                    cumulative_confirmed_deaths = [copy(cum) for _ in 1:ND],
                    k_confirmed_deaths = fill(1.0e8, ND),
                    province_death_shares = [copy(death_shares) for _ in 1:ND],
                )
            )
        end
        return base
    end

    ## A chain carrying both confirmed streams, and a national forecast the
    ## province summaries must not read.
    full_chain() = projection_chain(
        full_run(EPS);
        conf_daily = 1.0e3, shares = [0.7 0.6; 0.2 0.3; 0.1 0.1],
        deaths_daily = 1.0e2, death_shares = [0.5 0.5; 0.3 0.3; 0.2 0.2]
    )
    national() = DataFrame(
        confirmed_new = fill(1.0e6, ND),
        confirmed_deaths_new = fill(1.0e6, ND)
    )
end

@testitem "forecast_provinces continues the patch renewal" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces, PROVINCE_LABELS

    I = full_run(EPS)
    fc = forecast_provinces(
        projection_chain(I);
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
end

@testitem "forecast_provinces runs uncoupled patches" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces

    ## A chain sampled without the importation kernel carries no intensity,
    ## and the patches then renew on their own.
    I = full_run(zeros(NP))
    fc = forecast_provinces(
        projection_chain(I; eps = nothing);
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
        projection_chain(I; delta, halflife = 7.0);
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
        projection_chain(I; delta, drift_sd = [0.3, 0.2, 0.1]);
        horizon = H, n_patches = NP
    )
    for d in 1:ND
        change = [
            log(only(fc[(fc.patch .== p) .& (fc.draw .== d), :rt_forecast]) / RT[p])
                for p in 1:NP
        ]
        @test abs(sum(change)) < 1.0e-10
    end

    ## A chain carrying the fitted loading matrix draws its innovations
    ## through it. This one moves only the first two provinces, in opposite
    ## directions, so the third keeps its reverted value.
    F = [0.3 0.0; -0.3 0.0; 0.0 0.0]
    fc = forecast_provinces(
        projection_chain(
            I; delta, halflife = 7.0, drift_sd = [0.25, 0.25, 0.0],
            drift_factor = F
        );
        horizon = H, n_patches = NP
    )
    for d in 1:ND
        rt(p) = only(fc[(fc.patch .== p) .& (fc.draw .== d), :rt_forecast])
        change = [log(rt(p) / RT[p]) + 0.5 * delta[p] for p in 1:NP]
        @test abs(change[3]) < 1.0e-12
        @test change[1] ≈ -change[2]
        @test abs(change[1]) > 1.0e-6
    end
    ## Its rows sum to zero only over every patch, so asking for fewer is
    ## an error rather than a projection that no longer sums to zero.
    @test_throws ArgumentError forecast_provinces(
        projection_chain(
            I; delta, drift_sd = [0.25, 0.25, 0.0], drift_factor = F
        );
        horizon = H, n_patches = NP - 1
    )
end

@testitem "forecast_provinces grows confirmed counts at each patch's rate" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces

    I = full_run(EPS)
    shares = [0.7 0.6; 0.2 0.3; 0.1 0.1]
    conf_daily = 1.0e7
    fc = forecast_provinces(
        projection_chain(I; conf_daily, shares);
        horizon = H, n_patches = NP
    )
    @test !(:confirmed_deaths_new in propertynames(fc))
    for p in 1:NP
        ## The cut-off level is the national daily rate times the share at
        ## the last vintage, then each day grows with the patch's own
        ## infections. A near-Poisson replicate at this size sits within a
        ## fraction of a percent of its mean.
        growth = sum(I[p, (N + 1):(N + H)]) / I[p, N]
        expected = conf_daily * shares[p, end] * growth
        got = fc[fc.patch .== p, :confirmed_new]
        @test all(abs.(got .- expected) .< 0.01 * expected)
    end
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
    proj = forecast_provinces(chn; n_patches = NP)
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
    @test all(==("projection"), arch.method)
    @test nrow(arch) == NP * 2 * 2 * ND
    @test Set(arch.province) == Set(PROVINCE_NAMES[1:NP])
    @test all(arch.target_date .== made .+ Day.(arch.horizon))
    ## Each horizon's values are that horizon's projection, not the national
    ## forecast passed in.
    for h in (7, 14)
        proj = forecast_provinces(
            chn; horizon = h, n_patches = NP, patch_labels = PROVINCE_NAMES
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
    proj = forecast_provinces(chn; n_patches = NP)
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
    proj = forecast_provinces(chn; horizon = 14, n_patches = NP)
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
        projection_chain(full_run(EPS); sigma_rw = 0.3);
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
