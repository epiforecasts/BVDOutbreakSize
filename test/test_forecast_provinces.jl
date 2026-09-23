## Tests for the per-province forward projection, `forecast_provinces`.

@testsnippet ProvinceProjection begin
    using BVDOutbreakSize: patch_infections, province_importation_kernel,
        PROVINCE_POPULATIONS, cdf_nmax
    using Distributions: Gamma

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
            conf_daily = nothing, shares = nothing
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
        return base
    end
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
end

@testitem "forecast_provinces grows confirmed counts at each patch's rate" setup = [
    ProvinceProjection,
] begin
    using BVDOutbreakSize: forecast_provinces

    I = full_run(EPS)
    shares = [0.7 0.6; 0.2 0.3; 0.1 0.1]
    conf_daily = 1.0e5
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
