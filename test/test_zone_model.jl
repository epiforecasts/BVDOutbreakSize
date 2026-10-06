## Unit tests for the health-zone model: the pure share-renewal blocks in
## src/models/zone.jl, the Dirichlet-multinomial composition, the fixed
## inputs built from a parent chain, the bvd_zone model, and the chain
## post-processing in src/zone.jl.

## A synthetic two-patch zone dataset generated from a known share process,
## plus a stand-in parent chain (a `Dict` of iteration-by-chain matrices,
## which is all the input builder reads) so every item runs without a real
## joint fit. The truth is kept alongside so the recovery items can compare.
@testsnippet ZoneSynthetic begin
    using BVDOutbreakSize
    using BVDOutbreakSize: zone_initial_shares, deviation_knots,
        deviation_knot_dims, zone_fixed_terms, zone_binned_operator,
        zone_binned_increments, zone_forward, discretise_censored,
        lognormal_meansd, convolve_pmf, knot_days,
        zone_interpolation_weights, future_knot_days, zone_delay_pmfs,
        zone_delay_terms
    using Dates: Date, Day
    using Distributions: Gamma, Multinomial
    using Random: Xoshiro

    function zone_synthetic(;
            n = 120, first_vintage = 60, step = 3,
            lead_days = 42, seed = 11, total_scale = 1.0
        )
        rng = Xoshiro(seed)
        np = 2
        ## Fixed patch infections: growth then a plateau, two very different
        ## sizes so the small patch exercises the near-zero force paths.
        I_bar = zeros(Float64, np, n)
        for t in 1:n
            I_bar[1, t] = 30.0 * exp(0.03 * min(t, 80)) * total_scale
            I_bar[2, t] = 5.0 * exp(0.02 * min(t, 80)) * total_scale
        end
        ## The delay parameters of every parent draw, in the order of
        ## `_ZONE_DELAY_KEYS`, at the centres of the parent's priors.
        delay_params = [
            2.71 * 5.65, sqrt(2.71) * 5.65, 6.3, 3.5, 1.178, 3.694,
            4.5, 4.0, 1.178, 3.694, 2.151, 3.906,
        ]
        pmfs = zone_delay_pmfs(delay_params)
        g, f = pmfs.g, pmfs.f
        zones = [
            ("a", ["z1", "z2", "z3", "z4", "z5"]),
            ("b", ["y1", "y2", "y3"]),
        ]
        patch_ranges = [1:5, 6:8]
        nz = 8
        days = collect(first_vintage:step:n)
        days[end] == n || push!(days, n)
        t0 = clamp(days[1] - lead_days, 1, n)
        knots = knot_days(n; week = 7, start = t0)
        K = length(knots)
        ## The truth: uneven initial shares, a level contrast, and one zone
        ## in the big patch walking upward while another walks down. Every
        ## zone clears the walking threshold, so all carry innovations and
        ## the truth is representable in the fitted model.
        z_w = [1.2, 0.3, -0.2, -0.6, -0.7, 0.8, -0.3, -0.5]
        w0 = zone_initial_shares(z_w, patch_ranges, 2.0)
        walking = fill(true, nz)
        σ_level = 0.25
        ## One drift scale per patch, as the model samples it.
        σ_δ = fill(0.08, np)
        φ = exp2(-7 / 42)
        ## `n_p - 1` level draws per patch on the sum-to-zero basis.
        z_level = [0.6, -0.4, 0.2, -0.3, 0.5, -0.2]
        ## Per knot, the first basis direction of the big patch contrasts
        ## its first two zones: `+0.8` and `-0.8` in the draw's units.
        nd = deviation_knot_dims(patch_ranges, walking).drift
        z_drift = zeros(nd * (K - 1))
        for k in 2:K
            z_drift[(k - 2) * nd + 1] = 0.8 * sqrt(2)
        end
        patch_of_zone = [1, 1, 1, 1, 1, 2, 2, 2]
        δ_knots = deviation_knots(
            z_level, z_drift, σ_level, σ_δ, φ,
            patch_ranges, Matrix{Float64}[], Matrix{Float64}[], walking, K
        )
        terms = zone_delay_terms(I_bar, pmfs, t0, days, patch_of_zone)
        zd = (;
            I_bar, g, f, patch_ranges, patch_of_zone, knots, t0, n, days,
            terms.force_pre, terms.infections_pre, mixing = nothing,
            interp = zone_interpolation_weights(knots, t0, n),
            terms.report_bin, terms.death_bin,
        )
        fw = zone_forward(zd, δ_knots, w0, nothing)
        ## Observed counts: each vintage's allocated total is the rounded
        ## modelled patch total, split multinomially at the modelled shares.
        nv = length(days)
        counts = zeros(Int, nz, nv)
        for v in 1:nv, (p, zs) in enumerate(patch_ranges)

            C = fw.increments[zs, v]
            N = round(Int, sum(C))
            N > 0 || continue
            counts[zs, v] = rand(rng, Multinomial(N, C ./ sum(C)))
        end
        ## The manifest shape: cumulative per zone on the shared vintage days.
        hist = Dict{String, Dict{String, NamedTuple}}()
        z = 0
        for (prov, names) in zones
            hist[prov] = Dict{String, NamedTuple}()
            for nm in names
                z += 1
                hist[prov][nm] = (;
                    days = copy(days),
                    counts = cumsum(counts[z, :]),
                )
            end
        end
        ## Allocated zone deaths, a tenth of the confirmed, so the death
        ## composition the model always fits has something to score.
        death_hist = Dict{String, Dict{String, NamedTuple}}()
        for (prov, zones) in hist
            death_hist[prov] = Dict{String, NamedTuple}()
            for (nm, h) in zones
                death_hist[prov][nm] = (;
                    days = copy(h.days), counts = cld.(h.counts, 10),
                )
            end
        end
        seeding = Date("2026-02-13")
        obs = (;
            n, seeding, cutoff = seeding + Day(n - 1),
            zone_confirmed_history = hist,
            zone_death_history = death_hist,
        )
        ## Stand-in parent chain: four identical draws of the truth.
        ndraw = 4
        chain = Dict{Symbol, Any}(
            :infections_patch => reshape(
                [vec(I_bar) for _ in 1:ndraw],
                ndraw, 1
            ),
            Symbol("gi_state.gi_mean") => fill(2.71 * 5.65, ndraw, 1),
            Symbol("gi_state.gi_sd") => fill(sqrt(2.71) * 5.65, ndraw, 1),
            Symbol("inc_state.delay_mean") => fill(6.3, ndraw, 1),
            Symbol("inc_state.delay_sd") => fill(3.5, ndraw, 1),
            Symbol("cases_state.report_state.α") => fill(1.178, ndraw, 1),
            Symbol("cases_state.report_state.θ") => fill(3.694, ndraw, 1),
            Symbol("confirmed_state.receipt_state.d.delay_mean") =>
                fill(4.5, ndraw, 1),
            Symbol("confirmed_state.receipt_state.d.delay_sd") =>
                fill(4.0, ndraw, 1),
            Symbol("deaths_state.od_state.oa.α") => fill(1.178, ndraw, 1),
            Symbol("deaths_state.od_state.oa.θ") => fill(3.694, ndraw, 1),
            Symbol("deaths_state.od_state.ad.α") => fill(2.151, ndraw, 1),
            Symbol("deaths_state.od_state.ad.θ") => fill(3.906, ndraw, 1),
            :C_T => fill(sum(I_bar), ndraw, 1),
            ## The province model's own between-patch movement: a per-origin
            ## intensity and the arrivals into each patch, a twentieth of
            ## the second patch's infections and none into the first.
            :importation_epsilon_patch => reshape(
                [[0.02, 0.01] for _ in 1:ndraw], ndraw, 1
            ),
            :importation_patch => reshape(
                [
                    vec(
                        [
                            p == 2 ? 0.05 * I_bar[p, t] : 0.0
                                for p in 1:np, t in 1:n
                        ]
                    )
                        for _ in 1:ndraw
                ],
                ndraw, 1
            )
        )
        truth = (;
            z_w, w0, δ_knots, z_level, z_drift, σ_level, σ_δ, φ,
            δ_daily = zd.interp * transpose(δ_knots),
            walking, shares = fw.shares, increments = fw.increments,
            infections = fw.infections, forces = fw.forces,
        )
        return (;
            obs, chain, truth, I_bar, g, f, pmfs, delay_params,
            patch_ranges, days, t0, knots, counts, nz,
        )
    end

    ## Reference forms of the delay convolution, from the formula: the
    ## lower-triangular Toeplitz matrix with `F[j, j − s] = f_s` over the
    ## grid days, and each zone's patch's pre-`t0` term on those days.
    function reference_delay_matrix(f, nd)
        F = zeros(nd, nd)
        for j in 1:nd, s in 0:min(length(f) - 1, j - 1)

            F[j, j - s] = f[s + 1]
        end
        return F
    end
    reference_pre_rows(report_pre, patch_of_zone, t0, n) =
        [report_pre[p, t] for t in t0:n, p in patch_of_zone]

    function zone_inputs(
            syn; zones = nothing, walk_threshold = 30,
            kwargs...
        )
        return zone_fit_inputs(
            syn.chain, syn.obs; zones, walk_threshold,
            patch_names = ["a", "b"], patch_labels = ["A", "B"], kwargs...
        )
    end

    ## Stand-in parent forecast for the chain of `zone_synthetic`: each draw
    ## holds its patch infections at the cut-off level, a little apart from
    ## draw to draw, for `H` days, and predicts `totals` confirmed cases per
    ## patch in every forecast week. The arrivals into the second patch grow
    ## from a twentieth of its infections by a fiftieth a day.
    function zone_parent_forecast(syn; H = 28, totals = [80, 20])
        np, n = size(syn.I_bar)
        nw = length(future_knot_days(0, H))
        ndraw = size(syn.chain[:infections_patch], 1)
        inf = [
            vec(repeat(syn.I_bar[:, n], 1, H) .* (1 + 0.01 * i))
                for i in 1:ndraw
        ]
        conf = [vec(repeat(totals, 1, nw)) for _ in 1:ndraw]
        imp = [
            vec(
                [
                    p == 2 ? 0.05 * (1 + 0.02 * d) * reshape(v, np, H)[p, d] :
                        0.0 for p in 1:np, d in 1:H
                ]
            )
                for v in inf
        ]
        return Dict{Symbol, Any}(
            :forecast_infections_patch => reshape(inf, ndraw, 1),
            :forecast_importation_patch => reshape(imp, ndraw, 1),
            :forecast_province_confirmed => reshape(conf, ndraw, 1),
        )
    end

    ## The chain of `zone_synthetic` with the between-patch movement spread
    ## over the draws: each draw scales both origins' intensity, and the
    ## arrivals into the second patch with the first origin's.
    function zone_spread_chain(syn; scales = [0.5, 1.0, 2.0, 1.5])
        np, n = size(syn.I_bar)
        chain = copy(syn.chain)
        chain[:importation_epsilon_patch] = reshape(
            [[0.02 * s, 0.01 * sqrt(s)] for s in scales], :, 1
        )
        chain[:importation_patch] = reshape(
            [
                vec(
                    [
                        p == 2 ? 0.05 * s * syn.I_bar[p, t] : 0.0
                            for p in 1:np, t in 1:n
                    ]
                )
                    for s in scales
            ],
            :, 1
        )
        return chain
    end

    ## Health-zone metadata rows covering every synthetic zone, for the
    ## mixing kernel.
    function zone_metadata(syn)
        pairs = [("a", "z$i") for i in 1:5]
        append!(pairs, [("b", "y$i") for i in 1:3])
        return [
            (;
                zone = nm, label = uppercase(nm), province = prov,
                population = 10_000 * i, lat = 0.1 * i, lon = 0.2 * i,
                zscode = string(i),
            )
                for (i, (prov, nm)) in enumerate(pairs)
        ]
    end
end

@testitem "dirichlet_multinomial_logpdf: matches Distributions" begin
    using BVDOutbreakSize: dirichlet_multinomial_logpdf
    using Distributions: DirichletMultinomial, logpdf

    y = [7, 0, 3, 12]
    α = [2.0, 0.5, 1.5, 6.0]
    @test dirichlet_multinomial_logpdf(y, α) ≈
        logpdf(DirichletMultinomial(sum(y), α), y)
    ## A single-category split is certain.
    @test dirichlet_multinomial_logpdf([5], [3.0]) ≈ 0 atol = 1.0e-12
    ## Zero trials contribute nothing.
    @test dirichlet_multinomial_logpdf([0, 0], [1.0, 2.0]) ≈ 0 atol = 1.0e-12
end

@testitem "zone composition: one sum of Dirichlet-multinomials" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: zone_composition_logpdf,
        dirichlet_multinomial_logpdf, _zone_kappa
    using Distributions: DirichletMultinomial, logpdf

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    ρ = 0.05
    κ = (1 - ρ) / ρ
    C = syn.truth.increments
    lp = zone_composition_logpdf(
        zd.counts, C, zd.cell_patch, zd.cell_vintage,
        zd.cell_total, zd.cell_const, zd.patch_ranges, κ
    )
    ## The same sum written cell by cell through the full mass.
    cell_lp(c) = begin
        v = zd.cell_vintage[c]
        zs = zd.patch_ranges[zd.cell_patch[c]]
        π = C[zs, v] ./ sum(C[zs, v])
        logpdf(
            DirichletMultinomial(zd.cell_total[c], κ .* π),
            zd.counts[zs, v]
        )
    end
    expected = sum(cell_lp(c) for c in eachindex(zd.cell_patch))
    @test lp ≈ expected rtol = 1.0e-10
    ## Scale-free in the modelled level: only the split is scored.
    lp2 = zone_composition_logpdf(
        zd.counts, 100 .* C, zd.cell_patch,
        zd.cell_vintage, zd.cell_total, zd.cell_const, zd.patch_ranges, κ
    )
    @test lp2 ≈ lp rtol = 1.0e-10
    ## The concentration stays finite and positive at both ends of `ρ`'s
    ## support, so a saturated proposal never scores `NaN`.
    for ρ_end in (0.0, 1.0)
        κ_end = _zone_kappa(ρ_end)
        @test isfinite(κ_end) && κ_end > 0
        @test isfinite(
            zone_composition_logpdf(
                zd.counts, C, zd.cell_patch,
                zd.cell_vintage, zd.cell_total, zd.cell_const, zd.patch_ranges,
                κ_end
            )
        )
    end
    @test _zone_kappa(0.05) ≈ 19 rtol = 1.0e-12
    ## Every scored cell has a positive allocated total.
    @test all(>(0), zd.cell_total)
    @test all(
        c -> sum(
            zd.counts[
                zd.patch_ranges[zd.cell_patch[c]],
                zd.cell_vintage[c],
            ]
        ) == zd.cell_total[c],
        eachindex(zd.cell_total)
    )
end

@testitem "zone share renewal: conservation and the pre-t0 shortcut" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: zone_share_renewal, zone_fixed_terms

    syn = zone_synthetic()
    I_bar, g = syn.I_bar, syn.g
    np, n = size(I_bar)
    t0 = syn.t0
    nd = n - t0 + 1
    w0 = syn.truth.w0
    δ_daily = syn.truth.δ_daily
    fixed = zone_fixed_terms(I_bar, g, syn.f, t0)
    st = zone_share_renewal(
        I_bar, g, δ_daily, w0, syn.patch_ranges, t0,
        fixed.force_pre
    )
    ## Shares sum to one within every patch on every day, and the zone
    ## infections sum to the patch infections.
    for j in 1:nd, (p, zs) in enumerate(syn.patch_ranges)

        @test sum(st.shares[j, zs]) ≈ 1 rtol = 1.0e-12
        @test sum(st.infections[j, zs]) ≈ I_bar[p, t0 + j - 1] rtol = 1.0e-12
    end
    ## The pre-t0 constant force equals the explicit convolution of a
    ## history held at the initial share.
    L = length(g)
    for (p, zs) in enumerate(syn.patch_ranges), z in zs, j in 1:5:nd
        t = t0 + j - 1
        explicit = 0.0
        for s in 1:min(L, t - 1)
            I_past = t - s < t0 ? w0[z] * I_bar[p, t - s] :
                st.infections[t - s - t0 + 1, z]
            explicit += g[s] * I_past
        end
        @test st.forces[j, z] ≈ explicit rtol = 1.0e-10
    end
end

@testitem "zone operators: matrix forms equal the loops" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: zone_share_renewal, interpolate_knots

    syn = zone_synthetic()
    zd = zone_inputs(syn).model_data
    δk = syn.truth.δ_knots
    t0, n = zd.t0, zd.n
    nd = n - t0 + 1
    ## The interpolation weights reproduce the per-zone knot interpolation.
    δd = zd.interp * transpose(δk)
    for z in 1:syn.nz
        ref = interpolate_knots(δk[z, :], zd.knots, n)[t0:n]
        @test δd[:, z] ≈ ref rtol = 1.0e-12
    end
    ## The delay operator plus the pre-t0 rows reproduce the convolution,
    ## with the days before t0 held at the initial share.
    w0 = syn.truth.w0
    st = zone_share_renewal(
        zd.I_bar, zd.g, δd, w0, zd.patch_ranges, t0,
        zd.force_pre
    )
    fixed = zone_fixed_terms(zd.I_bar, zd.g, zd.f, t0)
    r_op = reference_delay_matrix(zd.f, nd) * st.infections .+
        reference_pre_rows(fixed.report_pre, zd.patch_of_zone, t0, n) .*
        transpose(w0)
    for (p, zs) in enumerate(zd.patch_ranges), z in zs, j in 1:7:nd
        t = t0 + j - 1
        explicit = 0.0
        for s in 0:min(length(zd.f) - 1, t - 1)
            I_past = t - s < t0 ? w0[z] * zd.I_bar[p, t - s] :
                st.infections[j - s, z]
            explicit += zd.f[s + 1] * I_past
        end
        @test r_op[j, z] ≈ explicit rtol = 1.0e-10
    end
    ## The model's forward pass matches the generator's increments.
    fw = zone_forward(zd, δk, w0, nothing)
    @test fw.increments ≈ syn.truth.increments rtol = 1.0e-10
    ## The patch index of every zone, for data that carries the ranges
    ## alone.
    @test BVDOutbreakSize._zone_patch_of_zone(zd.patch_ranges) ==
        zd.patch_of_zone
    @test isempty(BVDOutbreakSize._zone_patch_of_zone(UnitRange{Int}[]))
end

@testitem "zone Rt: the force-weighted mean is the implied patch Rt" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: zone_share_renewal, zone_fixed_terms,
        implied_national_Rt_at

    syn = zone_synthetic()
    I_bar, g = syn.I_bar, syn.g
    np, n = size(I_bar)
    t0 = syn.t0
    δ_daily = syn.truth.δ_daily
    fixed = zone_fixed_terms(I_bar, g, syn.f, t0)
    st = zone_share_renewal(
        I_bar, g, δ_daily, syn.truth.w0,
        syn.patch_ranges, t0, fixed.force_pre
    )
    for (p, zs) in enumerate(syn.patch_ranges), t in t0:7:n

        j = t - t0 + 1
        R = [st.infections[j, z] / st.forces[j, z] for z in zs]
        Λ = [st.forces[j, z] for z in zs]
        @test sum(R .* Λ) / sum(Λ) ≈
            implied_national_Rt_at(vec(I_bar[p, :]), g, t) rtol = 1.0e-10
    end
end

@testitem "zone deviations: centred, and level-only zones decay" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: deviation_knots, zone_initial_shares

    syn = zone_synthetic()
    K = size(syn.truth.δ_knots, 2)
    φ = syn.truth.φ
    ## Three walking zones in the big patch, none in the small one: two
    ## draws per knot.
    walking = [true, true, true, false, false, false, false, false]
    z_drift = randn(Xoshiro(3), 2 * (K - 1))
    δ = deviation_knots(
        syn.truth.z_level, z_drift, 0.25, fill(0.08, 2), φ,
        syn.patch_ranges, Matrix{Float64}[], Matrix{Float64}[], walking, K
    )
    for (p, zs) in enumerate(syn.patch_ranges), k in 1:K

        @test abs(sum(δ[zs, k])) < 1.0e-12
    end
    ## A level-only zone follows the AR mean path from its level.
    for z in findall(!, walking), k in 1:K

        @test δ[z, k] ≈ φ^(k - 1) * δ[z, 1] rtol = 1.0e-12
    end
    ## A walking zone leaves the mean path.
    @test !isapprox(δ[1, K], φ^(K - 1) * δ[1, 1]; rtol = 1.0e-3)
    ## The synthetic truth's own knots are centred too.
    for (p, zs) in enumerate(syn.patch_ranges), k in 1:K

        @test abs(sum(syn.truth.δ_knots[zs, k])) < 1.0e-12
    end
    ## The initial shares sum to one within each patch and are invariant to
    ## a common shift of the draws.
    w = zone_initial_shares(syn.truth.z_w, syn.patch_ranges, 2.0)
    w_shift = zone_initial_shares(syn.truth.z_w .+ 3.0, syn.patch_ranges, 2.0)
    for zs in syn.patch_ranges
        @test sum(w[zs]) ≈ 1 rtol = 1.0e-12
    end
    @test w ≈ w_shift rtol = 1.0e-12
end

@testitem "zone deviations: the basis keeps the centred covariance" begin
    using BVDOutbreakSize: deviation_knots, deviation_knot_dims,
        zone_correlation_factors
    using LinearAlgebra: I
    using Random: Xoshiro

    groups = [1:5, 6:6, 7:10]
    walking = Bool[1, 0, 1, 1, 0, 0, 1, 1, 0, 0]
    x = [0.0, 0.3, 1.1, 2.0, 2.4, 0.0, 0.5, 0.9, 1.7, 3.0]
    dist(zs) = [abs(x[a] - x[b]) for a in zs, b in zs]
    walkers(us) = [u for u in us if walking[u]]
    ℓ = 0.8
    A = zone_correlation_factors([dist(us) for us in groups], ℓ)
    Aw = zone_correlation_factors([dist(walkers(us)) for us in groups], ℓ)
    dims = deviation_knot_dims(groups, walking)
    @test dims == (; level = 4 + 0 + 3, drift = 2 + 0 + 1)
    σ_level, σ_δ, φ, K = 0.3, [0.1, 0.2, 0.05], 0.8, 3
    knots(zl, zδ) = deviation_knots(
        zl, zδ, σ_level, σ_δ, φ, groups, A, Aw, walking, K
    )
    unit(i, m) = [j == i ? 1.0 : 0.0 for j in 1:m]
    nδ = dims.drift * (K - 1)
    ## Both knots are linear in their draws, so the covariance is `J Jᵀ`.
    level(i) = knots(unit(i, dims.level), zeros(nδ))[:, 1]
    drift(i) = knots(zeros(dims.level), unit(i, nδ))[:, 2]
    Jl = reduce(hcat, [level(i) for i in 1:dims.level])
    Jδ = reduce(hcat, [drift(i) for i in 1:dims.drift])
    ## The covariance of the centred construction `σ P L z` with `L Lᵀ = C`:
    ## `σ² P C P` over each group, and over the walking zones for the
    ## innovations.
    function centred_cov(σs, sets)
        Σ = zeros(length(walking), length(walking))
        for (σ, zs) in zip(σs, sets)
            m = length(zs)
            m == 0 && continue
            C = exp.(.-dist(zs) ./ ℓ) + 1.0e-6 * I
            P = I - fill(1 / m, m, m)
            Σ[zs, zs] = σ^2 .* (P * C * P)
        end
        return Σ
    end
    @test Jl * Jl' ≈ centred_cov(fill(σ_level, 3), groups) atol = 1.0e-12
    @test Jδ * Jδ' ≈ centred_cov(σ_δ, walkers.(groups)) atol = 1.0e-12
    ## Every group sums to zero at every knot, and a level-only zone decays
    ## along the mean path.
    δ = knots(randn(Xoshiro(1), dims.level), randn(Xoshiro(2), nδ))
    for us in groups, k in 1:K
        @test abs(sum(δ[us, k])) < 1.0e-12
    end
    for z in findall(!, walking), k in 1:K
        @test δ[z, k] ≈ φ^(k - 1) * δ[z, 1] atol = 1.0e-14
    end
    ## One knot takes no innovation draw and is the level alone.
    δ1 = deviation_knots(
        ones(dims.level), Float64[], σ_level, σ_δ, φ, groups, A, Aw,
        walking, 1
    )
    @test δ1 ≈ knots(ones(dims.level), zeros(nδ))[:, 1:1]
    @test_throws DimensionMismatch deviation_knots(
        ones(dims.level), ones(dims.drift), σ_level, σ_δ, φ, groups, A, Aw,
        walking, 1
    )
    @test_throws DimensionMismatch knots(zeros(10), zeros(nδ))
    @test_throws DimensionMismatch knots(zeros(dims.level), zeros(10))
end

@testitem "zone_fit_inputs: units, cells and the walking set" setup = [
    ZoneSynthetic,
] begin
    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    @test inputs.zone_keys ==
        ["a.z1", "a.z2", "a.z3", "a.z4", "a.z5", "b.y1", "b.y2", "b.y3"]
    @test inputs.patch_ranges == [1:5, 6:8]
    @test inputs.patch_of_zone == [1, 1, 1, 1, 1, 2, 2, 2]
    @test inputs.days == syn.days
    @test inputs.t0 == syn.t0
    @test inputs.knots == syn.knots
    @test zd.counts == syn.counts
    @test size(zd.I_bar) == size(syn.I_bar)
    @test zd.I_bar ≈ syn.I_bar rtol = 1.0e-10
    @test zd.g ≈ syn.g
    @test zd.f ≈ syn.f
    ## The multiplier basis is built once, over the patches.
    @test zd.multiplier_basis ==
        BVDOutbreakSize.relative_multiplier_basis(inputs.patch_ranges)
    ## The walking set needs the threshold and at least two zones per patch.
    @test all(
        z -> inputs.walking[z] == (inputs.cumulative[z] >= 30) &&
            (
            !inputs.walking[z] ||
                count(
                inputs.walking[
                    inputs.patch_ranges[
                        inputs.patch_of_zone[z],
                    ],
                ]
            ) >= 2
        ),
        1:syn.nz
    )
    @test zd.n_walking == count(inputs.walking)
    ## Cumulative counts come from the manifest's last vintage.
    @test inputs.cumulative == vec(sum(syn.counts; dims = 2))
end

@testitem "bvd_zone: the truth scores above a perturbed truth" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, relative_multiplier_dims
    using Turing: DynamicPPL

    syn = zone_synthetic()
    ## The synthetic deaths are a rounded tenth of the cumulative cases, not
    ## a draw from the death composition, so the likelihood is scored on
    ## the cases alone: an empty death table leaves no death cell.
    obs_cases = merge(
        syn.obs, (; zone_death_history = Dict{String, Dict{String, NamedTuple}}())
    )
    inputs = zone_fit_inputs(
        syn.chain, obs_cases; zones = nothing, walk_threshold = 30,
        patch_names = ["a", "b"], patch_labels = ["A", "B"]
    )
    zd = inputs.model_data
    @test isempty(zd.death_cell_patch)
    m = bvd_zone(zd)
    K = length(zd.knots)
    truth = syn.truth
    nz = syn.nz
    nc = relative_multiplier_dims(zd.patch_ranges)
    ## Only the zones the inputs mark as walking carry innovations; the
    ## synthetic truth walks every zone, which is the marked set whenever
    ## the counts clear the threshold.
    @test inputs.walking == truth.walking
    ## The counts are multinomial, so the composition is evaluated near its
    ## multinomial limit (`ρ → 0`): at a coarser `ρ` the
    ## Dirichlet-multinomial expects more dispersion than the data carry and
    ## its likelihood is not maximised at the generating shares.
    ## Every other block is pinned at its neutral value, so the only thing
    ## that moves between evaluations is the perturbation.
    at(z_w, z_level, z_drift) = (;
        z_w, σ_level = truth.σ_level, z_level,
        δ_halflife = 42.0, σ_δ = truth.σ_δ, z_drift, ρ = 1.0e-4,
        ρ_death = 1.0e-4, σ_ascertainment = 1.0e-3,
        z_ascertainment = zeros(nc), σ_severity = 1.0e-3,
        z_severity = zeros(nc), η = zeros(zd.meld_d),
    )
    loglik(params) = begin
        vi = DynamicPPL.VarInfo(
            Xoshiro(1), m,
            DynamicPPL.InitFromParams(params)
        )
        DynamicPPL.loglikelihood(m, vi)
    end
    base = loglik(at(truth.z_w, truth.z_level, truth.z_drift))
    @test isfinite(base)
    ## Move one initial share, one level, and the innovations: each costs.
    z_w2 = copy(truth.z_w)
    z_w2[1] -= 0.5
    @test loglik(at(z_w2, truth.z_level, truth.z_drift)) < base
    z_l2 = copy(truth.z_level)
    z_l2[2] += 1.5
    @test loglik(at(truth.z_w, z_l2, truth.z_drift)) < base
    @test loglik(at(truth.z_w, truth.z_level, zero(truth.z_drift))) < base
    ## Profile of the first zone's initial share: the likelihood peaks at
    ## the truth to within the grid.
    grid = -0.6:0.1:0.6
    prof = [
        loglik(
            at(
                truth.z_w .+ [d, 0, 0, 0, 0, 0, 0, 0],
                truth.z_level, truth.z_drift
            )
        ) for d in grid
    ]
    @test abs(grid[argmax(prof)]) <= 0.1
end

@testitem "bvd_zone: prior draws carry the deterministics" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone
    using Turing: sample, Prior
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    nz = syn.nz
    K = length(zd.knots)
    chn = sample(
        bvd_zone(zd), Prior(), 20; chain_type = FlexiChains.VNChain,
        progress = false
    )
    for (q, len) in (
            (:delta_knots_zone, nz * K), (:share_knots_zone, nz * K),
            (:delta_T_zone, nz), (:share_T_zone, nz), (:share_start_zone, nz),
            (:R_T_zone, nz),
            (:region_drift_sd_zone, length(inputs.patch_ranges)),
        )
        @test all(v -> length(v) == len, vec(collect(chn[q])))
    end
    for q in (
            :region_sd_zone, :region_halflife_zone,
            :composition_rho_zone, :composition_rho_death_zone,
        )
        @test all(isfinite, vec(Array(chn[q])))
    end
    ## Shares sum to one within each patch in every draw.
    for v in vec(collect(chn[:share_T_zone])), zs in inputs.patch_ranges

        @test sum(v[zs]) ≈ 1 rtol = 1.0e-10
    end
    ## Deviations sum to zero within each patch at the cut-off.
    for v in vec(collect(chn[:delta_T_zone])), zs in inputs.patch_ranges

        @test abs(sum(v[zs])) < 1.0e-10
    end
    ## Reconstruction reproduces the stored cut-off values.
    shares = reconstruct_zone_shares(chn, inputs)
    rt = reconstruct_zone_rt(chn, inputs)
    sT = [collect(v) for v in vec(collect(chn[:share_T_zone]))]
    rT = [collect(v) for v in vec(collect(chn[:R_T_zone]))]
    w0 = [collect(v) for v in vec(collect(chn[:share_start_zone]))]
    for i in 1:20, z in 1:nz

        @test shares[z][i, inputs.n] ≈ sT[i][z] rtol = 1.0e-10
        @test shares[z][i, 1] ≈ w0[i][z] rtol = 1.0e-10
        ## The chain floors per draw; the reconstruction floors per zone at
        ## the median draw, so only draws finite under both rules compare.
        if isfinite(rT[i][z]) && isfinite(rt[z][i, inputs.n])
            @test rt[z][i, inputs.n] ≈ rT[i][z] rtol = 1.0e-10
        end
    end
    ## Zone infections sum to the draw's own patch trajectory, which under
    ## the stand-in parent is the truth to round-off.
    inf = zone_infections(chn, inputs)
    for i in 1:20, t in 1:inputs.n, (p, zs) in enumerate(inputs.patch_ranges)
        @test sum(inf[z][i, t] for z in zs) ≈ syn.I_bar[p, t] rtol = 1.0e-8
    end
    ## The zone stage's patch trajectories are the ones the zones split.
    pinf = zone_patch_infections(chn, inputs)
    @test length(pinf) == length(inputs.patch_ranges)
    for i in 1:20, t in 1:inputs.n, (p, zs) in enumerate(inputs.patch_ranges)
        @test pinf[p][i, t] ≈ sum(inf[z][i, t] for z in zs) rtol = 1.0e-8
    end
end

@testitem "zone forecast: the fitted model is unchanged and the split is proper" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_forecast
    using Turing: sample, Prior, DynamicPPL
    using Random: Xoshiro
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    inputs_f = zone_inputs(syn; parent_forecast = zone_parent_forecast(syn))
    zd, zdf = inputs.model_data, inputs_f.model_data
    zf = zdf.forecast
    np, n = size(zd.I_bar)
    d, K = zd.meld_d, length(zd.knots)
    nd = n - zd.t0 + 1
    ## The forecast inputs keep every fitted piece as it is.
    @test zf.horizon == 7
    @test zf.meld_L[1:d, 1:d] == zd.meld_L
    @test zf.meld_weights[1:(np * n), 1:d] == zd.meld_weights
    @test all(iszero, zf.meld_weights[1:(np * n), (d + 1):end])
    @test zf.interp[1:nd, 1:K] ≈ zd.interp
    @test all(iszero, zf.interp[1:nd, (K + 1):end])
    @test zf.I_bar[:, 1:n] == zd.I_bar
    @test size(zf.I_bar, 2) == n + 7
    @test all(r -> r == [80, 20], eachrow(zf.totals))
    ## So the fitted model's log density does not change.
    for seed in 1:3
        vi = DynamicPPL.VarInfo(Xoshiro(seed), bvd_zone(zd))
        @test DynamicPPL.logjoint(bvd_zone(zdf), vi) ==
            DynamicPPL.logjoint(bvd_zone(zd), vi)
    end
    ## A forecast of another horizon is refused.
    @test_throws ErrorException zone_inputs(
        syn; parent_forecast = zone_parent_forecast(syn), horizon = 5
    )
    chn = sample(
        bvd_zone(zd), Prior(), 5; chain_type = FlexiChains.VNChain,
        progress = false
    )
    @test_throws ArgumentError zone_forecast(chn, inputs)
    fc = zone_forecast(chn, inputs_f)
    draws = zone_forecast_draws(fc, inputs_f)
    @test size(draws.shares) == (5, syn.nz)
    @test all(isfinite, draws.shares) && all(>=(0), draws.shares)
    @test all(>(0), draws.kappa)
    for (p, zs) in enumerate(inputs_f.patch_ranges), i in 1:5
        @test sum(draws.shares[i, zs]) ≈ 1 rtol = 1.0e-10
        @test sum(draws.zones[z][i] for z in zs) == draws.patches[p][i]
        @test draws.patches[p][i] == (p == 1 ? 80 : 20)
    end
    @test all(v -> all(>=(0), v), draws.zones)
end

@testitem "zone tables: overview, forecast, truth and scores" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_forecast
    using Turing: sample, Prior
    using DataFrames: DataFrame, nrow, names
    using Dates: Day
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    chn = sample(
        bvd_zone(zd), Prior(), 8; chain_type = FlexiChains.VNChain,
        progress = false
    )
    ov = zone_overview_table(chn, inputs)
    @test nrow(ov) == syn.nz
    @test issorted(
        [
            (w ? 2.0 : 0.0) + (isnan(x) ? -1.0 : x)
                for (w, x) in zip(ov.walking, ov.p_R_above_1)
        ];
        rev = true
    )
    @test all(
        z -> isnan(ov.rt_median[z]) ||
            ov.rt_lo90[z] <= ov.rt_median[z] <= ov.rt_hi90[z], 1:syn.nz
    )
    @test Set(ov.zone) == Set(inputs.zone_labels)
    ## The zone forecast drawn from the model over a stand-in parent forecast.
    fcin = zone_inputs(syn; parent_forecast = zone_parent_forecast(syn))
    fc = zone_forecast(chn, fcin)
    ft = zone_forecast_table(fc, fcin)
    @test nrow(ft) == syn.nz + 2
    @test all(ft.lower_90 .<= ft.median .<= ft.upper_90)
    fd = zone_forecast_draws(fc, fcin)
    @test length(fd.zones) == syn.nz && length(fd.patches) == 2
    @test all(sum(fd.zones[z] for z in 1:5) .== fd.patches[1])
    ## Zone forecasts sum to the paired patch total draw by draw, so the
    ## patch rows carry the parent's totals exactly.
    tot = ft[ft.zone .== "Patch total", :]
    @test tot.upper_90 ≈ [80.0, 20.0]
    ## Truth from a later manifest: one more vintage a week on.
    later = deepcopy(syn.obs.zone_confirmed_history)
    made = inputs.cutoff
    for (prov, zones) in later, (nm, h) in zones

        push!(h.days, inputs.n + 7)
        push!(h.counts, h.counts[end] + 3)
    end
    obs2 = merge(syn.obs, (; zone_confirmed_history = later))
    truth = zone_forecast_truth(obs2, inputs; made_date = made)
    @test all(==(3), truth)
    ## Without a vintage at the target the truth is missing.
    truth0 = zone_forecast_truth(syn.obs, inputs; made_date = made)
    @test all(ismissing, truth0)
    vt = zone_forecast_vs_truth(fc, fcin; truth)
    @test nrow(vt) == syn.nz + 2
    @test all(vt.observed[vt.zone .!= "Patch total"] .== 3)
    @test all(vt.lower_90 .<= vt.lower_50 .<= vt.upper_50 .<= vt.upper_90)
    sc = zone_forecast_scores(fc, fcin; truth)
    @test Set(sc.method) ==
        Set(["zone model", "share persistence", "naive persistence"])
    @test all(isfinite, sc.log_score)
    @test nrow(sc) == 6
    ## The patch total carries every `score_draws` column. Share
    ## persistence forecasts no total, so its columns are missing.
    for c in (
            :crps, :log_crps, :dispersion, :overprediction,
            :underprediction, :coverage_50, :coverage_90, :bias, :n,
        )
        @test c in propertynames(sc)
    end
    modelled = sc[sc.method .== "zone model", :]
    @test all(isfinite, modelled.crps) && all(isfinite, modelled.bias)
    @test all(>(0), modelled.n)
    @test all(x -> x isa Bool, modelled.coverage_90)
    @test all(ismissing, sc[sc.method .== "share persistence", :crps])
    @test all(ismissing, sc[sc.method .== "share persistence", :coverage_90])
    cd = zone_composition_draws(chn, inputs)
    @test length(cd.expected) == syn.nz
    @test size(cd.expected[1]) == (8, length(inputs.days))
    for v in eachindex(inputs.days), zs in inputs.patch_ranges

        N = sum(syn.counts[zs, v])
        N > 0 || continue
        @test all(sum(cd.predictive[z][i, v] for z in zs) == N for i in 1:8)
        @test all(sum(cd.expected[z][i, v] for z in zs) ≈ 1 for i in 1:8)
    end
    cdc = zone_composition_draws(chn, inputs; cumulative = true)
    for zs in inputs.patch_ranges, i in 1:8

        @test cdc.predictive[zs[1]][i, :] == cumsum(cd.predictive[zs[1]][i, :])
        v = length(inputs.days)
        @test sum(cdc.predictive_share[z][i, v] for z in zs) ≈ 1
        @test sum(cdc.expected[z][i, v] for z in zs) ≈ 1
        @test sum(cdc.observed[z, v] for z in zs) ≈ 1
    end
    cal = zone_composition_calibration(chn, inputs)
    @test names(cal) == [
        "Stream", "Vintages", "Bias", "50% coverage",
        "90% coverage",
    ]
    @test nrow(cal) == 3 && cal.Stream[end] == "All zones"
    @test all(0 .<= cal[!, "90% coverage"] .<= 1)
    ## At the cut-off the rebuilt reproduction number is the chain's own
    ## wherever the chain reports it in every draw. The floor is decided
    ## per zone, so a reported zone carries every draw.
    rt = reconstruct_zone_rt(chn, inputs)
    R_chain = [collect(v) for v in vec(collect(chn[:R_T_zone]))]
    for z in 1:syn.nz
        col = rt[z][:, end]
        @test all(isnan, col) || !any(isnan, col)
        r = [R_chain[i][z] for i in 1:8]
        any(isnan, r) && continue
        @test col ≈ r rtol = 1.0e-8
    end
    ## Both compositions, each read on its own observed counts.
    for stream in (:cases, :deaths)
        ppc = zone_composition_ppc(
            chn, inputs; top = 3, prior_chain = chn, stream
        )
        @test nrow(ppc) > 0
        @test all(0 .<= ppc.observed_share .<= 1)
        @test all(ppc.lower_90 .<= ppc.median .<= ppc.upper_90)
        @test ppc.prior_median == ppc.median
        @test length(unique(ppc.zone)) <= 6
        figs = plot_zone_composition_ppc(
            chn, inputs; prior_chain = chn, top = 2, stream
        )
        @test length(figs) == 2
    end
    ## The deaths check is not a relabelling of the cases one: the fixture
    ## deaths are a tenth of the cases, so the counts differ.
    @test zone_composition_ppc(chn, inputs; top = 3).observed !=
        zone_composition_ppc(chn, inputs; top = 3, stream = :deaths).observed
    arch = zone_forecast_archive(fc, fcin; made_date = made, thin = 5)
    @test Set(propertynames(arch)) == Set(
        [
            :made_date, :horizon, :target_date,
            :province, :zone, :stream, :draw, :value, :method,
        ]
    )
    @test all(arch.method .== ZONE_FORECAST_METHOD)
    @test nrow(arch) == syn.nz * length(1:5:8)
    @test all(arch.target_date .== made + Day(7))
end

@testitem "zone_binned_increments: the explicit window sums" setup = [
    ZoneSynthetic,
] begin
    syn = zone_synthetic()
    t0, nz = syn.t0, syn.nz
    n = size(syn.I_bar, 2)
    nd = n - t0 + 1
    fixed = zone_fixed_terms(syn.I_bar, syn.g, syn.f, t0)
    infections = rand(Xoshiro(2), nd, nz)
    w0 = syn.truth.w0
    patch_of_zone = [p for (p, zs) in enumerate(syn.patch_ranges) for _ in zs]
    F = reference_delay_matrix(syn.f, nd)
    pre_rows = reference_pre_rows(fixed.report_pre, patch_of_zone, t0, n)
    reports = F * infections .+ pre_rows .* transpose(w0)
    ## Each window `(d_{v−1}, d_v]` sums the pre-`t0` patch term at the
    ## initial share and the zone's own daily reports from `t0`.
    explicit(days, v) = begin
        lo = v == 1 ? 1 : days[v - 1] + 1
        hi = days[v]
        [
            let p = patch_of_zone[z]
                w0[z] * sum(
                    (fixed.report_pre[p, t] for t in lo:hi if t < t0);
                    init = 0.0
                ) +
                    sum(
                    (reports[t - t0 + 1, z] for t in lo:hi if t >= t0);
                    init = 0.0
                )
            end
                for z in 1:nz
        ]
    end
    ## The fit's windows (the first opens before `t0`), windows entirely
    ## before `t0`, and windows that open after it.
    for days in (syn.days, [t0 - 10, t0 - 3, t0 + 5, n], [t0 + 3, t0 + 20, n])
        op = zone_binned_operator(syn.f, fixed, days, t0, n, patch_of_zone)
        C = zone_binned_increments(op, infections, w0)
        @test size(C) == (nz, length(days))
        for v in eachindex(days)
            @test C[:, v] ≈ explicit(days, v) rtol = 1.0e-12
        end
    end
end

@testitem "zone_parent_inputs: the mean over the parent draws" setup = [
    ZoneSynthetic,
] begin
    using Statistics: mean

    syn = zone_synthetic()
    ## Four draws at different scales of the truth.
    scales = [0.5, 1.0, 2.0, 1.5]
    chain = copy(syn.chain)
    chain[:infections_patch] = reshape(
        [s .* vec(syn.I_bar) for s in scales],
        4, 1
    )
    chain[:C_T] = reshape(scales .* sum(syn.I_bar), 4, 1)
    ## The patch trajectory the zone stage conditions on is the geometric
    ## mean over the parent draws; their spread reaches the fit through the
    ## shared quantity instead.
    avg = zone_parent_inputs(chain)
    @test exp.(avg.log_infections) ≈
        exp(mean(log.(scales))) .* vec(syn.I_bar) rtol = 1.0e-10
    inputs = zone_fit_inputs(
        chain, syn.obs;
        zones = nothing, patch_names = ["a", "b"], patch_labels = ["A", "B"]
    )
    @test inputs.model_data.I_bar ≈
        exp(mean(log.(scales))) .* syn.I_bar rtol = 1.0e-10
    ## A parent carrying the scales: each scale prior is the log-normal
    ## matched to its draws.
    with_priors = copy(chain)
    with_priors[:province_ascertainment_sd] = reshape(
        [0.1, 0.2, 0.3, 0.4], 4, 1
    )
    with_priors[:region_drift_sd] = reshape([0.04, 0.05, 0.06, 0.07], 4, 1)
    p = zone_parent_inputs(with_priors).priors
    @test p.ascertainment_sd[1] ≈ mean(log.([0.1, 0.2, 0.3, 0.4]))
    @test p.drift_sd[1] ≈ mean(log.([0.04, 0.05, 0.06, 0.07]))
    @test all(>(0), (p.ascertainment_sd[2], p.drift_sd[2]))
    ## The movement is centred as the infections are: the geometric mean of
    ## the intensity and of the odds that an infection was imported.
    spread = zone_spread_chain(syn)
    mv = zone_parent_inputs(spread)
    eps_draws = [collect(v) for v in vec(spread[:importation_epsilon_patch])]
    @test mv.origin_epsilon ≈ exp.(mean(e -> log.(e), eps_draws))
    odds = [0.05 * s / (1 - 0.05 * s) for s in [0.5, 1.0, 2.0, 1.5]]
    lo = reshape(mv.import_log_odds, 2, :)
    @test all(lo[2, :] .≈ mean(log.(odds)))
    @test all(lo[1, :] .< -700)
    ## The delays are the mean of the per-draw PMFs, and each draw's log
    ## parameters are carried for the meld.
    gi_means = [12.0, 14.0, 15.0, 17.0]
    varied = copy(chain)
    varied[Symbol("gi_state.gi_mean")] = reshape(gi_means, 4, 1)
    vd = zone_parent_inputs(varied)
    per_draw = [
        zone_delay_pmfs(
            [m; syn.delay_params[2:end]]
        ).g for m in gi_means
    ]
    @test vd.g ≈ mean(per_draw) rtol = 1.0e-12
    @test vd.delay_log_draws[:, 1] ≈ log.(gi_means)
    @test vd.delay_log_draws[1, 2:end] ≈ log.(syn.delay_params[2:end])
    @test size(vd.delay_log_draws) == (4, 12)
end

@testitem "zone_delay_pmfs: the joint's delays, onset to report included" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: ZONE_GI_NMAX, ZONE_INCUBATION_NMAX,
        ZONE_REPORT_NMAX, ZONE_RECEIPT_NMAX, ZONE_ONSET_DEATH_NMAX
    using Distributions: Gamma

    syn = zone_synthetic()
    p = syn.delay_params
    pm = zone_delay_pmfs(p)
    disc(d, nmax) = discretise_censored(d, nmax)
    inc = disc(lognormal_meansd(p[3], p[4]), ZONE_INCUBATION_NMAX)
    report = disc(Gamma(p[5], p[6]), ZONE_REPORT_NMAX)
    receipt = disc(lognormal_meansd(p[7], p[8]), ZONE_RECEIPT_NMAX)
    ## The case delay is the joint's onset-to-confirmation kernel after the
    ## incubation period: incubation, onset to report, then receipt.
    @test pm.f ≈ convolve_pmf(convolve_pmf(inc, report), receipt)
    ## The death delay is onset to death, the two atomic delays truncated
    ## and renormalised, then receipt.
    od = convolve_pmf(
        disc(Gamma(p[9], p[10]), ZONE_ONSET_DEATH_NMAX),
        disc(Gamma(p[11], p[12]), ZONE_ONSET_DEATH_NMAX)
    )[1:(ZONE_ONSET_DEATH_NMAX + 1)]
    @test pm.death_pmf ≈
        convolve_pmf(convolve_pmf(inc, od ./ sum(od)), receipt)
    ## The generation interval drops lag 0 and renormalises.
    gi = disc(Gamma((p[1] / p[2])^2, p[2]^2 / p[1]), ZONE_GI_NMAX)
    @test pm.g ≈ gi[2:end] ./ sum(gi[2:end])
    @test sum(pm.g) ≈ 1
    @test_throws DimensionMismatch zone_delay_pmfs(p[1:11])
end

@testitem "zone_fit_inputs: refuses inconsistent inputs" setup = [
    ZoneSynthetic,
] begin
    syn = zone_synthetic()
    kw = (;
        zones = nothing, patch_names = ["a", "b"],
        patch_labels = ["A", "B"],
    )
    ## Observations without zone tables.
    bare = (;
        n = syn.obs.n, seeding = syn.obs.seeding,
        cutoff = syn.obs.cutoff,
    )
    @test_throws ErrorException zone_fit_inputs(syn.chain, bare; kw...)
    ## A parent without the patch structure or a delay it needs.
    for key in (:infections_patch, Symbol("gi_state.gi_mean"))
        chain = copy(syn.chain)
        delete!(chain, key)
        @test_throws ErrorException zone_fit_inputs(chain, syn.obs; kw...)
    end
    ## A parent fitted to a shorter cut-off.
    chain = copy(syn.chain)
    chain[:infections_patch] = reshape(
        [vec(syn.I_bar[:, 1:(end - 5)]) for _ in 1:4], 4, 1
    )
    @test_throws ErrorException zone_fit_inputs(chain, syn.obs; kw...)
    ## Patches on different vintage days. `zone_increment_matrix` holds one
    ## day vector over every patch, so it is the refusal the caller sees
    ## rather than the per-patch guard that follows it.
    hist = deepcopy(syn.obs.zone_confirmed_history)
    for (nm, h) in hist["b"]
        hist["b"][nm] = (; days = h.days .- 1, counts = h.counts)
    end
    obs = merge(syn.obs, (; zone_confirmed_history = hist))
    @test_throws ErrorException zone_fit_inputs(syn.chain, obs; kw...)
    days_err = try
        zone_fit_inputs(syn.chain, obs; kw...)
    catch e
        e
    end
    @test occursin("same vintages", days_err.msg)
    @test occursin("days to `a.", days_err.msg)
    ## No zones at all.
    empty = merge(
        syn.obs,
        (; zone_confirmed_history = Dict{String, Dict{String, NamedTuple}}())
    )
    @test_throws ErrorException zone_fit_inputs(syn.chain, empty; kw...)
end

@testitem "zone_fit_inputs: a reattribution vintage is left out" setup = [
    ZoneSynthetic,
] begin
    syn = zone_synthetic()
    hist = deepcopy(syn.obs.zone_confirmed_history)
    nv = length(syn.days)
    ## Patch a carries 20 unallocated cases over the first three vintages,
    ## attributed to z1 at the fourth.
    un = vcat(fill(20, 3), zeros(Int, nv - 3))
    hist["a"]["unallocated"] = (; days = copy(syn.days), counts = un)
    z1 = hist["a"]["z1"]
    hist["a"]["z1"] = (; days = z1.days, counts = z1.counts .+ (20 .- un))
    obs = merge(syn.obs, (; zone_confirmed_history = hist))
    inputs = zone_fit_inputs(
        syn.chain, obs; zones = nothing,
        patch_names = ["a", "b"], patch_labels = ["A", "B"]
    )
    zd = inputs.model_data
    day4 = syn.days[4]
    @test zone_reattribution_days(hist) == Dict("a" => [day4])
    @test inputs.excluded ==
        [(; patch = "a", date = obs.seeding + Day(day4 - 1))]
    ## The fourth column is zeroed for patch a only, so its cell is not
    ## scored, while every other column is the plain difference.
    @test all(iszero, inputs.counts[1:5, 4])
    @test inputs.counts[6:8, 4] == syn.counts[6:8, 4]
    rest = setdiff(1:nv, 4)
    @test inputs.counts[:, rest] == syn.counts[:, rest]
    @test !any(
        c -> zd.cell_patch[c] == 1 && zd.cell_vintage[c] == 4,
        eachindex(zd.cell_patch)
    )
    @test any(
        c -> zd.cell_patch[c] == 2 && zd.cell_vintage[c] == 4,
        eachindex(zd.cell_patch)
    )
    ## The cumulative at the cut-off still carries the reattributed cases.
    @test inputs.cumulative[1] == sum(syn.counts[1, :]) + 20
    @test isempty(zone_inputs(syn).excluded)
    ## A death moved out of a named zone leaves the unallocated row flat,
    ## so the death composition drops that vintage on the widened rule
    ## while the case composition keeps it.
    deaths = deepcopy(syn.obs.zone_death_history)
    z2 = deaths["a"]["z2"]
    rising = collect(0:2:(2 * (nv - 1)))
    rising[5:end] .-= 6
    deaths["a"]["z2"] = (; z2.days, counts = rising)
    obs2 = merge(syn.obs, (; zone_death_history = deaths))
    zd2 = zone_fit_inputs(
        syn.chain, obs2; zones = nothing,
        patch_names = ["a", "b"], patch_labels = ["A", "B"]
    ).model_data
    @test any(!iszero, zone_inputs(syn).model_data.death_counts[1:5, 5])
    @test all(iszero, zd2.death_counts[1:5, 5])
    @test !any(
        c -> zd2.death_cell_patch[c] == 1 && zd2.death_cell_vintage[c] == 5,
        eachindex(zd2.death_cell_patch)
    )
    @test zd2.counts == zone_inputs(syn).counts
end

@testitem "zone mixing: the blocks, and what they conserve" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_importation_blocks, _zone_states,
        zone_deformation, zone_draw_mixing, destination_weighted_kernel,
        mobility_importation_kernel
    using Turing: sample, Prior
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn; zones = zone_metadata(syn))
    zd = inputs.model_data
    mix = zd.mixing
    @test mix !== nothing
    poz = inputs.patch_of_zone
    ## The within block is column-stochastic inside the origin's patch and
    ## empty across patches; the between block is the other way round.
    for q in 1:syn.nz
        zs = findall(==(poz[q]), poz)
        @test sum(mix.within[zs, q]) ≈ 1 rtol = 1.0e-12
        @test mix.within[q, q] == 0
        @test all(iszero, mix.within[setdiff(1:syn.nz, zs), q])
        @test all(iszero, mix.between[zs, q])
    end
    ## Summed over a destination patch's zones the between block is the
    ## province model's own patch-to-patch flow, so the same movement is
    ## not counted at both levels.
    parent_kernel = mobility_importation_kernel(PROVINCE_POPULATIONS[1:2])
    for q in 1:syn.nz, p in 1:2

        p == poz[q] && continue
        zs = findall(==(p), poz)
        @test sum(mix.between[zs, q]) ≈ parent_kernel[p, poz[q]] rtol = 1.0e-12
    end
    ## Arrivals stay a proper fraction of a patch's own infections.
    @test all(0 .<= zone_draw_mixing(mix, nothing, nothing).import_fraction .< 1)
    ## A parent's log weights reach the blocks through
    ## `destination_weighted_kernel`: a malformed one throws there rather
    ## than being skipped.
    function tilted(flow)
        chain = copy(syn.chain)
        chain[:importation_flow_effect] = reshape([flow for _ in 1:4], 4, 1)
        return zone_inputs(
            merge(syn, (; chain)); zones = zone_metadata(syn)
        ).model_data.mixing
    end
    @test_throws DimensionMismatch tilted(zeros(3))
    ## With two patches each origin has one destination, so any weighting
    ## cancels and the blocks are the unweighted ones.
    @test tilted([0.0, 0.4, -0.7, 0.0]).between ≈ mix.between
    ## The parent kernel is the posterior mean of the weighted kernels over
    ## draws, not the kernel at the mean log weight.
    Ws = [
        [0.0 1.5 -0.4; 2.0 0.0 0.3; -1.0 0.2 0.0],
        [0.0 -1.2 0.8; -0.5 0.0 1.1; 0.9 -0.7 0.0],
        [0.0 0.3 0.1; 0.4 0.0 -2.0; 1.6 0.5 0.0],
        [0.0 -0.2 -1.5; 1.0 0.0 0.0; -0.3 1.8 0.0],
    ]
    chain3 = copy(syn.chain)
    chain3[:importation_flow_effect] = reshape(vec.(Ws), 4, 1)
    K3p = province_importation_kernel(PROVINCE_POPULATIONS[1:3])
    mean_kernel = sum(destination_weighted_kernel(K3p, W) for W in Ws) ./ 4
    got = BVDOutbreakSize._mean_parent_kernel(chain3, :importation_flow_effect)
    @test got ≈ mean_kernel
    @test !(got ≈ destination_weighted_kernel(K3p, sum(Ws) ./ 4))
    ## With three patches a non-zero weighting changes the between block,
    ## and its column over a destination patch's zones is the weighted
    ## kernel's entry.
    K3 = province_importation_kernel(PROVINCE_POPULATIONS[1:3])
    η3 = [0.7, -0.2, -0.5]
    poz3 = [1, 1, 2, 2, 3, 3]
    pops3 = [5.0e5, 2.0e5, 3.0e5, 1.0e5, 4.0e5, 1.5e5]
    coords3 = [
        (-4.3, 15.3), (-4.6, 15.8), (-5.9, 22.4),
        (-6.2, 23.6), (0.5, 25.2), (1.2, 24.8),
    ]
    Kw3 = destination_weighted_kernel(K3, η3)
    plain3 = zone_importation_blocks(pops3, coords3, poz3, K3).between
    tilt3 = zone_importation_blocks(pops3, coords3, poz3, Kw3).between
    @test !(tilt3 ≈ plain3)
    for q in eachindex(poz3), p in 1:3

        p == poz3[q] && continue
        zs = findall(==(p), poz3)
        @test sum(tilt3[zs, q]) ≈ Kw3[p, poz3[q]] rtol = 1.0e-12
    end
    chn = sample(
        bvd_zone(zd), Prior(), 6;
        chain_type = FlexiChains.VNChain, progress = false
    )
    eps = [collect(v) for v in vec(collect(chn[:mixing_epsilon_zone]))]
    @test all(v -> length(v) == syn.nz && all(0 .< v .< 1), eps)
    ## The states read the fractions, and the mixed shares still sum to one
    ## within each patch but differ from the unmixed ones.
    states = _zone_states(chn, inputs)
    for (i, st) in enumerate(states)
        @test st.ε == eps[i]
        def = zone_deformation(zd, nothing)
        fw = zone_forward(zd, st.δ_knots, st.w0, st.ε, def)
        plain = zone_forward(zd, st.δ_knots, st.w0, nothing, def).shares
        for zs in inputs.patch_ranges
            @test all(sum(fw.shares[:, zs]; dims = 2) .≈ 1)
        end
        @test !(fw.shares ≈ plain)
        ## Every zone's infections still sum to the patch total the shared
        ## quantity gave, so mixing adds no infection and loses none.
        for (p, zs) in enumerate(inputs.patch_ranges), j in 1:size(fw.infections, 1)

            @test sum(fw.infections[j, zs]) ≈
                def.I_bar[p, inputs.t0 + j - 1] rtol = 1.0e-10
        end
        ## Some infection crosses a patch boundary.
        @test sum(fw.imports) > 0
    end
    ## Without metadata for every zone there is no mixing structure and the
    ## model carries no mixing parameters.
    zd0 = zone_inputs(syn).model_data
    @test zd0.mixing === nothing
    chn0 = sample(
        bvd_zone(zd0), Prior(), 2;
        chain_type = FlexiChains.VNChain, progress = false
    )
    @test !BVDOutbreakSize._has_key(chn0, :mixing_epsilon_zone)
    ## Metadata that misses one zone leaves both the mixing structure and
    ## the distance blocks off, since neither can be built for part of a
    ## patch.
    partial = zone_inputs(syn; zones = zone_metadata(syn)[2:end]).model_data
    @test partial.mixing === nothing
    @test isempty(partial.zone_distances)
    @test isempty(partial.zone_walk_distances)
    ## Metadata alone is not enough: the parent chain must also carry the
    ## province model's between-patch movement.
    still = copy(syn.chain)
    delete!(still, :importation_epsilon_patch)
    @test zone_inputs(
        merge(syn, (; chain = still)); zones = zone_metadata(syn)
    ).model_data.mixing === nothing
    flat = copy(syn.chain)
    delete!(flat, :importation_patch)
    @test zone_inputs(
        merge(syn, (; chain = flat)); zones = zone_metadata(syn)
    ).model_data.mixing === nothing
end

@testitem "zone mixing: the sampled decay and destination weights" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_gravity_blocks, gravity_pull,
        province_distance_matrix, _zone_states
    using Turing: sample, Prior
    using Statistics: median, quantile
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn; zones = zone_metadata(syn))
    zd = inputs.model_data
    mix = zd.mixing
    poz = inputs.patch_of_zone
    nz = syn.nz
    meta = zone_metadata(syn)
    pops = Float64[r.population for r in meta]
    coords = [(r.lat, r.lon) for r in meta]
    K = province_importation_kernel(PROVINCE_POPULATIONS[1:2])
    ## The gravity blocks written out from the pull at decay `γ`, each
    ## destination weighted by `exp(ω_z)`.
    function reference(γ, ω = zeros(nz))
        pull = gravity_pull(
            pops; distances = province_distance_matrix(coords), decay = γ
        ) .* exp.(ω)
        within = zeros(nz, nz)
        between = zeros(nz, nz)
        for q in 1:nz, z in 1:nz

            zs = findall(==(poz[z]), poz)
            s = sum(pull[zs, q])
            s > 0 || continue
            if poz[z] == poz[q]
                within[z, q] = pull[z, q] / s
            else
                between[z, q] = K[poz[z], poz[q]] * pull[z, q] / s
            end
        end
        return (; within, between)
    end
    ## Unit decay and no destination weight recover the fixed blocks.
    blocks(γ, ω) = zone_gravity_blocks(mix.gravity, γ, ω)
    b1 = blocks(1.0, zeros(nz))
    @test b1.within ≈ mix.within rtol = 1.0e-12
    @test b1.between ≈ mix.between rtol = 1.0e-12
    @test b1.within ≈ reference(1.0).within rtol = 1.0e-12
    @test b1.between ≈ reference(1.0).between rtol = 1.0e-12
    ## Another decay is the gravity pull at that decay.
    @test blocks(2.5, zeros(nz)).within ≈ reference(2.5).within rtol = 1.0e-12
    ## A weight shared within a patch cancels.
    ω_flat = [fill(0.7, 5); fill(-0.7, 3)]
    @test blocks(1.0, ω_flat).within ≈ mix.within
    @test blocks(1.0, ω_flat).between ≈ mix.between
    ## At any decay and weight the within block stays column-stochastic in
    ## the origin's patch and the between block's column over a destination
    ## patch is still the parent's flow; only the split moves.
    ω = [1.2, -0.4, 0.3, -0.9, -0.2, 0.8, -1.1, 0.3]
    for γ in (0.05, 0.6, 4.0)
        b = blocks(γ, ω)
        @test b.within ≈ reference(γ, ω).within rtol = 1.0e-10
        @test b.between ≈ reference(γ, ω).between rtol = 1.0e-10
        for q in 1:nz, p in 1:2

            zs = findall(==(p), poz)
            if p == poz[q]
                @test sum(b.within[zs, q]) ≈ 1 rtol = 1.0e-12
                @test b.within[q, q] == 0
                @test all(iszero, b.between[zs, q])
            else
                @test sum(b.between[zs, q]) ≈ K[p, poz[q]] rtol = 1.0e-12
                @test all(iszero, b.within[zs, q])
            end
        end
        @test !(b.within ≈ mix.within)
    end
    ## The model records the decay, the weight scale and the log weights,
    ## which sum to zero within every patch, and the states rebuild each
    ## draw's blocks from them.
    chn = sample(
        bvd_zone(zd), Prior(), 4;
        chain_type = FlexiChains.VNChain, progress = false
    )
    γs = vec(collect(chn[:mixing_decay_zone]))
    ## The decay prior keeps its median at the conventional 1 and its 95%
    ## point above 3.7, the upper end of the decay a 63-zone fit estimated.
    decay_prior = bvd_zone(zd).defaults.mixing_decay_prior
    @test median(decay_prior) ≈ 1
    @test quantile(decay_prior, 0.95) > 3.7
    ωs = [collect(v) for v in vec(collect(chn[:mixing_destination_zone]))]
    @test all(>(0), γs)
    @test all(>(0), vec(collect(chn[:mixing_destination_sd_zone])))
    for v in ωs, zs in inputs.patch_ranges
        @test abs(sum(v[zs])) < 1.0e-10
    end
    states = _zone_states(chn, inputs)
    for (i, st) in enumerate(states)
        b = blocks(γs[i], ωs[i])
        @test st.def.mixing.within ≈ b.within
        @test st.def.mixing.between ≈ b.between
    end
end

@testitem "bvd_zone: the two compositions sum to the likelihood" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_composition_logpdf, _zone_kappa,
        relative_multiplier_dims
    using Turing: DynamicPPL

    syn = zone_synthetic()
    ## Cumulative allocated deaths per zone, a tenth of the confirmed.
    deaths = Dict{String, Dict{String, NamedTuple}}()
    for (prov, zones) in syn.obs.zone_confirmed_history
        deaths[prov] = Dict{String, NamedTuple}()
        for (nm, h) in zones
            deaths[prov][nm] = (;
                days = copy(h.days),
                counts = cld.(h.counts, 10),
            )
        end
    end
    obs = merge(syn.obs, (; zone_death_history = deaths))
    inputs = zone_fit_inputs(
        syn.chain, obs; zones = nothing,
        patch_names = ["a", "b"], patch_labels = ["A", "B"]
    )
    zd = inputs.model_data
    ## Deaths are scored per vintage, on the same grid as the cases, and
    ## the rows sum to the cumulative the table ends on.
    @test size(zd.death_counts) == size(zd.counts)
    @test size(zd.death_bin.weights, 2) == length(zd.days)
    @test vec(sum(zd.death_counts; dims = 2)) == cld.(inputs.cumulative, 10)
    @test length(zd.death_cell_patch) == length(zd.death_cell_vintage)
    @test length(zd.death_cell_patch) > 2
    @test all(>(0), zd.death_cell_total)
    truth = syn.truth
    nc = relative_multiplier_dims(zd.patch_ranges)
    ## Zero contrasts leave both relative multipliers at one, so the
    ## likelihood is the two compositions at the modelled increments alone.
    params = (;
        z_w = truth.z_w, σ_level = truth.σ_level,
        z_level = truth.z_level, δ_halflife = 42.0, σ_δ = truth.σ_δ,
        z_drift = truth.z_drift, ρ = 0.05, ρ_death = 0.05,
        z_ascertainment = zeros(nc), z_severity = zeros(nc),
    )
    m = bvd_zone(zd)
    total = DynamicPPL.loglikelihood(
        m,
        DynamicPPL.VarInfo(Xoshiro(1), m, DynamicPPL.InitFromParams(params))
    )
    @test isfinite(total)
    fw = zone_forward(zd, truth.δ_knots, truth.w0, nothing)
    cases = zone_composition_logpdf(
        zd.counts, fw.increments, zd.cell_patch,
        zd.cell_vintage, zd.cell_total, zd.cell_const,
        zd.patch_ranges, _zone_kappa(0.05)
    )
    D = zone_binned_increments(zd.death_bin, fw.infections, truth.w0)
    extra = zone_composition_logpdf(
        zd.death_counts, D, zd.death_cell_patch,
        zd.death_cell_vintage, zd.death_cell_total, zd.death_cell_const,
        zd.patch_ranges, _zone_kappa(0.05)
    )
    @test total ≈ cases + extra rtol = 1.0e-8
    ## Without allocated zone deaths the fit refuses to run.
    bare = merge(
        syn.obs,
        (; zone_death_history = Dict{String, Dict{String, NamedTuple}}())
    )
    @test_throws ErrorException fit_zone(
        syn.chain, bare;
        zones = nothing, patch_names = ["a", "b"], patch_labels = ["A", "B"]
    )
end

@testitem "bvd_zone: the sampled dimension and the optional blocks" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, relative_multiplier_dims,
        deviation_knot_dims
    using Turing: DynamicPPL
    import FlexiChains
    using Turing: sample, Prior

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    K = length(zd.knots)
    np = length(inputs.patch_ranges)
    nc = relative_multiplier_dims(zd.patch_ranges)
    dimension(m) = length(DynamicPPL.link(DynamicPPL.VarInfo(m), m)[:])
    ## The initial-share draws over the zones, the level and the two
    ## multipliers' contrasts within each patch, the innovations on the
    ## walking zones' contrasts, the six scalar scales, one drift scale per
    ## patch and the shared draw.
    nd = deviation_knot_dims(zd.patch_ranges, zd.walking).drift
    dim = syn.nz + 3 * nc + nd * (K - 1) + 6 + np + zd.meld_d
    @test dimension(bvd_zone(zd)) == dim
    ## Mixing adds the within-patch intensity, a departure scale and one
    ## offset per zone, the zone decay, the destination weight scale and its
    ## contrasts within each patch, where the inputs carry the kernel, and one
    ## origin intensity per patch to the shared draw. The same metadata gives
    ## the zones centroids, so the correlation block switches on too and adds
    ## its reference correlation.
    zdm = zone_inputs(syn; zones = zone_metadata(syn)).model_data
    @test dimension(bvd_zone(zdm)) == dim + 5 + syn.nz + nc + np
    ## With no kept meld cell the shared draw is absent.
    zdc = merge(
        zd, (;
            meld_d = 0, meld_weights = zeros(0, 0), meld_L = zeros(0, 0),
        )
    )
    @test dimension(bvd_zone(zdc)) == dim - zd.meld_d
    ## With no walking zone the innovation block is absent from the model.
    inputs0 = zone_inputs(syn; walk_threshold = 10^6)
    zd0 = inputs0.model_data
    @test zd0.n_walking == 0
    model0 = bvd_zone(zd0)
    @test dimension(model0) == syn.nz + 3 * nc + 6 + np + zd0.meld_d
    chn0 = sample(
        model0, Prior(), 3; chain_type = FlexiChains.VNChain,
        progress = false
    )
    @test !any(p -> string(p) == "z_drift", FlexiChains.parameters(chn0))
    @test all(v -> length(v) == syn.nz, vec(collect(chn0[:delta_T_zone])))
end

@testitem "zone diagnostics: tables from a short chain" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone
    using Turing: sample, Prior, MCMCSerial
    using DataFrames: nrow, names
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    ## This seed leaves zones on both sides of `rt_floor`.
    chn = sample(
        Xoshiro(2), bvd_zone(inputs.model_data), Prior(), MCMCSerial(), 10, 2;
        chain_type = FlexiChains.VNChain, progress = false
    )
    diag = zone_diagnostics_table(chn, inputs)
    @test nrow(diag) == syn.nz
    @test diag.zone == inputs.zone_labels
    @test diag.walking == inputs.walking
    for stem in ("R_T", "share_T", "delta_T"), stat in (
                "rhat", "ess_bulk",
                "ess_tail",
            )

        @test "$(stat)_$(stem)" in names(diag)
    end
    @test all(isfinite, diag.rhat_share_T)
    @test all(>(0), diag.ess_bulk_share_T)
    ## Without inputs the zones are numbered and the walking column absent.
    bare = zone_diagnostics_table(chn)
    @test bare.zone == string.(1:syn.nz)
    @test !("walking" in names(bare))
    @test bare.rhat_share_T == diag.rhat_share_T
    ## A zone undefined in some draw is marked in the table, not warned
    ## about once per zone, and summarising `R_T_zone` apart leaves every
    ## value as the whole-chain summary gives it.
    @test any(v -> any(!isfinite, v), vec(collect(chn[:R_T_zone])))
    @test_logs min_level = Base.CoreLogging.Warn zone_diagnostics_table(
        chn, inputs
    )
    @test_logs min_level = Base.CoreLogging.Warn zone_sampler_diagnostics(
        chn, inputs
    )
    whole = FlexiChains.rhat(chn; warn = false)
    whole_rt = BVDOutbreakSize._zone_summary_vector(whole, :R_T_zone, syn.nz)
    defined = .!isnan.(diag.rhat_R_T)
    @test any(defined)
    @test diag.rhat_R_T[defined] ≈ whole_rt[defined]
    @test diag.rhat_share_T ≈
        BVDOutbreakSize._zone_summary_vector(whole, :share_T_zone, syn.nz)
    ## Sampler statistics are absent from a prior chain, so the per-chain
    ## fields are empty and the divergences zero.
    sd = zone_sampler_diagnostics(chn, inputs)
    @test isfinite(sd.max_rhat)
    @test sd.min_ess_bulk > 0 && sd.min_ess_tail > 0
    @test sd.n_divergent == 0
    @test isempty(sd.depth_cap_fraction) && isempty(sd.ebfmi) &&
        isempty(sd.step_size)
    @test sd.max_rhat_R_T_walking isa Float64
    @test isnan(zone_sampler_diagnostics(chn).max_rhat_R_T_walking)
end

@testitem "zone_meld_check: one row per shared cell" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone
    using Turing: sample, Prior
    using DataFrames: nrow
    using Dates: Day
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    meld = inputs.meld
    chn = sample(
        bvd_zone(zd), Prior(), 16; chain_type = FlexiChains.VNChain,
        progress = false
    )
    mc = zone_meld_check(chn, inputs)
    @test mc.dimension == meld.d
    ## The infection cells have a row each; the delay cells have none.
    n_cells = length(meld.cells_patch)
    @test nrow(mc.table) == n_cells
    @test mc.table.patch ==
        [inputs.patch_labels[p] for p in meld.cells_patch]
    @test mc.table.midpoint == [meld.midpoints[k] for k in meld.cells_week]
    @test mc.table.date ==
        [inputs.seeding + Day(m - 1) for m in mc.table.midpoint]
    @test mc.table.parent_sd ≈
        sqrt.(vec(sum(abs2, meld.L[1:n_cells, :]; dims = 2)))
    @test all(isfinite, mc.table.eta_mean)
    @test all(>(0), mc.table.eta_sd)
    @test isfinite(mc.mean_norm_sq) && mc.mean_norm_sq > 0
    ## With no kept cell the draw is empty, and a chain that carries no
    ## shared draw at all reads the same way.
    zdc = merge(
        zd, (;
            meld_d = 0, meld_weights = zeros(0, 0), meld_L = zeros(0, 0),
        )
    )
    inputs0 = merge(
        inputs,
        (; meld = merge(meld, (; d = 0, L = zeros(0, 0))))
    )
    chn0 = sample(
        bvd_zone(zdc), Prior(), 4; chain_type = FlexiChains.VNChain,
        progress = false
    )
    for empty_ in (
            zone_meld_check(chn0, inputs0),
            zone_meld_check(Dict{Symbol, Any}(), inputs0),
        )

        @test nrow(empty_.table) == 0
        @test empty_.dimension == 0
        @test isnan(empty_.mean_norm_sq)
    end
end

@testitem "zone_forecast_truth: needs a vintage on the target date" setup = [
    ZoneSynthetic,
] begin
    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    made = inputs.cutoff
    with_days(extra) = begin
        hist = deepcopy(syn.obs.zone_confirmed_history)
        for (prov, zones) in hist, (nm, h) in zones

            for d in extra
                push!(h.days, inputs.n + d)
                push!(h.counts, h.counts[end] + 1)
            end
        end
        merge(syn.obs, (; zone_confirmed_history = hist))
    end
    ## A vintage past the target alone, or one either side of it, would
    ## count a different window, so the truth is missing.
    @test all(
        ismissing, zone_forecast_truth(
            with_days([9]), inputs;
            made_date = made
        )
    )
    @test all(
        ismissing, zone_forecast_truth(
            with_days([5, 9]), inputs;
            made_date = made
        )
    )
    ## A vintage on the target date gives the count since the cut-off.
    truth = zone_forecast_truth(
        with_days([5, 7, 9]), inputs;
        made_date = made
    )
    @test all(==(2), truth)
    @test all(
        ==(1), zone_forecast_truth(
            with_days([5, 7, 9]), inputs;
            made_date = made, horizon = 5
        )
    )
    ## A zone absent from the later history is missing; the rest score.
    obs = with_days([7])
    delete!(obs.zone_confirmed_history["b"], "y3")
    truth = zone_forecast_truth(obs, inputs; made_date = made)
    @test ismissing(truth[end]) && all(==(1), truth[1:(end - 1)])
end

@testitem "fit_zone: refuses observations with no zone deaths" setup = [
    ZoneSynthetic,
] begin
    using Logging: NullLogger, with_logger

    syn = zone_synthetic()
    obs = merge(
        syn.obs,
        (; zone_death_history = Dict{String, Dict{String, NamedTuple}}())
    )
    ## The unmixed warning fires first, so the logger is silenced to keep
    ## the refusal the only thing the item reports.
    with_logger(NullLogger()) do
        @test_throws ErrorException fit_zone(
            syn.chain, obs; zones = nothing,
            patch_names = ["a", "b"], patch_labels = ["A", "B"]
        )
        ## The same through the `top_zones` subset.
        @test_throws ErrorException fit_zone(
            syn.chain, obs; zones = zone_metadata(syn), top_zones = 2,
            patch_names = ["a", "b"], patch_labels = ["A", "B"]
        )
    end
end

@testitem "AD gradient: bvd_zone differentiates (Mooncake)" tags = [:ad] setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, default_adtype
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    seed!(20260518)
    model = bvd_zone(inputs.model_data)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype()
    )
    logp, grad = logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
end

@testitem "fit_zone: a short NUTS fit recovers the synthetic shares" tags = [
    :slow,
] setup = [ZoneSynthetic] begin
    using BVDOutbreakSize: bvd_zone
    using Statistics: median
    using DataFrames: nrow, names

    syn = zone_synthetic()
    chn = fit_zone(
        syn.chain, syn.obs; samples = 150, chains = 2,
        n_adapts = 150, zones = nothing, patch_names = ["a", "b"],
        patch_labels = ["A", "B"]
    )
    inputs = zone_inputs(syn)
    nz = syn.nz
    sT = [collect(v) for v in vec(collect(chn[:share_T_zone]))]
    nd = inputs.n - inputs.t0 + 1
    truth = syn.truth.shares[nd, :]
    ## The big patch's cut-off shares sit within 0.05 absolute of the truth;
    ## its zones carry hundreds of cases per vintage.
    for z in 1:5
        med = median([s[z] for s in sT])
        @test abs(med - truth[z]) < 0.05
    end
    diag = zone_diagnostics_table(chn, inputs)
    @test nrow(diag) == nz
    @test "walking" in names(diag)
    sd = zone_sampler_diagnostics(chn, inputs; max_depth = 8)
    @test 0 <= sd.depth_cap_fraction[1] <= 1
    @test length(sd.ebfmi) == 2
    @test isfinite(sd.max_rhat_R_T_walking)
end

@testitem "AD gradient: bvd_zone with mixing, correlation and the meld" tags = [
    :ad,
] setup = [ZoneSynthetic] begin
    using BVDOutbreakSize: bvd_zone, default_adtype, relative_multiplier_dims
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!, Xoshiro

    ## The configuration the registry fits: metadata for every zone, so the
    ## mixing blocks and the distance correlation are on, and the shared
    ## draw is sampled.
    syn = zone_synthetic()
    inputs = zone_inputs(syn; zones = zone_metadata(syn))
    zd = inputs.model_data
    @test zd.mixing !== nothing && !isempty(zd.zone_distances) && zd.meld_d > 0
    seed!(20260518)
    model = bvd_zone(zd)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype()
    )
    logp, grad = logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
    ## The gradient stays finite at the extremes of the zone decay and the
    ## destination weights, and the decay reaches the likelihood.
    nc = relative_multiplier_dims(zd.patch_ranges)
    lps = map((0.02, 6.0)) do γ
        params = (;
            γ_zone = γ, σ_destination = 4.0,
            z_destination = [(-1)^k * 3.0 for k in 1:nc],
        )
        vix = DynamicPPL.link(
            DynamicPPL.VarInfo(
                Xoshiro(3), model, DynamicPPL.InitFromParams(params)
            ), model
        )
        lp, gr = logdensity_and_gradient(ldf, collect(vix[:]))
        @test isfinite(lp)
        @test all(isfinite, gr)
        lp
    end
    @test lps[1] != lps[2]
end

@testitem "fit_zone: a short NUTS fit with mixing and correlation on" tags = [
    :slow,
] setup = [ZoneSynthetic] begin
    using BVDOutbreakSize: fit_diagnostics

    syn = zone_synthetic()
    chn = fit_zone(
        syn.chain, syn.obs; samples = 60, chains = 1, n_adapts = 60,
        zones = zone_metadata(syn), patch_names = ["a", "b"],
        patch_labels = ["A", "B"]
    )
    @test size(chn, 1) == 60
    @test BVDOutbreakSize._has_key(chn, :mixing_epsilon_zone)
    @test BVDOutbreakSize._has_key(chn, :correlation_length_zone)
    d = fit_diagnostics(chn)
    @test isfinite(d.max_rhat) || isnan(d.max_rhat)
end

@testitem "zone forecast probabilities: thresholds and the table" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_forecast
    using Turing: sample, Prior
    using DataFrames: names
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    chn = sample(
        bvd_zone(inputs.model_data), Prior(), 6;
        chain_type = FlexiChains.VNChain, progress = false
    )
    fcin = zone_inputs(
        syn; parent_forecast = zone_parent_forecast(syn; totals = [32, 8])
    )
    fc = zone_forecast(chn, fcin)
    P = zone_forecast_probabilities(fc, fcin; thresholds = (1, 5, 10))
    @test size(P) == (syn.nz, 3)
    @test all(0 .<= P .<= 1)
    ## Non-increasing in the threshold, zone by zone.
    for z in 1:syn.nz
        @test P[z, 1] >= P[z, 2] >= P[z, 3]
    end
    ## A patch total below the threshold gives zero for every zone in it.
    fcin0 = zone_inputs(
        syn; parent_forecast = zone_parent_forecast(syn; totals = [0, 0])
    )
    P0 = zone_forecast_probabilities(
        zone_forecast(chn, fcin0), fcin0; thresholds = (1,)
    )
    @test all(iszero, P0)
    ## The table carries one column per threshold, missing on patch rows.
    t = zone_forecast_table(fc, fcin; thresholds = (1, 10))
    @test all(c -> c in names(t), ["p_ge_1", "p_ge_10"])
    @test all(ismissing, t[t.zone .== "Patch total", :p_ge_1])
    @test all(x -> ismissing(x) || 0 <= x <= 1, t.p_ge_1)
    ## Recent cases and last case dates follow the count matrix.
    recent = zone_recent_cases(inputs; window = 7)
    @test length(recent) == syn.nz && all(recent .>= 0)
    @test sum(recent) <= sum(inputs.counts)
    last = zone_last_case_dates(inputs)
    @test length(last) == syn.nz
    @test all(d -> ismissing(d) || d <= inputs.cutoff, last)
end

@testitem "parent priors: Beta shapes are floored above one" begin
    using BVDOutbreakSize: _zone_parent_beta, BETA_SHAPE_FLOOR

    ## Draws piled against zero would match a Beta with a shape below one,
    ## whose density diverges at the boundary.
    chn = Dict{Symbol, Any}(:corr => vcat(fill(0.0, 40), fill(0.1, 10)))
    a, b = _zone_parent_beta(chn, :corr, (1.0, 3.0); clamp_unit = true)
    @test a >= BETA_SHAPE_FLOOR && b >= BETA_SHAPE_FLOOR
    ## A well-spread posterior keeps its moment match.
    chn2 = Dict{Symbol, Any}(:rho => [0.02, 0.03, 0.025, 0.035, 0.028])
    a2, b2 = _zone_parent_beta(chn2, :rho, (1.0, 24.0))
    @test a2 / (a2 + b2) ≈ 0.0276 atol = 1.0e-3
end

@testitem "fit_diagnostics leaves out quantities undefined in some draw" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, fit_diagnostics, _nonfinite_keys
    using Turing: sample, Prior
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn; walk_threshold = 10^6)
    ## This seed leaves a zone below `rt_floor` in some draw.
    chn = sample(
        Xoshiro(6), bvd_zone(inputs.model_data), Prior(), 8;
        chain_type = FlexiChains.VNChain, progress = false
    )
    rT = [collect(v) for v in vec(collect(chn[:R_T_zone]))]
    keys_ = _nonfinite_keys(chn)
    @test "R_T_zone" in keys_
    @test ("R_T_zone" in keys_) == any(v -> any(isnan, v), rT)
    d = fit_diagnostics(chn)
    @test isfinite(d.max_rhat) || isnan(d.max_rhat)
end

@testitem "zone parent extract stands in for the chain" setup = [ZoneSynthetic] begin
    syn = zone_synthetic()
    parent = merge(
        syn.chain,
        Dict(
            :province_shares => reshape(
                [[0.8 0.8; 0.2 0.2] for _ in 1:4], 4, 1
            )
        )
    )
    ex = zone_parent_extract(parent; source = "test", n_patches = 2)
    @test ex.n == syn.obs.n
    @test :infections_patch in ex.keys && :province_shares in ex.keys
    @test ex.source == "test"
    kw = (;
        zones = nothing, walk_threshold = 30,
        patch_names = ["a", "b"], patch_labels = ["A", "B"],
    )
    a = zone_fit_inputs(parent, syn.obs; kw...)
    b = zone_fit_inputs(ex, syn.obs; kw...)
    @test a.model_data.I_bar == b.model_data.I_bar
    @test a.model_data.g == b.model_data.g
    @test a.model_data.meld_L == b.model_data.meld_L
    @test a.model_data.parent_priors == b.model_data.parent_priors
    ## A chain without the patch structure is refused.
    bare = Dict{Symbol, Any}(:C_T => fill(1.0, 4, 1))
    @test_throws ErrorException zone_parent_extract(bare)
end

@testitem "zone_subset keeps the largest zones and pools the rest" setup = [
    ZoneSynthetic,
] begin
    syn = zone_synthetic()
    meta = zone_metadata(syn)
    sub = zone_subset(syn.obs; top_zones = 2, zones = meta)
    conf = sub.obs.zone_confirmed_history
    ## Two kept plus the pooled rest in the five-zone patch; two kept and
    ## the single remaining zone under its own name in the three-zone patch.
    @test length(sub.kept["a"]) == 2 && haskey(conf["a"], "rest")
    @test length(conf["a"]) == 3
    @test length(sub.kept["b"]) == 3 && !haskey(conf["b"], "rest")
    ## The pooled series is the sum of what it pooled, so every patch total
    ## is unchanged at every vintage.
    for prov in ("a", "b"), v in eachindex(syn.days)
        before = sum(h.counts[v] for h in values(syn.obs.zone_confirmed_history[prov]))
        after = sum(h.counts[v] for h in values(conf[prov]))
        @test before == after
    end
    deaths = sub.obs.zone_death_history
    @test haskey(deaths["a"], "rest")
    ## Metadata: the kept rows, plus one pooled row with the summed
    ## population and a centroid inside the pooled zones' box.
    rows = sub.zones
    rest = only(r for r in rows if r.zone == "rest")
    pooled = [r for r in meta if r.province == "a" && !(r.zone in sub.kept["a"])]
    @test rest.population == sum(r.population for r in pooled)
    @test minimum(r.lat for r in pooled) <= rest.lat <= maximum(r.lat for r in pooled)
    @test count(r -> r.province == "a", rows) == 3
    ## The subset fits: mixing and correlation stay on with the pooled row.
    inputs = zone_fit_inputs(
        syn.chain, sub.obs; zones = rows, walk_threshold = 30,
        patch_names = ["a", "b"], patch_labels = ["A", "B"]
    )
    @test length(inputs.zone_keys) == 6
    @test inputs.model_data.mixing !== nothing
    @test !isempty(inputs.model_data.zone_distances)
    @test_throws ArgumentError zone_subset(syn.obs; top_zones = 0, zones = meta)
end

@testitem "zone diagnostics leave out a reproduction number undefined in some draw" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, _zone_rt_defined
    using Turing: sample, Prior
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn; walk_threshold = 10^6)
    ## This seed leaves a zone below `rt_floor` in some draw.
    chn = sample(
        Xoshiro(6), bvd_zone(inputs.model_data), Prior(), 8;
        chain_type = FlexiChains.VNChain, progress = false
    )
    defined = _zone_rt_defined(chn, syn.nz)
    @test any(defined) && !all(defined)
    rT = [collect(v) for v in vec(collect(chn[:R_T_zone]))]
    for z in 1:syn.nz
        @test defined[z] == all(r -> isfinite(r[z]), rT)
    end
    df = zone_diagnostics_table(chn, inputs)
    for z in 1:syn.nz
        defined[z] || @test isnan(df.rhat_R_T[z]) && isnan(df.ess_bulk_R_T[z])
    end
    sd = zone_sampler_diagnostics(chn, inputs)
    @test isnan(sd.max_rhat_R_T_walking) || isfinite(sd.max_rhat_R_T_walking)
end

@testitem "zone composition draws survive a share that underflows" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone
    using Turing: sample, Prior
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    chn = sample(
        bvd_zone(inputs.model_data), Prior(), 4;
        chain_type = FlexiChains.VNChain, progress = false
    )
    ## Push one zone's knots far below its patch so its share underflows.
    knots = [collect(v) for v in vec(collect(chn[:delta_knots_zone]))]
    K = length(inputs.knots)
    for v in knots, k in 1:K
        v[1 + (k - 1) * syn.nz] = -800.0
    end
    d = zone_composition_draws(chn, inputs)
    @test all(m -> all(isfinite, m), d.expected)
    @test all(m -> all(x -> x >= 0, m), d.predictive)
end

@testitem "zone multipliers carry the province level for reporting" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone
    using Turing: sample, Prior
    import FlexiChains

    ## A parent whose two patches have opposite ascertainment contrasts.
    syn = zone_synthetic()
    parent = merge(
        syn.chain,
        Dict(
            :province_ascertainment => reshape(
                [[1.5, 1 / 1.5] for _ in 1:4], 4, 1
            ),
            :province_cfr_relative => reshape(
                [[0.8, 1.25] for _ in 1:4], 4, 1
            )
        )
    )
    kw = (;
        zones = nothing, walk_threshold = 30,
        patch_names = ["a", "b"], patch_labels = ["A", "B"],
    )
    inputs = zone_fit_inputs(parent, syn.obs; kw...)
    zd = inputs.model_data
    ## One factor per zone, taken from its own patch.
    @test zd.province_ascertainment[1:5] == fill(1.5, 5)
    @test zd.province_ascertainment[6:8] ≈ fill(1 / 1.5, 3)
    @test zd.province_severity[1:5] == fill(0.8, 5)
    chn = sample(
        bvd_zone(zd), Prior(), 4;
        chain_type = FlexiChains.VNChain, progress = false
    )
    rel = [collect(v) for v in vec(collect(chn[:zone_ascertainment_relative]))]
    nat = [collect(v) for v in vec(collect(chn[:zone_ascertainment_national]))]
    sev = [collect(v) for v in vec(collect(chn[:zone_severity_relative]))]
    sevn = [collect(v) for v in vec(collect(chn[:zone_severity_national]))]
    for i in eachindex(rel)
        @test nat[i] ≈ zd.province_ascertainment .* rel[i]
        @test sevn[i] ≈ zd.province_severity .* sev[i]
        ## The within-patch contrast still has geometric mean one per patch,
        ## and the national one carries its province's level.
        for zs in inputs.patch_ranges
            @test prod(rel[i][zs])^(1 / length(zs)) ≈ 1 rtol = 1.0e-8
        end
    end
    ## The fatality scale is the model's own tight prior, not the parent's,
    ## so the death composition pins the incidence split.
    @test !haskey(zd.parent_priors, :severity_sd)
    sd = [only(v) for v in vec(collect(chn[:zone_severity_sd]))]
    @test all(>=(0), sd) && all(<(1), sd)

    ## A parent without the contrasts leaves the factors at one, so the two
    ## reported multipliers agree.
    zd0 = zone_fit_inputs(syn.chain, syn.obs; kw...).model_data
    @test all(isone, zd0.province_ascertainment)
    chn0 = sample(
        bvd_zone(zd0), Prior(), 2;
        chain_type = FlexiChains.VNChain, progress = false
    )
    r0 = [collect(v) for v in vec(collect(chn0[:zone_ascertainment_relative]))]
    n0 = [collect(v) for v in vec(collect(chn0[:zone_ascertainment_national]))]
    @test all(i -> n0[i] ≈ r0[i], eachindex(r0))
end

@testitem "fit_zone: model_kwargs reach the model" tags = [:slow] setup = [
    ZoneSynthetic,
] begin
    using Distributions: Normal, truncated

    syn = zone_synthetic()
    ## A fatality scale forced far from its default prior: the draws can
    ## only sit there if the keyword reached `bvd_zone` rather than the
    ## sampler, which drops what it does not know.
    chn = fit_zone(
        syn.chain, syn.obs; samples = 20, chains = 1, n_adapts = 20,
        zones = nothing, patch_names = ["a", "b"], patch_labels = ["A", "B"],
        model_kwargs = (;
            severity_sd_prior = truncated(Normal(0.5, 0.01); lower = 0),
        )
    )
    sd = [only(v) for v in vec(collect(chn[:zone_severity_sd]))]
    @test all(>(0.3), sd)
end

@testitem "zone compositions: the cases and the deaths streams" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_composition_draws,
        zone_composition_calibration, zone_binned_increments, safe_rate,
        _zone_composition, _zone_observed_counts,
        _zone_allocated_increments, _zone_multiplier_draws,
        _zone_modelled_shares, _zone_states
    using Turing: sample, Prior
    using DataFrames: eachrow, nrow, names
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    chn = sample(
        bvd_zone(zd), Prior(), 6; chain_type = FlexiChains.VNChain,
        progress = false
    )
    nv = length(inputs.days)

    ## Each stream reads its own observed counts, and the fixture's deaths
    ## are a tenth of the cases, so the two are genuinely different tables.
    @test _zone_observed_counts(inputs, :cases) === zd.counts
    @test _zone_observed_counts(inputs, :deaths) === zd.death_counts
    @test zd.death_counts != zd.counts
    @test sum(zd.death_counts) > 0
    @test _zone_composition(:cases).rho === :composition_rho_zone
    @test _zone_composition(:deaths).rho === :composition_rho_death_zone
    @test_throws ErrorException _zone_composition(:oops)

    ## Both streams normalise within a patch, draw a split of the observed
    ## allocated total and score in the stream-calibration columns.
    for (stream, counts) in (
            (:cases, zd.counts), (:deaths, zd.death_counts),
        )
        shares = _zone_modelled_shares(chn, inputs; stream)
        @test size(shares) == (6, syn.nz, nv)
        @test all(>=(0), shares)
        d = zone_composition_draws(chn, inputs; stream)
        for v in 1:nv, zs in inputs.patch_ranges

            N = sum(@view counts[zs, v])
            N > 0 || continue
            for i in 1:6
                @test sum(shares[i, zs, v]) ≈ 1
                @test sum(d.predictive[z][i, v] for z in zs) == N
                @test sum(d.expected[z][i, v] for z in zs) ≈ 1
            end
            @test sum(d.observed[z, v] for z in zs) ≈ 1
        end
        cal = zone_composition_calibration(chn, inputs; stream)
        @test names(cal) == [
            "Stream", "Vintages", "Bias", "50% coverage", "90% coverage",
        ]
        @test nrow(cal) == 3 && cal.Stream[end] == "All zones"
        ## A patch with no scored cell reports no coverage.
        @test all(
            r -> r.Vintages == 0 || 0 <= r["90% coverage"] <= 1,
            eachrow(cal)
        )
        @test cal.Vintages[end] > 0
    end

    ## The deaths run the same infections through the
    ## infection-to-confirmed-death delay rather than the case delay, so
    ## their allocated rates are the death operator's, not the cases'.
    st = first(_zone_states(chn, inputs))
    case_inc = _zone_allocated_increments(zd, st, :cases)
    death_inc = _zone_allocated_increments(zd, st, :deaths)
    @test size(death_inc) == size(case_inc)
    @test death_inc != case_inc
    fw = BVDOutbreakSize.zone_forward(zd, st.δ_knots, st.w0, st.ε, st.def)
    expected_death = zone_binned_increments(
        zd.death_bin, fw.infections, st.def.zone_pre .* st.w0
    )
    @test death_inc ≈ expected_death
    @test case_inc ≈ fw.increments

    ## Each stream's share carries its own within-patch relative multiplier,
    ## which is centred rather than constant and so does not cancel.
    mult = Dict(
        s => _zone_multiplier_draws(chn, s) for s in (:cases, :deaths)
    )
    @test mult[:cases] !== nothing && mult[:deaths] !== nothing
    @test mult[:cases][1] != mult[:deaths][1]
    shares = _zone_modelled_shares(chn, inputs; stream = :deaths)
    for zs in inputs.patch_ranges, v in 1:nv

        rates = [
            safe_rate(death_inc[z, v] * mult[:deaths][1][z]) for z in zs
        ]
        sum(rates) > eps() || continue
        @test shares[1, zs, v] ≈ rates ./ sum(rates)
    end
end

@testitem "the forecast inputs carry the mixing to the horizon" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_forecast, zone_share_renewal,
        zone_deformation, zone_draw_mixing
    using Turing: sample, Prior
    import FlexiChains
    syn = zone_synthetic()
    inputs = zone_inputs(
        syn; zones = zone_metadata(syn),
        parent_forecast = zone_parent_forecast(syn)
    )
    zd = inputs.model_data
    zf = zd.forecast
    nz = size(zd.counts, 1)
    horizon = zf.horizon
    @test zd.mixing !== nothing

    ## Every per-day term reaches the last day the renewal indexes. Past
    ## the cut-off the mean import odds are the parent forecast's own, and
    ## the arrival shares follow its mean infections.
    @test size(zf.mixing.import_log_odds, 2) == zd.n + horizon
    @test size(zf.mixing.arrival_shares, 3) == zd.n + horizon
    @test zf.mixing.import_log_odds[:, 1:zd.n] == zd.mixing.import_log_odds
    for d in 1:horizon
        f = 0.05 * (1 + 0.02 * d)
        @test zf.mixing.import_log_odds[2, zd.n + d] ≈ log(f / (1 - f))
        @test all(zf.mixing.import_log_odds[1, zd.n + d] .< -700)
        for p in 1:2
            @test sum(zf.mixing.arrival_shares[p, :, zd.n + d]) ≈ 1
        end
    end
    ## A parent forecast without the arrivals is refused when the zones mix.
    no_imports = zone_parent_forecast(syn)
    delete!(no_imports, :forecast_importation_patch)
    @test_throws ErrorException zone_inputs(
        syn; zones = zone_metadata(syn), parent_forecast = no_imports
    )
    ## So is one whose arrivals stop short of the horizon.
    short_imports = zone_parent_forecast(syn)
    short_imports[:forecast_importation_patch] = map(
        v -> v[1:(size(syn.I_bar, 1) * horizon - 1)],
        short_imports[:forecast_importation_patch]
    )
    @test_throws ErrorException zone_inputs(
        syn; zones = zone_metadata(syn), parent_forecast = short_imports
    )
    @test size(zf.I_bar, 2) == zd.n + horizon

    ## The forecast week's operator sums the delayed reports over
    ## `(n, n + H]`, and the extended operators bin the fitted vintages as
    ## the fitted ones do.
    nd = zd.n + horizon - zd.t0 + 1
    fx = zone_fixed_terms(zf.I_bar, zd.g, zd.f, zd.t0)
    I_ext = rand(Xoshiro(3), nd, nz)
    w = rand(Xoshiro(4), nz)
    daily = reference_delay_matrix(zd.f, nd) * I_ext .+
        reference_pre_rows(
        fx.report_pre, zd.patch_of_zone, zd.t0, zd.n + horizon
    ) .* transpose(w)
    week = (zd.n - zd.t0 + 2):nd
    @test vec(zone_binned_increments(zf.report_future, I_ext, w)) ≈
        vec(sum(daily[week, :]; dims = 1)) rtol = 1.0e-12
    fitted = 1:(zd.n - zd.t0 + 1)
    for s in (:report_bin, :death_bin)
        @test zone_binned_increments(getproperty(zf, s), I_ext, w) ≈
            zone_binned_increments(
            getproperty(zd, s), I_ext[fitted, :], w
        ) rtol = 1.0e-12
    end

    ## A term shorter than the grid is refused.
    def = zone_deformation(merge(zd, zf), nothing)
    short = zone_draw_mixing(zd.mixing, nothing, nothing)
    @test_throws DimensionMismatch zone_share_renewal(
        def.I_bar, zd.g, zeros(nd, nz), fill(1 / nz, nz),
        zd.patch_ranges, zd.t0, def.force_pre;
        mix = short, ε = fill(0.01, nz)
    )

    ## The mixed model forecasts a proper split.
    chn = sample(
        bvd_zone(zd), Prior(), 3; chain_type = FlexiChains.VNChain,
        progress = false
    )
    draws = zone_forecast_draws(zone_forecast(chn, inputs), inputs)
    for (p, zs) in enumerate(zd.patch_ranges), i in 1:3
        @test sum(draws.zones[z][i] for z in zs) == draws.patches[p][i]
    end
end

@testitem "zone meld: the origin intensities and delays join the shared quantity" setup = [
    ZoneSynthetic,
] begin
    syn = zone_synthetic()
    chain = zone_spread_chain(syn)
    spread = merge(syn, (; chain))
    plain = zone_inputs(spread)
    mixed = zone_inputs(spread; zones = zone_metadata(syn))
    np = 2
    nθ = length(syn.delay_params)
    d0 = length(plain.meld.cells_patch)
    ## After the infection cells, one log intensity per origin patch where
    ## the zones mix, then the log delay parameters.
    @test plain.meld.d == d0 + nθ
    @test mixed.meld.d == d0 + np + nθ
    @test mixed.meld.extra_cells == (d0 + 1):(d0 + np + nθ)
    @test plain.meld.extra_cells == (d0 + 1):(d0 + nθ)
    eps_cells = d0 .+ (1:np)
    @test mixed.meld.log_sums[:, eps_cells] ≈ reduce(
        vcat, [transpose(log.(v)) for v in vec(chain[:importation_epsilon_patch])]
    )
    ## The infection block of the shared quantity is unchanged, and the
    ## extra cells do not deform the patch trajectories directly.
    @test mixed.meld.L[1:d0, 1:d0] ≈ plain.meld.L[1:d0, 1:d0] rtol = 1.0e-12
    @test mixed.meld.weights[:, 1:d0] == plain.meld.weights[:, 1:d0]
    @test all(iszero, mixed.meld.weights[:, mixed.meld.extra_cells])
    zd = mixed.model_data
    @test zd.meld_d == d0 + np + nθ
    @test zd.meld_epsilon_rows == mixed.meld.L[eps_cells, :]
    @test zd.meld_delay_rows == mixed.meld.L[(d0 + np) .+ (1:nθ), :]
    @test zd.meld_delay_mean ≈ log.(syn.delay_params)
    @test size(plain.model_data.meld_epsilon_rows, 1) == 0
    @test plain.model_data.meld_delay_rows == plain.meld.L[d0 .+ (1:nθ), :]
end

@testitem "bvd_zone: each draw carries its own generation interval and delays" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_delay_parameters, _zone_states,
        _draw_vectors
    using Statistics: std
    using Turing: sample, Prior
    import FlexiChains

    syn = zone_synthetic()
    ## Parent draws that disagree on the generation-interval mean and the
    ## incubation mean.
    chain = copy(syn.chain)
    chain[Symbol("gi_state.gi_mean")] = reshape([12.0, 14.0, 15.0, 17.0], 4, 1)
    chain[Symbol("inc_state.delay_mean")] = reshape([5.5, 6.0, 6.5, 7.5], 4, 1)
    inputs = zone_inputs(merge(syn, (; chain)))
    zd = inputs.model_data
    chn = sample(
        bvd_zone(zd), Prior(), 8; chain_type = FlexiChains.VNChain,
        progress = false
    )
    θ = _draw_vectors(chn, :delay_parameters_zone)
    η = _draw_vectors(chn, :parent_eta_zone)
    ## The drawn parameters are the shared draw's delay cells, and they move
    ## from draw to draw as the parent's do.
    for i in eachindex(θ)
        @test θ[i] ≈ zone_delay_parameters(zd, η[i])
    end
    @test std([t[1] for t in θ]) > 0.1
    @test std([t[3] for t in θ]) > 0.05
    ## The render rebuilds each draw's generation interval and delays from
    ## the same draw.
    states = _zone_states(chn, inputs)
    for i in (1, 5)
        pm = zone_delay_pmfs(θ[i])
        @test states[i].def.g ≈ pm.g
        @test states[i].def.report_bin.weights ≈ zone_binned_operator(
            pm.f, BVDOutbreakSize.zone_fixed_terms(zd.I_bar, pm.g, pm.f, zd.t0),
            zd.days, zd.t0, zd.n, zd.patch_of_zone
        ).weights
    end
    @test !(states[1].def.g ≈ states[5].def.g)
end

@testitem "zone mixing: each draw's arrivals follow its origin intensities" setup = [
    ZoneSynthetic,
] begin
    using BVDOutbreakSize: bvd_zone, zone_draw_mixing, zone_parent_scale,
        zone_parent_epsilon, zone_deformation, _zone_states
    using StatsFuns: logit
    using Turing: sample, Prior
    using DataFrames: nrow
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(
        merge(syn, (; chain = zone_spread_chain(syn)));
        zones = zone_metadata(syn)
    )
    zd = inputs.model_data
    mix = zd.mixing
    np, n = size(zd.I_bar)
    ## At the parent's centre a draw reads the mean movement.
    base = zone_draw_mixing(mix, nothing, nothing)
    @test base.origin_weight == mix.origin_weight
    @test logit.(base.import_fraction[2, :]) ≈ mix.import_log_odds[2, :]
    ## The second patch receives only from the first, so its odds of an
    ## imported infection move with the first origin's intensity and the
    ## ratio of the two deformed trajectories.
    η = zeros(zd.meld_d)
    η[length(inputs.meld.cells_patch) + 1] = 1.5
    η[1] = 0.7
    scale = zone_parent_scale(zd.meld_weights, zd.meld_L, η, np, zd.n)
    e = zone_parent_epsilon(zd.meld_epsilon_rows, η)
    @test e[1] > 1.5
    draw = zone_draw_mixing(mix, scale, e)
    @test draw.origin_weight ≈ mix.origin_weight .* e[inputs.patch_of_zone]
    for t in zd.t0:zd.n
        @test logit(draw.import_fraction[2, t]) ≈ mix.import_log_odds[2, t] +
            log(e[1]) + log(scale[1, t]) - log(scale[2, t])
    end
    @test all(0 .<= draw.import_fraction .< 1)
    @test zone_deformation(zd, scale, e).mixing.import_fraction ==
        draw.import_fraction
    ## The rebuilt draws read the same movement the model sampled with.
    chn = sample(
        bvd_zone(zd), Prior(), 4; chain_type = FlexiChains.VNChain,
        progress = false
    )
    nd = zd.n - zd.t0 + 1
    shares = [collect(v) for v in vec(collect(chn[:import_share_T_zone]))]
    for (i, st) in enumerate(_zone_states(chn, inputs))
        fw = zone_forward(zd, st.δ_knots, st.w0, st.ε, st.def)
        @test fw.imports[nd, :] ./ fw.infections[nd, :] ≈ shares[i] rtol = 1.0e-8
    end
    ## The check reports the patch-week cells and counts every dimension.
    mc = zone_meld_check(chn, inputs)
    @test nrow(mc.table) == length(inputs.meld.cells_patch)
    @test mc.dimension == inputs.meld.d
end
