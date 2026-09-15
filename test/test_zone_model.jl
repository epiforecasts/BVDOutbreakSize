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
    using BVDOutbreakSize: zone_initial_shares, zone_deviation_knots,
                           zone_fixed_terms, zone_report_increments,
                           zone_forward, discretise_censored,
                           lognormal_meansd, convolve_pmf, knot_days,
                           zone_interpolation_weights, zone_delay_operator,
                           zone_report_pre_rows
    using Dates: Date, Day
    using Distributions: Gamma, Multinomial
    using Random: Xoshiro

    function zone_synthetic(; n = 120, first_vintage = 60, step = 3,
            lead_days = 42, seed = 11, total_scale = 1.0)
        rng = Xoshiro(seed)
        np = 2
        ## Fixed patch infections: growth then a plateau, two very different
        ## sizes so the small patch exercises the near-zero force paths.
        I_bar = zeros(Float64, np, n)
        for t in 1:n
            I_bar[1, t] = 30.0 * exp(0.03 * min(t, 80)) * total_scale
            I_bar[2, t] = 5.0 * exp(0.02 * min(t, 80)) * total_scale
        end
        g = let pmf = discretise_censored(Gamma(2.71, 5.65), 40)
            pmf[2:end] ./ sum(pmf[2:end])
        end
        f = convolve_pmf(
            discretise_censored(lognormal_meansd(6.3, 3.5), 16),
            discretise_censored(lognormal_meansd(4.5, 4.0), 17))
        zones = [("a", ["z1", "z2", "z3", "z4", "z5"]),
            ("b", ["y1", "y2", "y3"])]
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
        walk_index = collect(1:nz)
        n_walking = nz
        σ_level = 0.25
        σ_δ = 0.08
        φ = exp2(-7 / 42)
        z_level = [0.6, -0.4, 0.2, -0.3, -0.1, 0.5, -0.2, -0.3]
        z_drift = zeros(n_walking * (K - 1))
        for k in 2:K
            z_drift[(k - 2) * n_walking + 1] = 0.8
            z_drift[(k - 2) * n_walking + 2] = -0.8
        end
        δ_knots = zone_deviation_knots(z_level, z_drift, σ_level, σ_δ, φ,
            patch_ranges, walking, walk_index, n_walking, K)
        fixed = zone_fixed_terms(I_bar, g, f, t0)
        patch_of_zone = [1, 1, 1, 1, 1, 2, 2, 2]
        zd = (; I_bar, g, f, patch_ranges, knots, t0, n, days,
            fixed.force_pre, fixed.report_pre, fixed.report_pre_cum,
            fixed.infections_pre, mixing_kernel = zeros(nz, nz),
            interp = zone_interpolation_weights(knots, t0, n),
            report_matrix = zone_delay_operator(f, n - t0 + 1),
            report_pre_rows = zone_report_pre_rows(fixed.report_pre,
                patch_of_zone, t0, n))
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
                hist[prov][nm] = (; days = copy(days),
                    counts = cumsum(counts[z, :]))
            end
        end
        seeding = Date("2026-02-13")
        obs = (; n, seeding, cutoff = seeding + Day(n - 1),
            zone_confirmed_history = hist)
        ## Stand-in parent chain: four identical draws of the truth.
        ndraw = 4
        chain = Dict{Symbol, Any}(
            :infections_patch => reshape([vec(I_bar) for _ in 1:ndraw],
                ndraw, 1),
            Symbol("gi_state.α") => fill(2.71, ndraw, 1),
            Symbol("gi_state.θ") => fill(5.65, ndraw, 1),
            Symbol("inc_state.delay_mean") => fill(6.3, ndraw, 1),
            Symbol("inc_state.delay_sd") => fill(3.5, ndraw, 1),
            Symbol("confirmed_state.receipt_state.d.delay_mean") =>
                fill(4.5, ndraw, 1),
            Symbol("confirmed_state.receipt_state.d.delay_sd") =>
                fill(4.0, ndraw, 1),
            :C_T => fill(sum(I_bar), ndraw, 1))
        truth = (; z_w, w0, δ_knots, z_level, z_drift, σ_level, σ_δ, φ,
            δ_daily = zd.interp * transpose(δ_knots),
            walking, shares = fw.shares, increments = fw.increments,
            infections = fw.infections, forces = fw.forces)
        return (; obs, chain, truth, I_bar, g, f, patch_ranges, days, t0,
            knots, counts, nz)
    end

    zone_inputs(syn; kwargs...) = zone_fit_inputs(syn.chain, syn.obs;
        zones = nothing, patch_names = ["a", "b"], patch_labels = ["A", "B"],
        walk_threshold = 30, kwargs...)
end

@testitem "dirichlet_multinomial_logpdf: matches Distributions" begin
    using BVDOutbreakSize: dirichlet_multinomial_logpdf
    using Distributions: DirichletMultinomial, logpdf

    y = [7, 0, 3, 12]
    α = [2.0, 0.5, 1.5, 6.0]
    @test dirichlet_multinomial_logpdf(y, α) ≈
          logpdf(DirichletMultinomial(sum(y), α), y)
    ## A single-category split is certain.
    @test dirichlet_multinomial_logpdf([5], [3.0]) ≈ 0 atol = 1e-12
    ## Zero trials contribute nothing.
    @test dirichlet_multinomial_logpdf([0, 0], [1.0, 2.0]) ≈ 0 atol = 1e-12
end

@testitem "zone composition: one sum of Dirichlet-multinomials" setup=[
    ZoneSynthetic
] begin
    using BVDOutbreakSize: zone_composition_logpdf,
                           dirichlet_multinomial_logpdf
    using Distributions: DirichletMultinomial, logpdf

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    ρ = 0.05
    κ = (1 - ρ) / ρ
    C = syn.truth.increments
    lp = zone_composition_logpdf(zd.counts, C, zd.cell_patch, zd.cell_vintage,
        zd.cell_total, zd.cell_const, zd.patch_ranges, κ)
    ## The same sum written cell by cell through the full mass.
    cell_lp(c) = begin
        v = zd.cell_vintage[c]
        zs = zd.patch_ranges[zd.cell_patch[c]]
        π = C[zs, v] ./ sum(C[zs, v])
        logpdf(DirichletMultinomial(zd.cell_total[c], κ .* π),
            zd.counts[zs, v])
    end
    expected = sum(cell_lp(c) for c in eachindex(zd.cell_patch))
    @test lp ≈ expected rtol = 1e-10
    ## Scale-free in the modelled level: only the split is scored.
    lp2 = zone_composition_logpdf(zd.counts, 100 .* C, zd.cell_patch,
        zd.cell_vintage, zd.cell_total, zd.cell_const, zd.patch_ranges, κ)
    @test lp2 ≈ lp rtol = 1e-10
    ## Every scored cell has a positive allocated total.
    @test all(>(0), zd.cell_total)
    @test all(
        c -> sum(zd.counts[zd.patch_ranges[zd.cell_patch[c]],
            zd.cell_vintage[c]]) == zd.cell_total[c],
        eachindex(zd.cell_total))
end

@testitem "zone share renewal: conservation and the pre-t0 shortcut" setup=[
    ZoneSynthetic
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
    st = zone_share_renewal(I_bar, g, δ_daily, w0, syn.patch_ranges, t0,
        fixed.force_pre)
    ## Shares sum to one within every patch on every day, and the zone
    ## infections sum to the patch infections.
    for j in 1:nd, (p, zs) in enumerate(syn.patch_ranges)

        @test sum(st.shares[j, zs]) ≈ 1 rtol = 1e-12
        @test sum(st.infections[j, zs]) ≈ I_bar[p, t0 + j - 1] rtol = 1e-12
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
        @test st.forces[j, z] ≈ explicit rtol = 1e-10
    end
end

@testitem "zone operators: matrix forms equal the loops" setup=[
    ZoneSynthetic
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
        @test δd[:, z] ≈ ref rtol = 1e-12
    end
    ## The delay operator plus the pre-t0 rows reproduce the convolution,
    ## with the days before t0 held at the initial share.
    w0 = syn.truth.w0
    st = zone_share_renewal(zd.I_bar, zd.g, δd, w0, zd.patch_ranges, t0,
        zd.force_pre)
    r_op = zd.report_matrix * st.infections .+
           zd.report_pre_rows .* transpose(w0)
    for (p, zs) in enumerate(zd.patch_ranges), z in zs, j in 1:7:nd
        t = t0 + j - 1
        explicit = 0.0
        for s in 0:min(length(zd.f) - 1, t - 1)
            I_past = t - s < t0 ? w0[z] * zd.I_bar[p, t - s] :
                     st.infections[j - s, z]
            explicit += zd.f[s + 1] * I_past
        end
        @test r_op[j, z] ≈ explicit rtol = 1e-10
    end
    ## The model's forward pass matches the generator's increments.
    fw = zone_forward(zd, δk, w0, nothing)
    @test fw.increments ≈ syn.truth.increments rtol = 1e-10
end

@testitem "zone Rt: the force-weighted mean is the implied patch Rt" setup=[
    ZoneSynthetic
] begin
    using BVDOutbreakSize: zone_share_renewal, zone_fixed_terms,
                           implied_national_Rt_at

    syn = zone_synthetic()
    I_bar, g = syn.I_bar, syn.g
    np, n = size(I_bar)
    t0 = syn.t0
    δ_daily = syn.truth.δ_daily
    fixed = zone_fixed_terms(I_bar, g, syn.f, t0)
    st = zone_share_renewal(I_bar, g, δ_daily, syn.truth.w0,
        syn.patch_ranges, t0, fixed.force_pre)
    for (p, zs) in enumerate(syn.patch_ranges), t in t0:7:n

        j = t - t0 + 1
        R = [st.infections[j, z] / st.forces[j, z] for z in zs]
        Λ = [st.forces[j, z] for z in zs]
        @test sum(R .* Λ) / sum(Λ) ≈
              implied_national_Rt_at(vec(I_bar[p, :]), g, t) rtol = 1e-10
    end
end

@testitem "zone deviations: centred, and level-only zones decay" setup=[
    ZoneSynthetic
] begin
    using BVDOutbreakSize: zone_deviation_knots, zone_initial_shares

    syn = zone_synthetic()
    K = size(syn.truth.δ_knots, 2)
    φ = syn.truth.φ
    ## Three walking zones in the big patch, none in the small one.
    walking = [true, true, true, false, false, false, false, false]
    walk_index = [1, 2, 3, 0, 0, 0, 0, 0]
    z_drift = randn(Xoshiro(3), 3 * (K - 1))
    δ = zone_deviation_knots(syn.truth.z_level, z_drift, 0.25, 0.08, φ,
        syn.patch_ranges, walking, walk_index, 3, K)
    for (p, zs) in enumerate(syn.patch_ranges), k in 1:K

        @test abs(sum(δ[zs, k])) < 1e-12
    end
    ## A level-only zone follows the AR mean path from its level.
    for z in findall(!, walking), k in 1:K

        @test δ[z, k] ≈ φ^(k - 1) * δ[z, 1] rtol = 1e-12
    end
    ## A walking zone leaves the mean path.
    @test !isapprox(δ[1, K], φ^(K - 1) * δ[1, 1]; rtol = 1e-3)
    ## The synthetic truth's own knots are centred too.
    for (p, zs) in enumerate(syn.patch_ranges), k in 1:K

        @test abs(sum(syn.truth.δ_knots[zs, k])) < 1e-12
    end
    ## The initial shares sum to one within each patch and are invariant to
    ## a common shift of the draws.
    w = zone_initial_shares(syn.truth.z_w, syn.patch_ranges, 2.0)
    w_shift = zone_initial_shares(syn.truth.z_w .+ 3.0, syn.patch_ranges, 2.0)
    for zs in syn.patch_ranges
        @test sum(w[zs]) ≈ 1 rtol = 1e-12
    end
    @test w ≈ w_shift rtol = 1e-12
end

@testitem "zone_fit_inputs: units, cells, walking set and start values" setup=[
    ZoneSynthetic
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
    @test zd.I_bar ≈ syn.I_bar rtol = 1e-10
    @test zd.g ≈ syn.g
    @test zd.f ≈ syn.f
    ## The walking set needs the threshold and at least two zones per patch.
    @test all(
        z -> inputs.walking[z] == (inputs.cumulative[z] >= 30) &&
             (!inputs.walking[z] ||
              count(inputs.walking[inputs.patch_ranges[
            inputs.patch_of_zone[z]]]) >= 2),
        1:syn.nz)
    @test zd.n_walking == count(inputs.walking)
    @test all(z -> (zd.walk_index[z] > 0) == inputs.walking[z], 1:syn.nz)
    ## The initial-share start is the centred log observed first-vintage
    ## cumulative, halved.
    for zs in inputs.patch_ranges
        lg = log.(syn.counts[zs, 1] .+ 0.5)
        @test inputs.z_w_start[zs] ≈ (lg .- sum(lg) / length(lg)) ./ 2
    end
    ## Cumulative counts come from the manifest's last vintage.
    @test inputs.cumulative == vec(sum(syn.counts; dims = 2))
end

@testitem "bvd_zone: the truth scores above a perturbed truth" setup=[
    ZoneSynthetic
] begin
    using BVDOutbreakSize: bvd_zone
    using Turing: DynamicPPL

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    m = bvd_zone(zd)
    K = length(zd.knots)
    truth = syn.truth
    ## Only the zones the inputs mark as walking carry innovations; the
    ## synthetic truth walks the first three zones, which is the marked set
    ## whenever the counts clear the threshold.
    @test inputs.walking == truth.walking
    ## The counts are multinomial, so the composition is evaluated near its
    ## multinomial limit (`ρ → 0`): at a coarser `ρ` the
    ## Dirichlet-multinomial expects more dispersion than the data carry and
    ## its likelihood is not maximised at the generating shares.
    at(z_w, z_level, z_drift) = (; z_w, σ_level = truth.σ_level, z_level,
        δ_halflife = 42.0, σ_δ = truth.σ_δ, z_drift, ρ = 1e-4)
    loglik(params) = begin
        vi = DynamicPPL.VarInfo(Xoshiro(1), m, DynamicPPL.InitFromParams(params))
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
    prof = [loglik(at(truth.z_w .+ [d, 0, 0, 0, 0, 0, 0, 0],
                truth.z_level, truth.z_drift)) for d in grid]
    @test abs(grid[argmax(prof)]) <= 0.1
end

@testitem "bvd_zone: prior draws carry the deterministics" setup=[
    ZoneSynthetic
] begin
    using BVDOutbreakSize: bvd_zone
    using Turing: sample, Prior
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    nz = syn.nz
    K = length(zd.knots)
    chn = sample(bvd_zone(zd), Prior(), 20; chain_type = FlexiChains.VNChain,
        progress = false)
    for (q, len) in ((:delta_knots_zone, nz * K), (:share_knots_zone, nz * K),
        (:delta_T_zone, nz), (:share_T_zone, nz), (:share_start_zone, nz),
        (:R_T_zone, nz))
        @test all(v -> length(v) == len, vec(collect(chn[q])))
    end
    for q in (:region_sd_zone, :region_drift_sd_zone, :region_halflife_zone,
        :composition_rho_zone)
        @test all(isfinite, vec(Array(chn[q])))
    end
    ## Shares sum to one within each patch in every draw.
    for v in vec(collect(chn[:share_T_zone])), zs in inputs.patch_ranges

        @test sum(v[zs]) ≈ 1 rtol = 1e-10
    end
    ## Deviations sum to zero within each patch at the cut-off.
    for v in vec(collect(chn[:delta_T_zone])), zs in inputs.patch_ranges

        @test abs(sum(v[zs])) < 1e-10
    end
    ## Reconstruction reproduces the stored cut-off values.
    shares = reconstruct_zone_shares(chn, inputs)
    rt = reconstruct_zone_rt(chn, inputs)
    sT = [collect(v) for v in vec(collect(chn[:share_T_zone]))]
    rT = [collect(v) for v in vec(collect(chn[:R_T_zone]))]
    w0 = [collect(v) for v in vec(collect(chn[:share_start_zone]))]
    for i in 1:20, z in 1:nz

        @test shares[z][i, inputs.n] ≈ sT[i][z] rtol = 1e-10
        @test shares[z][i, 1] ≈ w0[i][z] rtol = 1e-10
        if isfinite(rT[i][z])
            @test rt[z][i, inputs.n] ≈ rT[i][z] rtol = 1e-10
        else
            @test isnan(rt[z][i, inputs.n])
        end
    end
    ## Zone infections paired with the parent sum to the parent's patches.
    inf = zone_infections(chn, syn.chain, inputs)
    for i in 1:20, t in 1:inputs.n, (p, zs) in enumerate(inputs.patch_ranges)
        @test sum(inf[z][i, t] for z in zs) ≈ syn.I_bar[p, t] rtol = 1e-8
    end
end

@testitem "zone forecast: continuous with the fit and a proper split" setup=[
    ZoneSynthetic
] begin
    using BVDOutbreakSize: bvd_zone, _zone_states, _zone_extended_data,
                           _zone_extended_deviations, zone_forward_daily
    using Turing: sample, Prior
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    chn = sample(bvd_zone(zd), Prior(), 5; chain_type = FlexiChains.VNChain,
        progress = false)
    H = 7
    fc = zone_forecast_shares(chn, inputs; horizon = H)
    @test size(fc) == (5, syn.nz, H)
    for i in 1:5, d in 1:H, zs in inputs.patch_ranges
        @test sum(fc[i, zs, d]) ≈ 1 rtol = 1e-10
    end
    ## The extended recursion reproduces the fitted days exactly and the
    ## deviations continue on the AR mean path.
    st = first(_zone_states(chn, inputs))
    fitted = zone_forward(zd, st.δ_knots, st.w0, st.ε)
    zd_ext = _zone_extended_data(inputs; horizon = H)
    δ_ext = _zone_extended_deviations(st, inputs; horizon = H)
    ext = zone_forward_daily(zd_ext, δ_ext, st.w0, st.ε)
    nd = inputs.n - inputs.t0 + 1
    @test ext.shares[1:nd, :] ≈ fitted.shares rtol = 1e-10
    @test ext.infections[1:nd, :] ≈ fitted.infections rtol = 1e-10
    for d in 1:H, z in 1:syn.nz

        @test δ_ext[nd + d, z] ≈ st.φ^(d / 7) * δ_ext[nd, z] rtol = 1e-12
    end
    ## Patch infections continue at the cut-off weekly growth.
    for p in 1:2
        ratio = syn.I_bar[p, inputs.n] / syn.I_bar[p, inputs.n - 7]
        @test zd_ext.I_bar[p, inputs.n + 7] ≈ syn.I_bar[p, inputs.n] * ratio
    end
end

@testitem "zone tables: overview, forecast, truth and scores" setup=[
    ZoneSynthetic
] begin
    using BVDOutbreakSize: bvd_zone
    using Turing: sample, Prior
    using DataFrames: DataFrame, nrow, names
    using Dates: Day
    import FlexiChains

    syn = zone_synthetic()
    inputs = zone_inputs(syn)
    zd = inputs.model_data
    chn = sample(bvd_zone(zd), Prior(), 8; chain_type = FlexiChains.VNChain,
        progress = false)
    ov = zone_overview_table(chn, inputs)
    @test nrow(ov) == syn.nz
    @test issorted(
        [(w ? 2.0 : 0.0) + (isnan(x) ? -1.0 : x)
         for (w, x) in zip(ov.walking, ov.p_R_above_1)];
        rev = true)
    @test all(
        z -> isnan(ov.rt_median[z]) ||
             ov.rt_lo90[z] <= ov.rt_median[z] <= ov.rt_hi90[z], 1:syn.nz)
    @test Set(ov.zone) == Set(inputs.zone_labels)
    ## A stand-in national forecast and parent with a province split.
    nd1 = 30
    fc = DataFrame(confirmed_new = fill(100.0, nd1))
    parent = merge(syn.chain,
        Dict(:province_shares => reshape(
            [[0.8 0.8; 0.2 0.2] for _ in 1:nd1], nd1, 1)))
    ft = zone_forecast_table(chn, parent, fc, inputs)
    @test nrow(ft) == syn.nz + 2
    @test all(ft.lower_90 .<= ft.median .<= ft.upper_90)
    fd = zone_forecast_draws(chn, parent, fc, inputs)
    @test length(fd.zones) == syn.nz && length(fd.patches) == 2
    @test all(sum(fd.zones[z] for z in 1:5) .≈ fd.patches[1])
    ## Zone forecasts sum to the patch total draw by draw, so the patch rows
    ## carry the national split exactly.
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
    vt = zone_forecast_vs_truth(chn, parent, fc, inputs; truth)
    @test nrow(vt) == syn.nz + 2
    @test all(vt.observed[vt.zone .!= "Patch total"] .== 3)
    @test all(vt.lower_90 .<= vt.lower_50 .<= vt.upper_50 .<= vt.upper_90)
    sc = zone_forecast_scores(chn, parent, fc, inputs; truth)
    @test Set(sc.method) ==
          Set(["zone model", "share persistence", "naive persistence"])
    @test all(isfinite, sc.log_score)
    @test nrow(sc) == 6
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
    @test names(cal) == ["Stream", "Vintages", "Bias", "50% coverage",
        "90% coverage"]
    @test nrow(cal) == 3 && cal.Stream[end] == "All zones"
    @test all(0 .<= cal[!, "90% coverage"] .<= 1)
    ## Paired with a parent whose draws all equal the cut's mean trajectory
    ## (to round-off), the reproduction number is the unpaired one, and at
    ## the cut-off the chain's own wherever the chain reports it in every
    ## draw. The floor is decided per zone, so a reported zone carries every
    ## draw.
    rt = reconstruct_zone_rt(chn, inputs)
    rt_paired = reconstruct_zone_rt(chn, inputs; parent_chain = syn.chain)
    R_chain = [collect(v) for v in vec(collect(chn[:R_T_zone]))]
    same(a, b) = all(isequal.(isnan.(a), isnan.(b))) &&
                 isapprox(filter(!isnan, a), filter(!isnan, b); rtol = 1e-8)
    for z in 1:syn.nz
        @test same(rt[z], rt_paired[z])
        col = rt[z][:, end]
        @test all(isnan, col) || !any(isnan, col)
        r = [R_chain[i][z] for i in 1:8]
        any(isnan, r) && continue
        @test col ≈ r rtol = 1e-8
    end
    ov_paired = zone_overview_table(chn, inputs; parent_chain = syn.chain)
    @test same(ov_paired.rt_median, ov.rt_median)
    ppc = zone_composition_ppc(chn, inputs; top = 3, prior_chain = chn)
    @test all(0 .<= ppc.observed_share .<= 1)
    @test all(ppc.lower_90 .<= ppc.median .<= ppc.upper_90)
    @test ppc.prior_median == ppc.median
    @test length(unique(ppc.zone)) <= 6
    figs = plot_zone_composition_ppc(chn, inputs; prior_chain = chn, top = 2)
    @test length(figs) == 2
    arch = zone_forecast_archive(chn, parent, [(7, fc)], inputs;
        made_date = made, thin = 5)
    @test Set(propertynames(arch)) == Set([:made_date, :horizon, :target_date,
        :province, :zone, :stream, :draw, :value])
    @test nrow(arch) == syn.nz * length(1:5:nd1)
    @test all(arch.target_date .== made + Day(7))
end

@testitem "AD gradient: bvd_zone differentiates (Mooncake)" tags=[:ad] setup=[
    ZoneSynthetic
] begin
    using BVDOutbreakSize: bvd_zone, default_adtype
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!

    syn=zone_synthetic()
    inputs=zone_inputs(syn)
    seed!(20260518)
    model=bvd_zone(inputs.model_data)
    vi=DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0=collect(vi[:])
    ldf=DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype())
    logp, grad=logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
end

@testitem "fit_zone: a short NUTS fit recovers the synthetic shares" tags=[
    :slow
] setup=[ZoneSynthetic] begin
    using BVDOutbreakSize: bvd_zone
    using Statistics: median
    using DataFrames: nrow, names

    syn=zone_synthetic()
    chn=fit_zone(syn.chain, syn.obs; samples = 150, chains = 2,
        n_adapts = 150, zones = nothing, patch_names = ["a", "b"],
        patch_labels = ["A", "B"])
    inputs=zone_inputs(syn)
    nz=syn.nz
    sT=[collect(v) for v in vec(collect(chn[:share_T_zone]))]
    nd=inputs.n-inputs.t0+1
    truth=syn.truth.shares[nd, :]
    ## Every zone's cut-off share sits near the truth: within 0.05 absolute
    ## for the big patch's zones, which carry hundreds of cases per vintage.
    for z in 1:5
        med=median([s[z] for s in sT])
        @test abs(med - truth[z]) < 0.05
    end
    diag=zone_diagnostics_table(chn, inputs)
    @test nrow(diag) == nz
    @test "walking" in names(diag)
    sd=zone_sampler_diagnostics(chn, inputs; max_depth = 8)
    @test 0 <= sd.depth_cap_fraction[1] <= 1
    @test length(sd.ebfmi) == 2
    @test isfinite(sd.max_rhat_R_T_walking)
end
