"""
    ZoneFixtures

A synthetic health-zone dataset generated from a known share process, with
a stand-in parent chain, so the zone model runs without a real joint fit.

This file is included by the `ZoneSynthetic` snippet in
`test/test_zone_model.jl` and by `test/ad_fixtures.jl`, so the zone tests
and the zone benchmarks build their inputs one way. Its imports come from
whichever environment includes it.
"""
module ZoneFixtures

using BVDOutbreakSize: zone_initial_shares, deviation_knots,
    deviation_knot_dims, knot_days, zone_interpolation_weights,
    zone_delay_pmfs, zone_delay_terms, zone_forward
using Dates: Date, Day
using Distributions: Multinomial
using Random: Xoshiro, randn

export zone_synthetic, zone_metadata

"""
    zone_synthetic(; n = 120, first_vintage = 60, step = 3, lead_days = 42,
                   seed = 11, total_scale = 1.0, zones_per_patch = [5, 3])

A zone dataset of `sum(zones_per_patch)` zones over `length(zones_per_patch)`
patches and `n` days, with a vintage every `step` days from
`first_vintage` to `n`. Patch `p` is named by the `p`-th letter from `a`
and its zones by the `p`-th letter back from `z`, numbered from one.

The patch infections grow then plateau, the first patch largest and the
last smallest, so the small patches exercise the near-zero force paths.
The default two patches of five and three zones carry a hand-set truth that
the recovery items are written against. Any other shape draws its initial
shares and level contrasts at `seed + 1`.

Returns the observations, the stand-in parent chain (a `Dict` of
iteration-by-chain matrices, which is all the input builder reads) and the
truth alongside.
"""
function zone_synthetic(;
        n = 120, first_vintage = 60, step = 3,
        lead_days = 42, seed = 11, total_scale = 1.0,
        zones_per_patch = [5, 3]
    )
    rng = Xoshiro(seed)
    np = length(zones_per_patch)
    ## Fixed patch infections: growth then a plateau, the sizes spread from
    ## the first patch to the last.
    size_p = collect(range(30.0, 5.0; length = np))
    rate_p = collect(range(0.03, 0.02; length = np))
    I_bar = zeros(Float64, np, n)
    for t in 1:n, p in 1:np
        I_bar[p, t] = size_p[p] * exp(rate_p[p] * min(t, 80)) * total_scale
    end
    ## The delay parameters of every parent draw, in the order of
    ## `_ZONE_DELAY_KEYS`, at the centres of the parent's priors.
    delay_params = [
        2.71 * 5.65, sqrt(2.71) * 5.65, 6.3, 3.5, 1.178, 3.694,
        4.5, 4.0, 1.178, 3.694, 2.151, 3.906,
    ]
    pmfs = zone_delay_pmfs(delay_params)
    g, f = pmfs.g, pmfs.f
    patch_names = [string('a' + (p - 1)) for p in 1:np]
    ## Zero-padded past nine, so the sorted names keep the zone order.
    zones = [
        (
            patch_names[p],
            [
                string('z' - (p - 1), lpad(i, ndigits(m), '0'))
                    for i in 1:m
            ],
        )
            for (p, m) in enumerate(zones_per_patch)
    ]
    ends = cumsum(zones_per_patch)
    patch_ranges = [(e - m + 1):e for (e, m) in zip(ends, zones_per_patch)]
    nz = ends[end]
    patch_of_zone = reduce(
        vcat, [fill(p, m) for (p, m) in enumerate(zones_per_patch)]
    )
    days = collect(first_vintage:step:n)
    days[end] == n || push!(days, n)
    t0 = clamp(days[1] - lead_days, 1, n)
    knots = knot_days(n; week = 7, start = t0)
    K = length(knots)
    ## Every zone clears the walking threshold, so all carry innovations and
    ## the truth is representable in the fitted model.
    walking = fill(true, nz)
    nd = deviation_knot_dims(patch_ranges, walking)
    if zones_per_patch == [5, 3]
        ## Uneven initial shares, a level contrast, and one zone in the big
        ## patch walking upward while another walks down.
        z_w = [1.2, 0.3, -0.2, -0.6, -0.7, 0.8, -0.3, -0.5]
        z_level = [0.6, -0.4, 0.2, -0.3, 0.5, -0.2]
    else
        truth_rng = Xoshiro(seed + 1)
        z_w = 0.6 .* randn(truth_rng, nz)
        z_level = 0.4 .* randn(truth_rng, nd.level)
    end
    w0 = zone_initial_shares(z_w, patch_ranges, 2.0)
    σ_level = 0.25
    ## One drift scale per patch, as the model samples it.
    σ_δ = fill(0.08, np)
    φ = exp2(-7 / 42)
    ## Per knot, the first basis direction of the first patch contrasts its
    ## first two zones: `+0.8` and `-0.8` in the draw's units.
    z_drift = zeros(nd.drift * (K - 1))
    for k in 2:K
        z_drift[(k - 2) * nd.drift + 1] = 0.8 * sqrt(2)
    end
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
        ## intensity and the arrivals into each patch, a twentieth of the
        ## infections of every patch after the first and none into it.
        :importation_epsilon_patch => reshape(
            [[0.02 / p for p in 1:np] for _ in 1:ndraw], ndraw, 1
        ),
        :importation_patch => reshape(
            [
                vec(
                    [
                        p >= 2 ? 0.05 * I_bar[p, t] : 0.0
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
        patch_ranges, days, t0, knots, counts, nz, zones, patch_names,
    )
end

"""
    zone_metadata(syn)

Health-zone metadata rows covering every zone of `syn`, for the mixing
kernel and the zone correlation.
"""
function zone_metadata(syn)
    pairs = [(prov, nm) for (prov, names) in syn.zones for nm in names]
    return [
        (;
            zone = nm, label = uppercase(nm), province = prov,
            population = 10_000 * i, lat = 0.1 * i, lon = 0.2 * i,
            zscode = string(i),
        )
            for (i, (prov, nm)) in enumerate(pairs)
    ]
end

end # module
