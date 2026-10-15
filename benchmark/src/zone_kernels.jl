# The health-zone model's inner kernels, each timed on its own.
#
# The whole-model gradient sits under "AD gradients"/"Zone"; these are the
# pieces it spends that time in. They run on the inputs of the `bvd_zone`
# scenario in `test/ad_fixtures.jl`, at one prior draw, so a kernel's time
# here is its share of that scenario's gradient. The patch quantities are
# the cut, the province model's posterior mean, rather than one draw of the
# shared quantity: the kernels' cost does not depend on which.
#
# The share renewal is timed forward and through Mooncake, whose gradient
# runs the hand-written rule in `src/mooncake_rules.jl`. The gradient is of
# a fixed weighting of the zone infections, so every output feeds back.

using BVDOutbreakSize: zone_share_renewal_kernel, zone_composition_logpdf,
    deviation_knots, zone_binned_increments, zone_initial_shares,
    zone_correlation_factors, zone_deformation, _zone_kappa
using Mooncake: Mooncake
using Random: Xoshiro
using Turing: DynamicPPL

## The prior draw's quantities the kernels take, built as `bvd_zone` builds
## them.
function zone_kernel_inputs(model)
    zd = model.args.zd
    draw = rand(Xoshiro(ADFixtures.SEED), model)
    names = (
        :z_w, :σ_level, :z_level, :δ_halflife, :σ_δ, :ρ, :ρ_death, :z_drift,
        :ε_within, :τ_mix, :z_mix, :ρ_corr,
    )
    θ = NamedTuple{names}(Tuple(draw[DynamicPPL.VarName{s}()] for s in names))
    K = length(zd.knots)
    ℓ = -zd.correlation_distance / log(θ.ρ_corr)
    level_factors = zone_correlation_factors(zd.zone_distances, ℓ)
    drift_factors = zone_correlation_factors(zd.zone_walk_distances, ℓ)
    knot_args = (
        θ.z_level, θ.z_drift, θ.σ_level, θ.σ_δ, exp2(-zd.week / θ.δ_halflife),
        zd.patch_ranges, level_factors, drift_factors, zd.walking, K,
    )
    δ_knots = deviation_knots(knot_args...)
    w0 = zone_initial_shares(θ.z_w, zd.patch_ranges, zd.share_scale)
    nz = length(w0)
    logit_ε = log(θ.ε_within / (1 - θ.ε_within)) .+
        θ.τ_mix .* (θ.z_mix .- sum(θ.z_mix) / nz)
    ε = @. 1 / (1 + exp(-logit_ε))
    def = zone_deformation(zd, nothing)
    δ_daily = zd.interp * transpose(δ_knots)
    st = zone_share_renewal_kernel(
        def.I_bar, def.g, δ_daily, w0, zd.patch_ranges, zd.t0,
        def.force_pre, def.mixing, ε
    )
    pre = def.zone_pre .* w0
    cases = zone_binned_increments(def.report_bin, st.infections, pre)
    deaths = zone_binned_increments(def.death_bin, st.infections, pre)
    return (;
        zd, θ, knot_args, w0, ε, def, δ_daily, st, pre, cases, deaths,
    )
end

## Scalar summary of the renewal for the gradient: a fixed weighting of the
## zone infections, each day and zone weighted differently.
function zone_renewal_loss(I_bar, g, δ_daily, w0, force_pre, mix, ε, a)
    st = zone_share_renewal_kernel(
        I_bar, g, δ_daily, w0, a.patch_ranges, a.t0, force_pre, mix, ε
    )
    return sum(a.weights .* st.infections)
end

let scen = only(filter(s -> s.group == "Zone", SCENARIOS))
    k = zone_kernel_inputs(scen.model)
    zd, def = k.zd, k.def
    grp = SUITE["Zone kernels"] = BenchmarkGroup()

    renewal = grp["share renewal"] = BenchmarkGroup()
    renewal["forward"] = @benchmarkable zone_share_renewal_kernel(
        $(def.I_bar), $(def.g), $(k.δ_daily), $(k.w0), $(zd.patch_ranges),
        $(zd.t0), $(def.force_pre), $(def.mixing), $(k.ε)
    )
    a = (;
        zd.patch_ranges, zd.t0,
        weights = reshape(
            sin.(1:length(k.st.infections)) .+ 1.5, size(k.st.infections)
        ),
    )
    args = (
        zone_renewal_loss, def.I_bar, def.g, k.δ_daily, k.w0, def.force_pre,
        def.mixing, k.ε, a,
    )
    cache = Mooncake.prepare_gradient_cache(args...)
    renewal["Mooncake gradient"] = @benchmarkable(
        Mooncake.value_and_gradient!!($cache, $(args)...)
    )

    composition = grp["composition logpdf"] = BenchmarkGroup()
    composition["cases"] = @benchmarkable zone_composition_logpdf(
        $(zd.counts), $(k.cases), $(zd.cell_patch), $(zd.cell_vintage),
        $(zd.cell_total), $(zd.cell_const), $(zd.patch_ranges),
        $(_zone_kappa(k.θ.ρ))
    )
    composition["deaths"] = @benchmarkable zone_composition_logpdf(
        $(zd.death_counts), $(k.deaths), $(zd.death_cell_patch),
        $(zd.death_cell_vintage), $(zd.death_cell_total),
        $(zd.death_cell_const), $(zd.patch_ranges),
        $(_zone_kappa(k.θ.ρ_death))
    )

    grp["deviation knots"] = @benchmarkable deviation_knots($(k.knot_args)...)

    binned = grp["binned increments"] = BenchmarkGroup()
    binned["cases"] = @benchmarkable zone_binned_increments(
        $(def.report_bin), $(k.st.infections), $(k.pre)
    )
    binned["deaths"] = @benchmarkable zone_binned_increments(
        $(def.death_bin), $(k.st.infections), $(k.pre)
    )
end
