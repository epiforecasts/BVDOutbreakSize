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
    deviation_knots, zone_binned_increments, zone_deformation, zone_forward,
    _zone_draw_terms, _zone_kappa
using Mooncake: Mooncake
using Random: Xoshiro
using Turing: DynamicPPL

## The prior draw's quantities the kernels take, built by the model's own
## `_zone_draw_terms` and `zone_forward` at the cut.
function zone_kernel_inputs(model)
    zd = model.args.zd
    draw = rand(Xoshiro(ADFixtures.SEED), model)
    names = (
        :z_w, :σ_level, :z_level, :δ_halflife, :σ_δ, :ρ, :ρ_death, :z_drift,
        :ε_within, :τ_mix, :z_mix, :ρ_corr,
    )
    θ = NamedTuple{names}(Tuple(draw[DynamicPPL.VarName{s}()] for s in names))
    terms = _zone_draw_terms(zd, θ, length(zd.knots))
    def = zone_deformation(zd, nothing)
    fw = zone_forward(zd, terms.δ_knots, terms.w0, terms.ε, def)
    deaths = zone_binned_increments(def.death_bin, fw.infections, fw.pre)
    return (; zd, θ, terms.knot_args, terms.w0, terms.ε, def, fw, deaths)
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
        $(def.I_bar), $(def.g), $(k.fw.δ_daily), $(k.w0), $(zd.patch_ranges),
        $(zd.t0), $(def.force_pre), $(def.mixing), $(k.ε)
    )
    a = (;
        zd.patch_ranges, zd.t0,
        weights = reshape(
            sin.(1:length(k.fw.infections)) .+ 1.5, size(k.fw.infections)
        ),
    )
    args = (
        zone_renewal_loss, def.I_bar, def.g, k.fw.δ_daily, k.w0, def.force_pre,
        def.mixing, k.ε, a,
    )
    cache = Mooncake.prepare_gradient_cache(args...)
    renewal["Mooncake gradient"] = @benchmarkable(
        Mooncake.value_and_gradient!!($cache, $(args)...)
    )

    composition = grp["composition logpdf"] = BenchmarkGroup()
    composition["cases"] = @benchmarkable zone_composition_logpdf(
        $(zd.counts), $(k.fw.increments), $(zd.cell_patch), $(zd.cell_vintage),
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
        $(def.report_bin), $(k.fw.infections), $(k.fw.pre)
    )
    binned["deaths"] = @benchmarkable zone_binned_increments(
        $(def.death_bin), $(k.fw.infections), $(k.fw.pre)
    )
end
