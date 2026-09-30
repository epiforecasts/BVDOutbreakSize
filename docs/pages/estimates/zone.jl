# # Health-zone estimates
#
# Each patch of the headline joint fit split across its health zones: the
# reproduction number, the share of the patch, the one-week forecast and the
# probability of at least a few cases, zone by zone.
# The [health-zone model](@ref "Health-zone model") on the methods page gives the maths.
# This page carries its results and the checks of the fit, and the interactive map is on the [dashboard](@ref "Dashboard").
# The change in each zone's estimate over the past week has its own section below.
# The one-week zone forecast is on the [health-zone forecasts](@ref "Health-zone forecasts") page and its scores in the [health-zone forecast evaluation](@ref "Health-zone forecast evaluation").

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations and the fit registry. See
## `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))
#-
## The fits this page reads, loaded from the cache here: the headline joint,
## the zone fit melded from it, and the frozen zone fit the week-on-week
## comparison reads.
chn_joint = load_fit("joint");
chn_local = load_fit("local");
frozen_local = load_fit("local_frozen_validation");

## The zone stage's fixed inputs with the joint's forecast, and the frozen
## fit's inputs for the comparison below, both from the shared setup, so the
## health-zone forecast page draws the same forecast from the same inputs.
zone_inputs = zone_stage_inputs(; forecast = true);
zone_patch = zone_inputs.patch_of_zone;
frozen_zone_inputs = frozen_zone_stage_inputs();

#md # ```@raw html
#md # </details>
#md # ```

include(joinpath(pkgdir(BVDOutbreakSize), "docs", "front_matter.jl")) #hide
MarkdownTable(report_dates(obs.cutoff)) #hide

# ## Estimates by zone
#
# The maps below show, for the [health-zone model](@ref "Health-zone model"), the reproduction number at the cut-off, the bounds of the 90% interval on the forecast confirmed cases over the coming week, and the confirmed cases to date, zone by zone.
# The forecast is mapped as its two bounds rather than a single number, so a zone's colour reads as a range.
# On the reproduction-number map a zone whose 90% interval straddles one is washed towards white.
# Zones with no confirmed case, or too few infections for a reproduction number, are grey.
# The four maps share the zone boundaries.
# A zone's forecast can be read against its reproduction number and its cases to date.
# The interactive map on the [dashboard](@ref "Dashboard") adds the probability of at least $K$ cases at a chosen $K$ and a filter for zones with or without a case over the past one, two or four weeks.

#md # ```@raw html
#md # <details><summary>Health-zone post-processing</summary>
#md # ```

## The geojson keys a zone without the province prefix the manifest carries.
zone_map_keys = [
    String(last(split(k, "."; limit = 2)))
        for k in zone_inputs.zone_keys
];
## Daily zone reproduction numbers over the zone grid, each zone draw on
## its own draw of the patch trajectory so the patch uncertainty is
## carried, and the cut-off values the map and ranking read. A zone's
## reproduction number is reported from the day its cumulative infections
## reach the floor in the median draw.
zone_rt_traj = reconstruct_zone_rt(chn_local, zone_inputs);
zone_RT_finite = [filter(isfinite, m[:, obs.n]) for m in zone_rt_traj];
zone_RT_reported = findall(!isempty, zone_RT_finite);
zone_share_T = let vs = vec(collect(chn_local[:share_T_zone]))
    [Float64[v[z] for v in vs] for z in eachindex(zone_map_keys)]
end;
_zq(v, p) = quantile(v, p)
## The one-week zone forecast drawn from the zone model, and the
## probability of at least K cases per zone from the same draws.
zone_fc = zone_forecast(chn_local, zone_inputs);
zone_fc_draws = zone_forecast_draws(zone_fc, zone_inputs);
ZONE_THRESHOLDS = (1, 5, 10, 20)
zone_fc_probs = zone_forecast_probabilities(
    zone_fc, zone_inputs; thresholds = ZONE_THRESHOLDS, draws = zone_fc_draws
);
zone_overview = zone_overview_table(chn_local, zone_inputs);
## Cases allocated to each zone over the past one, two and four weeks,
## and the last vintage on which each zone's count rose.
zone_recent = Dict(
    w => zone_recent_cases(zone_inputs; window = w) for w in (7, 14, 28)
);
zone_last_case = zone_last_case_dates(zone_inputs);

zone_map_fig = plot_zone_map_panels(
    [
        (;
            values = [median(zone_RT_finite[z]) for z in zone_RT_reported],
            zones = zone_map_keys[zone_RT_reported],
            lower = [_zq(zone_RT_finite[z], 0.05) for z in zone_RT_reported],
            upper = [_zq(zone_RT_finite[z], 0.95) for z in zone_RT_reported],
            diverging_at = 1.0, scale = log10,
            title = "Reproduction number at the cut-off",
            colorbar_label = "R",
        ),
        (;
            values = [_zq(v, 0.05) for v in zone_fc_draws.zones],
            zones = zone_map_keys, scale = CairoMakie.Makie.pseudolog10,
            title = "Confirmed cases over the coming week (lower 90%)",
            colorbar_label = "cases",
        ),
        (;
            values = [_zq(v, 0.95) for v in zone_fc_draws.zones],
            zones = zone_map_keys, scale = CairoMakie.Makie.pseudolog10,
            title = "Confirmed cases over the coming week (upper 90%)",
            colorbar_label = "cases",
        ),
        (;
            values = Float64.(zone_inputs.cumulative), zones = zone_map_keys,
            scale = CairoMakie.Makie.pseudolog10,
            title = "Confirmed cases to date", colorbar_label = "cases",
        ),
    ];
    ncols = 2, title = "Health zones at the cut-off"
);

#md # ```@raw html
#md # </details>
#md # ```

zone_map_fig #hide

# The panels below trace the reproduction number of the twelve zones with the most infections over the past week as coloured bands, each against its whole patch's implied reproduction number as a dark line with a grey band.
# Where a zone's bands depart from the patch line, the gap is the zone's fitted deviation from its patch.
# The patch line averages all the patch's zones weighted by their recent infections, so the zones shown tend to move it most.

#md # ```@raw html
#md # <details><summary>Zone reproduction-number trajectories</summary>
#md # ```

## Each patch's implied reproduction number from the joint draws with the
## same generation interval the zone stage fixes, so the patch line is close
## to the quantity the zone values average to.
zone_grid = zone_inputs.t0:obs.n
_patch_infection_draws = vec(collect(chn_joint[:infections_patch]));
patch_implied_rt = [
    let m = Matrix{Float64}(
            undef,
            length(_patch_infection_draws), obs.n
        )
        for (i, v) in enumerate(_patch_infection_draws)
            I = reshape(Float64.(v), N_PATCHES, obs.n)
            m[i, :] .= implied_national_Rt(I[p, :], zone_inputs.g)
        end
        m
    end
        for p in 1:N_PATCHES
];
## Each zone's infections over the past week, the median over draws, which
## chooses the panels.
zone_week_infections = [
    median(vec(sum(m[:, (obs.n - 6):obs.n]; dims = 2)))
        for m in zone_infections(chn_local, zone_inputs)
];
zone_rt_fig = plot_rt_zones(
    [replace(m[:, zone_grid], NaN => missing) for m in zone_rt_traj],
    zone_inputs.zone_labels, zone_patch;
    patch_labels = zone_inputs.patch_labels,
    dates = grid_date.(zone_grid), as_of_date = obs.cutoff,
    ranking = zone_week_infections,
    ranking_label = "infections over the past week", top = 12,
    modelled = zone_inputs.walking,
    patch_rt = [m[:, zone_grid] for m in patch_implied_rt]
);

#md # ```@raw html
#md # </details>
#md # ```

zone_rt_fig #hide

# The ranking below orders the zones by the posterior probability that their reproduction number exceeds one.
# A zone whose reproduction number is from its province, not modelled separately, is drawn hollow in grey.

#md # ```@raw html
#md # <details><summary>Zone ranking</summary>
#md # ```

zone_ranking_fig = plot_zone_ranking(
    zone_overview;
    patch_labels = zone_inputs.patch_labels
);
## The overview as displayed: the interval strings, without the numeric
## columns the figure reads.
zone_overview_display = let d = zone_overview[
        :,
        [:zone, :patch, :cases, :share, :R_T, :p_R_above_1, :delta_T],
    ]
    d[!, "R modelled separately"] = [
        w ? "yes" : "no" for w in zone_overview.walking
    ]
    d
end;

#md # ```@raw html
#md # </details>
#md # ```

zone_ranking_fig #hide

# The table gives the twenty highest-ranked zones: the confirmed cases to date, the zone's share of its patch's infections at the cut-off in percent, its reproduction number, the probability that it exceeds one, its log-transmission deviation at the cut-off and whether its reproduction number is modelled separately.
# A zone whose reproduction number is not modelled separately takes it from its province.
# The share, the reproduction number and the deviation are each a median with a 90% interval.
# Every zone is listed in the fold below it.

MarkdownTable(first(zone_overview_display, 20)) #hide

#md # ```@raw html
#md # <details><summary>All health zones</summary>
#md # ```

MarkdownTable(zone_overview_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Composition checks
#
# Whether the model reproduces each zone's observed share of its patch's confirmed cases and deaths is on the [in-sample checks](@ref zone-compositions) page.

# ## [Data currency](@id zone-data-currency)
#
# The zone blocks are the per-zone confirmed-case and confirmed-death tables of the situation reports.
# Each is read here against the cut-off, so a block that stops being picked up shows as a date before it rather than as a flat series.

#md # ```@raw html
#md # <details><summary>Zone block currency</summary>
#md # ```

## Every zone block whose last vintage falls before the cut-off. The grace
## the national table allows is not applied here: one vintage behind is
## already worth reading.
zone_currency = let
    status = stream_report_status(obs; stratum = :zone)
    behind = status[[ismissing(d) || d > 0 for d in status.days_since], :]
    isempty(behind) ?
        Markdown.parse(
            "Every zone block reports to the cut-off, $(obs.cutoff)."
        ) :
        MarkdownTable(
            DataFrame(
                "Block" => behind.label,
                "Last reported" => [
                    ismissing(d) ? "never" : string(d)
                    for d in behind.last_date
                ],
                "Days before cut-off" => [
                    ismissing(d) ? "-" : string(d)
                    for d in behind.days_since
                ]
            )
        )
end;

#md # ```@raw html
#md # </details>
#md # ```

zone_currency #hide

# ## Health-zone fit diagnostics
#
# The table gives the sampler diagnostics of the two zone fits: the worst R-hat and smallest effective sample sizes over every stored quantity but the zone reproduction number, and the divergences.
# Per chain it gives the fraction of iterations at the tree-depth cap, the energy fraction of missing information and the adapted step size.
# The per-zone R-hat and effective sample sizes of the cut-off reproduction number, share and deviation are in the fold.
# The rows whose reproduction number is modelled separately are the ones to read.

#md # ```@raw html
#md # <details><summary>Zone fit diagnostics</summary>
#md # ```

_zone_per_chain(v) = join(string.(round.(v; sigdigits = 3)), " / ")
function _zone_sampler_row(label, chn)
    d = zone_sampler_diagnostics(chn)
    return (
        fit = label, max_rhat = round(d.max_rhat; digits = 3),
        min_ess_bulk = round(d.min_ess_bulk; digits = 0),
        min_ess_tail = round(d.min_ess_tail; digits = 0),
        divergences = d.n_divergent,
        depth_cap = _zone_per_chain(d.depth_cap_fraction),
        ebfmi = _zone_per_chain(d.ebfmi),
        step_size = _zone_per_chain(d.step_size),
    )
end
zone_sampler_table = DataFrame(
    [
        _zone_sampler_row("health zones", chn_local),
        _zone_sampler_row("health zones (frozen)", frozen_local.chn),
    ]
);
zone_diagnostics = let d = zone_diagnostics_table(chn_local, zone_inputs)
    for c in names(d)[3:end]
        d[!, c] = round.(d[!, c]; digits = startswith(c, "rhat") ? 3 : 0)
    end
    DataFrame(
        [
            n == "walking" ?
                "R modelled separately" => [w ? "yes" : "no" for w in d[!, n]] :
                n => d[!, n]
                for n in names(d)
        ]
    )
end;

#md # ```@raw html
#md # </details>
#md # ```

MarkdownTable(zone_sampler_table) #hide

#md # ```@raw html
#md # <details><summary>Per-zone R-hat and effective sample sizes</summary>
#md # ```

MarkdownTable(zone_diagnostics) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The figure below sets the reproduction number implied by the zone stage's own patch trajectories against the one implied by the headline joint fit, nationally and for each patch.
# Both are computed from infections with the generation interval the zone stage fixes.
# Agreement says the melding stage has kept the joint's patch trajectories rather than moved them to fit the zone data.

#md # ```@raw html
#md # <details><summary>Reproduction number from the zone stage and the joint fit</summary>
#md # ```

## The implied reproduction number of each row of `draws` (draws × days).
function _implied_rt_matrix(draws::AbstractMatrix)
    m = similar(draws, Float64)
    for i in axes(draws, 1)
        m[i, :] .= implied_national_Rt(view(draws, i, :), zone_inputs.g)
    end
    return m
end
## The zone stage's deformed patch trajectories on one side and the joint's
## draws on the other, nationally then per patch. The joint's per-patch
## values are the grey references of the trajectory panels above.
zone_stage_rt = let I = zone_patch_infections(chn_local, zone_inputs)
    vcat([_implied_rt_matrix(sum(I))], _implied_rt_matrix.(I))
end;
joint_stage_rt = let draws = _patch_infection_draws
    national = Float64[
        sum(reshape(Float64.(v), N_PATCHES, obs.n)[:, t])
            for v in draws, t in 1:obs.n
    ]
    vcat([_implied_rt_matrix(national)], patch_implied_rt)
end;
zone_meld_rt_fig = plot_rt_zones(
    [m[:, zone_grid] for m in zone_stage_rt],
    vcat(["National"], zone_inputs.patch_labels),
    vcat([N_PATCHES + 1], 1:N_PATCHES);
    patch_labels = vcat(zone_inputs.patch_labels, ["National"]),
    patch_colours = [:firebrick, :steelblue, :seagreen, :darkorange, :black],
    dates = grid_date.(zone_grid), as_of_date = obs.cutoff,
    top = N_PATCHES + 1, ncols = 3, unit = "patch",
    reference_rt = [m[:, zone_grid] for m in joint_stage_rt],
    reference_label = "Headline joint fit",
    title = "Reproduction number from the zone stage and the joint fit"
);

#md # ```@raw html
#md # </details>
#md # ```

zone_meld_rt_fig #hide

# ## Health-zone parameters against their priors
#
# The table sets the posterior of each zone hyperparameter against its prior, with the ratio of their standard deviations.
# A ratio near one says the zone data add little to the prior.
# A ratio above one says the posterior is wider than the prior, which happens when the data move a parameter into the prior's wider tail.
# The drift scale is one per patch.
# The pair plot overlays the prior on the posterior of the scalar hyperparameters.

#md # ```@raw html
#md # <details><summary>Draw from the zone model's prior</summary>
#md # ```

prior_chn_zone = zone_prior_draws(zone_inputs);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Compute the zone prior and posterior table</summary>
#md # ```

## Each scalar zone hyperparameter with its table label and pair-plot axis
## label, kept where both chains carry it: the mixing and correlation
## blocks are sampled only when their inputs are on.
zone_hyper = [
    h for h in (
            (:region_sd_zone, "Level spread σ_level", "σ_level"),
            (:region_halflife_zone, "Deviation half-life (days)", "half-life (days)"),
            (:correlation_reference_zone, "Correlation ρ_corr", "ρ_corr"),
            (
                :zone_ascertainment_sd, "Ascertainment spread σ_ascertainment",
                "σ_ascertainment",
            ),
            (:zone_severity_sd, "Severity spread σ_severity", "σ_severity"),
            (:mixing_within_zone, "Within-patch mixing ε_within", "ε_within"),
            (:mixing_departure_zone, "Mixing departure τ_mix", "τ_mix"),
            (:composition_rho_zone, "Case composition ρ", "ρ"),
            (:composition_rho_death_zone, "Death composition ρ_death", "ρ_death"),
        )
        if BVDOutbreakSize._has_key(chn_local, h[1]) &&
        BVDOutbreakSize._has_key(prior_chn_zone, h[1])
]
zone_hyper_keys = first.(zone_hyper)
_hyper_draws(chn, k) = Float64.(vec(collect(chn[k])))
_drift_draws(chn, p) = Float64[
    v[p] for v in vec(collect(chn[:region_drift_sd_zone]))
]
function _prior_posterior_row(label, post, prior)
    f(x) = string(round(x; sigdigits = 3))
    ci(v) = string(
        f(median(v)), " (", f(quantile(v, 0.05)), "–",
        f(quantile(v, 0.95)), ")"
    )
    return (
        parameter = label, posterior = ci(post), prior = ci(prior),
        sd_ratio = round(std(post) / std(prior); digits = 2),
    )
end
zone_prior_table = DataFrame(
    vcat(
        [
            _prior_posterior_row(
                label, _hyper_draws(chn_local, k),
                _hyper_draws(prior_chn_zone, k)
            )
                for (k, label, _) in zone_hyper
        ],
        [
            _prior_posterior_row(
                "Drift scale σ_δ ($(zone_inputs.patch_labels[p]))",
                _drift_draws(chn_local, p), _drift_draws(prior_chn_zone, p)
            )
                for p in eachindex(zone_inputs.patch_labels)
        ]
    )
);
## The pair plot's axes take the symbols alone.
zone_hyper_pair_fig = plot_pair(
    chn_local, zone_hyper_keys;
    prior = prior_chn_zone,
    labels = Dict(k => sym for (k, _, sym) in zone_hyper)
);

#md # ```@raw html
#md # </details>
#md # ```

MarkdownTable(zone_prior_table) #hide

#md # ```@raw html
#md # <details><summary>Zone hyperparameter pair plot (prior overlaid)</summary>
#md # ```

zone_hyper_pair_fig #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Change over the past week
#
# The [health-zone model](@ref "Health-zone model") is melded onto the headline fit's patch infections one way, so the zone data do not update the national and province estimates.
# The comparison below reads the reproduction number of every zone walking in both fits at the frozen and live cut-offs, matched by key.
# The dot plot shows the fifteen zones the frozen fit ranks highest, the trajectories the twelve with most confirmed cases, and the table the ten with most confirmed cases.

#md # ```@raw html
#md # <details><summary>Zone cut-off summaries shared by the comparisons</summary>
#md # ```

## One row per zone in `zs` with the median and the 50% and 90% intervals
## of its draws (one vector per zone), in the schema the zone dot plots
## read, plus the zone's key to match variants by and its confirmed cases
## to rank by. A zone with no finite draw is left out.
function _zone_cutoff_summary(draws, inputs, zs = eachindex(draws))
    keep = [i for (i, v) in enumerate(draws) if any(isfinite, v)]
    z = zs[keep]
    tbl = zone_summary_table(
        [filter(isfinite, draws[i]) for i in keep],
        inputs.zone_labels[z], inputs.patch_of_zone[z]
    )
    tbl.key = inputs.zone_keys[z]
    tbl.cases = inputs.cumulative[z]
    return tbl
end

## The `top` zones by confirmed cases in the first variant, one column per
## variant for the reproduction number. Each cell is a median with its 90%
## interval. Variants are matched by zone key.
function _zone_comparison_table(variants; top::Integer = 10)
    function cell(t, key, d, scale)
        i = findfirst(==(key), t.key)
        i === nothing && return ""
        f(x) = string(round(scale * x; digits = d))
        return string(f(t.median[i]), " (", f(t.lo90[i]), "–", f(t.hi90[i]), ")")
    end
    base = first(last(first(variants)))
    order = sortperm(base.cases; rev = true)[1:min(top, size(base, 1))]
    df = DataFrame(
        "Zone" => base.label[order],
        "Province" => [PROVINCE_LABELS[p] for p in base.patch[order]],
        "Cases" => base.cases[order]
    )
    for (q, name, d, scale) in ((:R, "R", 2, 1),),
            (label, s) in variants

        haskey(s, q) || continue
        df[!, "$name ($label)"] = [
            cell(s[q], k, d, scale)
                for k in base.key[order]
        ]
    end
    return df
end

#md # ```@raw html
#md # </details>
#md # ```

#
# The frozen zone fit and the live one condition on different weeks of data and on different parent fits.
# The comparison reads each zone's reproduction number on the frozen cut-off day from both fits, alongside the live estimate at the current cut-off.
# A zone carries its own walk only once it has reported 30 confirmed cases, and the comparison is restricted to zones walking in both fits.
# The trajectory panels draw the frozen fit behind the live one for the twelve such zones with most confirmed cases.

#md # ```@raw html
#md # <details><summary>Zone reproduction numbers from the frozen and live fits</summary>
#md # ```

## Daily reproduction numbers rebuilt from both fits, and the zones walking
## in both as live index => frozen index, matched by key. The summaries
## take their labels and cases from the live inputs.
zone_rt_live = reconstruct_zone_rt(chn_local, zone_inputs)
zone_rt_frozen = reconstruct_zone_rt(frozen_local.chn, frozen_zone_inputs)
zone_week_pairs = [
    z => j
        for (z, k) in enumerate(zone_inputs.zone_keys)
        for j in (findfirst(==(k), frozen_zone_inputs.zone_keys),)
        if j !== nothing && zone_inputs.walking[z] &&
        frozen_zone_inputs.walking[j]
]
zone_week_variants = let zs = first.(zone_week_pairs), js = last.(zone_week_pairs),
        n_f = frozen_zone_inputs.n

    [
        "frozen fit at its cut-off" => (;
            R = _zone_cutoff_summary(
                [zone_rt_frozen[j][:, n_f] for j in js],
                zone_inputs, zs
            ),
        ),
        "live fit on the same day" => (;
            R = _zone_cutoff_summary(
                [zone_rt_live[z][:, n_f] for z in zs],
                zone_inputs, zs
            ),
        ),
        "live fit at its cut-off" => (;
            R = _zone_cutoff_summary(
                [zone_rt_live[z][:, end] for z in zs],
                zone_inputs, zs
            ),
        ),
    ]
end
zone_week_table = _zone_comparison_table(zone_week_variants);
zone_week_fig = plot_zone_comparison(
    [l => s.R for (l, s) in zone_week_variants];
    xlabel = "Reproduction number", reference_line = 1.0,
    title = "Zone reproduction number from the frozen and live fits"
);

## The frozen trajectories padded onto the live grid, undefined past the
## frozen cut-off, behind the live ones over the zone grid. The plot reads
## an undefined day as `missing`, where the reconstruction writes `NaN`.
zone_week_rt_fig = let grid = zone_inputs.t0:obs.n, n_f = frozen_zone_inputs.n
    asmissing(m) = replace(m, NaN => missing)
    frozen = map(zone_week_pairs) do (z, j)
        m = fill(NaN, size(zone_rt_frozen[j], 1), obs.n)
        m[:, 1:n_f] .= zone_rt_frozen[j]
        asmissing(m[:, grid])
    end
    zs = first.(zone_week_pairs)
    plot_rt_zones(
        [asmissing(zone_rt_live[z][:, grid]) for z in zs],
        zone_inputs.zone_labels[zs], zone_inputs.patch_of_zone[zs];
        patch_labels = zone_inputs.patch_labels,
        dates = grid_date.(grid), as_of_date = obs.cutoff,
        ranking = zone_inputs.cumulative[zs],
        ranking_label = "confirmed cases to date", top = 12,
        reference_rt = frozen, reference_label = "Frozen fit",
        title = "Zone reproduction number from the live fit, " *
            "with the frozen fit behind"
    )
end

#md # ```@raw html
#md # </details>
#md # ```

MarkdownTable(zone_week_table) #hide

#-

zone_week_fig #hide

#-

zone_week_rt_fig #hide

# ## Saving zone assets
#
# The [dashboard](@ref "Dashboard") map reads the per-zone estimates and time series written here, and the summary dashboard's health-zone table is built from them.
# The zone forecast figure is written by the [health-zone forecasts](@ref "Health-zone forecasts") page and the frozen zone forecast and its scores by the [health-zone forecast evaluation](@ref "Health-zone forecast evaluation") page.

#md # ```@raw html
#md # <details><summary>Write the zone dashboard and release assets</summary>
#md # ```

dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
## The health-zone maps and the one-week zone forecast for the summary
## dashboard, and the per-zone estimates the interactive map reads: one row per zone keyed as the
## geojson keys it, with the cases and deaths to date, the reproduction
## number and the chance it exceeds one, the one-week forecast and the
## share of the patch's infections, whether the zone's reproduction number
## is modelled separately (`walking`) and the data cut-off. A zone below the
## reporting floor carries no reproduction number.
CairoMakie.save(joinpath(dashboard_dir, "zone_rt_map.png"), zone_map_fig)
_zone_deaths = [
    let h = obs.zone_death_history
        haskey(h, prov) && haskey(h[prov], z) &&
            !isempty(h[prov][z].counts) ? Int(h[prov][z].counts[end]) :
            0
    end
        for (prov, z) in zip(
            zone_inputs.zone_province,
            zone_inputs.zone_names
        )
]
_zone_rt_stat(f) = [isempty(r) ? missing : f(r) for r in zone_RT_finite]
zone_estimates = DataFrame(
    zone = zone_map_keys,
    label = zone_inputs.zone_labels, province = zone_inputs.zone_province,
    patch = zone_inputs.patch_labels[zone_patch],
    cases = zone_inputs.cumulative, deaths = _zone_deaths,
    R_T_median = _zone_rt_stat(median),
    p_rt_above_one = [
        isempty(r) ? missing : mean(r .> 1) for r in zone_RT_finite
    ],
    R_T_lower = _zone_rt_stat(r -> quantile(r, 0.05)),
    R_T_upper = _zone_rt_stat(r -> quantile(r, 0.95)),
    forecast_median = [median(v) for v in zone_fc_draws.zones],
    forecast_lower = [quantile(v, 0.05) for v in zone_fc_draws.zones],
    forecast_upper = [quantile(v, 0.95) for v in zone_fc_draws.zones],
    share_median = [median(v) for v in zone_share_T],
    share_lower = [quantile(v, 0.05) for v in zone_share_T],
    share_upper = [quantile(v, 0.95) for v in zone_share_T],
    p_ge_1 = zone_fc_probs[:, 1], p_ge_5 = zone_fc_probs[:, 2],
    p_ge_10 = zone_fc_probs[:, 3], p_ge_20 = zone_fc_probs[:, 4],
    cases_last_7 = zone_recent[7], cases_last_14 = zone_recent[14],
    cases_last_28 = zone_recent[28],
    last_case_date = [ismissing(d) ? "" : string(d) for d in zone_last_case],
    walking = Int.(zone_inputs.walking),
    as_of = fill(string(obs.cutoff), length(zone_map_keys))
)
CSV.write(joinpath(dashboard_dir, "zone_estimates.csv"), zone_estimates)
## Each zone's daily reproduction number and weekly allocated confirmed
## cases, for the dashboard map's detail column.
_with_series(t, s) = (t[!, :series] .= s; t)
CSV.write(
    joinpath(dashboard_dir, "zone_timeseries.csv"),
    vcat(
        _with_series(
            rt_quantile_table(
                zone_rt_traj, zone_map_keys; cutoff = obs.cutoff,
                n = obs.n, from = zone_inputs.t0
            ), "rt"
        ),
        _with_series(
            weekly_count_table(
                zone_inputs.days, zone_inputs.counts, zone_map_keys;
                cutoff = obs.cutoff, n = zone_inputs.n
            ), "cases"
        );
        cols = :union
    )
)
open(joinpath(dashboard_dir, "zone_headline.md"), "w") do io
    print(io, zone_headline(zone_estimates; top = 10))
end
## The same frame into the release outputs.
output_dir = get(
    ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output")
)
mkpath(output_dir)
CSV.write(joinpath(output_dir, "zone_estimates.csv"), zone_estimates)

#md # ```@raw html
#md # </details>
#md # ```

# ---
#
# The full analysis code, data and model definitions are in the
# [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
# repository.
