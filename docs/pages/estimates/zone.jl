# # Health-zone estimates
#
# This page splits each patch of the joint fit across its health zones, giving each zone's reproduction number, share of the patch, one-week forecast and probability of at least a few cases.
# The [health-zone model](@ref "Health-zone model") on the methods page gives the maths, and the interactive map is on the [dashboard](@ref "Dashboard").
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
## the zone fit melded from it, and the frozen zone fit the sampler
## diagnostics read.
chn_joint = load_fit("joint");
chn_local = load_fit("local");
frozen_local = load_fit("local_frozen_validation");

## The zone stage's fixed inputs with the joint's forecast, from the shared
## setup, so the health-zone forecast page draws the same forecast from the
## same inputs.
zone_inputs = zone_stage_inputs(; forecast = true);
zone_patch = zone_inputs.patch_of_zone;

#md # ```@raw html
#md # </details>
#md # ```

include(joinpath(pkgdir(BVDOutbreakSize), "docs", "front_matter.jl")) #hide
MarkdownTable(report_dates(obs.cutoff)) #hide

# ## Estimates by zone
#
# The interactive map on the [dashboard](@ref "Dashboard") adds the probability of at least $K$ cases at a chosen $K$ and a filter by recent cases.

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
## carried, and the cut-off values the map and table read. A zone's
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
#md # <figure class="bvd-figure">
#md # ```

zone_map_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The reproduction number at the cut-off, the bounds of the 90% interval on next week's confirmed cases, and the confirmed cases to date, by health zone.
#md # A zone whose reproduction number interval straddles one is washed towards white, and zones with no confirmed case or too few infections for a reproduction number are grey.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>Compute the weekly cases and reproduction number by zone</summary>
#md # ```

zone_cases_rt = cases_rt_table(
    zone_rt_traj, zone_inputs.days, zone_inputs.counts;
    cutoff = obs.cutoff, n = obs.n, areas = zone_inputs.zone_labels,
    patch = zone_patch
);
zone_cases_rt_fig = plot_cases_rt(
    zone_cases_rt; patch_labels = zone_inputs.patch_labels, top = 12,
    unit = "health zone"
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

zone_cases_rt_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The twelve zones with the most confirmed cases over the past two weeks, each against its reproduction number, with the week before.
#md # Faint bars span the 90% credible interval of the reproduction number.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>Zone overview table</summary>
#md # ```

## The overview as displayed: the interval strings, without the numeric
## columns.
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
#md # <figure class="bvd-figure">
#md # ```


MarkdownTable(first(zone_overview_display, 20)) #hide

#md # ```@raw html
#md # <figcaption>
#md # The twenty highest-ranked zones, with every zone in the fold below.
#md # The share is of the patch's infections at the cut-off, in percent, and the share, reproduction number and deviation are medians with 90% intervals.
#md # A zone whose reproduction number is not modelled separately takes it from its province.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>All health zones</summary>
#md # ```

MarkdownTable(zone_overview_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## [Point prevalence](@id zone-prevalence)
#
# The zone estimates split the infections in the community at the cut-off in the same way as the [province estimates](@ref province-prevalence).
# The [point prevalence](@ref methods-prevalence) section of the Methods gives the calculation.
# The figure and table show the twenty zones with the highest rate in the community per 100,000 residents.
# Every zone is listed in the fold below the table.

#md # ```@raw html
#md # <details><summary>Compute point prevalence by zone</summary>
#md # ```

zone_population = Dict(
    string(z.province, ".", z.zone) => z.population for z in load_health_zones()
);
zone_prevalence_overview = prevalence_overview(
    zone_prevalence(chn_local, zone_inputs, chn_joint),
    [zone_population[k] for k in zone_inputs.zone_keys];
    labels = zone_inputs.zone_labels, patch = zone_patch
);
zone_prevalence_fig = plot_prevalence_ranking(
    zone_prevalence_overview;
    patch_labels = zone_inputs.patch_labels, max_rows = 20,
    title = "Point prevalence by health zone"
);

#md # ```@raw html
#md # </details>
#md # ```

zone_prevalence_fig #hide

#-

MarkdownTable(
    prevalence_table(zone_prevalence_overview; max_rows = 20, area = "Zone")
) #hide

#md # ```@raw html
#md # <details><summary>Point prevalence in every health zone</summary>
#md # ```

MarkdownTable(prevalence_table(zone_prevalence_overview; area = "Zone")) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Composition checks
#
# Whether the model reproduces each zone's observed share of its patch's confirmed cases and deaths is on the [in-sample checks](@ref zone-compositions) page.

# ## [Data currency](@id zone-data-currency)
#
# The zone blocks are the per-zone confirmed-case and confirmed-death tables of the situation reports.
# A block that stops being picked up shows its last date before the cut-off rather than a flat series.

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
#md # <figure class="bvd-figure">
#md # ```

MarkdownTable(zone_sampler_table) #hide

#md # ```@raw html
#md # <figcaption>
#md # Sampler diagnostics of the two zone fits.
#md # The worst R-hat and smallest effective sample sizes cover every stored quantity except the zone reproduction number.
#md # Per chain, the table also gives the fraction of iterations at the tree-depth cap, the energy fraction of missing information and the adapted step size.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Per-zone R-hat and effective sample sizes</summary>
#md # <figure class="bvd-figure">
#md # ```

MarkdownTable(zone_diagnostics) #hide

#md # ```@raw html
#md # <figcaption>
#md # Per-zone R-hat and effective sample sizes of the cut-off reproduction number, share and deviation.
#md # The rows whose reproduction number is modelled separately are the ones to read.
#md # </figcaption>
#md # </figure>
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Reproduction number from the zone stage and the joint fit</summary>
#md # ```

## Each patch's implied reproduction number from the joint draws with the
## joint posterior mean generation interval, so the patch line is close to
## the quantity the zone values average to.
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
## The implied reproduction number of each row of `draws` (draws × days).
function _implied_rt_matrix(draws::AbstractMatrix)
    m = similar(draws, Float64)
    for i in axes(draws, 1)
        m[i, :] .= implied_national_Rt(view(draws, i, :), zone_inputs.g)
    end
    return m
end
## The zone stage's deformed patch trajectories on one side and the joint's
## draws on the other, nationally then per patch.
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
    patch_colours = vcat(BVDOutbreakSize._ZONE_PATCH_COLOURS, [:black]),
    dates = grid_date.(zone_grid), as_of_date = obs.cutoff,
    top = N_PATCHES + 1, ncols = 3, unit = "patch",
    reference_rt = [m[:, zone_grid] for m in joint_stage_rt],
    reference_label = "Headline joint fit",
    title = "Reproduction number from the zone stage and the joint fit"
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

zone_meld_rt_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The reproduction number implied by the zone stage's patch trajectories against that implied by the joint fit, nationally and for each patch, both using the joint posterior mean generation interval.
#md # Agreement means the melding has kept the joint's patch trajectories rather than moved them to fit the zone data.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Health-zone parameters against their priors

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
            (
                :zone_severity_sd, "Severity spread σ_severity (logit)",
                "σ_severity",
            ),
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
#md # <figure class="bvd-figure">
#md # ```

MarkdownTable(zone_prior_table) #hide

#md # ```@raw html
#md # <figcaption>
#md # Posterior against prior for each zone hyperparameter, with the ratio of their standard deviations, and the drift scale given per patch.
#md # A ratio near one means the zone data add little to the prior, and a ratio above one means the data move the parameter into the prior's wider tail.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Zone hyperparameter pair plot (prior overlaid)</summary>
#md # <figure class="bvd-figure">
#md # ```

zone_hyper_pair_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Joint posterior of the scalar zone hyperparameters, with the prior overlaid.
#md # </figcaption>
#md # </figure>
#md # </details>
#md # ```

# ## Change over the past week
#
# How each zone's reproduction number has moved between the fit a week earlier and the current one is on the [in-sample checks](@ref zone-change-past-week) page.

# ## Saving zone assets
#
# The [dashboard](@ref "Dashboard") map reads the per-zone estimates and time series written here, and the summary dashboard's health-zone table is built from them.
# The zone forecast assets are written by the [health-zone forecasts](@ref "Health-zone forecasts") and [health-zone forecast evaluation](@ref "Health-zone forecast evaluation") pages.

#md # ```@raw html
#md # <details><summary>Write the zone dashboard and release assets</summary>
#md # ```

dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
## The health-zone maps, the weekly cases against R and the one-week zone forecast for the summary
## dashboard, and the per-zone estimates the interactive map reads: one row per zone keyed as the
## geojson keys it, with the cases and deaths to date, the reproduction
## number and the chance it exceeds one, the one-week forecast and the
## share of the patch's infections, whether the zone's reproduction number
## is modelled separately (`walking`) and the data cut-off. A zone below the
## reporting floor carries no reproduction number.
CairoMakie.save(joinpath(dashboard_dir, "zone_rt_map.png"), zone_map_fig)
CairoMakie.save(
    joinpath(dashboard_dir, "cases_rt_zones.png"),
    plot_cases_rt(
        zone_cases_rt; patch_labels = zone_inputs.patch_labels, top = 12,
        unit = "health zone", size = (680, 560), fontsize = 12
    )
)
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
## Each zone's daily reproduction number, weekly allocated confirmed cases
## and weekly confirmed deaths, for the dashboard map's detail column. A
## zone's deaths are its own cumulative series, so each zone is binned on
## its own vintages.
_with_series(t, s) = (t[!, :series] .= s; t)
_zone_weekly_deaths = let h = obs.zone_death_history
    out = DataFrame(area = String[], date = Date[], count = Int[])
    for (prov, z, key) in zip(
            zone_inputs.zone_province, zone_inputs.zone_names, zone_map_keys
        )
        haskey(h, prov) && haskey(h[prov], z) &&
            !isempty(h[prov][z].counts) || continue
        c = h[prov][z].counts
        append!(
            out,
            weekly_count_table(
                h[prov][z].days, reshape([c[1]; max.(diff(c), 0)], 1, :),
                [key]; cutoff = obs.cutoff, n = obs.n
            )
        )
    end
    out
end
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
        ),
        _with_series(_zone_weekly_deaths, "deaths");
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
# The full analysis code, data and model definitions are in the [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize) repository.
