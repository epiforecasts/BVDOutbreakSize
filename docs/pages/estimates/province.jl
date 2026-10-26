# # Province estimates
#
# This page gives the province estimates from the joint model, whose methods are on the [Methods](@ref "Methods") page.
# The [national estimates](@ref "National estimates") are the sum of the provinces here.
# The per-province forecast is on the [province forecasts](@ref "Province forecasts") page and its scoring is in the [forecast by province](@ref "Forecast by province").
#
# This page is generated from
# [`docs/pages/estimates/province.jl`](https://github.com/epiforecasts/BVDOutbreakSize/blob/main/docs/pages/estimates/province.jl).
# The model code it calls is in
# [`src/`](https://github.com/epiforecasts/BVDOutbreakSize/tree/main/src).
# See [aim and origins](@ref "Aim and origins") and [limitations](@ref "Limitations").

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations and the fit registry. See
## `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))
#-
## The fits this page reads, loaded from the cache here.
chn_joint = load_fit("joint");

#md # ```@raw html
#md # </details>
#md # ```

include(joinpath(pkgdir(BVDOutbreakSize), "docs", "front_matter.jl")) #hide
MarkdownTable(report_dates(obs.cutoff)) #hide

# ## Summary
#
#md # ```@raw html
#md # <details><summary>Compute the weekly cases and reproduction number by province</summary>
#md # ```

## The provincial Rt trajectories, rebuilt on the same grid as the
## reproduction number by province figure further down.
province_summary_rt = reconstruct_patch_rt(
    chn_joint;
    n = obs.n, n_patches = N_PATCHES, rt_start = _rt_start_plot,
    rt_walk_start = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n)
);
## The week-ahead province forecast, read from the same draws as the
## province forecasts page.
province_summary_forecast = forecast_provinces(
    fit_forecast("joint"); horizon = 7, n_patches = N_PATCHES
);
province_cases_rt = cases_rt_table(
    province_summary_rt, province_cases.days, province_cases.increments;
    cutoff = obs.cutoff, n = obs.n, forecast = province_summary_forecast
);
province_cases_rt_fig = plot_cases_rt(province_cases_rt);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

province_cases_rt_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Each province's confirmed cases over the most recent week against its reproduction number, with the week before and the forecast week ahead.
#md # Faint bars span the 90% credible interval of the reproduction number, and the forecast also carries a bar across its 90% interval for cases.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>Compute the province comparison</summary>
#md # ```

province_headline_md = patch_headline(chn_joint, N_PATCHES)
province_headline = Markdown.parse(province_headline_md);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

province_headline #hide

#md # ```@raw html
#md # <figcaption>
#md # The provinces compared from the joint posterior, as equal-tailed 90% credible intervals, with shares and probabilities computed draw by draw.
#md # The reproduction number and the relative case ascertainment are identified only as a product, and the per-province deaths break the tie.
#md # </figcaption>
#md # </figure>
#md # ```

# ### Maps

#md # ```@raw html
#md # <details><summary>Maps of infections, reproduction number and ascertainment</summary>
#md # ```

## The infections panel is also drawn alone for the summary page.
province_infections_panel = (;
    values = province_map_summary(chn_joint, :C_T_patch, N_PATCHES).values,
    title = "Infections to date", scale = log10,
    colorbar_label = "Infections (median)",
)
province_map_fig = plot_province_map(
    [
        province_infections_panel,
        (;
            province_map_summary(chn_joint, :R_T_patch, N_PATCHES)...,
            title = "Reproduction number at the cut-off",
            diverging_at = 1.0, colorbar_label = "R (median)",
        ),
        (;
            province_map_summary(
                chn_joint, :province_ascertainment, N_PATCHES
            )...,
            title = "Relative case ascertainment", diverging_at = 1.0,
            scale = log10, colorbar_label = "Relative to national (median)",
        ),
    ]
);

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_map_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Each province is shaded by its posterior median.
#md # The reproduction number at the cut-off and the relative case ascertainment are centred on one, and a province whose 90% credible interval spans one is washed out.
#md # </figcaption>
#md # </figure>
#md # </details>
#md # ```

# ### Detail by province
#
# Each province's case-fatality ratio is in the [case-fatality ratio by province](@ref "Case-fatality ratio by province").

#md # ```@raw html
#md # <details><summary>Compute the per-province tables and figure</summary>
#md # ```

## The tables below and the pair-plot dropdowns further down name the
## provinces in this order.
@assert PROVINCE_LABELS[1:N_PATCHES] == [
    "Ituri", "Nord-Kivu", "Haut-Uele",
    "Sud-Kivu", "Tshopo", "Bas-Uele", "Sud-Ubangi", "Nord-Ubangi",
]
province_detail_tables = [
    patch_summary_table(chn_joint, N_PATCHES; patch = p) for p in 1:N_PATCHES
];
province_detail_fig = plot_patch_summary(chn_joint, N_PATCHES);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

province_detail_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Each province's estimates as medians with 30%, 60% and 90% credible intervals, repeated in the tables below.
#md # The log-Rt deviation is from the national trend, the drift is the deviation's walk scale and the contrast is against Ituri.
#md # </figcaption>
#md # </figure>
#md # ```

# #### Ituri

province_detail_tables[1] #hide

# #### Nord-Kivu

province_detail_tables[2] #hide

# #### Haut-Uele

province_detail_tables[3] #hide

# #### Sud-Kivu

province_detail_tables[4] #hide

# #### Tshopo

province_detail_tables[5] #hide

# #### Bas-Uele

province_detail_tables[6] #hide

# #### Sud-Ubangi

province_detail_tables[7] #hide

# #### Nord-Ubangi

province_detail_tables[8] #hide

# #### Isolation beds

#md # ```@raw html
#md # <details><summary>Province bed table</summary>
#md # ```

province_bed_overview = province_bed_table(chn_joint, N_PATCHES);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

province_bed_overview #hide

#md # ```@raw html
#md # <figcaption>
#md # Beds, bed demand, occupied beds, utilisation and patients above the beds by province at the cut-off, from the occupancy and bed counts in the situation reports.
#md # Each province's beds are its modelled share of the national capacity, and its occupied beds are its demand share of the national modelled occupancy.
#md # Patients above the beds are the occupied beds in excess of them, and the provinces sum to the national figures.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Size and infections

#md # ```@raw html
#md # <details><summary>Modelled infections by province</summary>
#md # ```

province_infections_fig = plot_infections_patches(
    chn_joint;
    n = obs.n, seeding = obs.seeding, n_patches = N_PATCHES
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

province_infections_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Modelled infections by province, daily (top) and cumulative (bottom).
#md # </figcaption>
#md # </figure>
#md # ```

# The provinces are coupled by the gravity kernel in the [mixing and importation](@ref "Mixing and importation") Methods section.
# The distances between the patches' population centres run from 191 km between Sud-Ubangi and Nord-Ubangi to 1 220 km between Sud-Kivu and Sud-Ubangi.
# Ituri is 324 km from Nord-Kivu and 216 km from Haut-Uele.
# About half of what leaves Nord-Kivu lands in Sud-Kivu, 260 km away and the second most populous of the eight.
# Nord-Ubangi sends about a third of its exports to Sud-Ubangi, its nearest neighbour.

#md # ```@raw html
#md # <details><summary>Importation intensity and imports by province</summary>
#md # ```

importation_table = summary_table(
    chn_joint, [:importation_epsilon];
    digits = 4,
    labels = Dict(:importation_epsilon => "Importation intensity")
);

province_imports_fig = plot_imports_patches(
    chn_joint;
    n = obs.n, seeding = obs.seeding, n_patches = N_PATCHES
);

#md # ```@raw html
#md # </details>
#md # ```

importation_table #hide

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_imports_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Imported infections by province.
#md # Every arrival is debited from its origin on the same day, so the figure shows where infection occurred rather than extra infection.
#md # </figcaption>
#md # </figure>
#md # ```

# ## [Point prevalence](@id province-prevalence)
#
# Point prevalence is the number of people infected at the cut-off who have not yet died, recovered or been detected as a suspected case.
# We split it into people infected but not yet symptomatic and people with symptoms who have not yet been detected.
# Together they are the infections in the community.
# The [point prevalence](@ref methods-prevalence) section of the Methods gives the calculation.
# A province's undetected count depends on its suspected-case detection probability.
# The province compositions see that probability only as a product with case confirmation, so the province occupancy splits and its prior set it.
# The table gives each rate per 100,000 residents and the number of people in the community, as medians with 90% credible intervals.

#md # ```@raw html
#md # <details><summary>Compute point prevalence by province</summary>
#md # ```

province_prevalence = prevalence_overview(
    patch_prevalence(chn_joint, N_PATCHES),
    PROVINCE_POPULATIONS[1:N_PATCHES];
    labels = PROVINCE_LABELS[1:N_PATCHES], patch = collect(1:N_PATCHES)
);
province_prevalence_fig = plot_prevalence_ranking(
    province_prevalence; title = "Point prevalence by province"
);

#md # ```@raw html
#md # </details>
#md # ```

MarkdownTable(prevalence_table(province_prevalence; area = "Province")) #hide

#-

province_prevalence_fig #hide

# A province rate is an average over zones with different levels of transmission.
# The [health-zone estimates](@ref "Health-zone estimates") give the same rates zone by zone.

# ## Reproduction number by province

# The national reproduction number is in [reproduction number over time](@ref "Reproduction number over time").

#md # ```@raw html
#md # <details><summary>Reproduction number by province</summary>
#md # ```

province_rt_fig = plot_rt_patches(
    chn_joint;
    n = obs.n, breakpoint = _BREAKPOINT,
    n_patches = N_PATCHES,
    rt_start = _rt_start_plot,
    rt_walk_start = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n),
    display_start = _rt_start_plot,
    as_of_date = string(obs.cutoff), seeding = obs.seeding,
    ramp = RT_INTERVENTION_RAMP
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

province_rt_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The reproduction number of each province, with the national trajectory in grey.
#md # The deviations sum to zero, so the grey band is the incidence-weighted middle of the panels rather than any one province.
#md # Sud-Kivu, Tshopo, Bas-Uele, Sud-Ubangi and Nord-Ubangi hold few confirmed cases, so their panels are carried mostly by the deviation prior and their widths are not measurements.
#md # </figcaption>
#md # </figure>
#md # ```

# The spread of the provinces' log-Rt deviations is the spatial diagnostic.
# Its prior gives a 26% probability that the Ituri to Nord-Kivu ratio moves by more than 25% over the window.
# With eight patches and five of them carrying little signal, the cross-province correlation is not identified and tracks its prior.

#md # ```@raw html
#md # <details><summary>Spatial hyperparameter summary table</summary>
#md # ```

## Labels for the spatial hyperparameters, shared by this table and the pair
## plot against their prior below.
spatial_labels = Dict(
    :region_sd => "Rt deviation spread",
    :region_halflife => "Rt deviation half-life (days)",
    :region_corr_primary_secondary => "Ituri-N.Kivu Rt correlation",
    :province_suspected_ascertainment_sd => "Detection spread (logit)",
    :province_ascertainment_sd => "Confirmation spread",
    :importation_epsilon => "Importation intensity",
    :importation_destination_sd => "Importation destination spread",
    :importation_origin_destination_correlation =>
        "Origin-destination correlation",
    :importation_flow_sd => "Importation flow spread",
    :importation_reciprocity => "Importation reciprocity",
    :province_cfr_sd => "Lethality spread (logit)",
    :province_death_ascertainment_sd => "Death-confirmation spread"
)
spatial_hyper_table = summary_table(
    chn_joint,
    [
        :region_sd, :region_halflife, :region_corr_primary_secondary,
        :province_suspected_ascertainment_sd, :province_ascertainment_sd,
    ];
    digits = 3, labels = spatial_labels
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

spatial_hyper_table #hide

#md # ```@raw html
#md # <figcaption>
#md # The spatial hyperparameters of the province deviations.
#md # The correlation is between deviation innovations that sum to zero, so it leans negative, and with equal spreads each province's correlations with the other seven average −1/7.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Case-fatality ratio by province
#
# The case-fatality ratio of each province is set against the national ratios in the [confirmed case-fatality ratio](@ref "Confirmed case-fatality ratio").
# Only the product of a province's lethality and death confirmation is identified, and the [province compositions](@ref methods-province-compositions) Methods section says how the priors split it.

#md # ```@raw html
#md # <details><summary>Province case-fatality spread</summary>
#md # ```

province_cfr_spread = summary_table(
    chn_joint,
    [:province_cfr_sd, :province_death_ascertainment_sd];
    digits = 3, labels = spatial_labels
);

#md # ```@raw html
#md # </details>
#md # ```

province_cfr_spread #hide

#md # ```@raw html
#md # <details><summary>Province case-fatality table</summary>
#md # ```

confirmed_cfr = delay_corrected_confirmed_cfr(
    chn_joint;
    obs_confirmed = obs.confirmed_cases,
    obs_confirmed_deaths = obs.confirmed_deaths
);
province_cfr = province_cfr_table(
    chn_joint, confirmed_cfr;
    province_cases = vec(sum(province_cases.increments; dims = 2)),
    province_deaths = vec(sum(province_deaths.increments; dims = 2)),
    n_patches = N_PATCHES
);

#md # ```@raw html
#md # </details>
#md # ```

province_cfr #hide

# ## Province parameters against their priors
#
# The pair plots below set the province parameters against their prior, drawn from the province model.

#md # ```@raw html
#md # <details><summary>Draw from the patch model's prior</summary>
#md # ```

prior_patch_chn = patch_prior_draws(obs);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Spatial hyperparameter pair plot (prior overlaid)</summary>
#md # ```

spatial_pair_fig = plot_pair(
    chn_joint,
    [
        :region_sd, :region_halflife, :region_corr_primary_secondary,
        :province_suspected_ascertainment_sd, :province_ascertainment_sd,
        :importation_epsilon, :importation_destination_sd,
        :importation_origin_destination_correlation, :importation_flow_sd,
        :importation_reciprocity, :province_cfr_sd,
        :province_death_ascertainment_sd,
    ];
    prior = prior_patch_chn, labels = spatial_labels
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

spatial_pair_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Joint posterior of the spatial hyperparameters, with the prior overlaid.
#md # These are the spread, half-life and correlation of the Rt deviations, the spreads of suspected-case detection and case confirmation, the importation intensity, its destination and flow spreads and their correlations, and the spreads of lethality and death confirmation.
#md # </figcaption>
#md # </figure>
#md # ```

# Each province's pair plot sets the reproduction number at the cut-off against the relative case ascertainment, the case-fatality ratio against the relative death confirmation, and the importation intensity against the export weight relative to Ituri.
# Each pair is identified only as a product, so a ridge between the two is expected and the prior sets the position along it.

#md # ```@raw html
#md # <details><summary>Compute the per-province pair plots</summary>
#md # ```

province_pair_labels = Dict(
    :R_T_patch => "Reproduction number",
    :province_ascertainment => "Case ascertainment",
    :CFR_patch => "Case-fatality ratio",
    :province_death_ascertainment => "Death confirmation",
    :importation_epsilon_patch => "Importation intensity",
    :export_weight => "Export weight"
)
## Ituri is the export reference, its weight fixed at one, so its pair plot
## leaves the weight out.
province_pair_figs = [
    plot_pair(
        chn_joint,
        [
            :R_T_patch, :province_ascertainment,
            :CFR_patch, :province_death_ascertainment,
            :importation_epsilon_patch,
            (p == 1 ? () : (:export_weight,))...,
        ];
        patch = p, prior = prior_patch_chn, labels = province_pair_labels
    )
        for p in 1:N_PATCHES
];

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Ituri pair plot (prior overlaid)</summary>
#md # ```

province_pair_figs[1] #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Nord-Kivu pair plot (prior overlaid)</summary>
#md # ```

province_pair_figs[2] #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Haut-Uele pair plot (prior overlaid)</summary>
#md # ```

province_pair_figs[3] #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Sud-Kivu pair plot (prior overlaid)</summary>
#md # ```

province_pair_figs[4] #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Tshopo pair plot (prior overlaid)</summary>
#md # ```

province_pair_figs[5] #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Bas-Uele pair plot (prior overlaid)</summary>
#md # ```

province_pair_figs[6] #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Sud-Ubangi pair plot (prior overlaid)</summary>
#md # ```

province_pair_figs[7] #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Nord-Ubangi pair plot (prior overlaid)</summary>
#md # ```

province_pair_figs[8] #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Composition checks
#
# Whether the model reproduces each province's observed share of the national total is on the [in-sample checks](@ref province-compositions) page.

# ## Data currency
#
# The province blocks are the provincial tables of the situation reports: confirmed cases, confirmed deaths and the 24-hour laboratory volumes.
# A block that stops updating shows its last date before the cut-off rather than a flat series.

#md # ```@raw html
#md # <details><summary>Province block currency</summary>
#md # ```

## Every province block whose last vintage falls before the cut-off, from
## the shared stream registry rather than a per-page list of dates. The
## grace the national table allows is not applied here: one vintage behind
## is already worth reading.
province_currency = let
    status = stream_report_status(obs; stratum = :province)
    behind = status[[ismissing(d) || d > 0 for d in status.days_since], :]
    isempty(behind) ?
        Markdown.parse(
            "Every province block reports to the cut-off, $(obs.cutoff)."
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

province_currency #hide

# ## Province fit diagnostics
#
# The national [fit diagnostics](@ref national-fit-diagnostics) summarise the same joint fit as a whole, and its [breakdown by parameter](@ref "Fit diagnostics by parameter") covers every parameter.

#md # ```@raw html
#md # <details><summary>Province-level parameter diagnostics</summary>
#md # ```

## The quantities the sections above hold one per province: the
## reproduction number and infections to date from the renewal equation,
## the case-fatality ratio and the two ascertainment ratios the province
## compositions identify, and the importation intensity and export weight
## the per-province pair plots show.
province_diagnostic_quantities = [
    :R_T_patch, :C_T_patch, :CFR_patch, :province_ascertainment,
    :province_death_ascertainment, :importation_epsilon_patch, :export_weight,
]
province_diagnostics_labels = merge(
    province_pair_labels, Dict(:C_T_patch => "Infections to date")
)
province_diagnostics_frame = let df = parameter_diagnostics(chn_joint)
    df[in.(Symbol.(df.parameter), Ref(province_diagnostic_quantities)), :]
end
province_diagnostics_table = MarkdownTable(
    family_diagnostics_table(
        province_diagnostics_frame;
        n = length(province_diagnostic_quantities),
        labels = province_diagnostics_labels
    )
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

province_diagnostics_table #hide

#md # ```@raw html
#md # <figcaption>
#md # Sampler diagnostics for the province-level parameters of the joint fit: the reproduction number, infections to date, case-fatality ratio, case and death ascertainment ratios, importation intensity and export weight.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Saving province assets
#
# The summary dashboard shows the province comparison, the province map, the reproduction number by province and the weekly cases against the reproduction number, so they are written here.

#md # ```@raw html
#md # <details><summary>Write the province dashboard assets</summary>
#md # ```

dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
CairoMakie.save(joinpath(dashboard_dir, "rt_provinces.png"), province_rt_fig)
## The summary's cases against R figure is drawn at the page width with
## larger labels.
CairoMakie.save(
    joinpath(dashboard_dir, "cases_rt_provinces.png"),
    plot_cases_rt(province_cases_rt; size = (680, 560), fontsize = 13)
)
## The summary's province map: infections to date alone.
CairoMakie.save(
    joinpath(dashboard_dir, "province_map.png"),
    plot_province_map(
        [province_infections_panel]; ncols = 1, panel_size = (680, 500)
    )
)
open(joinpath(dashboard_dir, "provinces.md"), "w") do io
    print(io, province_headline_md)
end
## The summary's compact province table: infections to date, R at the
## cut-off and the confirmed cases over the latest week of reports, the
## same week the cases against R figure ends on.
province_recent_cases = let last_day = min(last(province_cases.days), obs.n)
    weekly_count_table(
        province_cases.days, province_cases.increments[1:N_PATCHES, :],
        PROVINCE_LABELS[1:N_PATCHES];
        cutoff = obs.cutoff - Day(obs.n - last_day), n = last_day, weeks = 1
    ).count
end
open(joinpath(dashboard_dir, "provinces_summary.md"), "w") do io
    print(
        io,
        province_summary_table(
            chn_joint, N_PATCHES; recent_cases = province_recent_cases
        )
    )
end

#md # ```@raw html
#md # </details>
#md # ```

# ---
#
# The full analysis code, data and model definitions are in the [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize) repository.
