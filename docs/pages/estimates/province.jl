# # Province estimates
#
# This page gives the province estimates from the joint model.
# The methods for this model are on the [Methods](@ref "Methods") page.
#
# The model runs one renewal equation per province and fits the national streams against the summed provinces, so the [national estimates](@ref "National estimates") are the sum of the provinces here.
# The per-province case-fatality ratio sits with the national one in the [confirmed case-fatality ratio](@ref "Confirmed case-fatality ratio").
# The per-province forecast is on the [province forecasts](@ref "Province forecasts") page and its scoring is in the [forecast by province](@ref "Forecast by province").
# Each needs a quantity the page it sits on already computes.
#
# This page is generated from
# [`docs/pages/estimates/province.jl`](https://github.com/epiforecasts/BVDOutbreakSize/blob/main/docs/pages/estimates/province.jl).
# The model code it calls is in
# [`src/`](https://github.com/epiforecasts/BVDOutbreakSize/tree/main/src).
# See [aim and origins](@ref "Aim and origins") and [limitations](@ref "Limitations").

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations, the fit registry and every model fit
## (loaded from the content-addressed cache). See `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))

#md # ```@raw html
#md # </details>
#md # ```

include(joinpath(pkgdir(BVDOutbreakSize), "docs", "front_matter.jl")) #hide
MarkdownTable(report_dates(obs.cutoff)) #hide

# ## Summary
#
# The bullets below compare the provinces, from the joint posterior.
# Each range is an equal-tailed 90% credible interval, and the shares and probabilities are computed draw by draw.
# The reproduction number and the relative case ascertainment are identified only as a product, and the per-province deaths break the tie.

#md # ```@raw html
#md # <details><summary>Compute the province comparison</summary>
#md # ```

province_headline_md = patch_headline(chn_joint, N_PATCHES)
province_headline = Markdown.parse(province_headline_md);

#md # ```@raw html
#md # </details>
#md # ```

province_headline #hide

# ### Detail by province
#
# Each province's own estimates are below, as equal-tailed 30%, 60% and 90% credible intervals.

#md # ```@raw html
#md # <details><summary>Compute the per-province ranges</summary>
#md # ```

province_detail_headline = Markdown.parse(
    patch_detail_headline(chn_joint, N_PATCHES)
);

#md # ```@raw html
#md # </details>
#md # ```

province_detail_headline #hide

# The table below gives the same intervals for each province, including its log-Rt deviation from the trend, the deviation's walk scale and the contrast against Ituri.

#md # ```@raw html
#md # <details><summary>Per-province summary table</summary>
#md # ```

province_detail_table = patch_summary_table(chn_joint, N_PATCHES);

province_detail_table #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Per-province summary figure</summary>
#md # ```

province_detail_fig = plot_patch_summary(chn_joint, N_PATCHES);

province_detail_fig #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Size and infections

# The national outbreak size in the [joint model estimates](@ref "Joint model estimates") is the sum of the four patches' renewal equations.
# The figure below shows the modelled infections behind those totals, daily on the top row and cumulative on the bottom.

#md # ```@raw html
#md # <details><summary>Modelled infections by province</summary>
#md # ```

province_infections_fig = plot_infections_patches(
    chn_joint;
    n = obs.n, seeding = obs.seeding, n_patches = N_PATCHES
);

#md # ```@raw html
#md # </details>
#md # ```

province_infections_fig #hide

# The provinces are coupled by a gravity kernel weighted by destination population, described in the [mixing and importation](@ref "Mixing and importation") Methods section, with its intensity estimated.
# Every arrival is debited from its origin the same day, so the figure reads as where infection occurred rather than as extra infection.
# The distances between the patch capitals are 379 km from Bunia to Goma, 322 km from Bunia to Isiro and 206 km from Goma to the pooled patch's centre, so most of what leaves Nord-Kivu lands in the pooled patch.

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

#-

province_imports_fig #hide

# ## Reproduction number by province

# The national reproduction number is in [reproduction number over time](@ref "Reproduction number over time").
# Below it is split by province, one panel per province with the national trajectory in grey behind it.
# The deviations sum to zero, so the grey band is the incidence-weighted middle of the panels rather than any one province.
# A panel tracking grey says that province moves with the national trajectory.
# The pooled patch holds almost no confirmed cases, so its panel is carried by the deviation prior and its width is not a measurement.

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
#md # ```

province_rt_fig #hide

# The spread of the provinces' log-Rt deviations is the spatial diagnostic.
# The prior admits real divergence, with a 31% prior probability that the Ituri to Nord-Kivu ratio moves by more than 25% over the window, so a shrunken posterior is a finding rather than an artefact of the prior.
# With four patches and the pooled one carrying almost no signal the cross-province correlation is not identified, and it tracks its prior.

#md # ```@raw html
#md # <details><summary>Spatial hyperparameter summary table</summary>
#md # ```

## Labels for the spatial hyperparameters, shared by this table and the pair
## plot against their prior below.
spatial_labels = Dict(
    :region_sd => "Rt deviation spread",
    :region_halflife => "Rt deviation half-life (days)",
    :region_corr_primary_secondary => "Ituri-N.Kivu Rt correlation",
    :province_ascertainment_sd => "Ascertainment spread",
    :province_testing_coefficient => "Testing effect on ascertainment",
    :importation_epsilon => "Importation intensity",
    :province_cfr_sd => "Lethality spread",
    :province_death_ascertainment_sd => "Death-confirmation spread"
)
spatial_hyper_table = summary_table(
    chn_joint,
    [
        :region_sd, :region_halflife, :region_corr_primary_secondary,
        :province_ascertainment_sd, :province_testing_coefficient,
    ];
    digits = 3, labels = spatial_labels
);

#md # ```@raw html
#md # </details>
#md # ```

spatial_hyper_table #hide

# ## Province parameters against their priors
#
# The pair plots below set the posterior of the province parameters against their prior.
# The prior is the four-patch model's, as the shared prior draws carry no province parameters.

#md # ```@raw html
#md # <details><summary>Draw from the patch model's prior</summary>
#md # ```

prior_patch_chn = patch_prior_draws(obs);

#md # ```@raw html
#md # </details>
#md # ```

# The first pair plot covers the spatial hyperparameters: the spread, half-life and correlation of the Rt deviations, the spread of case ascertainment and its testing effect, the importation intensity, and the spreads of lethality and death confirmation.

#md # ```@raw html
#md # <details><summary>Spatial hyperparameter pair plot (prior overlaid)</summary>
#md # ```

spatial_pair_fig = plot_pair(
    chn_joint,
    [
        :region_sd, :region_halflife, :region_corr_primary_secondary,
        :province_ascertainment_sd, :province_testing_coefficient,
        :importation_epsilon, :province_cfr_sd,
        :province_death_ascertainment_sd,
    ];
    prior = prior_patch_chn, labels = spatial_labels
);

#md # ```@raw html
#md # </details>
#md # ```

spatial_pair_fig #hide

# The pair plots below take one province at a time.
# Each sets the reproduction number at the cut-off against the relative case ascertainment, and the case-fatality ratio against the relative death confirmation.
# Each pair is identified only as a product, so a ridge between the two is expected and its position along the ridge is set by the prior.

#md # ```@raw html
#md # <details><summary>Compute the per-province pair plots</summary>
#md # ```

## The dropdowns below name the provinces in this order.
@assert PROVINCE_LABELS[1:N_PATCHES] ==
    ["Ituri", "Nord-Kivu", "Haut-Uele", "Other provinces"]
province_pair_labels = Dict(
    :R_T_patch => "Reproduction number",
    :province_ascertainment => "Case ascertainment",
    :CFR_patch => "Case-fatality ratio",
    :province_death_ascertainment => "Death confirmation"
)
province_pair_figs = [
    plot_pair(
        chn_joint,
        [
            :R_T_patch, :province_ascertainment,
            :CFR_patch, :province_death_ascertainment,
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
#md # <details><summary>Other provinces pair plot (prior overlaid)</summary>
#md # ```

province_pair_figs[4] #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Composition checks
#
# Whether the model reproduces each province's observed share of the national total is on the [in-sample checks](@ref province-compositions) page.

# ## Saving province assets
#
# The summary dashboard shows the province comparison and the reproduction
# number by province, so they are written here rather than on the National
# page.

#md # ```@raw html
#md # <details><summary>Write the province dashboard assets</summary>
#md # ```

dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
CairoMakie.save(joinpath(dashboard_dir, "rt_provinces.png"), province_rt_fig)
open(joinpath(dashboard_dir, "provinces.md"), "w") do io
    print(io, province_headline_md)
end

#md # ```@raw html
#md # </details>
#md # ```

# ---
#
# The full analysis code, data and model definitions are in the
# [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
# repository.
