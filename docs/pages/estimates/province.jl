# # Province estimates
#
# The outbreak by province.
# The model runs one renewal equation per province and fits the national
# streams against the summed provinces, so the national estimates on the
# [analysis](@ref "Results") page are the sum of the provinces here.
#
# The per-province case-fatality ratio sits with the national one in the
# [confirmed case-fatality ratio](@ref "Confirmed case-fatality ratio").
# The per-province forecast is in the
# [one-week-ahead forecast results](@ref "One-week-ahead forecast results")
# and its scoring in the [forecast by province](@ref "Forecast by province").
# Each needs a quantity the page it sits on already computes.

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

# ## Size and infections

# The national outbreak size on the [analysis](@ref "Joint model estimates") page is the sum of the four patches' renewal equations.
# The reproduction number and the relative case ascertainment are identified only as a product, and the per-province deaths break the tie.

#md # ```@raw html
#md # <details><summary>Cross-province overview table</summary>
#md # ```

province_overview_table = patch_overview_table(chn_joint, N_PATCHES);

#md # ```@raw html
#md # </details>
#md # ```

province_overview_table #hide

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

province_imports_fig #hide

# ## Reproduction number by province

# The national reproduction number is on the [analysis](@ref "Reproduction number over time") page.
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

# The figure below gives each province's estimates, including its log-Rt deviation from the trend, the deviation's walk scale and the contrast against Ituri.

#md # ```@raw html
#md # <details><summary>Per-province summary figure</summary>
#md # ```

province_detail_fig = plot_patch_summary(chn_joint, N_PATCHES);

#md # ```@raw html
#md # </details>
#md # ```

province_detail_fig #hide

#md # ```@raw html
#md # <details><summary>Per-province summary table</summary>
#md # ```

province_detail_table = patch_summary_table(chn_joint, N_PATCHES);

province_detail_table #hide

#md # ```@raw html
#md # </details>
#md # ```

# The spread of those deviations is the spatial diagnostic.
# The prior admits real divergence, with a 31% prior probability that the Ituri to Nord-Kivu ratio moves by more than 25% over the window, so a shrunken posterior is a finding rather than an artefact of the prior.
# With four patches and the pooled one carrying almost no signal the cross-province correlation is not identified, and it tracks its prior.

#md # ```@raw html
#md # <details><summary>Spatial hyperparameter summary table</summary>
#md # ```

spatial_hyper_table = summary_table(
    chn_joint,
    [
        :region_sd, :region_halflife, :region_corr_primary_secondary,
        :province_ascertainment_sd, :province_testing_coefficient,
    ];
    digits = 3,
    labels = Dict(
        :region_sd => "Rt deviation spread",
        :region_halflife => "Rt deviation half-life (days)",
        :region_corr_primary_secondary => "Ituri-N.Kivu Rt correlation",
        :province_ascertainment_sd => "Ascertainment spread",
        :province_testing_coefficient => "Testing effect on ascertainment"
    )
);

#md # ```@raw html
#md # </details>
#md # ```

spatial_hyper_table #hide

# ## Composition checks
#
# Whether the model reproduces each province's observed share of the national total is on the [in-sample checks](@ref "Province compositions") page.

# ## Saving province assets
#
# The summary dashboard shows the per-province table and the reproduction
# number by province, so they are written here rather than on the analysis
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
    print(io, markdown_table(province_overview_table))
end

#md # ```@raw html
#md # </details>
#md # ```

# ---
#
# The full analysis code, data and model definitions are in the
# [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
# repository.
