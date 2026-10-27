# # Province in-sample checks
#
# This page checks whether the joint model reproduces the per-province data it was fitted to, follows each province's published estimates across releases, and refits it as a single population to check the spatial structure.
# The national checks are on the [national in-sample checks](@ref "In-sample checks") page.

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations and the fit registry. See
## `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))
#-
## The fits this page reads, loaded from the cache here.
chn_joint = load_fit("joint")
## `sens_no_patches` is the headline with `n_patches = 1`, the check on the
## spatial structure.
chn_no_patches = load_fit("sens_no_patches")
posterior_C_joint = vec(Array(chn_joint[:C_T]))
posterior_C_no_patches = vec(Array(chn_no_patches[:C_T]));

#md # ```@raw html
#md # </details>
#md # ```

# ## Summary
#
# Whether each province's cases and deaths are reproduced, from the checks below.
# Bias and coverage are defined as on the [national in-sample checks](@ref "In-sample checks") page, with coverage of the 90% predictive interval.

#md # ```@eval
#md # using Markdown, BVDOutbreakSize
#md # dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
#md # Markdown.parse(read(joinpath(dir, "evaluation_insample_province.md"), String))
#md # ```

# ## Province prior predictive check
#
# Whether the province model's prior, before any data are fitted, brackets each province's observed share of the confirmed cases and deaths.

#md # ```@raw html
#md # <details><summary>Draw from the province model's prior</summary>
#md # ```

prior_patch_chn = patch_prior_draws(obs);

prior_province_table = patch_overview_table(prior_patch_chn, N_PATCHES);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Show prior province summary table</summary>
#md # ```

MarkdownTable(prior_province_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Province composition prior predictive checks</summary>
#md # ```

prior_case_ppc_fig = plot_province_composition_ppc(
    prior_patch_chn;
    share_key = :province_shares,
    obs_increments = province_cases.increments,
    days = province_cases.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed case share by province, prior"
);

prior_death_ppc_fig = plot_province_composition_ppc(
    prior_patch_chn;
    share_key = :province_death_shares,
    obs_increments = province_deaths.increments,
    days = province_deaths.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed death share by province, prior"
);

#md # ```@raw html
#md # </details>
#md # ```

prior_case_ppc_fig #hide

#-

prior_death_ppc_fig #hide

#md # ```@raw html
#md # <details><summary>Province prior and posterior pair plot</summary>
#md # ```

## Per-province entries of a vector deterministic, one draw vector each,
## keyed `<name>_<patch>` so they can sit side by side in one table.
function _patch_columns(chn, key::Symbol, name::AbstractString; f = identity)
    draws = [collect(v) for v in vec(collect(chn[key]))]
    return [
        Symbol(name, "_", p) => Float64[f(d[p]) for d in draws]
            for p in 1:N_PATCHES
    ]
end
_province_pair_draws(chn) = NamedTuple(
    vcat(
        _patch_columns(chn, :R_T_patch, "R_T"),
        _patch_columns(chn, :C_T_patch, "log10_C_T"; f = x -> log10(x + 1))
    )
)
_province_pair_labels = Dict(
    Symbol(name, "_", p) => string(label, ", ", PROVINCE_LABELS[p])
        for p in 1:N_PATCHES
        for (name, label) in (("R_T", "R_T"), ("log10_C_T", "log10 C_T"))
)
province_pair_fig = plot_pair(
    _province_pair_draws(chn_joint);
    prior = _province_pair_draws(prior_patch_chn),
    labels = _province_pair_labels
);

#md # ```@raw html
#md # </details>
#md # ```

province_pair_fig #hide

# ## [Province compositions](@id province-compositions)
#
# Whether each province's modelled share of the national confirmed cases and deaths reproduces the observed share at every spatial vintage (see the [province compositions](@ref methods-province-compositions) Methods section).

#md # ```@raw html
#md # <details><summary>Province composition posterior predictive checks</summary>
#md # ```

province_case_ppc_fig = plot_province_composition_ppc(
    chn_joint;
    share_key = :province_shares,
    obs_increments = province_cases.increments,
    days = province_cases.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed case share by province"
);

province_death_ppc_fig = plot_province_composition_ppc(
    chn_joint;
    share_key = :province_death_shares,
    obs_increments = province_deaths.increments,
    days = province_deaths.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed death share by province"
);

#md # ```@raw html
#md # </details>
#md # ```

province_case_ppc_fig #hide

#-

province_death_ppc_fig #hide

# The four province terms for the isolation and laboratory data are checked the same way.
# A gap in a panel is a day on which that province printed nothing.

#md # ```@raw html
#md # <details><summary>Province laboratory, occupancy, bed and admission split posterior predictive checks</summary>
#md # ```

province_lab_bin_days = [
    maximum(province_lab.days[province_lab.bins .== b])
        for b in 1:maximum(province_lab.bins)
];

province_lab_ppc_fig = plot_province_composition_ppc(
    chn_joint;
    share_key = :province_lab_shares,
    obs_increments = province_lab.increments,
    days = province_lab_bin_days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Analysed specimen share by province, by week"
);

province_occupancy_ppc_fig = plot_province_split_ppc(
    chn_joint;
    share_key = :province_occupancy_share, rows = province_isolation,
    seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Isolation occupancy share by province"
);

province_beds_ppc_fig = plot_province_split_ppc(
    chn_joint;
    share_key = :province_capacity_share, rows = province_capacity,
    seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Isolation bed share by province"
);

province_admissions_ppc_fig = plot_province_split_ppc(
    chn_joint;
    share_key = :province_admissions_share, rows = province_admissions,
    seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Admission share by province"
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

province_lab_ppc_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The split of each calendar week's analysed specimens among the provinces, which identifies the background split.
#md # </figcaption>
#md # </figure>
#md # ```

#-

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_occupancy_ppc_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The split of the patients in isolation among the provinces printed that day, on the weekly days the fit scores, against the per-patch bed demand.
#md # </figcaption>
#md # </figure>
#md # ```

#-

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_beds_ppc_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The split of the beds among the provinces printed that day, on the days a count changed, against each patch's modelled share of the national capacity.
#md # </figcaption>
#md # </figure>
#md # ```

#-

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_admissions_ppc_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The split of each day's 24h admissions among the provinces printed that day, against each patch's modelled admissions.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Province stream calibration
#
# Whether each province's share predictions are calibrated at the observed national total, scored as in the national [stream calibration](@ref "Stream calibration").
# The analysed-specimen rows score the weekly laboratory split the same way, at each week's observed national analysed total.

#md # ```@raw html
#md # <details><summary>Build the province calibration panels</summary>
#md # ```

province_case_panels = province_composition_panels(
    chn_joint;
    share_key = :province_shares,
    obs_increments = province_cases.increments,
    stream = "Confirmed cases", n_patches = N_PATCHES
);
province_death_panels = province_composition_panels(
    chn_joint;
    share_key = :province_death_shares,
    obs_increments = province_deaths.increments,
    stream = "Confirmed deaths", n_patches = N_PATCHES
);
province_lab_panels = province_composition_panels(
    chn_joint;
    share_key = :province_lab_shares,
    obs_increments = province_lab.increments,
    stream = "Analysed specimens", n_patches = N_PATCHES
);
province_panels = vcat(
    province_case_panels, province_death_panels, province_lab_panels
);
province_calibration_table = stream_calibration(province_panels);
province_calibration_fig = plot_stream_calibration(
    province_calibration_table
);

#md # ```@raw html
#md # </details>
#md # ```

province_calibration_fig #hide

#md # ```@raw html
#md # <details><summary>Province calibration table</summary>
#md # ```

MarkdownTable(province_calibration_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Province counts
#
# Whether each province's cases and deaths are reproduced as counts once the national total is predicted as well.

#md # ```@raw html
#md # <details><summary>Build the province count panels</summary>
#md # ```

## Each province's confirmed cases and deaths split from the national
## posterior predictive by the fitted shares, shared with the dashboard's
## province series (see `joint_count_replicates` in `docs/pages/_setup.jl`).
province_count_replicates = joint_count_replicates();
province_case_count_panels = province_count_replicates.cases.panels;
province_death_count_panels = province_count_replicates.deaths.panels;
province_count_calibration_table = stream_calibration(
    vcat(province_case_count_panels, province_death_count_panels)
);
province_count_calibration_fig = plot_stream_calibration(
    province_count_calibration_table
);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Province count calibration plot</summary>
#md # ```

province_count_calibration_fig #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Province count calibration table</summary>
#md # ```

MarkdownTable(province_count_calibration_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Province correlations and totals
#
# Which province quantities trade off against each other, and whether each province's case and death totals agree with the observed ones.

#md # ```@raw html
#md # <details><summary>Province posterior correlation heatmap</summary>
#md # ```

## LaTeX label for one province's entry of a quantity. Spaces in a province
## name need escaping inside `\mathrm`.
_province_tex(sym, p) = string(
    sym, raw"\ (\mathrm{", replace(PROVINCE_LABELS[p], " " => raw"\ "), "})"
)
province_correlation_draws = NamedTuple(
    vcat(
        [:C_T => vec(Array(chn_joint[:C_T]))],
        _patch_columns(chn_joint, :R_T_patch, "R_T"),
        _patch_columns(chn_joint, :province_ascertainment, "asc"),
        _patch_columns(chn_joint, :CFR_patch, "CFR")
    )
);
province_correlation_labels = merge(
    Dict(:C_T => raw"C_T"),
    Dict(
        Symbol(name, "_", p) => _province_tex(tex, p)
            for p in 1:N_PATCHES
            for (name, tex) in (
                ("R_T", raw"R_T"), ("asc", raw"p_\mathrm{asc}"),
                ("CFR", raw"\mathrm{CFR}"),
            )
    )
);
province_correlation_fig = plot_correlation_heatmap(
    province_correlation_draws; labels = province_correlation_labels
);

#md # ```@raw html
#md # </details>
#md # ```

province_correlation_fig #hide

#md # ```@raw html
#md # <details><summary>Province totals against observed</summary>
#md # ```

## Keyed by panel title, which names the stream and the province. A total
## that is the same in every draw (a patch with no deaths at any vintage) has
## no spread to plot and gives the corner plot a zero-width axis, so it is
## left out.
_varying_panels = filter(
    p -> length(unique(sum(r) for r in p.replicates)) > 1, province_panels
);
province_totals = NamedTuple(
    Symbol(p.title) => [sum(Float64.(r)) for r in p.replicates]
        for p in _varying_panels
);
province_observed = NamedTuple(
    Symbol(p.title) => Float64(sum(p.observed)) for p in _varying_panels
);
province_pairs_fig = plot_stream_pairs(province_totals, province_observed);

#md # ```@raw html
#md # </details>
#md # ```

province_pairs_fig #hide

# ## Estimate evolution across releases
#
# How each province's estimates have moved as releases accrued, read from the estimates each release published at its own cut-off.
# Releases began publishing per-province estimates recently, so each panel fills in as releases accumulate.

#md # ```@raw html
#md # <details><summary>Load the per-province estimates of each release</summary>
#md # ```

## Written by `scripts/score_releases.jl` from each release's
## `province_estimates.csv`, one row per release, province and quantity. A
## missing file reads as an empty table, which the figures report as nothing
## collected.
province_release_df = _release_data(
    joinpath("province", "released_estimates.csv"),
    (;
        release = String, date = Date, area = String, patch = String,
        pooled = Int, quantity = String, median = Float64, lo90 = Float64,
        hi90 = Float64,
    )
)
## One panel per province, the patches first in their usual order. A
## province pooled into a shared patch is marked, since it carries its
## patch's reproduction number and case-fatality ratio.
_release_provinces = unique(vcat(PROVINCE_NAMES, province_release_df.area))
_release_pooled = Set(
    r.area for r in eachrow(province_release_df)
        if coalesce(r.pooled, 0) == 1
)
_release_titles = Dict(
    a => string(
        get(Dict(zip(PROVINCE_NAMES, PROVINCE_LABELS)), a, a),
        a in _release_pooled ? " (pooled)" : ""
    ) for a in _release_provinces
)
_province_release_groups(q) = area_evolution_groups(
    province_release_df, q;
    areas = _release_provinces, labels = _release_titles
)
_no_release_note = "No releases scored yet. This fills in as releases accumulate."
province_rt_release_fig = plot_evolution_by_group(
    _province_release_groups("R_T");
    ylabel = "Reproduction number",
    title = "Province reproduction number by release",
    released_label = "Released estimate", marks_label = "median, 90% bar",
    refline = 1.0, ymax = 3.0, ncols = 4,
    empty_note = _no_release_note
);
province_cases_release_fig = plot_evolution_by_group(
    _province_release_groups("cases");
    ylabel = "Confirmed cases",
    title = "Province confirmed cases by release",
    released_label = "Reported by the release cut-off", marks_label = "",
    shared_yrange = false, ncols = 4, empty_note = _no_release_note
);
province_cfr_release_fig = plot_evolution_by_group(
    _province_release_groups("CFR");
    ylabel = "Case-fatality ratio",
    title = "Province case-fatality ratio by release",
    released_label = "Released estimate", marks_label = "median, 90% bar",
    ymax = 1.0, ncols = 4,
    empty_note = _no_release_note
);

#md # ```@raw html
#md # </details>
#md # ```

# ### Reproduction number by release

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_rt_release_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Each province's reproduction number at the cut-off of each release, as a median with its 90% interval bar.
#md # The axis is fixed at three, and an interval past it is clamped and marked with an open triangle.
#md # </figcaption>
#md # </figure>
#md # ```

# ### Confirmed cases by release

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_cases_release_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The confirmed cases each province had reported by the cut-off of each release, each panel on its own scale.
#md # </figcaption>
#md # </figure>
#md # ```

# ### Case-fatality ratio by release

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_cfr_release_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Each province's case-fatality ratio at each release, as a median with its 90% interval bar.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Province sensitivity to assumptions
#
# ### Spatial structure sensitivity
#
# We refit the joint model with its eight patches reduced to one well-mixed population.
# Splitting the country into provinces adds no national data, so the national outbreak size should not move far either way.
# The two still differ, because the country grows at the force-weighted mean of the provinces' reproduction numbers, which sits above the central trend they pool toward.
# That gap is first order in the deviation scale, so a large gap between the two posteriors would point at the deviation priors rather than at the data.

spatial_sensitivity_table = streams_table(
    "Meta-population (headline)" => posterior_C_joint,
    "Single population (n_patches = 1)" => posterior_C_no_patches
);
spatial_sensitivity_table

#md # ```@raw html
#md # <details><summary>Spatial-structure density overlay</summary>
#md # ```

spatial_sensitivity_fig = plot_density_overlay(
    "Meta-population (headline)" => posterior_C_joint,
    "Single population" => posterior_C_no_patches;
    xlabel = "Cumulative infections"
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

spatial_sensitivity_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Posterior densities of the cumulative infections under the headline province model and the single-population fit.
#md # </figcaption>
#md # </figure>
#md # ```

# The spatial structure could also move other national quantities that both fits estimate, so their posteriors should sit on top of each other.

#md # ```@raw html
#md # <details><summary>National quantities under both structures</summary>
#md # ```

spatial_rt_fig = plot_rt_streams(
    [
        (;
            label = "Single population (n_patches = 1)",
            chn = chn_no_patches, rt_start = _rt_start_plot,
            rt_walk_start = clamp(
                _BREAKPOINT - RT_WALK_LEAD, _rt_start_plot,
                obs.n
            ), colour = :steelblue,
        ),
    ];
    joint = (;
        chn = chn_joint, rt_start = _rt_start_plot,
        rt_walk_start = clamp(
            _BREAKPOINT - RT_WALK_LEAD, _rt_start_plot,
            obs.n
        ),
    ),
    n = obs.n, breakpoint = _BREAKPOINT,
    as_of_date = string(obs.cutoff), seeding = obs.seeding,
    display_start = _rt_start_plot, ncols = 1,
    title = "National reproduction number under both structures",
    reference_label = "the meta-population headline",
    panel_label = "the single-population fit"
);

spatial_cfr_fig = plot_density_overlay(
    "Meta-population (headline)" => vec(Array(chn_joint[:CFR])),
    "Single population" => vec(Array(chn_no_patches[:CFR]));
    xlabel = "Case-fatality ratio"
);

spatial_rt_density_fig = plot_density_overlay(
    "Meta-population (headline)" => vec(Array(chn_joint[:R_T])),
    "Single population" => vec(Array(chn_no_patches[:R_T]));
    xlabel = "Reproduction number at the cut-off"
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

spatial_rt_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The national reproduction number from the single-population fit (coloured) against the headline province model (grey), which is hidden where the two agree.
#md # </figcaption>
#md # </figure>
#md # ```

#-

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

spatial_cfr_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Posterior densities of the case-fatality ratio under both structures.
#md # </figcaption>
#md # </figure>
#md # ```

#-

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

spatial_rt_density_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Posterior densities of the reproduction number at the cut-off under both structures.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>National quantities under both structures, as a table</summary>
#md # ```

## One row per national quantity, one column per structure, each cell a
## median with a 90% credible interval. Built here rather than by stacking
## two `summary_table` calls, so the two structures sit side by side and the
## reader compares along a row.
spatial_quantities_table = let
    ## A count rounded to zero decimals still prints a trailing ".0", so
    ## whole-number quantities go through `Int`.
    fmt(x, d) = d <= 0 ? string(round(Int, x)) : string(round(x; digits = d))
    cell(v, d) = string(
        fmt(quantile(v, 0.5), d), " (",
        fmt(quantile(v, 0.05), d), "–", fmt(quantile(v, 0.95), d), ")"
    )
    rows = [
        ("Cumulative infections", :C_T, 0),
        ("Reproduction number at the cut-off", :R_T, 2),
        ("Case-fatality ratio", :CFR, 2),
        ("Outbreak age (days)", :T, 0),
        ("Latest growth rate (per day)", :r, 3),
    ]
    DataFrame(
        "Quantity" => [r[1] for r in rows],
        "Meta-population (headline)" => [
            cell(
                vec(Array(chn_joint[r[2]])), r[3]
            )
                for r in rows
        ],
        "Single population" => [
            cell(vec(Array(chn_no_patches[r[2]])), r[3])
                for r in rows
        ]
    )
end;

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

spatial_quantities_table #hide

#md # ```@raw html
#md # <figcaption>
#md # The same three quantities as credible intervals under both structures, with the outbreak start date each fit implies.
#md # This is where the agreement is read rather than from the figures.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Province parameter recovery
#
# Whether the model recovers each province's values when fitted to data it simulated itself, from the same runs as the national [parameter recovery](@ref "Parameter recovery").

#md # ```@raw html
#md # <details><summary>Recovery figure</summary>
#md # ```

province_recovery_results = recovery_results()
_recovery_short = [
    "C_T_patch" => "C_T", "R_T_patch" => "R_T", "CFR_patch" => "CFR",
    "province_ascertainment" => "ascertainment",
    "region_drift_sd" => "drift SD",
]
province_recovery_quantities = [
    "$(k)[$(p)]" for (k, _) in _recovery_short for p in PROVINCE_LABELS
]
province_recovery_fig = isempty(province_recovery_results.params) ? nothing :
    plot_recovery(
        province_recovery_results.params, province_recovery_results.draws,
        province_recovery_results.prior;
        quantities = province_recovery_quantities,
        labels = Dict(
            "$(k)[$(p)]" => "$(v), $(p)" for (k, v) in _recovery_short
            for p in PROVINCE_LABELS
        ),
        panel_labels = Dict(
            "$(k)[$(p)]" => p for (k, _) in _recovery_short for p in PROVINCE_LABELS
        ),
        row_labels = last.(_recovery_short),
        log_x = ["C_T_patch[$(p)]" for p in PROVINCE_LABELS]
    );

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

isempty(province_recovery_results.params) ? Markdown.parse("No parameter-recovery run is available for this build.") : province_recovery_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Each seed's posterior median with its 50% and 90% intervals divided by its true value (top).
#md # Below, each quantity by province is on its own scale, with the prior in grey, each seed's posterior in its colour and its true value dashed.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>Summary table across seeds</summary>
#md # ```

province_recovery_summary = isempty(province_recovery_results.params) ?
    DataFrame() :
    recovery_summary_table(province_recovery_results.params; province = true);

province_recovery_summary_display = isempty(province_recovery_summary) ?
    Markdown.parse("No parameter-recovery run is available for this build.") :
    MarkdownTable(province_recovery_summary);

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_recovery_summary_display #hide

#md # ```@raw html
#md # <figcaption>
#md # The error of each seed's posterior median relative to the truth, and the z-score of the truth, summarised across seeds.
#md # </figcaption>
#md # </figure>
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Each seed's recovered province values</summary>
#md # ```

province_recovery = let r = province_recovery_results.params
    isempty(r) ? DataFrame() :
        r[
            occursin.("[", r.quantity), [
                :seed, :quantity, :truth, :median, :lower_90, :upper_90,
                :covered_90,
            ],
        ]
end;

province_recovery_display = isempty(province_recovery) ?
    Markdown.parse("No parameter-recovery run is available for this build.") : MarkdownTable(province_recovery);
province_recovery_display #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Saving province in-sample outputs

#md # ```@raw html
#md # <details><summary>Write the summary bullets</summary>
#md # ```

## The bullets under the summary heading at the top of the page. They read
## tables built further down, so they are written here and read back when
## the site is assembled.
evaluation_insample_province_summary = let
    fmt(x) = ismissing(x) || !isfinite(x) ? "n/a" :
        string(round(x; digits = 2))
    fitted = filter(r -> isfinite(r["Bias"]), province_calibration_table)
    worst = first(
        sort(fitted, "Bias"; by = abs, rev = true), min(3, size(fitted, 1))
    )
    function coverage_bullet(heading, tbl, scale)
        low = first(sort(tbl, "90% coverage"), min(1, size(tbl, 1)))
        return string(
            "- **", heading, ":** ", count(>=(0.8), tbl[!, "90% coverage"]),
            " of ", size(tbl, 1), " province streams have 90% coverage of ",
            "at least 0.8", scale,
            size(low, 1) == 0 ? "." :
                string(
                    "; the lowest is ", low[1, "Stream"], " at ",
                    fmt(low[1, "90% coverage"]), "."
                )
        )
    end
    overall = [
        string(
            "- **Least well reproduced:** ",
            join(
                [
                    string(
                        r["Stream"], " (bias ", fmt(r["Bias"]),
                        ", 90% coverage ", fmt(r["90% coverage"]), ")"
                    )
                        for r in eachrow(worst)
                ], "; "
            ), "."
        ),
        coverage_bullet("Coverage", fitted, ""),
        coverage_bullet(
            "Count-scale coverage",
            filter(
                r -> isfinite(r["Bias"]), province_count_calibration_table
            ),
            ""
        ),
    ]
    rows(tbl) = Dict(r["Stream"] => r for r in eachrow(tbl))
    cal = rows(province_calibration_table)
    count_cal = rows(province_count_calibration_table)
    function calibration(tbl, kind, p, label)
        r = get(tbl, string(kind, ", ", PROVINCE_LABELS[p]), nothing)
        r === nothing && return string("- ", label, ": not scored.")
        return string(
            "- ", label, ": bias ", fmt(r["Bias"]), " and 90% coverage ",
            fmt(r["90% coverage"]), " over ", r["Vintages"], " vintages."
        )
    end
    detail = [
        join(
            [
                string("**", PROVINCE_LABELS[p], "**"), "",
                calibration(cal, "Confirmed cases", p, "Confirmed cases"),
                calibration(cal, "Confirmed deaths", p, "Confirmed deaths"),
                calibration(
                    count_cal, "Confirmed cases", p,
                    "Confirmed cases as counts"
                ),
                calibration(
                    count_cal, "Confirmed deaths", p,
                    "Confirmed deaths as counts"
                ),
            ], "\n"
        )
            for p in 1:N_PATCHES
    ]
    join(vcat([join(overall, "\n")], detail), "\n\n")
end
dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
write(
    joinpath(dashboard_dir, "evaluation_insample_province.md"),
    evaluation_insample_province_summary
);

#md # ```@raw html
#md # </details>
#md # ```
