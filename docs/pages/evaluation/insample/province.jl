# # Province in-sample checks
#
# Whether the fitted joint model reproduces the per-province data it was fitted to.
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
chn_joint = load_fit("joint");

#md # ```@raw html
#md # </details>
#md # ```

# ## Summary
#
# Whether each province's cases and deaths are reproduced, from the checks further down this page.
# Bias runs from −1 to 1 and is negative when the model under-predicts, and coverage is the fraction of vintages inside the 90% predictive interval.

#md # ```@eval
#md # using Markdown, BVDOutbreakSize
#md # dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
#md # Markdown.parse(read(joinpath(dir, "evaluation_insample_province.md"), String))
#md # ```

# ## Province prior predictive check
#
# Whether the four-patch prior, before any data are fitted, brackets each province's observed share of the confirmed cases and deaths.

#md # ```@raw html
#md # <details><summary>Draw from the four-patch prior</summary>
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

province_death_ppc_fig #hide

# The three province terms added with the isolation and laboratory data are checked the same way.
# The laboratory panel is the split of each calendar week's analysed specimens, which the background split identifies.
# The occupancy panel is the split of the patients in isolation among the provinces printed that day, on the weekly days the fit scores, over the per-patch bed demand.
# The bed panel is the split of the beds among the provinces printed that day, on the days a count changed, over each patch's static share of the national capacity.
# A gap in a panel is a day on which that province printed nothing.

#md # ```@raw html
#md # <details><summary>Province laboratory, occupancy and bed split posterior predictive checks</summary>
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

#md # ```@raw html
#md # </details>
#md # ```

province_lab_ppc_fig #hide

province_occupancy_ppc_fig #hide

province_beds_ppc_fig #hide

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

## The joint posterior predictive, shared with the national page (see
## `joint_posterior_predictive` in `docs/pages/_setup.jl`). Its draws pair
## one to one with `chn_joint`, whose shares split them.
pp_joint = joint_posterior_predictive();
_pp_draws(vn) = collect(pp_joint[BVDOutbreakSize.FlexiChains.Prefixed(vn)]);

## The national confirmed cases are replicated per laboratory window, oldest
## first, as on the national page. The first confirmed vintage is an
## unreplicated baseline.
_conf_windows = BVDOutbreakSize.confirmed_positivity_windows(
    obs.confirmed_history, obs.lab_history, obs.lab_daily_history
);
_conf_obs = collect(
    first(
        pp_joint[k] for k in keys(pp_joint)
            if occursin(
                "confirmed_state.confirmed_positives.positives", string(k)
            )
    )
);
_conf_replicates = [
    vcat(collect(e), collect(p), collect(l))
        for (e, p, l) in zip(
            vec(_pp_draws(@varname(early_increments.increments))),
            vec(_conf_obs),
            vec(_pp_draws(@varname(late_increments.increments)))
        )
];
province_case_count_panels = province_count_panels(
    chn_joint;
    share_key = :province_shares,
    obs_increments = province_cases.increments,
    province_days = province_cases.days,
    national_days = vcat(
        _conf_windows.early_days, _conf_windows.obs_days,
        _conf_windows.late_days
    ),
    national_replicates = _conf_replicates,
    baseline = Int(first(obs.confirmed_history.counts)),
    stream = "Confirmed cases", n_patches = N_PATCHES
);
## The national confirmed deaths are replicated from zero at every vintage.
province_death_count_panels = province_count_panels(
    chn_joint;
    share_key = :province_death_shares,
    obs_increments = province_deaths.increments,
    province_days = province_deaths.days,
    national_days = obs.confirmed_deaths_history.days,
    national_replicates = _pp_draws(@varname(cdeath_increments.increments)),
    stream = "Confirmed deaths", n_patches = N_PATCHES
);
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
# ## Province parameter recovery
#
# Whether the model recovers each province's values when fitted to data it simulated itself, from the same runs as the national [parameter recovery](@ref "Parameter recovery").
# The top panel shows each seed's posterior median with its 50% and 90% intervals divided by that seed's true value, and below each quantity by province is on its own scale with the prior in grey, each seed's posterior in its colour and its true value as a dashed line.

#md # ```@raw html
#md # <details><summary>Recovery figure</summary>
#md # ```

province_recovery_results = recovery_results()
_recovery_short = [
    "C_T_patch" => "C_T", "R_T_patch" => "R_T", "CFR_patch" => "CFR",
    "province_ascertainment" => "ascertainment",
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
#md # ```

isempty(province_recovery_results.params) ? Markdown.parse("No parameter-recovery run is available for this build.") : province_recovery_fig #hide

# The error of each seed's posterior median relative to the truth, and the z-score of the truth, summarised across seeds.

#md # ```@raw html
#md # <details><summary>Summary across seeds</summary>
#md # ```

province_recovery_summary = isempty(province_recovery_results.params) ?
    DataFrame() :
    recovery_summary_table(province_recovery_results.params; province = true);

province_recovery_summary_display = isempty(province_recovery_summary) ?
    Markdown.parse("No parameter-recovery run is available for this build.") :
    MarkdownTable(province_recovery_summary);

#md # ```@raw html
#md # </details>
#md # ```

province_recovery_summary_display #hide

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
