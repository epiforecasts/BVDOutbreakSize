# # Province forecasts
#
# This page gives the one-week-ahead forecast for each province from the joint model, as defined in the [province forecast](@ref "Province forecast") Methods section.
# The national forecast is on the [forecasts](@ref "Forecasts") page.
# How these forecasts have scored against what each province went on to report is in the [forecast by province across releases](@ref "Forecast by province across releases") evaluation.
# How well the model reproduces each province's share is on the [in-sample checks](@ref province-compositions) page.

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
# The table gives each province's forecast for the week to $T + 7$ as a 90% predictive interval, and the bullets compare the provinces draw by draw.

#md # ```@raw html
#md # <details><summary>Project each province a week ahead</summary>
#md # ```

## Read from the same draws as the national forecast page, so the provinces
## add up to the national forecast shown there.
province_draws = fit_forecast("joint");
province_projection = forecast_provinces(
    province_draws; horizon = 7, n_patches = N_PATCHES
);
province_forecast = province_forecast_table(
    province_draws, province_projection;
    n_patches = N_PATCHES
);
province_forecast_fig = plot_province_forecast(
    province_draws, province_projection;
    n_patches = N_PATCHES
);
province_week_end = obs.cutoff + Day(7);
province_forecast_md = province_forecast_headline(
    province_projection; n_patches = N_PATCHES
);

#md # ```@raw html
#md # </details>
#md # ```

Markdown.parse(province_forecast_md) #hide

# ## One-week-ahead forecast by province

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

province_forecast_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Each panel is one forecast target for the week after the cut-off, with the provinces side by side.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>Build the forecast maps</summary>
#md # ```

province_forecast_draws(col) = [
    float.(province_projection[province_projection.patch .== p, col])
        for p in 1:N_PATCHES
]
province_forecast_map = plot_province_map(
    filter(
        !isnothing,
        [
            (;
                province_map_summary(province_forecast_draws(:rt_forecast))...,
                title = "Reproduction number at T+7", diverging_at = 1.0,
                colorbar_label = "R (median)",
            ),
            (;
                values = province_map_summary(
                    province_forecast_draws(:confirmed_new)
                ).values,
                title = "New confirmed cases by T+7",
                scale = CairoMakie.Makie.pseudolog10,
                colorbar_label = "Confirmed cases (median)",
            ),
            :isolation_level in propertynames(province_projection) ? (;
                    values = province_map_summary(
                        province_forecast_draws(:isolation_level)
                    ).values,
                    title = "Patients in isolation at T+7",
                    scale = CairoMakie.Makie.pseudolog10,
                    colorbar_label = "Patients (median)",
                ) : nothing,
        ]
    )
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

province_forecast_map #hide

#md # ```@raw html
#md # <figcaption>
#md # Each province is shaded by its forecast.
#md # The reproduction number is centred on one, and a province whose 90% predictive interval spans one is washed out.
#md # </figcaption>
#md # </figure>
#md # ```

MarkdownTable(province_forecast) #hide

# ## Forecast for each province
#
# Each province's figure has one histogram per target with its 90% predictive interval shaded.
# The dashed rule on the case and death panels is the count the province reported over its most recent week.

#md # ```@raw html
#md # <details><summary>Build the per-province figures</summary>
#md # ```

## The spatial tables' most recent week in each province, clamped as the
## compositions are, so a downward revision reads as no new cases.
recent_cases = province_recent_counts(
    obs.province_confirmed_history, PROVINCE_NAMES, N_PATCHES
)
recent_deaths = province_recent_counts(
    obs.province_death_history, PROVINCE_NAMES, N_PATCHES
)
function forecast_province_observed(p)
    o = (;)
    recent_cases === nothing ||
        (o = merge(o, (; confirmed_new = recent_cases.counts[p])))
    recent_deaths === nothing ||
        (o = merge(o, (; confirmed_deaths_new = recent_deaths.counts[p])))
    return o
end
forecast_province_fig(p) = plot_province_forecast_detail(
    province_draws, province_projection;
    province = p, n_patches = N_PATCHES,
    observed = forecast_province_observed(p)
)
## The province blocks below are written out one per patch.
@assert N_PATCHES == 8 && PROVINCE_LABELS[1:8] == [
    "Ituri", "Nord-Kivu", "Haut-Uele",
    "Sud-Kivu", "Tshopo", "Bas-Uele", "Sud-Ubangi", "Nord-Ubangi",
]
## Cases and deaths come from separate tables, so each stream's observed
## week is stated on its own.
function forecast_province_recent_week(r, stream)
    r === nothing &&
        return "The spatial tables carry no recent week of $(stream)."
    return string(
        "The observed week of $(stream) runs from ",
        grid_date(r.start_day), " to ", grid_date(r.last_day), "."
    )
end
recent_note = join(
    [
        forecast_province_recent_week(recent_cases, "confirmed cases"),
        forecast_province_recent_week(recent_deaths, "confirmed deaths"),
        "Each ends at its last spatial vintage, which can be earlier than " *
            "the cut-off.",
    ], " "
);

#md # ```@raw html
#md # </details>
#md # ```

Markdown.parse(recent_note) #hide

# ### Ituri

forecast_province_fig(1) #hide

# ### Nord-Kivu

forecast_province_fig(2) #hide

# ### Haut-Uele

forecast_province_fig(3) #hide

# ### Sud-Kivu

forecast_province_fig(4) #hide

# ### Tshopo

forecast_province_fig(5) #hide

# ### Bas-Uele

forecast_province_fig(6) #hide

# ### Sud-Ubangi

forecast_province_fig(7) #hide

# ### Nord-Ubangi

forecast_province_fig(8) #hide

# ## Saving province map assets
#
# The [dashboard](@ref "Dashboard") map's province layer and its detail column read the per-province and national estimates and time series written here.

#md # ```@raw html
#md # <details><summary>Write the province map estimates and time series</summary>
#md # ```

## The national row reads the joint's national reproduction number and
## case-fatality ratio and the provinces' summed forecasts, patients in
## isolation and beds.
province_map_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(province_map_dir)
has_death_forecast = :confirmed_deaths_new in propertynames(province_projection)
has_care_forecast = :isolation_level in propertynames(province_projection) &&
    :bed_capacity in propertynames(province_projection)
has_care = BVDOutbreakSize._has_key(chn_joint, :province_bed_capacity)
## Patients in isolation over beds, draw by draw.
care_ratio(iso, beds) = [i ./ b for (i, b) in zip(iso, beds)]
province_care_forecast = has_care_forecast ? (;
        isolation_forecast = province_forecast_draws(:isolation_level),
        beds_forecast = province_forecast_draws(:bed_capacity),
        bed_use_forecast = care_ratio(
            province_forecast_draws(:isolation_level),
            province_forecast_draws(:bed_capacity)
        ),
    ) : (;)
CSV.write(
    joinpath(province_map_dir, "province_estimates.csv"),
    province_map_estimates(
        chn_joint, province_forecast_draws(:confirmed_new);
        n_patches = N_PATCHES,
        deaths_forecast = has_death_forecast ?
            province_forecast_draws(:confirmed_deaths_new) : nothing,
        confirmed_history = obs.province_confirmed_history,
        death_history = obs.province_death_history, cutoff = obs.cutoff,
        isolation_history = obs.province_isolation_history,
        bed_history = obs.province_bed_capacity_history, n = obs.n,
        province_care_forecast...
    )
);
national_forecast_draws(col) = let d = province_projection
    [float(sum(d[d.draw .== i, col])) for i in sort(unique(d.draw))]
end
## The provinces' cut-off draws summed into one national draw vector.
national_care_draws(key) = let v = BVDOutbreakSize._per_patch(
        chn_joint, key, N_PATCHES
    )
    [sum(x[i] for x in v) for i in eachindex(first(v))]
end
national_care = has_care ? (;
        isolation = [national_care_draws(:province_expected_isolation)],
        beds = [national_care_draws(:province_bed_capacity)],
        bed_use = care_ratio(
            [national_care_draws(:province_expected_isolation)],
            [national_care_draws(:province_bed_capacity)]
        ),
    ) : (;)
national_care_forecast = has_care_forecast ? (;
        isolation_forecast = [national_forecast_draws(:isolation_level)],
        beds_forecast = [national_forecast_draws(:bed_capacity)],
        bed_use_forecast = care_ratio(
            [national_forecast_draws(:isolation_level)],
            [national_forecast_draws(:bed_capacity)]
        ),
    ) : (;)
CSV.write(
    joinpath(province_map_dir, "national_estimates.csv"),
    province_map_estimates(
        [vec(Array(chn_joint[:R_T]))],
        [national_forecast_draws(:confirmed_new)];
        deaths_forecast = has_death_forecast ?
            [national_forecast_draws(:confirmed_deaths_new)] : nothing,
        cfr = [vec(Array(chn_joint[:CFR]))],
        confirmed_history = Dict("national" => obs.confirmed_history),
        death_history = Dict("national" => obs.confirmed_deaths_history),
        cutoff = obs.cutoff, patch_names = ["national"],
        patch_labels = ["National"],
        members = Dict("national" => ["national"]),
        isolation_history = Dict("national" => obs.isolation_history),
        bed_history = Dict("national" => obs.bed_capacity_history),
        n = obs.n, national_care..., national_care_forecast...
    )
);
## The daily reproduction number, the weekly confirmed cases and deaths and
## their posterior predictive, of each patch and of the country, keyed by
## the patch label, and the
## reported patients in isolation and beds over the same eight weeks, keyed
## by source province.
province_ts_rt_args = (;
    n = obs.n, breakpoint = _BREAKPOINT, rt_start = _rt_start_plot,
    rt_walk_start = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n),
    ramp = RT_INTERVENTION_RAMP,
)
province_ts_inc = province_increment_matrix(
    obs.province_confirmed_history, PROVINCE_NAMES, N_PATCHES
)
province_ts_deaths = province_increment_matrix(
    obs.province_death_history, PROVINCE_NAMES, N_PATCHES
)
national_increments(c) = isempty(c) ? zeros(Int, 1, 0) :
    reshape([c[1]; max.(diff(c), 0)], 1, :)
## Each patch's and the country's weekly counts from their increments.
function province_weekly(inc, national)
    return vcat(
        weekly_count_table(
            inc.days, inc.increments, PROVINCE_LABELS[1:N_PATCHES];
            cutoff = obs.cutoff, n = obs.n
        ),
        weekly_count_table(
            national.days, national_increments(national.counts),
            ["National"]; cutoff = obs.cutoff, n = obs.n
        )
    )
end
## The joint's posterior predictive of the same weekly counts, binned into
## the same windows (see `joint_count_replicates` in `docs/pages/_setup.jl`).
## The country is left out, with a warning, when its replicates sit on
## other vintage days than its observed series over those weeks, since they
## would then be binned differently.
function province_weekly_fit(stream::Symbol, inc, national_days)
    r = getproperty(joint_count_replicates(), stream)
    as_matrix(reps) = permutedims(reduce(hcat, reps))
    t = weekly_quantile_table(
        inc.days, [as_matrix(p.replicates) for p in r.panels],
        PROVINCE_LABELS[1:N_PATCHES]; cutoff = obs.cutoff, n = obs.n
    )
    recent(d) = filter(>(obs.n - 7 * 8), d)
    if recent(r.national_days) != recent(national_days)
        @warn "National $(stream) fit left out of the province time series: " *
            "its replicates and observed counts are on different vintage days" replicated = recent(r.national_days) observed = recent(national_days)
        return t
    end
    return vcat(
        t, weekly_quantile_table(
            r.national_days, [as_matrix(r.national)], ["National"];
            cutoff = obs.cutoff, n = obs.n
        )
    )
end
## Tag a long table with the series it holds.
_with_series(t, s) = (t[!, :series] .= s; t)
CSV.write(
    joinpath(province_map_dir, "province_timeseries.csv"),
    vcat(
        _with_series(
            rt_quantile_table(
                [
                    reconstruct_patch_rt(
                        chn_joint; n_patches = N_PATCHES,
                        province_ts_rt_args...
                    );
                    [reconstruct_rt(chn_joint; province_ts_rt_args...)]
                ],
                [PROVINCE_LABELS[1:N_PATCHES]; "National"];
                cutoff = obs.cutoff, n = obs.n, from = _rt_start_plot
            ), "rt"
        ),
        _with_series(
            province_weekly(province_ts_inc, obs.confirmed_history), "cases"
        ),
        _with_series(
            province_weekly(province_ts_deaths, obs.confirmed_deaths_history),
            "deaths"
        ),
        _with_series(
            province_weekly_fit(
                :cases, province_ts_inc, obs.confirmed_history.days
            ), "cases_fit"
        ),
        _with_series(
            province_weekly_fit(
                :deaths, province_ts_deaths, obs.confirmed_deaths_history.days
            ), "deaths_fit"
        ),
        _with_series(
            reported_level_table(
                merge(
                    obs.province_isolation_history,
                    Dict("National" => obs.isolation_history)
                ); cutoff = obs.cutoff, n = obs.n, from = obs.n - 55
            ), "isolation"
        ),
        _with_series(
            reported_level_table(
                merge(
                    obs.province_bed_capacity_history,
                    Dict("National" => obs.bed_capacity_history)
                ); cutoff = obs.cutoff, n = obs.n, from = obs.n - 55
            ), "beds"
        );
        cols = :union
    )
);
## The same three tables into the release outputs, so each results release
## carries the dashboard map's province and national inputs at its cut-off.
output_dir = get(
    ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output")
)
mkpath(output_dir)
for f in (
        "province_estimates.csv", "national_estimates.csv",
        "province_timeseries.csv",
    )
    cp(joinpath(province_map_dir, f), joinpath(output_dir, f); force = true)
end

#md # ```@raw html
#md # </details>
#md # ```

# ---
#
# The full analysis code, data and model definitions are in the [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize) repository.
