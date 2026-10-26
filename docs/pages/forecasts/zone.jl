# # Health-zone forecasts
#
# This page splits the one-week-ahead confirmed-case forecast across the health zones, as defined in the [health-zone model](@ref "Health-zone model") Methods section.
# The national forecast is on the [forecasts](@ref "Forecasts") page and the split by province on the [province forecasts](@ref "Province forecasts") page.
# How these forecasts have scored against what each zone went on to report is in the [forecast by health zone](@ref "Forecast by health zone") evaluation.
# The zone estimates behind the split are on the [health-zone estimates](@ref "Health-zone estimates") page.

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations and the fit registry. See
## `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))
#-
## The zone fit this page reads, loaded from the cache here.
chn_local = load_fit("local");

#md # ```@raw html
#md # </details>
#md # ```

# ## Summary
#
# The bullets rank the zones, and the figure and table below give each zone's projection.

#md # ```@raw html
#md # <details><summary>Split the one-week-ahead forecast across the zones</summary>
#md # ```

## `zone_stage_inputs` is defined in the shared setup, so the zone
## estimates page draws the same forecast from the same fixed inputs.
zone_inputs = zone_stage_inputs(; forecast = true);
zone_patch = zone_inputs.patch_of_zone;
zone_week_end = obs.cutoff + Day(7);
ZONE_THRESHOLDS = (1, 5, 10, 20)
zone_fc = zone_forecast(chn_local, zone_inputs);
zone_fc_draws = zone_forecast_draws(zone_fc, zone_inputs);
zone_fc_probs = zone_forecast_probabilities(
    zone_fc, zone_inputs; thresholds = ZONE_THRESHOLDS, draws = zone_fc_draws
);
zone_forecast_summary = zone_forecast_table(
    zone_fc, zone_inputs; thresholds = ZONE_THRESHOLDS
);
## Cases allocated to each zone over the past one, two and four weeks, and
## the last vintage on which each zone's count rose.
zone_recent = Dict(
    w => zone_recent_cases(zone_inputs; window = w) for w in (7, 14, 28)
);
zone_last_case = zone_last_case_dates(zone_inputs);
zone_forecast_bullets = let n_zones = length(zone_inputs.zone_labels)
    pct(x) = string(round(Int, 100 * x), "%")
    lead = first(
        sortperm([median(v) for v in zone_fc_draws.zones]; rev = true), 3
    )
    ranked = join(
        [
            "$(zone_inputs.zone_labels[z]) " *
                "$(median_interval_text(zone_fc_draws.zones[z]))"
                for z in lead
        ], "; "
    )
    quiet = findall(iszero, zone_recent[14])
    quiet_text = isempty(quiet) ?
        "every zone has had a case allocated over the past two weeks." :
        let z = quiet[argmax(zone_fc_probs[quiet, 1])]
            "$(length(quiet)) zones have had no case allocated over the " *
            "past two weeks. $(zone_inputs.zone_labels[z]) is the " *
            "most likely of them to report one, at " *
            "$(pct(zone_fc_probs[z, 1]))."
    end
    join(
        [
            "Projected counts are for the week to $(zone_week_end).", "",
            "- **Most new confirmed cases:** $(ranked).",
            "- **Chance of a case:** " *
                "$(count(>=(0.5), zone_fc_probs[:, 1])) of $(n_zones) " *
                "zones have at least an even chance of reporting a " *
                "confirmed case, and " *
                "$(count(>=(0.5), zone_fc_probs[:, 3])) of reporting at " *
                "least ten.",
            "- **Quiet zones:** $(quiet_text)",
        ], "\n"
    )
end;

#md # ```@raw html
#md # </details>
#md # ```

Markdown.parse(zone_forecast_bullets) #hide

# ## One-week-ahead forecast by health zone
#
# The patch totals in the table are the projections on the [province forecasts](@ref "Province forecasts") page.

#md # ```@raw html
#md # <details><summary>One-week-ahead zone forecast</summary>
#md # ```

zone_forecast_fig = plot_zone_forecast(
    zone_fc_draws.zones,
    zone_inputs.zone_labels, zone_patch;
    patch_labels = zone_inputs.patch_labels, top = 15
);
## The fifteen largest zone forecasts, then every patch total.
zone_forecast_display = let t = zone_forecast_summary
    zones = sort(t[t.zone .!= "Patch total", :], :median; rev = true)
    vcat(first(zones, 15), t[t.zone .== "Patch total", :])
end;

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

zone_forecast_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The fifteen zones with the largest forecasts of confirmed cases over the week after the cut-off.
#md # </figcaption>
#md # </figure>
#md # <figure class="bvd-figure">
#md # ```

MarkdownTable(zone_forecast_display) #hide

#md # ```@raw html
#md # <figcaption>
#md # The fifteen zones with the largest forecasts, then each patch total.
#md # The last four columns are the probability that the zone reports at least 1, 5, 10 and 20 confirmed cases over the week.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Zones quiet for two weeks

#md # ```@raw html
#md # <details><summary>Quiet zones ranked by the chance of a case</summary>
#md # ```

zone_quiet_display = let
    quiet = findall(iszero, zone_recent[14])
    order = quiet[sortperm(zone_fc_probs[quiet, 1]; rev = true)]
    DataFrame(
        zone = zone_inputs.zone_labels[order],
        patch = zone_inputs.patch_labels[zone_patch[order]],
        cases = zone_inputs.cumulative[order],
        cases_past_4_weeks = zone_recent[28][order],
        last_case = [
            ismissing(d) ? "" : string(d) for d in zone_last_case[order]
        ],
        p_ge_1 = round.(zone_fc_probs[order, 1]; digits = 2),
        p_ge_5 = round.(zone_fc_probs[order, 2]; digits = 2),
    )
end;

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

MarkdownTable(zone_quiet_display) #hide

#md # ```@raw html
#md # <figcaption>
#md # Zones with no allocated case over the past two weeks, ranked by the probability that they report at least one, with their cases over the past four weeks and the date of their last case.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Saving zone forecast assets
#
# The summary dashboard shows the summary bullets.

#md # ```@raw html
#md # <details><summary>Write the zone forecast asset</summary>
#md # ```

dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
open(joinpath(dashboard_dir, "zone_forecast.md"), "w") do io
    print(io, zone_forecast_bullets)
end

#md # ```@raw html
#md # </details>
#md # ```

# ## Past zone forecasts against what was observed
#
# Each release archives its zone forecast, and later releases score it in the [forecast by health zone across releases](@ref "Forecast by health zone across releases") evaluation.

#md # ```@raw html
#md # <details><summary>Load the archived zone forecasts and their outcomes</summary>
#md # ```

## Written by `scripts/score_releases.jl` from each release's
## `zone_forecast.csv`, in the national overlay's schema. A missing file
## reads as an empty table, which the figure reports as nothing scored.
zone_overlay_df = _release_data(
    joinpath("zone", "forecast_overlay.csv"),
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, observed = Float64,
        median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
        hi60 = Float64, lo90 = Float64, hi90 = Float64,
    )
)
zone_past_keys = zone_inputs.zone_keys[
    first(sortperm([median(v) for v in zone_fc_draws.zones]; rev = true), 10),
]
zone_past_fig = plot_forecast_overlay(
    zone_score_rows(
        scored_overlay(zone_overlay_df), zone_past_keys,
        zone_inputs.zone_labels[indexin(zone_past_keys, zone_inputs.zone_keys)]
    );
    empty_message = "No zone forecast has been scored yet. Scores appear " *
        "from the first release that archives the zone forecast onwards."
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

zone_past_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Past forecasts for the ten zones with the largest forecasts in this build, by the cut-off each was made from.
#md # Each shows its median and 90% predictive interval beside the persistence baseline and the count the zone went on to report.
#md # A week is left out when the zone tables carry no vintage on its last day, when it holds a harmonisation-break day, or when its province moved unallocated cases into named zones during it.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Archiving the zone forecast
#
# Each release carries the zone forecast draws, for scoring once the week has been reported.

#md # ```@raw html
#md # <details><summary>Write the zone forecast archive</summary>
#md # ```

output_dir = get(
    ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output")
)
mkpath(output_dir)
CSV.write(
    joinpath(output_dir, "zone_forecast.csv"),
    zone_forecast_archive(
        zone_fc, zone_inputs; made_date = zone_inputs.cutoff, thin = 5
    )
);

#md # ```@raw html
#md # </details>
#md # ```

# ---
#
# The full analysis code, data and model definitions are in the [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize) repository.
