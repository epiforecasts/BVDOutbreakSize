# # Health-zone forecast evaluation
#
# How the health-zone split on the [health-zone forecasts](@ref "Health-zone forecasts") page has scored against what each zone went on to report.
# Scoring follows the national [forecast evaluation](@ref "Forecast evaluation"), against the same persistence baseline.
# The same scoring by province is on the [province forecast evaluation](@ref "Province forecast evaluation") page.
# The in-sample zone checks are on the [health-zone in-sample checks](@ref "Health-zone in-sample checks") page.

#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations and the fit registry. See
## `docs/pages/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "pages", "_setup.jl"))
#-
## The fits this page reads, loaded from the cache here: the frozen joint
## and the zone fit melded from it.
frozen_lastweek = load_fit("frozen_validation");
frozen_local = load_fit("local_frozen_validation");

#md # ```@raw html
#md # </details>
#md # ```

# ## Summary
#
# The overall bullets come first, then a block per province, from the scores further down this page.
# The log score of the split is higher when the model allocates the week's cases better than the persistence rules.

#md # ```@eval
#md # using Markdown, BVDOutbreakSize
#md # dir = joinpath(pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
#md # Markdown.parse(read(joinpath(dir, "evaluation_forecast_zone.md"), String))
#md # ```

# ## Forecast by health zone
#
# The one-week-ahead forecast split by health zone, scored against what each zone went on to report.
# Each zone's forecast is the frozen national projection times its province's modelled share times the zone's projected share of that province's confirmed reports over the week, from the [health-zone model](@ref "Health-zone model") melded from the frozen fit.
# The observed count is the change in each zone's cumulative confirmed cases between the frozen cut-off and the current data, clamped at zero.
# When the zone tables end before the current cut-off, the week scored is shortened to the days they cover.
# The scores and the two persistence rules they are set against are those of the [zone forecast scoring](@ref "Forecast scoring against a persistence baseline") in the analysis methods.
# The figure shows the fifteen zones with the largest forecast medians and the fold below the scores holds every zone.

#md # ```@raw html
#md # <details><summary>Zone forecast against observed</summary>
#md # ```

## `frozen_zone_stage_inputs` is defined in the shared setup, so the zone
## estimates page reads the same frozen inputs.
frozen_zone_inputs = frozen_zone_stage_inputs();
## The frozen fit's one-week-ahead national forecast, the same one the
## forecast evaluation page validates. `validation_forecast_from` is defined
## in the shared setup.
validation_forecast = validation_forecast_from("frozen_validation");
## The horizon is the days from the frozen cut-off to the last zone
## vintage, at most the week. A shorter week needs its own national
## forecast at that horizon, since the zone split scales its total; at the
## full week it is the forecast above. Zone tables that end at or before
## the frozen cut-off leave no day to score, and the section shows a note
## in place of its outputs.
zone_validation_horizon = min(
    7,
    value(frozen_zone_inputs.dates[end] - frozen_local.o.cutoff)
)
_zone_validation_missing = Markdown.parse(
    "_The zone tables end at or before the frozen cut-off, so the zone " *
        "forecast is not scored in this build._"
)
if zone_validation_horizon >= 1
    zone_validation_forecast = zone_validation_horizon == 7 ?
        validation_forecast :
        forecast_reported(
            fit_forecast("frozen_validation");
            horizon = zone_validation_horizon,
            obs_cases = frozen_lastweek.o.reported_cases,
            obs_deaths = frozen_lastweek.o.total_deaths,
            obs_confirmed = frozen_lastweek.o.confirmed_cases,
            obs_confirmed_deaths = frozen_lastweek.o.confirmed_deaths,
            obs_recovered = frozen_lastweek.o.recovered_cases
        )
    zone_truth = zone_forecast_truth(
        obs, frozen_zone_inputs;
        made_date = frozen_local.o.cutoff, horizon = zone_validation_horizon
    )
    zone_validation_table = zone_forecast_vs_truth(
        frozen_local.chn,
        frozen_lastweek.chn, zone_validation_forecast, frozen_zone_inputs;
        truth = zone_truth, horizon = zone_validation_horizon
    )
    zone_validation_scores = zone_forecast_scores(
        frozen_local.chn,
        frozen_lastweek.chn, zone_validation_forecast, frozen_zone_inputs;
        truth = zone_truth, horizon = zone_validation_horizon
    )
    zone_validation_fig = plot_zone_forecast(
        zone_validation_table;
        patch_labels = frozen_zone_inputs.patch_labels,
        xlabel = "New confirmed cases over $(zone_validation_horizon) days",
        title = "Zone forecast from $(frozen_local.o.cutoff) against observed"
    )
else
    ## The frames stay frames, since the release assets below write them.
    zone_validation_table = DataFrame()
    zone_validation_scores = DataFrame()
    zone_validation_fig = _zone_validation_missing
end;
## The scores rounded for display. Share persistence forecasts the split
## alone, so its total columns are blank. No patch is scored when every
## observed total is zero or incomplete, which leaves an empty frame
## without the columns.
zone_validation_scores_table = isempty(zone_validation_scores) ?
    _zone_validation_missing :
    let s = zone_validation_scores
        fmt(x, d) = ismissing(x) || isnan(x) ? "" :
        string(round(x; digits = d))
        DataFrame(
            "Province" => s.patch, "Rule" => s.method,
            "Observed total" => s.observed,
            "Log score of the split" => [fmt(x, 2) for x in s.log_score],
            "CRPS of the total" => [fmt(x, 1) for x in s.crps],
            "Log CRPS" => [fmt(x, 2) for x in s.log_crps],
            "Dispersion" => [fmt(x, 1) for x in s.dispersion],
            "Bias" => [fmt(x, 2) for x in s.bias],
            "Total within 90%" => [
                ismissing(x) ? "" : string(x)
                for x in s.coverage_90
            ]
        )
end;
zone_validation_table_display = isempty(zone_validation_table) ?
    _zone_validation_missing :
    zone_validation_table;

#md # ```@raw html
#md # </details>
#md # ```

zone_validation_fig #hide

MarkdownTable(zone_validation_scores_table) #hide

#md # ```@raw html
#md # <details><summary>Zone forecast against observed, every zone</summary>
#md # ```

MarkdownTable(zone_validation_table_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Saving zone forecast outputs
#
# The frozen zone fit's one-week-ahead forecast against the observed zone increments, and its scores against the two persistence rules, are written as release assets.

#md # ```@raw html
#md # <details><summary>Write zone forecast outputs</summary>
#md # ```

output_dir = get(
    ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output")
)
mkpath(output_dir)
CSV.write(
    joinpath(output_dir, "zone_forecast_validation.csv"),
    zone_validation_table
)
CSV.write(
    joinpath(output_dir, "zone_forecast_scores.csv"),
    zone_validation_scores
)

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Write the summary bullets</summary>
#md # ```

## The bullets under the summary heading at the top of the page. They read
## the scores built further down, so they are written here and read back
## when the site is assembled.
evaluation_forecast_zone_summary = let s = zone_validation_scores
    fmt(x) = ismissing(x) || (x isa Real && isnan(x)) ? "n/a" :
        string(round(x; digits = 2))
    if isempty(s)
        string(
            "- **Zone forecasts:** the zone tables end at or before the ",
            "frozen cut-off, so no zone forecast is scored in this build."
        )
    else
        patches = unique(s.patch)
        rule(p, m) = let r = s[(s.patch .== p) .& (s.method .== m), :]
            size(r, 1) == 0 ? nothing : first(eachrow(r))
        end
        beats(p) = let m = rule(p, "zone model")
            m === nothing ? false :
                all(
                    r -> r === nothing || m.log_score >= r.log_score,
                    [rule(p, "share persistence"), rule(p, "naive persistence")]
                )
        end
        observed_total = sum(s[s.method .== "zone model", :observed])
        overall = [
            string(
                "- **Split:** the zone model allocates the week better ",
                "than both persistence rules in ", count(beats, patches),
                " of ", length(patches), " provinces scored, over ",
                observed_total, " observed confirmed cases."
            ),
        ]
        function detail(p)
            m = rule(p, "zone model")
            m === nothing && return string("**", p, "**", "\n\n- Not scored.")
            return join(
                [
                    string("**", p, "**"), "",
                    string(
                        "- Log score of the split: ", fmt(m.log_score),
                        ", against ",
                        fmt(
                            rule(p, "share persistence") === nothing ?
                                missing : rule(p, "share persistence").log_score
                        ),
                        " for share persistence and ",
                        fmt(
                            rule(p, "naive persistence") === nothing ?
                                missing : rule(p, "naive persistence").log_score
                        ),
                        " for naive persistence."
                    ),
                    string(
                        "- Total over ", m.observed,
                        " observed cases: bias ", fmt(m.bias),
                        ", within the 90% interval: ",
                        ismissing(m.coverage_90) ? "n/a" :
                            string(m.coverage_90), "."
                    ),
                ], "\n"
            )
        end
        join(vcat([join(overall, "\n")], [detail(p) for p in patches]), "\n\n")
    end
end
dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
write(
    joinpath(dashboard_dir, "evaluation_forecast_zone.md"),
    evaluation_forecast_zone_summary
);

#md # ```@raw html
#md # </details>
#md # ```
