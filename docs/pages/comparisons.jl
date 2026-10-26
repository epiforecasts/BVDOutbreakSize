# # Comparisons
#
# Our estimates and projections set against those other groups have published for the same outbreak.

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
frozen_by_cutoff = frozen_fits_by_cutoff()
frozen_C(c) = vec(Array(frozen_by_cutoff[c].chn[:C_T]))
posterior_C_joint = vec(Array(chn_joint[:C_T]));

#md # ```@raw html
#md # </details>
#md # ```

# ## Comparison with McCabe et al.
#
# McCabe et al. published their estimates as scenarios at fixed situation-report cut-offs, each with a 95% confidence interval.
# Their back-calculation-from-deaths scenarios differ between reports, since the 18 May report used 88 reported deaths and the 20 May update 131 with corrected case-fatality ratios.
# Their scenarios estimate cumulative cases at the report dates, though the report does not say whether these are symptomatic cases or all infections.
# We take the like-for-like quantity to be our cumulative symptom onsets on the same dates, rather than the latent infections or our total at the current cut-off.

#md # ```@raw html
#md # <details><summary>McCabe scenarios with uncertainty against our estimates</summary>
#md # ```

function _ci90row(xs)
    return (
        round(Int, quantile(xs, 0.5)),
        round(Int, quantile(xs, 0.05)),
        round(Int, quantile(xs, 0.95)),
    )
end

## Our modelled cumulative symptom onsets on a McCabe report date, read off
## the joint fit's per-draw `cumulative_onsets` trajectory. The grid runs to
## the cut-off on day `n`, so the day-index for a date is `n` minus the days
## from that date back to the cut-off (`grid_day("2026-06-07") = n`,
## `"2026-05-20") = n - 18`, `"2026-05-18") = n - 20`).
_onset_trajs = let mat = chn_joint[:cumulative_onsets]
    [collect(v) for v in vec(collect(mat))]
end
## Inverse of `grid_date(day) = obs.cutoff - Day(obs.n - day)`: the day-index
## whose calendar date is `date`, using `value` (imported above) for the
## day count rather than the non-exported `Dates.date2epochdays`.
_grid_day(date) = obs.n - value(obs.cutoff - Date(date))
function _ours_on(date)
    d = _grid_day(date)
    return _ci90row(Float64[t[d] for t in _onset_trajs])
end

## Our matched cumulative-onset estimate for each report date, keyed by date so
## it lands beside that vintage's scenarios in its own panel.
mccabe_ours = Dict(
    "2026-05-18" => _ours_on("2026-05-18"),
    "2026-05-20" => _ours_on("2026-05-20"),
    "2026-05-27" => _ours_on("2026-05-27")
)

## One panel per report date; within a panel each method-and-family is one row,
## with the case-fatality / window / doubling-time sweep dodged onto that single
## line, so the ~40 scenarios keep their intervals without becoming ~40 rows.
matched_comparison_fig = plot_scenario_comparison(
    REPORT_SCENARIOS_CI;
    ours = mccabe_ours,
    date_titles = [
        "2026-05-18" => "18 May report",
        "2026-05-20" => "20 May update",
        "2026-05-27" => "27 May (Lancet)",
    ],
    xlabel = "Cumulative cases"
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

matched_comparison_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # McCabe et al.'s scenarios from the 18 May report, the 20 May update and the 27 May Lancet publication, one panel each, against our cumulative symptom onsets on each report date with their credible interval.
#md # Each method and scenario family is one line sweeping the case-fatality ratio, the geographic window and the doubling time.
#md # Their 95% confidence intervals come from exact negative-binomial counts for geographic spread, which uses exported cases and travel volume, and from a Poisson likelihood profile for back-calculation from deaths.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Frozen-fit C_T intervals (kept for the CSV export, not shown)</summary>
#md # ```

## The frozen-fit infection counts for the published
## `frozen_matched_cutoffs.csv` export, not shown on the page.
frozen_streams_table = streams_table(
    "frozen 20 May" => frozen_C("2026-05-20"),
    "frozen 23 May" => frozen_C("2026-05-23"),
    "frozen 27 May" => frozen_C("2026-05-27"),
    "frozen 8 June" => frozen_C(default_chamla_cutoff()),
    "current data" => posterior_C_joint
);

#md # ```@raw html
#md # </details>
#md # ```

# ## Comparison with Chamla et al.
#
# A second group, [chamla2026](@citet) at the World Health Organization Regional Office for Africa, published a stochastic compartmental model of the same outbreak on 25 June 2026.
# Their model is a discrete-time susceptible-exposed-infectious-recovered-dead ensemble, recalibrated by simulation filtering to the laboratory-confirmed cases and anchored on the 598 confirmed cases reported by 8 June.
# It projects the confirmed-case trajectory under a low, central and high transmissibility scenario.
#
# Their published quantity is the cumulative confirmed-case count with the reporting fraction held at one, so it excludes the unconfirmed and unascertained cases that this analysis and McCabe et al. estimate.
# The like-for-like comparison is therefore against our own confirmed-case projection.
# We take our fit frozen at 8 June, the date of their calibration anchor, and roll its confirmed cases forward to the dates Chamla report with the one-week-ahead forecast, rather than refitting to their assumptions.

#md # ```@raw html
#md # <details><summary>Project the 8 June fit forward and assemble the Chamla comparison</summary>
#md # ```

## The 8 June frozen joint fit matches Chamla's confirmed-case calibration
## anchor exactly and carries the confirmed-case testing history through then,
## so we forecast its confirmed-case stream to the dates Chamla report.
chamla_anchor = frozen_by_cutoff["2026-06-08"]

## Our projected cumulative confirmed cases at a horizon of `h` days past the
## 8 June cut-off, drawn from the 8 June model run past its cut-off and
## summarised as (median, 5%, 95%).
function _our_confirmed_h(h)
    fc = forecast_reported(
        fit_forecast("frozen_2026-06-08");
        horizon = h,
        obs_cases = chamla_anchor.o.reported_cases,
        obs_deaths = chamla_anchor.o.total_deaths,
        obs_confirmed = chamla_anchor.o.confirmed_cases,
        obs_confirmed_deaths = chamla_anchor.o.confirmed_deaths
    )
    return _ci90row(float.(fc.confirmed_cum))
end

## Our projection at Chamla's forward report dates (10 and 24 June, week 12):
## the anchor day is the fitted confirmed total at 8 June, each later date a
## forward forecast. Reused for the matched-date table and the week-12 figure.
chamla_fan = map(["2026-06-08", "2026-06-10", "2026-06-24"]) do d
    h = value(Date(d) - chamla_anchor.cutoff)
    row = h == 0 ?
        (
            chamla_anchor.o.confirmed_cases, chamla_anchor.o.confirmed_cases,
            chamla_anchor.o.confirmed_cases,
        ) : _our_confirmed_h(h)
    (d, row...)
end
_fan_at(date) =
let r = first(x for x in chamla_fan if x[1] == date)
    (r[2], r[3], r[4])
end
ours_10jun = _fan_at("2026-06-10")
ours_24jun = _fan_at("2026-06-24")

## Observed confirmed totals on Chamla's forward report dates, read from the
## confirmed-case history frozen at each date.
observed_10jun = freeze_observations("2026-06-10").confirmed_cases
observed_24jun = freeze_observations("2026-06-24").confirmed_cases

## Observed confirmed cases over the comparison window: the daily cumulative
## series read off the chain's grid from 18 May (Chamla's first projected point)
## to the cut-off.
chamla_obs_series = let
    ds = [grid_date(d) for d in obs.confirmed_history.days]
    cs = obs.confirmed_history.counts
    [(string(ds[i]), cs[i]) for i in eachindex(ds) if ds[i] >= Date("2026-05-18")]
end

## Chamla's central confirmed-case projection over the comparison window; their
## later, far-larger horizons are noted in the text rather than plotted so the
## window stays legible.
chamla_central_window = CHAMLA_CONFIRMED_CENTRAL[1:4]

chamla_projection_fig = plot_projection_comparison(;
    external = chamla_central_window,
    ours = chamla_fan,
    observed = chamla_obs_series,
    external_label = "Chamla et al. central (R₀=1.71)",
    ours_label = "Our projection (from 8 June)",
    observed_label = "Observed confirmed",
    title = "Confirmed-case projections versus observed, from mid-May"
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

chamla_projection_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Chamla et al.'s central confirmed-case projection and ours from the fit frozen at 8 June, against the observed confirmed cases from mid-May.
#md # </figcaption>
#md # </figure>
#md # ```

#md # ```@raw html
#md # <details><summary>Week-12 (24 June) scenario spread against ours and observed</summary>
#md # ```

chamla_w12_rows = vcat(
    [(label, m, lo, hi) for (label, m, lo, hi) in CHAMLA_CONFIRMED_W12],
    [("Our projection (from 8 June)", ours_24jun...)],
    [
        (
            "Observed by 24 June", observed_24jun,
            observed_24jun, observed_24jun,
        ),
    ]
)
chamla_w12_groups = vcat(
    fill("Chamla et al. scenarios", 3),
    ["Our projection"], ["Observed"]
)

chamla_w12_fig = plot_estimate_comparison(
    chamla_w12_rows;
    xlabel = "Cumulative confirmed cases by 24 June",
    groups = chamla_w12_groups,
    group_colours = [
        "Chamla et al. scenarios" => :steelblue,
        "Our projection" => :firebrick,
        "Observed" => :black,
    ]
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

chamla_w12_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # Chamla et al.'s scenario spread at 24 June against our projection and the observed count.
#md # Their central scenario projected just under a thousand confirmed cases, and their low and high scenarios roughly 870 to 1360.
#md # </figcaption>
#md # </figure>
#md # <details><summary>Matched-date projection numbers (10 and 24 June)</summary>
#md # ```

chamla_comparison_table = let
    fmt(t) = string(t[1], " (", t[2], "–", t[3], ")")
    central(date) =
    let r = first(
            x for x in CHAMLA_CONFIRMED_CENTRAL
                if x[1] == date
        )
        fmt((r[2], r[3], r[4]))
    end
    DataFrame(
        "Date" => ["10 June", "24 June"],
        "Chamla central (90% PI)" => [
            central("2026-06-10"),
            central("2026-06-24"),
        ],
        "Our projection (90% CrI)" => [fmt(ours_10jun), fmt(ours_24jun)],
        "Observed confirmed" => [
            string(observed_10jun),
            string(observed_24jun),
        ]
    )
end;

#md # ```@raw html
#md # <figure class="bvd-figure">
#md # ```

MarkdownTable(chamla_comparison_table) #hide

#md # ```@raw html
#md # <figcaption>
#md # Our projection and Chamla et al.'s scenarios on the matched dates of 10 and 24 June.
#md # </figcaption>
#md # </figure>
#md # </details>
#md # ```

# Beyond the comparison window their central scenario continues to roughly 8200 confirmed cases by mid-September, with the high scenario far higher.

# ## Reproduction number behind the projection

#md # ```@raw html
#md # <details><summary>Reproduction number as estimated by the 8 June fit</summary>
#md # ```

## Reconstruct the reproduction-number trajectory the 8 June fit estimated,
## mirroring the current-data R_t figure but with the frozen vintage's own grid,
## breakpoint and renewal start.
chamla_rt_obs = chamla_anchor.o
chamla_rt_breakpoint = chamla_rt_obs.n - chamla_rt_obs.who_first_sitrep_days
chamla_rt_start = clamp(
    chamla_rt_obs.n - round(Int, chamla_rt_obs.tmrca_days) + RENEWAL_START_LEAD,
    1, chamla_rt_obs.n
)
chamla_rt_fig = plot_rt(
    chamla_anchor.chn;
    n = chamla_rt_obs.n, breakpoint = chamla_rt_breakpoint,
    rt_start = chamla_rt_start,
    rt_walk_start = clamp(
        chamla_rt_breakpoint - RT_WALK_LEAD,
        chamla_rt_start, chamla_rt_obs.n
    ),
    as_of_date = string(chamla_rt_obs.cutoff),
    seeding = chamla_rt_obs.seeding
);

#md # ```@raw html
#md # </details>
#md # <figure class="bvd-figure">
#md # ```

chamla_rt_fig #hide

#md # ```@raw html
#md # <figcaption>
#md # The reproduction number estimated by the fit frozen at 8 June, which carries the forward projection above.
#md # Its decline over the weeks to the cut-off bends the projection away from sustained early growth.
#md # </figcaption>
#md # </figure>
#md # ```

# ## Saving comparison outputs
#
# The frozen-fit table is written to the shared output directory.

#md # ```@raw html
#md # <details><summary>Write comparison outputs</summary>
#md # ```

output_dir = get(
    ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output")
)
mkpath(output_dir)
CSV.write(
    joinpath(output_dir, "frozen_matched_cutoffs.csv"),
    frozen_streams_table
)

#md # ```@raw html
#md # </details>
#md # ```
