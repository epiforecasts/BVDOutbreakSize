# # National estimates
#
# This page gives the national estimates from the joint model.
# The methods for this model are on the [Methods](@ref "Methods") page.
#
# This page is generated from
# [`docs/pages/estimates/national.jl`](https://github.com/epiforecasts/BVDOutbreakSize/blob/main/docs/pages/estimates/national.jl).
# The model code it calls is in
# [`src/`](https://github.com/epiforecasts/BVDOutbreakSize/tree/main/src).
# See [aim and origins](@ref "Aim and origins") and [limitations](@ref "Limitations").
#
# **Offline copy.** A self-contained single-file HTML version of this report, built from the same run, is attached to each results release: [download the latest](https://github.com/epiforecasts/BVDOutbreakSize/releases/latest/download/analysis.html).
# It carries the front matter, the methods and the national results, and links to the other pages resolve against the hosted site.
#
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

# ## Results
#
# ### Summary
#
# The numbers below are our estimate of the underlying infections to date, reported and unreported, from the joint posterior.
# Each is given as equal-tailed 30%, 60% and 90% credible intervals.

#md # ```@raw html
#md # <details><summary>Compute the headline ranges</summary>
#md # ```

summary_ranges = let
    med(x) = quantile(x, 0.5)
    iqr(x) = quantile(x, 0.75) - quantile(x, 0.25)
    ## Posterior-minus-prior shift in units of the parameter's prior IQR,
    ## reusing the prior draws so nothing is respecified here.
    shift(post, prior) = round(
        (med(post) - med(prior)) / iqr(prior);
        digits = 2
    )

    C = posterior_C_joint
    Td = vec(Array(chn_joint[:T]))
    r0d = vec(Array(chn_joint[:r0]))
    rd = vec(Array(chn_joint[:r]))
    dt = vec(Array(chn_joint[:doubling_time]))
    R0d = vec(Array(chn_joint[:R0]))
    RTd = vec(Array(chn_joint[:R_T]))
    cfrd = vec(Array(chn_joint[:CFR]))
    sC = posterior_summary(C)
    sT = posterior_summary(Td)
    sr0 = posterior_summary(r0d)
    sr = posterior_summary(rd)
    ## Doubling time is unbounded at zero growth, so its own quantiles do
    ## not bound it once the growth rate spans zero. Map the growth rate's
    ## interval instead, as the summary tables do.
    sdt0 = map(doubling_time, posterior_summary(r0d))
    sdt = map(doubling_time, posterior_summary(rd))
    sR0 = posterior_summary(R0d)
    sRT = posterior_summary(RTd)
    scfr = posterior_summary(cfrd)

    ints_i(s) = string(
        "30% ", round(Int, s.lo30), "–", round(Int, s.hi30),
        ", 60% ", round(Int, s.lo60), "–", round(Int, s.hi60),
        ", 90% ", round(Int, s.lo90), "–", round(Int, s.hi90)
    )
    ## Per-province cumulative infections, read off the patch deterministic
    ## one draw at a time so the provinces stay coupled draw for draw.
    C_patch = [collect(v) for v in vec(collect(chn_joint[:C_T_patch]))]
    sprov = [
        posterior_summary([v[p] for v in C_patch])
            for p in 1:N_PATCHES
    ]
    ## The same draws give each province's reproduction number at the cut-off
    ## and its case-fatality ratio, so the three read coupled draw for draw.
    Rt_patch_draws = [collect(v) for v in vec(collect(chn_joint[:R_T_patch]))]
    sprov_rt = [
        posterior_summary([v[p] for v in Rt_patch_draws])
            for p in 1:N_PATCHES
    ]
    cfr_patch_draws = [
        collect(v)
            for v in vec(collect(chn_joint[:CFR_patch]))
    ]
    sprov_cfr = [
        posterior_summary([100 * v[p] for v in cfr_patch_draws])
            for p in 1:N_PATCHES
    ]
    ints_f(
        s,
        d
    ) = string(
        "30% ", round(s.lo30; digits = d), "–", round(s.hi30; digits = d),
        ", 60% ", round(s.lo60; digits = d), "–", round(s.hi60; digits = d),
        ", 90% ", round(s.lo90; digits = d), "–", round(s.hi90; digits = d)
    )
    start_from(t) = obs.cutoff - Day(round(Int, t))
    ints_d(s) = string(
        "30% ", start_from(s.hi30), "–", start_from(s.lo30),
        ", 60% ", start_from(s.hi60), "–", start_from(s.lo60),
        ", 90% ", start_from(s.hi90), "–", start_from(s.lo90)
    )
    ## One interval as a bare `lo–hi`, for a table cell that takes its level
    ## from the column header rather than repeating it in every cell.
    bound(s, lvl, d) = string(
        round(getproperty(s, Symbol("lo", lvl)); digits = d), "–",
        round(getproperty(s, Symbol("hi", lvl)); digits = d)
    )
    bound_i(s, lvl) = string(
        round(Int, getproperty(s, Symbol("lo", lvl))), "–",
        round(Int, getproperty(s, Symbol("hi", lvl)))
    )
    ## One province block of the per-province table: a row per province, a
    ## column per interval level. Three of these stacked read down each
    ## province in one pass.
    prov_rows(cell) = join(
        [
            "| $(PROVINCE_LABELS[p]) | $(cell(p, 30)) | $(cell(p, 60)) | " *
                "$(cell(p, 90)) |" for p in 1:N_PATCHES
        ], "\n"
    )
    f_lo = round(sC.lo90 / obs.confirmed_cases; digits = 1)
    f_hi = round(sC.hi90 / obs.confirmed_cases; digits = 1)

    ## How far the data has moved each estimate from its prior, in prior
    ## interquartile ranges, reusing the prior draws.
    moves = [
        "cumulative infection count" => shift(C, vec(Array(prior_chn[:C_T]))),
        "outbreak age" => shift(Td, vec(Array(prior_chn[:T]))),
        "doubling time" => shift(dt, vec(Array(prior_chn[:doubling_time]))),
    ]
    biggest = argmax(p -> abs(p.second), moves)

    Markdown.parse(
        """
        - **Cumulative infections:** the outbreak is estimated to have caused
          $(ints_i(sC)) infections to date, reported and unreported.
        - Against the $(obs.confirmed_cases) laboratory-confirmed cases by the
          cut-off that is roughly $(f_lo)–$(f_hi)× as many infections, so
          confirmed cases are estimated to capture only a small share of the
          outbreak.
        - **Outbreak start and age:** the outbreak is estimated to have begun on
          a start date of $(ints_d(sT)), an elapsed age to the cut-off of
          $(ints_i(sT)) days.
        - **Growth rate and doubling time:** the initial growth rate is
          estimated to have been $(ints_f(sr0, 3)) per day, an initial doubling
          time of $(ints_f(sdt0, 1)) days.
          The latest growth rate is estimated to be $(ints_f(sr, 3)) per day, a
          latest doubling time of $(ints_f(sdt, 1)) days.
        - **Reproduction number:** the initial reproduction number is estimated
          to have been $(ints_f(sR0, 2)) and the latest to be $(ints_f(sRT, 2)).
        - **Case-fatality ratio:** the case-fatality ratio is estimated to be
          $(ints_f(scfr, 2)).
        - **Shift from priors:** how far the data has moved each estimate from
          its prior, in prior interquartile ranges, where a value of one means
          the posterior median sits one prior interquartile range from the prior
          median, zero means unchanged, and the sign gives the direction.
          The fit moves the cumulative infection count by $(moves[1].second),
          the outbreak age by $(moves[2].second) and the doubling time by
          $(moves[3].second); the largest move is in the $(biggest.first).

        **By province.** Equal-tailed credible intervals at the cut-off.

        Infections to date:

        | Province | 30% | 60% | 90% |
        |---|---|---|---|
        $(prov_rows((p, l) -> bound_i(sprov[p], l)))

        Reproduction number:

        | Province | 30% | 60% | 90% |
        |---|---|---|---|
        $(prov_rows((p, l) -> bound(sprov_rt[p], l, 2)))

        Case-fatality ratio (%):

        | Province | 30% | 60% | 90% |
        |---|---|---|---|
        $(prov_rows((p, l) -> bound(sprov_cfr[p], l, 1)))
        """
    )
end;

#md # ```@raw html
#md # </details>
#md # ```

summary_ranges #hide

# #### Fit diagnostics
#
# Fit diagnostics for the joint fit and each individual fit.
# These indicate how reliable the results are from the perspective of the inference algorithm.
# The [breakdown by parameter](@ref "Fit diagnostics by parameter") can be used to further diagnose any issues.

#md # ```@raw html
#md # <details><summary>Build the fit diagnostics table</summary>
#md # ```

fit_diagnostics_table = diagnostics_table(
    "joint" => chn_joint,
    "joint, no patches" => chn_no_patches,
    "exports" => chn_exports,
    "deaths (DRC)" => chn_deaths,
    "cases (DRC)" => chn_cases,
    "confirmed (DRC)" => chn_confirmed,
    "confirmed deaths (DRC)" => chn_confirmed_deaths,
    "isolation (DRC)" => chn_treatment,
    "onsets (DRC)" => chn_onsets,
    "frozen (1wk back)" => frozen_lastweek.chn,
    (
        RUN_SENSITIVITY ?
            [
                "delay sensitivity" => chn_joint_community_delay,
                "clock sensitivity (ExpGrowth)" => chn_joint_exp_growth_clock,
            ] :
            []
    )...
);

#md # ```@raw html
#md # </details>
#md # ```

fit_diagnostics_table #hide

# #### Data currency
#
# The cut-off is the last date any stream reports.
# Streams that stopped before it are carried frozen.

#md # ```@raw html
#md # <details><summary>Streams that stop before the cut-off</summary>
#md # ```

stream_currency_table = let status = stream_report_status(obs),
        stale = status[.!status.reporting, :]
    DataFrame(
        "Stream" => stale.label,
        "Last reported" => [
            ismissing(d) ? "never" : string(d) for d in stale.last_date
        ],
        "Days before cut-off" => [
            ismissing(d) ? "-" : string(d) for d in stale.days_since
        ]
    )
end;

#md # ```@raw html
#md # </details>
#md # ```

MarkdownTable(stream_currency_table) #hide

# ### Joint model estimates
#

#md # ```@raw html
#md # <details><summary>Cumulative infection count summary table</summary>
#md # ```

cumulative_cases_summary = summary_table(
    chn_joint, [:C_T]; digits = 0
);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_cases_summary #hide

# The figure below shows the cumulative trajectories and current-cut-off densities for three latent quantities: infections, symptom onsets and deaths.
# The infection density is the headline outbreak size, a count of infections rather than reported cases.

#md # ```@raw html
#md # <details><summary>Cumulative infections, onsets and deaths figure</summary>
#md # ```

cumulative_traj_fig = plot_cumulative_trajectories(
    chn_joint;
    n = obs.n, seeding = obs.seeding
);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_traj_fig #hide

# The national count above is the sum of the four patches' renewal equations.
# The per-province sizes, the modelled infections behind them and the importation between provinces are on the [province estimates](@ref "Province estimates") page.

# The cumulative infection count is set by the reproduction number trajectory and the outbreak age.
# The left panel below shows the posterior for the outbreak start date.
# The right panel shows the joint posterior of the outbreak age and the early doubling time.

#md # ```@raw html
#md # <details><summary>Outbreak start date and seeding-time posterior</summary>
#md # ```

start_date_fig = plot_start_date_pair(
    chn_joint;
    as_of_date = string(obs.cutoff)
);

#md # ```@raw html
#md # </details>
#md # ```

start_date_fig #hide

# The table below reports credible intervals on the infection-process parameters: the growth rate and doubling time, the reproduction number, the outbreak age, the case-fatality ratio and the cumulative infection count.
# The pair plot beside it shows their joint posterior, with the prior overlaid so the data's contribution to each marginal is visible.

#md # ```@raw html
#md # <details><summary>Infection-parameter summary table</summary>
#md # ```

infection_summary = summary_table(
    chn_joint,
    [:r, :doubling_time, :T, :R_T, :CFR, :C_T]; digits = 2
);

#md # ```@raw html
#md # </details>
#md # ```

infection_summary #hide

#md # ```@raw html
#md # <details><summary>Infection-parameter pair plot (prior overlaid)</summary>
#md # ```

infection_pair_fig = plot_pair(
    chn_joint,
    [
        :R_T, :r, :T, :CFR,
        Symbol("rt_state.sigma_rw"), Symbol("rt_state.intervention_effect"),
    ];
    prior = prior_chn, labels = display_names
);

#md # ```@raw html
#md # </details>
#md # ```

infection_pair_fig #hide

# The infection model carries two delays: the generation interval and the incubation period.
# The table and pair plot below report their posterior means and standard deviations.

#md # ```@raw html
#md # <details><summary>Infection-delay summary table</summary>
#md # ```

infection_delay_summary = summary_table(
    chn_joint,
    [
        Symbol("gi_state.α"), Symbol("gi_state.θ"),
        Symbol("inc_state.delay_mean"), Symbol("inc_state.delay_sd"),
    ];
    digits = 2, labels = display_names
);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Show infection-delay summary table</summary>
#md # ```

infection_delay_summary #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Infection-delay pair plot (prior overlaid)</summary>
#md # ```

infection_delay_pair_fig = plot_pair(
    chn_joint,
    [
        Symbol("gi_state.α"), Symbol("gi_state.θ"),
        Symbol("inc_state.delay_mean"), Symbol("inc_state.delay_sd"),
    ];
    prior = prior_chn, labels = display_names
);

#md # ```@raw html
#md # </details>
#md # ```

infection_delay_pair_fig #hide

# ### Reproduction number over time
#
# The figure below shows the daily reproduction number over the established outbreak, from the genetic bound to the cut-off.
# The 30%, 60% and 90% credible ribbons are shown with about a hundred sampled trajectories, and the no-growth threshold at one as a grey dashed line.
# The first situation report on 18 May 2026 marks the start of the response scale-up (red dashed).
# The end of the three-week scale-up is the red dotted line, and the data cut-off is grey dashed.

#md # ```@raw html
#md # <details><summary>Reproduction-number trajectory</summary>
#md # ```

## `rt_start` is the renewal/established-window start the plot shows from;
## `rt_walk_start` is where the random walk's knots begin — `RT_WALK_LEAD`
## days (a month) before the first situation report, matching `bvd_joint`'s
## `rt_walk_lead` — so the chain reconstruction uses the same knot grid the
## model did, floored at the renewal start. R_t is flat at R0 between the two.
rt_fig = plot_rt(
    chn_joint;
    n = obs.n, breakpoint = _BREAKPOINT,
    rt_start = _rt_start_plot,
    rt_walk_start = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n),
    as_of_date = string(obs.cutoff), seeding = obs.seeding,
    ramp = RT_INTERVENTION_RAMP
);

#md # ```@raw html
#md # </details>
#md # ```

rt_fig #hide

# The same trajectory split by province, the per-province summary and the spatial hyperparameters are on the [province estimates](@ref "Reproduction number by province") page.

# The table reports the posterior of the response effect on the reproduction number as a multiplier, where a value below one is the factor by which the response lowers the reproduction number once the scale-up completes.

#md # ```@raw html
#md # <details><summary>Intervention-effect summary table</summary>
#md # ```

intervention_effect = vec(
    Array(
        chn_joint[Symbol("rt_state.intervention_effect")]
    )
);
intervention_table = streams_table(
    "Rt multiplier exp(effect)" => exp.(intervention_effect);
    digits = 2
);

#md # ```@raw html
#md # </details>
#md # ```

intervention_table #hide

# ### Observation delays
#
# The delays carry latent infections through to each observed event: reporting, death, hospitalisation abroad and laboratory receipt.
# The onset-to-report and onset-to-hospitalisation delays are the same line-list onset-to-admission delay, sampled on its natural Gamma shape and scale.
# The onset-to-death delay is the convolution of two atomic Gamma delays, onset-to-admission and admission-to-death, each with its own shape and scale.
# The report-to-receipt delay is sampled by its mean and standard deviation.
# The length-of-stay delays are also shown: the isolation-bed BVD treatment stay (prior: the line-list admission-to-death delay), the non-BVD rule-out stay (prior: the report-to-receipt turnaround), and the confirmation-to-recovery delay.

#md # ```@raw html
#md # <details><summary>Observation-delay summary table</summary>
#md # ```

obs_delay_summary = summary_table(
    chn_joint,
    [
        Symbol("cases_state.report_state.α"),
        Symbol("cases_state.report_state.θ"),
        Symbol("deaths_state.od_state.oa.α"),
        Symbol("deaths_state.od_state.oa.θ"),
        Symbol("deaths_state.od_state.ad.α"),
        Symbol("deaths_state.od_state.ad.θ"),
        Symbol("exports_state.detect_state.α"),
        Symbol("exports_state.detect_state.θ"),
        Symbol("confirmed_state.receipt_state.d.delay_mean"),
        Symbol("confirmed_state.receipt_state.d.delay_sd"),
        :isolation_bvd_los_mean,
        :isolation_ruleout_los_mean,
        :recovery_delay_mean,
    ];
    digits = 2, labels = display_names
);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Show observation-delay summary table</summary>
#md # ```

obs_delay_summary #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Observation-delay pair plot (prior overlaid)</summary>
#md # ```

obs_delay_pair_fig = plot_pair(
    chn_joint,
    [
        Symbol("cases_state.report_state.α"),
        Symbol("deaths_state.od_state.oa.α"),
        Symbol("exports_state.detect_state.α"),
        Symbol("confirmed_state.receipt_state.d.delay_mean"),
        :isolation_bvd_los_mean,
        :isolation_ruleout_los_mean,
        :recovery_delay_mean,
    ];
    prior = prior_chn, labels = display_names
);

#md # ```@raw html
#md # </details>
#md # ```

obs_delay_pair_fig #hide

# ### Surveillance parameters
#
# The surveillance-data parameters cover the reporting fractions for the DRC and Uganda, the surveillance dispersions, and the laboratory pipeline: the testing fraction and receipt delay, the specimens analysed per suspect sampled, the per-suspected and per-test positivity, the non-BVD background rate, and the death-confirmation probability.
# The six passive-surveillance count streams (suspected cases, suspected deaths, confirmed cases, confirmed deaths, isolation occupancy and recovered) each have their own negative-binomial dispersion, partially pooled from a shared population.
# $k$ is the population-level dispersion, $k_{\text{cases}}$, $k_{\text{deaths}}$, $k_{\text{confirmed}}$ and $k_{\text{confirmed deaths}}$ the per-stream values for the four DRC count streams, and a pooling spread completes the group.
# The isolation and recovered streams add the proportion of suspects admitted to a bed and the recovery probability among confirmed cases.
# Their dispersions ($k_{\text{iso}}$, $k_{\text{rec}}$) are drawn from the same pooled population as the length-of-stay delays above.

#md # ```@raw html
#md # <details><summary>Surveillance-parameter summary table</summary>
#md # ```

surveillance_summary = summary_table(
    chn_joint,
    [
        :p_drc, :p_uganda, :k, :k_cases, :k_deaths, :k_confirmed,
        :k_confirmed_deaths, :dispersion_sd, :tau_test,
        :specimens_per_suspect, :lambda_bg,
        :suspected_positivity, :test_positivity, :expected_confirmed_T,
        :expected_analysed_T, :death_ascertainment, :background_cfr,
        :tau_death, :death_composition,
        :death_confirmation, :expected_confirmed_deaths_T,
        :isolation_admission, :isolation_dispersion, :expected_isolation_T,
        :expected_bed_demand_T, :bed_capacity, :bed_shortfall_T,
        :incare_cfr, :incare_cfr_modifier, :incare_confirm_modifier,
        :isolation_death_los_mean,
        :isolation_recovery_los_mean, :abscond_fraction,
        :recovery_probability, :recovered_dispersion, :expected_recovered_T,
    ];
    digits = 3
);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Show surveillance-parameter summary table</summary>
#md # ```

surveillance_summary #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Surveillance-parameter pair plot (prior overlaid)</summary>
#md # ```

surveillance_pair_fig = plot_pair(
    chn_joint,
    [
        :p_drc, :p_uganda, :k, :tau_test, :lambda_bg, :test_positivity,
        :death_confirmation,
    ];
    prior = prior_chn
);

#md # ```@raw html
#md # </details>
#md # ```

surveillance_pair_fig #hide

# ### Model against observed
#
# Whether the fit reproduces the data it was fitted to, stream by stream, is on the [in-sample checks](@ref "In-sample checks") page.

# ### Symptom-onset reporting delay and ascertainment
#
# The table reports the onset-report hazard's hyperparameters and two derived quantities.
# These are the share of a representative onset date's eventual reports that arrive within 7 days, and the median modelled ascertainment over the onset dates the ascertainment walk spans.
# The first comes from the delay hazard and the second from the ascertainment level anchored on the confirmed pipeline, so they are separate estimates (see the [symptom-onset reporting delay](@ref "Symptom-onset reporting delay") Methods section).
# The ascertainment offset is the row to read first, since it is the triangle's departure from the confirmed pipeline's own ascertainment and its prior is centred on no departure at all.
# The noise-walk and per-scan level rows are diagnostics.

#md # ```@raw html
#md # <details><summary>Reconstruct the onset-report hazard and calendar walk</summary>
#md # ```

## The report-date grid the calendar walk spans, `_onset_hazard_grid_start`
## to `_onset_grid_end`, is a fixed function of the digitised triangle, so
## the shared setup builds it once from `obs.onset_curve_history`.
## `_onset_grid_start` (the ascertainment walk's own, always the earliest
## scored onset date) is generally earlier; see `onset_hazard_grid_start`.
## Every posterior draw's `logit_h0` (the baseline delay hazard) and `γ`
## (the report-date calendar walk), rebuilt from the non-centred
## innovations the chain stores. `reconstruct_onset_hazard` is the package
## function the onset forecast also uses, so the hazard plotted here and
## the one projected forward are the same object rather than two copies of
## the same reconstruction that could drift apart.
_onset_hazard = reconstruct_onset_hazard(
    chn_joint;
    grid_start = _onset_hazard_grid_start, grid_end = _onset_grid_end,
    alpha_grid_start = _onset_grid_start
)

## Cell indices grouped by their snapshot's report day, the daily onsets
## series per posterior draw (`diff` of the chain's `cumulative_onsets`
## trajectory; the model stores only the running sum), and the chain's
## per-vintage scan level and noise scale. Built here so the noise-scale
## table below and the predictive figures further down share one
## computation.
_onset_cells_by_report = Dict{Int, Vector{Int}}()
for (i, r) in enumerate(obs.onset_curve_history.report_days)
    push!(get!(_onset_cells_by_report, r, Int[]), i)
end
_onset_report_grid_days = sort(collect(keys(_onset_cells_by_report)))
_onset_daily_draws = [
    vcat(v[1], diff(v))
        for v in vec(collect(chn_joint[:cumulative_onsets]))
]
_onset_scan_level = [
    collect(v) for v in vec(collect(chn_joint[:onset_scan_level]))
]
_onset_noise_scale = [
    collect(v) for v in vec(collect(chn_joint[:onset_noise_scale]))
]

## The digitised, deduplicated, cut-off-filtered snapshot blocks
## `load_onset_curve` scores, kept here for their raw cumulative
## onset-date counts: the fitted stream only ever sees between-vintage
## increments, so the observed cumulative levels the figures below plot
## are read back from the source blocks directly rather than
## reconstructed from the fitted increments.
_onset_path = joinpath(
    pkgdir(BVDOutbreakSize), "data",
    "onset_curve_scanned.csv"
)
_onset_snaps = filter(
    b -> b.report_date <= obs.cutoff,
    BVDOutbreakSize._dedup_onset_blocks(
        BVDOutbreakSize._read_onset_curve_blocks(_onset_path)
    )
)
## Keyed by report day rather than kept in order: a snapshot whose printed
## extent misses the scored window contributes no cells, so the panels and
## the snapshot blocks are not guaranteed to line up positionally.
_onset_snap_by_day = Dict(
    obs.n - value(obs.cutoff - b.report_date) => b for b in _onset_snaps
)

## A representative onset day (the median scored onset date), so the 7-day
## fraction below reflects a typical, not an edge, calendar day.
_onset_u_ref = isempty(obs.onset_curve_history.onset_days) ?
    _onset_grid_start :
    round(Int, quantile(obs.onset_curve_history.onset_days, 0.5))
## The delay profile is the normalised delay CDF, which reaches one by
## construction, so the 7-day fraction is read straight off it.
_onset_7d_fraction = [
    onset_report_G(
        6, _onset_hazard.logit_h0[i],
        _onset_hazard.γ[i], _onset_u_ref,
        _onset_hazard_grid_start
    )
        for i in eachindex(_onset_hazard.logit_h0)
]
filter!(isfinite, _onset_7d_fraction)

## Per-draw median ascertainment over the onset dates the ascertainment
## walk spans, then summarised across draws. Ascertainment is its own
## level rather than the delay hazard's asymptote (see the [symptom-onset
## reporting delay](@ref "Symptom-onset reporting delay") Methods section).
_onset_ascertainment_draws = [
    quantile(_onset_hazard.alpha[i], 0.5)
        for i in eachindex(_onset_hazard.alpha)
]

_onset_labels = merge(
    display_names,
    Dict(
        Symbol("onset_report_state.η0") => "onset-report hazard baseline (logit)",
        Symbol("onset_report_state.σ_h0") => "onset-report hazard pooling SD",
        Symbol("onset_report_state.σ_γ") => "onset-report calendar-walk step size",
        Symbol("onset_report_state.β") => "onset ascertainment offset (logit)",
        Symbol("onset_report_state.σ_a") => "onset ascertainment walk step size",
        Symbol("onset_report_state.log_τ0") => "onset-report noise scale (log, first snapshot)",
        Symbol("onset_report_state.σ_τ") => "onset-report noise-walk step size",
        Symbol("onset_report_state.σ_scan") => "shared per-scan level error"
    )
);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Symptom-onset reporting-delay summary table</summary>
#md # ```

## Derived-quantity rows built with the same pre-prettify column schema
## `summary_table` uses, so they `vcat` cleanly onto its output.
onset_derived_raw = DataFrame(
    quantity = String[],
    lower_90 = Float64[], lower_60 = Float64[],
    lower_30 = Float64[], upper_30 = Float64[],
    upper_60 = Float64[], upper_90 = Float64[]
)
for (label, draws) in [
        (
            "share of eventual reports arriving within 7 days (median onset date)",
            _onset_7d_fraction,
        ),
        ("median modelled ascertainment", _onset_ascertainment_draws),
    ]
    s = posterior_summary(draws)
    push!(
        onset_derived_raw,
        (
            label, round(s.lo90; digits = 3), round(s.lo60; digits = 3),
            round(s.lo30; digits = 3), round(s.hi30; digits = 3),
            round(s.hi60; digits = 3), round(s.hi90; digits = 3),
        )
    )
end
onset_derived_table = BVDOutbreakSize._prettify(onset_derived_raw)

onset_summary = vcat(
    summary_table(
        chn_joint,
        [
            Symbol("onset_report_state.η0"), Symbol("onset_report_state.σ_h0"),
            Symbol("onset_report_state.σ_γ"),
            Symbol("onset_report_state.log_τ0"),
            Symbol("onset_report_state.σ_τ"),
            Symbol("onset_report_state.σ_scan"),
        ];
        digits = 3, labels = _onset_labels
    ),
    onset_derived_table
);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Show symptom-onset reporting-delay summary table</summary>
#md # ```

onset_summary #hide

#md # ```@raw html
#md # </details>
#md # ```

# The digitisation-noise scale `τ` is fitted per snapshot, reported here as its own table alongside `log_τ0`, `σ_τ` and `σ_scan` above: median and 90% credible interval of each surviving vintage's own `τ`, in report-date order.

#md # ```@raw html
#md # <details><summary>Symptom-onset digitisation-noise scale per snapshot</summary>
#md # ```

onset_tau_raw = DataFrame(
    report_date = String[], median = Float64[],
    lower_90 = Float64[], upper_90 = Float64[]
)
for (v_idx, R) in enumerate(_onset_report_grid_days)
    draws = [d[v_idx] for d in _onset_noise_scale]
    push!(
        onset_tau_raw,
        (
            string(grid_date(R)), round(median(draws); digits = 2),
            round(quantile(draws, 0.05); digits = 2),
            round(quantile(draws, 0.95); digits = 2),
        )
    )
end
onset_tau_table = BVDOutbreakSize._prettify(onset_tau_raw);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Show symptom-onset digitisation-noise scale per snapshot</summary>
#md # ```

onset_tau_table #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Symptom-onset reporting-delay pair plot (prior overlaid)</summary>
#md # ```

onset_pair_fig = plot_pair(
    chn_joint,
    [
        Symbol("onset_report_state.η0"), Symbol("onset_report_state.σ_h0"),
        Symbol("onset_report_state.σ_γ"),
        Symbol("onset_report_state.β"), Symbol("onset_report_state.σ_a"),
        Symbol("onset_report_state.log_τ0"),
        Symbol("onset_report_state.σ_τ"),
        Symbol("onset_report_state.σ_scan"),
    ];
    prior = prior_chn, labels = _onset_labels
);

#md # ```@raw html
#md # </details>
#md # ```

onset_pair_fig #hide

# The figure below reads the fitted hazard differently: the onset-to-report delay distribution implied for a report on a given calendar day, holding that day's calendar-walk level fixed across the delay axis (`onset_report_delay_pmf`).
# The top panel is the mean delay, the bottom the delay's SD, both with 50%/90% credible ribbons over report time.

#md # ```@raw html
#md # <details><summary>Symptom-onset reporting delay over report time</summary>
#md # ```

onset_delay_profile_fig = plot_onset_delay_profile(
    _onset_hazard;
    grid_start = _onset_hazard_grid_start, grid_end = _onset_grid_end,
    seeding = obs.seeding
);

#md # ```@raw html
#md # </details>
#md # ```

onset_delay_profile_fig #hide

# The figure below is the latest snapshot's printed curve against the model's fitted level for those same bars, right-truncated at its report day.

#md # ```@raw html
#md # <details><summary>Latest snapshot: complete curve against fitted level</summary>
#md # ```

_onset_last_R = last(_onset_report_grid_days)
_onset_last_snap = _onset_snap_by_day[_onset_last_R]
_onset_last_us = sort(
    [
        u for u in eachindex(_onset_daily_draws[1])
            if haskey(_onset_last_snap.onsets, grid_date(u))
    ]
)
_onset_last_observed = Float64[
    _onset_last_snap.onsets[grid_date(u)] for u in _onset_last_us
]
_onset_last_draws = [
    onset_level_predictive_draws(
        u, _onset_daily_draws, _onset_hazard, _onset_scan_level,
        _onset_noise_scale, length(_onset_report_grid_days);
        grid_start = _onset_hazard_grid_start,
        alpha_grid_start = _onset_grid_start,
        target_delay = _onset_last_R - u
    )
        for u in _onset_last_us
]
_onset_last_title = "Latest snapshot " *
    "($(string(_onset_last_snap.report_date))): " *
    "complete curve vs fitted level"
onset_last_snapshot_fig = plot_onset_level_band(
    grid_date.(_onset_last_us), _onset_last_observed, _onset_last_draws;
    title = _onset_last_title, band_colour = :mediumpurple
);

#md # ```@raw html
#md # </details>
#md # ```

onset_last_snapshot_fig #hide

# Each panel below is one digitised snapshot, plotted by onset date.
# The grey crosses are the counts that snapshot's own figure printed.
# The band is the model's posterior predictive of each onset date's eventual reported total: modelled onsets times ascertainment, sampled through `onset_increments_model`'s `missing` branch so it carries the fitted τ, `scan_level` and Student-t tail as the likelihood scores a level cell.
# The black points are the latest figure's own reading, shown for context rather than as the target the band is read against.
# A panel stops short at the snapshot's own last printed onset date.
# Every snapshot is fitted, but only the first and the eight most recent are shown here.

#md # ```@raw html
#md # <details><summary>Nowcasts of the digitised reporting-triangle snapshots</summary>
#md # ```

## Latest printed value for each onset date the digitised figures cover,
## and the report day that reading came off. Ordered by report date, so the
## last block carrying a date gives the current reading. That is not the
## newest snapshot for every date: the figures do not all print the same
## range of onset dates, so a date a later figure stops short of keeps its
## reading, and its shorter delay, from an earlier one. A date inside a
## block's printed extent but with no row is a zero-height bar and does
## count; a date outside that extent is not covered by that figure at all
## and is skipped (the same rule the loader applies, see the [Data](@ref
## methods-data) section).
_onset_last_printed = Dict{Int, Float64}()
_onset_last_report_day = Dict{Int, Int}()
for snap in _onset_snaps
    lo, hi = extrema(keys(snap.onsets))
    R = obs.n - value(obs.cutoff - snap.report_date)
    for d in lo:Day(1):hi
        u = obs.n - value(obs.cutoff - d)
        (1 <= u <= obs.n) || continue
        _onset_last_printed[u] = Float64(get(snap.onsets, d, 0))
        _onset_last_report_day[u] = R
    end
end

## The first snapshot plus the eight most recent, not every surviving
## vintage: up to 48 panels is too many, and the earliest and latest carry
## the most information (complete-curve level cells, most current
## nowcast). Every snapshot still enters the likelihood; this only trims
## which get their own panel.
_onset_selected_report_days = length(_onset_report_grid_days) <= 9 ?
    _onset_report_grid_days :
    sort(
        unique(
            vcat(
                first(_onset_report_grid_days),
                last(_onset_report_grid_days, 8)
            )
        )
    )

## One panel per selected snapshot: printed counts against the model's
## posterior predictive of the eventual reported total
## (`onset_level_predictive_draws`'s default `target_delay = nothing`).
## `v_idx` is the snapshot's position in the full
## `_onset_report_grid_days`, the index its `scan_level`/`τ` are stored
## under.
_onset_panels = map(_onset_selected_report_days) do R
    v_idx = findfirst(==(R), _onset_report_grid_days)
    snap = _onset_snap_by_day[R]
    us = sort(obs.onset_curve_history.onset_days[_onset_cells_by_report[R]])
    observed = Float64[get(snap.onsets, grid_date(u), 0) for u in us]
    nowcast = [
        onset_level_predictive_draws(
            u, _onset_daily_draws, _onset_hazard, _onset_scan_level,
            _onset_noise_scale, v_idx;
            grid_start = _onset_hazard_grid_start,
            alpha_grid_start = _onset_grid_start
        )
            for u in us
    ]
    (;
        title = string(snap.report_date), dates = grid_date.(us), observed,
        nowcast, latest = [_onset_last_printed[u] for u in us],
    )
end

onset_fit_fig = plot_onset_nowcast_grid(_onset_panels);

onset_fit_fig #hide

#md # ```@raw html
#md # </details>
#md # ```
#
# The posterior predictive below reads the same reporting triangle along the onset date instead of the report date.
# It compares the latest digitised bar for each onset date against the model's posterior predictive for that bar, read at the current cut-off's own delay, and against the modelled onsets themselves.
# The gap between the two bands is the part of the epidemic the latest figure does not carry, whether because it is never ascertained or because it has not been reported yet.

#md # ```@raw html
#md # <details><summary>Reconstruct symptom onsets by date of onset</summary>
#md # ```

_onset_by_date_days = sort(collect(keys(_onset_last_printed)))

## Modelled onsets on each of those days.
_onset_by_date_onsets = [
    [
        _onset_daily_draws[i][u]
            for i in eachindex(_onset_daily_draws)
    ]
        for u in _onset_by_date_days
]
## The count the latest figure should print for each date, read at the
## current cut-off's delay (`_onset_grid_end - u`) rather than the eventual
## total, through the same predictive measurement error a single digitised
## bar carries (`onset_level_predictive_draws`'s level-cell case). `v_idx`
## is the vintage that last printed each date, the lookup the snapshot
## panels use.
_onset_by_date_vidx = [
    findfirst(==(_onset_last_report_day[u]), _onset_report_grid_days)
        for u in _onset_by_date_days
]
_onset_by_date_reps = [
    onset_level_predictive_draws(
        u, _onset_daily_draws, _onset_hazard, _onset_scan_level,
        _onset_noise_scale, _onset_by_date_vidx[k];
        grid_start = _onset_hazard_grid_start,
        alpha_grid_start = _onset_grid_start,
        target_delay = _onset_grid_end - u
    )
        for (k, u) in enumerate(_onset_by_date_days)
]

onset_ppc_by_date_fig = let
    fig = CairoMakie.Figure(; size = (900, 380))
    ax = CairoMakie.Axis(
        fig[1, 1];
        title = "Symptom onsets by date of onset: modelled vs digitised",
        xlabel = "onset date", ylabel = "cases"
    )
    xs = Float64.(_onset_by_date_days)
    q(ds, p) = [quantile(d, p) for d in ds]
    CairoMakie.band!(
        ax, xs, q(_onset_by_date_onsets, 0.05),
        q(_onset_by_date_onsets, 0.95); color = (:seagreen, 0.2)
    )
    CairoMakie.lines!(
        ax, xs, q(_onset_by_date_onsets, 0.5);
        color = :seagreen, linewidth = 2
    )
    CairoMakie.band!(
        ax, xs, q(_onset_by_date_reps, 0.05),
        q(_onset_by_date_reps, 0.95); color = (:mediumpurple, 0.2)
    )
    CairoMakie.lines!(
        ax, xs, q(_onset_by_date_reps, 0.5);
        color = :mediumpurple, linewidth = 2
    )
    CairoMakie.scatter!(
        ax, xs,
        [_onset_last_printed[u] for u in _onset_by_date_days];
        color = :black, marker = :cross, markersize = 9
    )
    ## Calendar labels on a grid-day axis, at weekly ticks so they do not
    ## collide at this width.
    _ticks = _onset_by_date_days[1:7:end]
    ax.xticks = (Float64.(_ticks), string.(grid_date.(_ticks)))
    ax.xticklabelrotation = pi / 4
    CairoMakie.Legend(
        fig[2, 1],
        [
            CairoMakie.MarkerElement(color = :black, marker = :cross),
            CairoMakie.PolyElement(color = (:mediumpurple, 0.3)),
            CairoMakie.PolyElement(color = (:seagreen, 0.3)),
        ],
        [
            "digitised (latest)", "modelled posterior predictive",
            "modelled onsets",
        ];
        orientation = :horizontal, tellwidth = false, tellheight = true
    )
    fig
end;

#md # ```@raw html
#md # </details>
#md # ```

onset_ppc_by_date_fig #hide

# ### Counterfactual: lower bound under no further transmission
#
# $\Delta D$ is the committed future deaths if transmission stopped at the report date, defined in the [no-onward-transmission counterfactual](@ref "No-onward-transmission counterfactual") Methods section.

#md # ```@raw html
#md # <details><summary>Project no-onward deaths and summarise</summary>
#md # ```

no_onward = predict_no_onward_deaths(
    chn_joint; obs_deaths = obs.total_deaths
);

no_onward_table = streams_table(
    "no-onward total" => no_onward.total_projected;
    digits = 0
);

#md # ```@raw html
#md # </details>
#md # ```

no_onward_table #hide

# The left panel shows the *still expected* deaths $\Delta D$: future deaths in cases already infected by $T$, net of those already observed.
# The right panel shows the *projected total* $D(T) + \Delta D$, whose axis starts at the observed death count.

#md # ```@raw html
#md # <details><summary>No-onward projected-deaths plot</summary>
#md # ```

no_onward_fig = plot_no_onward_deaths(
    no_onward; obs_deaths = obs.total_deaths
);

#md # ```@raw html
#md # </details>
#md # ```

no_onward_fig #hide

# ### Confirmed case-fatality ratio
#
# The delay-corrected confirmed CFR, defined in the [delay-corrected confirmed CFR](@ref "Delay-corrected confirmed case-fatality ratio") Methods section, is set against the structural (infection-based) CFR and the naive confirmed ratio.
# The corrected ratio debiases the naive confirmed ratio for the real-time delay between a case being confirmed and a death being confirmed.
# The structural CFR is the onset-level estimate the joint model fits.

#md # ```@raw html
#md # <details><summary>Compute the confirmed-CFR comparison</summary>
#md # ```

confirmed_cfr = delay_corrected_confirmed_cfr(
    chn_joint;
    obs_confirmed = obs.confirmed_cases,
    obs_confirmed_deaths = obs.confirmed_deaths
);

confirmed_cfr_summary = confirmed_cfr_table(confirmed_cfr);

## Summary line carrying the data-anchored corrected estimate alongside the
## infection-based structural CFR, so both can be quoted together.
confirmed_cfr_line = let r = confirmed_cfr
    pct(x) = round(100 * x; digits = 1)
    corr = filter(isfinite, r.corrected)
    struc = filter(isfinite, r.structural)
    cs = posterior_summary(corr)
    ss = posterior_summary(struc)
    Markdown.parse(
        string(
            "**Delay-corrected confirmed CFR:** ",
            pct(quantile(corr, 0.5)), "% (90% CrI ",
            pct(cs.lo90), "–", pct(cs.hi90), "%), versus a naive confirmed ratio ",
            "of ", pct(r.naive_observed), "% and a structural (infection-based) ",
            "CFR of ", pct(quantile(struc, 0.5)), "% (90% CrI ",
            pct(ss.lo90), "–", pct(ss.hi90), "%)."
        )
    )
end;

#md # ```@raw html
#md # </details>
#md # ```

confirmed_cfr_line #hide

confirmed_cfr_summary #hide

# The figure below shows the posterior densities of the delay-corrected confirmed CFR and the structural CFR.
# The naive observed confirmed ratio is drawn as a solid vertical rule, and the median uncorrected modelled confirmed ratio as a dashed rule.
# The gap from the naive rule to the corrected density is the real-time delay debiasing.
# The gap to the structural density is the residual case/death ascertainment difference.

#md # ```@raw html
#md # <details><summary>Confirmed-CFR density plot</summary>
#md # ```

confirmed_cfr_fig = plot_confirmed_cfr(confirmed_cfr);

#md # ```@raw html
#md # </details>
#md # ```

confirmed_cfr_fig #hide

# The same three ratios by province are below.
# Each province has its own case-fatality ratio, partially pooled toward the national value, and its own death confirmation, pooled far more tightly.
# The delay-corrected ratio is the national corrected ratio scaled by a province's lethality and death confirmation over its case ascertainment, which is what varies by province once the delays are corrected for.
# The structural ratio is the national ratio times that province's lethality contrast.
#
# The death composition identifies only the product of the lethality and death-confirmation contrasts, so their split is set by their priors rather than by the data.
# The lethality prior is the looser of the two, so a provincial excess of deaths over cases is read first as lethality and only marginally as death-finding.
# The spread of the lethality contrast is reported below against its prior, so a posterior that has not moved can be read as the prior's rather than as a finding.
# The two spreads are on the same log scale, so their sizes are comparable directly.

#md # ```@raw html
#md # <details><summary>Province case-fatality spread</summary>
#md # ```

province_cfr_spread = summary_table(
    chn_joint,
    [:province_cfr_sd, :province_death_ascertainment_sd];
    digits = 3,
    labels = Dict(
        :province_cfr_sd => "Lethality spread",
        :province_death_ascertainment_sd => "Death-confirmation spread"
    )
);

#md # ```@raw html
#md # </details>
#md # ```

province_cfr_spread #hide

#md # ```@raw html
#md # <details><summary>Province case-fatality table</summary>
#md # ```

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

# ### Forecast results
#
# The one-week-ahead projection of every stream, the symptom-onset nowcast and the scoring against what arrived are on the [forecasts](@ref "Forecasts") page.

# ## Saving results
#
# The tables above are written to an output directory at the repo root so they can be archived and shared.
# On every push to the main branch a GitHub Actions workflow regenerates these files and publishes them as a GitHub Release, downloadable from the repository's releases page (<https://github.com/epiforecasts/BVDOutbreakSize/releases>).
# The release bundles the summary tables, a thinned set of posterior draws, the latent symptom-onset ("symptomatic cases") trajectory over time, the one- to four-week-ahead forecasts of the observed streams, and a copy of the input data manifest.

#md # ```@raw html
#md # <details><summary>Write outputs to output/</summary>
#md # ```

## Outputs default to `output/` in the package directory (where the
## docs build and Release workflow expect them). Set `BVD_OUTPUT_DIR`
## to redirect them, e.g. when running from a read-only package
## install.
output_dir = get(
    ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output")
)
mkpath(output_dir)

## Full parameter summary for the published CSV (infection, surveillance and
## export parameters together).
joint_summary = summary_table(
    chn_joint,
    [
        :r, :r0, :doubling_time, :T, :R_T, :CFR, :C_T,
        :p_drc, :p_uganda, :k, :tau_test, :lambda_bg,
        Symbol("exports_state.travel_state.daily_travellers"),
    ]; digits = 2
)
CSV.write(joinpath(output_dir, "posterior_summary.csv"), joint_summary)
CSV.write(
    joinpath(output_dir, "confirmed_cfr_summary.csv"),
    confirmed_cfr_summary
)

## Copy the input data so the release records what produced these
## results.
cp(
    joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml"),
    joinpath(output_dir, "observations.toml"); force = true
)

## Thinned posterior draws of the key joint parameters (every 10th
## draw) so downstream users can recompute their own summaries.
## `cumulative_onsets_T` is the cumulative symptom onsets by the cut-off,
## the latent "symptomatic cases" outcome (the onset analogue of `C_T`),
## read off the last day of each draw's `cumulative_onsets` trajectory.
_cum_onset_draws = vec(collect(chn_joint[:cumulative_onsets]))
cumulative_onsets_T = Float64[v[end] for v in _cum_onset_draws]
## Raw walk base `log_R0`, the renewal reproduction-number walk's starting
## point on the log scale, kept unexponentiated so downstream scoring takes
## its own exp. This is a distinct quantity from `r0`, the growth-clock
## initial rate, so both columns are published side by side.
log_R0_draws = vec(Array(chn_joint[Symbol("rt_state.log_R0")]))
posterior_draws = DataFrame(
    r = vec(Array(chn_joint[:r])),
    r0 = vec(Array(chn_joint[:r0])),
    doubling_time = vec(Array(chn_joint[:doubling_time])),
    T = vec(Array(chn_joint[:T])),
    R_T = vec(Array(chn_joint[:R_T])),
    CFR = vec(Array(chn_joint[:CFR])),
    p_drc = vec(Array(chn_joint[:p_drc])),
    p_uganda = vec(Array(chn_joint[:p_uganda])),
    C_T = vec(Array(chn_joint[:C_T])),
    cumulative_onsets_T = cumulative_onsets_T,
    confirmed_cfr_corrected = confirmed_cfr.corrected
)
posterior_draws[!, Symbol("rt_state.log_R0")] = log_R0_draws
posterior_draws = posterior_draws[1:10:end, :]
CSV.write(joinpath(output_dir, "posterior_draws.csv"), posterior_draws);

## One- to four-week-ahead forecasts of the observed streams, saved as a
## release asset so each release records the forecast it made and it can later
## be scored against what is observed. Only the incident and level quantities
## are archived (see `forecast_archive`), thinned to keep the asset compact.
forecast_horizons = (7, 14, 21, 28)
forecast_runs = [
    (
        h,
        forecast_reported(
            chn_joint; horizon = h,
            obs_cases = obs.reported_cases,
            obs_deaths = obs.total_deaths,
            obs_confirmed = obs.confirmed_cases,
            obs_confirmed_deaths = obs.confirmed_deaths,
            obs_recovered = obs.recovered_cases,
            grid_n = obs.n,
            onset_grid_start = _onset_hazard_grid_start,
            onset_grid_end = _onset_grid_end,
            onset_alpha_grid_start = _onset_grid_start
        ),
    )
        for h in forecast_horizons
]
CSV.write(
    joinpath(output_dir, "forecast.csv"),
    forecast_archive(forecast_runs; made_date = obs.cutoff, thin = 5)
);

## The per-province split of the same forecasts, in the `forecast.csv`
## schema plus the province each row is a share of, so a release records the
## provincial forecast it made alongside the national one.
CSV.write(
    joinpath(output_dir, "province_forecast.csv"),
    province_forecast_archive(
        chn_joint, forecast_runs;
        made_date = obs.cutoff, n_patches = N_PATCHES, thin = 5
    )
);

## The same one- to four-week-ahead forecast made from each FROZEN joint
## re-fit (the McCabe-matched cut-offs, the Chamla anchor and the one-week-back
## validation fit), each stamped with its OWN cut-off as the made date. This
## gives historical forecast evaluation using the current model at past data
## cut-offs, scored against what has since been observed, without
## reconstructing old release tags. Each frozen fit uses its own frozen
## observations for the cut-off counts. The May cut-offs predate the isolation
## and recovered streams, so those are simply absent for them; the per-stream
## guard in `forecast_archive` skips a stream a fit does not carry.
frozen_forecast_fits = unique(
    f -> f.o.cutoff,
    [frozen_results; frozen_by_cutoff[chamla_cutoff]; frozen_lastweek]
)
## The `fit` column tells the frozen joint and each frozen single-stream fit
## apart when scored. `score_release` falls back to one default where an
## archive carries no such column, so an older release still scores as the
## joint.
frozen_forecast_archive = DataFrame(
    made_date = Date[], horizon = Int[],
    target_date = Date[], stream = String[], draw = Int[], value = Float64[],
    fit = String[]
)
## The onset grid belongs to the triangle each frozen fit actually saw, not
## to the live one: the May cut-offs predate the digitised figure entirely,
## so their grid is empty and the onset block is simply absent for them.
## Returns `(alpha_grid_start, hazard_grid_start, grid_end)`; the two grid
## starts diverge exactly as for the live fit (see
## `onset_hazard_grid_start`).
function _frozen_onset_grid(o)
    isempty(o.onset_curve_history.onset_days) &&
        return (nothing, nothing, nothing)
    gs = minimum(o.onset_curve_history.onset_days)
    ge = max(maximum(o.onset_curve_history.report_days), gs)
    hgs = onset_hazard_grid_start(
        o.onset_curve_history.onset_days, o.onset_curve_history.report_days
    )
    return (gs, hgs, ge)
end
for f in frozen_forecast_fits
    _fgs, _fhgs, _fge = _frozen_onset_grid(f.o)
    runs = [
        (
            h,
            forecast_reported(
                f.chn; horizon = h,
                obs_cases = f.o.reported_cases,
                obs_deaths = f.o.total_deaths,
                obs_confirmed = f.o.confirmed_cases,
                obs_confirmed_deaths = f.o.confirmed_deaths,
                obs_recovered = f.o.recovered_cases,
                grid_n = f.o.n,
                onset_grid_start = _fhgs, onset_grid_end = _fge,
                onset_alpha_grid_start = _fgs
            ),
        )
            for h in forecast_horizons
    ]
    _rows = forecast_archive(runs; made_date = f.o.cutoff, thin = 5)
    _rows[!, :fit] = fill(FROZEN_FIT, size(_rows, 1))
    append!(frozen_forecast_archive, _rows)
end

## The frozen single-stream fits, forecast from their own chains as
## `stream_forecasts.csv` forecasts the live ones. They are registered only
## at the validation cut-off and only for still-reported streams.
_frozen_stream_of = Dict(
    "cases" => (:reported_cases, "reported cases"),
    "deaths" => (:suspected_deaths, "suspected deaths"),
    "confirmed" => (:confirmed_cases, "confirmed cases"),
    "confirmed_deaths" => (:confirmed_deaths, "confirmed deaths"),
    "treatment" => (:isolation_beds, "isolation beds")
)
for (_sid, _sf) in sort(collect(pairs(frozen_lastweek_streams)); by = first)
    _stream, _label = _frozen_stream_of[_sid]
    _o = _sf.o
    _bp = _o.n - _o.who_first_sitrep_days
    ## Each stream on its own cut-off count, the beds on their occupancy.
    _base = if _stream === :isolation_beds
        isempty(_o.isolation_history.counts) ? 0 :
            _o.isolation_history.counts[end]
    elseif _stream === :reported_cases
        _o.reported_cases
    elseif _stream === :suspected_deaths
        _o.total_deaths
    elseif _stream === :confirmed_cases
        _o.confirmed_cases
    else
        _o.confirmed_deaths
    end
    for h in forecast_horizons
        _vals = forecast_stream(
            _sf.chn, _stream; horizon = h,
            obs_value = _base, n = _o.n, breakpoint = _bp,
            rt_start = 1, rt_walk_start = 1
        )
        for (_d, _i) in enumerate(1:5:length(_vals))
            push!(
                frozen_forecast_archive,
                (
                    _o.cutoff, h, _o.cutoff + Day(h), _label, _d,
                    Float64(_vals[_i]), _sid,
                )
            )
        end
    end
end
CSV.write(
    joinpath(output_dir, "forecast_frozen.csv"),
    frozen_forecast_archive
);

## Per-fit release assets: the reproduction number, outbreak size and forecasts
## for every fit rather than the joint alone, so a release records what each
## dataset implies on its own and can later be scored against the joint. The
## single-stream fits walk Rt from day 1 while the joint walks from
## `RT_WALK_LEAD` days before the first situation report, so each fit carries
## the starts its own fit used. Each single-stream fit forecasts only the
## dataset it observes; the joint forecasts every shared stream. Recovered has
## no single-stream fit, so the joint is the only fit that carries it.
stream_thin = 5
_rt_walk_start_joint = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n)
## Observed bed occupancy at the cut-off, the level the isolation forecast
## anchors on.
_iso_at_cutoff = isempty(obs.isolation_history.counts) ? 0 :
    obs.isolation_history.counts[end]
## The reporting triangle's own cumulative total at the cut-off. It anchors
## the reported quantity rather than changing it: the onset forecast is the
## INCREMENT this total should add over the horizon, not the level (see the
## methods section on the nowcast and forecast).
_onset_at_cutoff = something(obs.onset_curve_history.last_total, 0)
## Cumulative recovered at the cut-off. The loader leaves it missing when the
## manifest carries no recovered vintages, and the forecast returns the
## increment rather than this base, so a zero stands in for that case.
_recovered_at_cutoff = coalesce(obs.recovered_cases, 0)
stream_fits = [
    (;
        fit = "joint", chn = chn_joint, rt_start = _rt_start_plot,
        rt_walk_start = _rt_walk_start_joint,
        streams = [
            (:reported_cases, "reported cases", obs.reported_cases),
            (:suspected_deaths, "suspected deaths", obs.total_deaths),
            (:confirmed_cases, "confirmed cases", obs.confirmed_cases),
            (:confirmed_deaths, "confirmed deaths", obs.confirmed_deaths),
            (:recovered, "recovered", _recovered_at_cutoff),
            (:isolation_beds, "isolation beds", _iso_at_cutoff),
            (:exports, "exports", obs.exported_cases),
            (:onset_reports, "onset reports", _onset_at_cutoff),
        ],
    ),
    (;
        fit = "cases", chn = chn_cases, rt_start = 1, rt_walk_start = 1,
        streams = [(:reported_cases, "reported cases", obs.reported_cases)],
    ),
    (;
        fit = "deaths", chn = chn_deaths, rt_start = 1, rt_walk_start = 1,
        streams = [(:suspected_deaths, "suspected deaths", obs.total_deaths)],
    ),
    (;
        fit = "confirmed", chn = chn_confirmed, rt_start = 1, rt_walk_start = 1,
        streams = [(:confirmed_cases, "confirmed cases", obs.confirmed_cases)],
    ),
    (;
        fit = "confirmed_deaths", chn = chn_confirmed_deaths, rt_start = 1,
        rt_walk_start = 1,
        streams = [
            (
                :confirmed_deaths, "confirmed deaths",
                obs.confirmed_deaths,
            ),
        ],
    ),
    (;
        fit = "treatment", chn = chn_treatment, rt_start = 1,
        rt_walk_start = 1,
        streams = [(:isolation_beds, "isolation beds", _iso_at_cutoff)],
    ),
    (;
        fit = "exports", chn = chn_exports, rt_start = 1, rt_walk_start = 1,
        streams = [(:exports, "exports", obs.exported_cases)],
    ),
    (;
        fit = "onsets", chn = chn_onsets, rt_start = 1, rt_walk_start = 1,
        streams = [(:onset_reports, "onset reports", _onset_at_cutoff)],
    ),
]

## Cut-off reproduction number per fit. The joint exposes it as `R_T`; the
## single-stream composers do not (the alias lives in `bvd_joint`, not the
## shared latent submodel), so theirs is rebuilt from the walk parameters every
## chain carries and read at the cut-off, the last day of the reconstructed
## path. `ramp` matches the model's 21-day intervention scale-up.
function _fit_rt_draws(f)
    f.fit == "joint" && return vec(Array(f.chn[:R_T]))
    rt = reconstruct_rt(
        f.chn; n = obs.n, breakpoint = _BREAKPOINT,
        rt_start = f.rt_start, rt_walk_start = f.rt_walk_start, ramp = RT_INTERVENTION_RAMP
    )
    return Float64[rt[i, obs.n] for i in axes(rt, 1)]
end

## Basic reproduction number per fit, the renewal walk's starting value, a
## single distribution per fit rather than a daily series. Every current
## single-stream composer shares the joint's latent submodel and so carries
## its own walk base, but `r0_walk_draws` probes rather than assumes, so a
## single-stream model built without its own renewal walk drops out of this
## quantity instead of breaking the release.
_stream_quantities = [
    (
        f.fit, _fit_rt_draws(f), vec(Array(f.chn[:C_T])),
        r0_walk_draws(f.chn),
    ) for f in stream_fits
]

## One row per fit and quantity, with the median and the 30/60/90% credible
## bounds the report's tables use.
function _stream_estimate_row(fit, quantity, draws)
    s = posterior_summary(draws)
    return (
        fit = fit, quantity = quantity, median = quantile(draws, 0.5),
        lo30 = s.lo30, hi30 = s.hi30, lo60 = s.lo60, hi60 = s.hi60,
        lo90 = s.lo90, hi90 = s.hi90,
    )
end
stream_estimates = DataFrame(
    [
        _stream_estimate_row(fit, q, d)
            for (fit, rt, ct, r0) in _stream_quantities
            for (q, d) in (
                ("R_T", rt), ("C_T", ct),
                ("R0", r0),
            )
            if !isnothing(d)
    ]
)
CSV.write(joinpath(output_dir, "stream_estimates.csv"), stream_estimates);

## Thinned reproduction-number and outbreak-size draws per fit, so downstream
## scoring can recompute its own summaries rather than reuse the intervals.
stream_draws = DataFrame(
    [
        (fit = fit, quantity = q, draw = d, value = v)
            for (fit, rt, ct, r0) in _stream_quantities
            for (q, vals) in (
                ("R_T", rt), ("C_T", ct),
                ("R0", r0),
            )
            if !isnothing(vals)
            for (d, v) in enumerate(vals[1:stream_thin:end])
    ]
)
CSV.write(joinpath(output_dir, "stream_draws.csv"), stream_draws);

## Per-fit forecasts of each fit's own observed stream, in the `forecast.csv`
## long schema plus the fit that made them. Rebuilding a single-stream fit's
## cut-off growth rate needs the grid length and the breakpoint, which are data
## rather than chain contents, so both are passed.
stream_forecasts = DataFrame(
    made_date = Date[], horizon = Int[],
    target_date = Date[], stream = String[], draw = Int[], value = Float64[],
    fit = String[]
)
for f in stream_fits, (stream, label, obs_value) in f.streams,
        h in forecast_horizons
    ## The onset grid is ignored by every other stream, so it is passed
    ## unconditionally rather than branching the loop on the stream name.
    _vals = forecast_stream(
        f.chn, stream; horizon = h,
        obs_value = obs_value, n = obs.n, breakpoint = _BREAKPOINT,
        rt_start = f.rt_start, rt_walk_start = f.rt_walk_start,
        onset_grid_start = _onset_hazard_grid_start,
        onset_grid_end = _onset_grid_end,
        onset_alpha_grid_start = _onset_grid_start
    )
    for (d, i) in enumerate(1:stream_thin:length(_vals))
        push!(
            stream_forecasts, (
                obs.cutoff, h, obs.cutoff + Day(h), label,
                d, Float64(_vals[i]), f.fit,
            )
        )
    end
end

## Confirmed/suspect ward-bed occupancy for the joint, read from the joint's
## own `forecast_reported` runs so the ward beds are scored on the same
## footing as the total occupancy in the preferred `stream_forecasts.csv`
## asset. `forecast_stream` cannot project these. The split is not a growable
## stream but a partition of the total.
##
## Dead as it stands. `forecast_reported` projects the total occupancy alone
## and emits neither column, so the guard below always skips and nothing
## reaches `stream_forecasts.csv` under these two labels. The joint does
## carry the cut-off split (`expected_confirmed_incare_T`), so what is
## missing is the partition of the projected level, not the quantity to
## partition it by. Kept, with the guard, so the loop starts writing the
## moment `forecast_reported` grows those columns.
for (h, fc) in forecast_runs,
        (col, label) in (
            (:suspect_occupancy, "isolation beds (suspected)"),
            (:confirmed_occupancy, "treatment beds"),
        )

    col in propertynames(fc) || continue
    _wvals = fc[!, col]
    for (d, i) in enumerate(1:stream_thin:length(_wvals))
        push!(
            stream_forecasts, (
                obs.cutoff, h, obs.cutoff + Day(h), label,
                d, Float64(_wvals[i]), "joint",
            )
        )
    end
end
CSV.write(joinpath(output_dir, "stream_forecasts.csv"), stream_forecasts);

## Latent symptom-onset trajectory over time, the "symptomatic cases" curve,
## showing outbreak growth: one row per grid day with the 30/60/90%
## credible intervals of both the daily new and cumulative onsets.
onsets_over_time_table = onsets_over_time(
    chn_joint;
    n = obs.n, seeding = obs.seeding
)
CSV.write(
    joinpath(output_dir, "onsets_over_time.csv"),
    onsets_over_time_table
);

#md # ```@raw html
#md # </details>
#md # ```

# ### Summary-page assets
#
# The one-page [Summary dashboard](@ref) reuses the results computed above rather than re-fitting.
# Here we save its headline text, headline tables and the figures it shows (the national reproduction number and infections over time) into `docs/src/summary_assets/`.
# The static dashboard page embeds them after this build step has run.

#md # ```@raw html
#md # <details><summary>Write the dashboard assets</summary>
#md # ```

dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)

## Figures: estimated R(t) nationally and latent infections over time, both
## produced in the Results sections above and written here at the dashboard
## size. The province and in-sample pages write the figures they own.
CairoMakie.save(joinpath(dashboard_dir, "rt.png"), rt_fig)
CairoMakie.save(
    joinpath(dashboard_dir, "infections.png"),
    cumulative_traj_fig
)

## Headline prose: the same bullet summary shown at the top of the Results
## section, serialised to markdown so the dashboard renders it verbatim.
open(joinpath(dashboard_dir, "headline.md"), "w") do io
    print(io, sprint(Markdown.plain, summary_ranges))
end

## Headline tables: outbreak size and timing as whole numbers, and the
## growth and severity parameters to two decimals, each with reader-friendly
## quantity names.
dashboard_counts = summary_table(
    chn_joint, [:C_T, :T]; digits = 0,
    labels = Dict(
        :C_T => "Cumulative infections",
        :T => "Outbreak age (days)"
    )
)
dashboard_rates = summary_table(
    chn_joint,
    [:R0, :R_T, :r, :doubling_time, :CFR]; digits = 2,
    labels = Dict(
        :R0 => "Initial reproduction number",
        :R_T => "Latest reproduction number",
        :r => "Latest growth rate (per day)",
        :doubling_time => "Latest doubling time (days)",
        :CFR => "Case-fatality ratio"
    )
)
open(joinpath(dashboard_dir, "headline_counts.md"), "w") do io
    print(io, markdown_table(dashboard_counts))
end
open(joinpath(dashboard_dir, "headline_rates.md"), "w") do io
    print(io, markdown_table(dashboard_rates))
end

## Fit diagnostics: the same table the Results section shows, so the
## dashboard reports how the fit behind its numbers sampled without
## building a second table.
open(joinpath(dashboard_dir, "diagnostics.md"), "w") do io
    print(io, markdown_table(fit_diagnostics_table))
end

## The data cut-off the dashboard reports as of, written as a plain date.
open(joinpath(dashboard_dir, "cutoff.md"), "w") do io
    print(io, string(obs.cutoff))
end

#md # ```@raw html
#md # </details>
#md # ```

# ---
#
# The full analysis code, data and model definitions are in the [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize) repository.
# Issues, corrections and suggestions are welcome there.
# Maintained by Sam Abbott, Kath Sherratt, Samuel Brand and Sebastian Funk.
