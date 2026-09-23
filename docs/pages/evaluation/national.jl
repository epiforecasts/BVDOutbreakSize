# # National evaluation
#
# How well the joint model reproduces the national data it was fitted to, and how its forecasts scored against the data that arrived afterwards.
# The same checks by province are on the [province evaluation](@ref "Province evaluation") page.

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

# ## In-sample checks
#
# Whether the fitted model reproduces the national data it was fitted to.

# ### Prior predictive check
#
# Before any observation is taken into account, what does the prior imply about replicated exports, deaths and reported cases?
# Draws from the prior over the unobserved data should bracket the observed counts.

# The draws come from the shared setup, with every observation withheld, so
# the national page overlays the same ones on each posterior.

#md # ```@raw html
#md # <details><summary>Summarise the joint prior</summary>
#md # ```

prior_C_table = summary_table(prior_chn, [:C_T]; digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Show prior summary table</summary>
#md # ```

prior_C_table #hide

#md # ```@raw html
#md # </details>
#md # ```

# Pair plot of the prior over the latent quantities.

#md # ```@raw html
#md # <details><summary>Prior pair plot</summary>
#md # ```

prior_pair_fig = plot_pair(
    prior_chn,
    [
        :C_T, :R_T, :r, :T, :CFR, :k,
        :p_drc, :p_uganda,
    ]
);

#md # ```@raw html
#md # </details>
#md # ```

prior_pair_fig #hide

# ### Posterior predictive checks
#
# A posterior predictive check draws replicated observations from the fitted joint model and compares them to the observed counts.
# The checks cover two groups: the dated DRC surveillance streams and the Uganda exports.
# The latent infection process is not checked here, as it carries no direct observation, and is shown instead as the estimated cumulative trajectories in the [joint model estimates](@ref "Joint model estimates") figure.
#
# The surveillance group is checked first, split by whether a stream was still being reported at the cut-off.
# A stream counts as still reporting when its last situation-report vintage falls within a week of the cut-off.
# Every panel runs over its own reporting dates with the observed series overlaid, and its date axis is labelled about once a week.

# #### Streams still reporting
#
# ##### Cumulative
#
# A cumulative panel is drawn as replicated cumulative trajectories.
# A daily panel (the isolation-bed occupancy, the 24h analysed volume) is drawn day by day, each day's replicated count against the observed count.

#md # ```@raw html
#md # <details><summary>Joint posterior predictive plot</summary>
#md # ```

## Drop the increment counts but keep each stream's vintage day grid, so
## `predict` resamples the per-vintage increments rather than holding them
## at the observed values. The confirmed-case windows and the per-window
## positivity random effect are defined by the confirmed and laboratory
## histories, so those are passed with their counts intact (only the
## cut-off scalars are set to `missing`) to keep the generator's latent
## dimensions identical to the fitted chain.
_days_only(h) = (; days = h.days, counts = Int[]);

pp_joint = predict(
    bvd_joint(
        obs.n, missing, missing, missing, missing, missing, missing;
        confirmed_deaths = missing,
        recovered_cases = missing,
        deaths_history = _days_only(obs.deaths_history),
        reported_history = _days_only(obs.reported_history),
        suspected_daily_history = _days_only(obs.suspected_daily_history),
        suspected_daily_deaths_history =
            _days_only(obs.suspected_daily_deaths_history),
        isolation_history = _days_only(obs.isolation_history),
        bed_capacity_history = _days_only(obs.bed_capacity_history),
        ## Kept so the generator's occupancy-break dimension matches the fitted
        ## chain (the offset step on the `[occupancy_break_dates]` days).
        occupancy_break_days = obs.occupancy_break_days,
        recovered_history = _days_only(obs.recovered_history),
        treatment_admissions_history =
            _days_only(obs.treatment_admissions_history),
        treatment_deaths_history = _days_only(obs.treatment_deaths_history),
        treatment_ruleout_history = _days_only(obs.treatment_ruleout_history),
        treatment_absconded_history =
            _days_only(obs.treatment_absconded_history),
        treatment_confirmed_incare_history =
            _days_only(obs.treatment_confirmed_incare_history),
        treatment_suspect_incare_history =
            _days_only(obs.treatment_suspect_incare_history),
        confirmed_history = obs.confirmed_history,
        ## Counts kept, like the confirmed cases above: the cut-off scalar
        ## (`confirmed_deaths = missing`) is this stream's generator gate, so
        ## `predict` still resamples the increments while the dated history
        ## supplies both the vintage grid and the published break discrepancy
        ## the step is centred on. Differencing an emptied history cannot
        ## recover that discrepancy, which would leave the harmonised vintage
        ## replicated as a day of real deaths.
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        ## Kept, like the occupancy break above, so the generator's confirmed
        ## break dimension matches the fitted chain (the level step and the
        ## de-anchored positivity denominator on the
        ## `[confirmed_break_dates]` days). Without them the harmonised
        ## vintage is replicated as though its whole increment were one day of
        ## incidence, so 22 July plots as a gross outlier against a chain that
        ## fitted it as mostly backlog, and the `confirmed_step` columns go
        ## unused.
        confirmed_break_days = obs.confirmed_break_days,
        confirmed_break_gross_cases = obs.confirmed_break_gross_cases,
        confirmed_break_gross_deaths = obs.confirmed_break_gross_deaths,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        ## Kept with its real cell grid (`onset_days`/`report_days`/
        ## `prev_report_days`) but `increments = missing`, so `predict`
        ## resamples the reporting-triangle increments over the actual
        ## scored cells rather than the default empty grid.
        onset_curve_history = (;
            onset_days = obs.onset_curve_history.onset_days,
            report_days = obs.onset_curve_history.report_days,
            prev_report_days = obs.onset_curve_history.prev_report_days,
            increments = missing,
        ),
        breakpoint = _BREAKPOINT,
        background_pooling = background_pooling_model,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days,
        ## The generator must be the model that was fitted. `n_patches`
        ## defaults to one, so leaving these out regenerates every stream
        ## from a single well-mixed population while the chain carries a
        ## four-patch fit: the draws still apply, the latent trajectory they
        ## are replayed through does not, and every stream driven by BVD
        ## cases comes out short by the difference. The province grids are
        ## kept with `missing` increments, like the onset triangle above, so
        ## `predict` resamples the compositions over the real cells.
        n_patches = N_PATCHES,
        province_increments = missing,
        province_days = province_cases.days,
        province_testing_covariate = province_testing,
        province_death_increments = missing,
        province_death_days = province_deaths.days
    ),
    chn_joint
);

## `predict` stores each stream's per-vintage increments as one
## vector-valued variable (`<stream>_increments.increments`); the slice is
## an iter×chain matrix of per-draw increment vectors, exactly the
## `replicates` shape `plot_vintage_conditional_ppc` grounds on each
## vintage's observed previous cumulative for the one-step-ahead
## predictive. Look it up by its VarName with FlexiChains' `Prefixed`, which
## matches a (submodel-prefixed) key by its varname tail: `Prefixed(@varname(
## reported_increments.increments))` finds `cases_state.reported_increments.
## increments` without hard-coding the `cases_state.` prefix, and matches by
## the varname tail rather than a loose substring, so it cannot be fooled by a
## scalar `expected_*_T` deterministic. `FlexiChains` is a package
## dependency (imported, not exported), so it is reached through the package
## namespace.
const _Prefixed = BVDOutbreakSize.FlexiChains.Prefixed;
_vintage_replicates(pp, vn) = collect(pp[_Prefixed(vn)]);

## Grid day-index → INSP situation-report date label.
_vintage_dates(days) = string.(obs.seeding .+ Day.(days .- 1));

reported_panel = (;
    id = :suspected_cases,
    title = "Suspected cases",
    dates = _vintage_dates(obs.reported_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(reported_increments.increments)
    ),
    observed = obs.reported_history.counts, colour = :steelblue,
);
## Daily new-suspect inflow: a per-day count (not cumulative), so the panel
## is drawn with `cumulative = false` — each replicate is its own daily
## count against the observed daily count rather than a running total. Its
## days pick up where the cumulative suspected panel freezes.
suspected_daily_panel = (;
    id = :suspected_daily,
    title = "New suspects/day",
    dates = _vintage_dates(obs.suspected_daily_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(suspected_daily.increments)
    ),
    observed = obs.suspected_daily_history.counts,
    colour = :slateblue, cumulative = false,
);
## Isolation/treatment-bed occupancy: a census stock, so the panel is drawn
## with `cumulative = false` — each replicate is the modelled bed count on a
## report day against the observed "Patients en isolement" count. It is a
## level, not a count of new events, so it carries its own `ylabel` rather
## than the "Daily count" and "New per vintage" defaults, which would read as
## an accumulating total on a series that rises through the outbreak. The count
## is the suspect inflow carried through a length-of-stay survival, so its
## level and lag reflect the admission proportion and the stays. The censored-
## occupancy likelihood stores its per-day predictive draws under the submodel
## `obs` variable (not `increments`), so the replicates are read from that key.
## The treatment model scores the total occupancy only on the days without a
## published confirmed/suspect split (a per-day total-or-split switch): on the
## split days the two sub-stock census panels carry the fit instead, so the
## `isolation.obs` predictive holds only the non-split days. Drop the split
## days from the panel's dates and observed counts to match that length.
_iso_split_days = Set(Int.(obs.treatment_confirmed_incare_history.days))
_iso_keep = [!(Int(d) in _iso_split_days) for d in obs.isolation_history.days]
isolation_panel = (;
    id = :isolation_beds,
    title = "Patients in isolation",
    dates = _vintage_dates(obs.isolation_history.days[_iso_keep]),
    replicates = _vintage_replicates(
        pp_joint, @varname(isolation.obs)
    ),
    observed = obs.isolation_history.counts[_iso_keep],
    colour = :darkorange, cumulative = false,
    ylabel = "Beds occupied",
);
deaths_panel = (;
    id = :suspected_deaths,
    title = "Suspected deaths",
    dates = _vintage_dates(obs.deaths_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(death_increments.increments)
    ),
    observed = obs.deaths_history.counts, colour = :firebrick,
);
## Daily new suspected deaths: a per-day count (not cumulative), so the panel
## is drawn with `cumulative = false` — each replicate is its own daily count
## against the observed daily count rather than a running total. Its days
## pick up where the cumulative suspected-death panel freezes, the deaths
## analogue of the new-suspects-per-day panel.
suspected_daily_deaths_panel = (;
    id = :suspected_daily_deaths,
    title = "New suspected deaths/day",
    dates = _vintage_dates(obs.suspected_daily_deaths_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(suspected_daily_deaths.increments)
    ),
    observed = obs.suspected_daily_deaths_history.counts,
    colour = :indianred, cumulative = false,
);
## Specimens analysed is the single modelled laboratory volume (the
## report-to-analysed delay and tested-fraction throughput), fit to the
## cumulative analysed series, so it gets the same cumulative conditional
## check as the suspected streams. This is the testing volume the
## confirmed-positivity denominator is built from.
tests_analysed_panel = (;
    id = :tests_analysed,
    title = "Specimens analysed (cumulative)",
    dates = _vintage_dates(obs.lab_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(analysed_increments.increments)
    ),
    observed = obs.lab_history.counts, colour = :seagreen,
);
## Post-cutoff 24h analysed volume: once the cumulative series stops, INSP
## reports a 24h analysed count on some days. These are fitted as per-day
## volumes (not cumulative), so the panel is a standalone daily check
## (`cumulative = false`): the modelled daily analysed volume against the
## observed 24h count on each reported day.
tests_analysed_daily_panel = (;
    id = :tests_analysed_daily,
    title = "Specimens analysed (24h)",
    dates = _vintage_dates(obs.lab_daily_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(analysed_daily_increments.increments)
    ),
    observed = obs.lab_daily_history.counts, colour = :teal,
    cumulative = false,
);

## Confirmed cases are scored over two groups of laboratory windows: the
## early confirmed vintages (no per-vintage analysed denominator, scored
## as counts against the modelled laboratory volume) and the observed
## windows (a Binomial of the observed analysed denominator). Both groups
## produce per-window replicate increments in `predict`, so concatenating
## them oldest-first gives the per-vintage cumulative confirmed-case
## trajectory, grounded on the observed cumulative confirmed at each window
## end-day. The 24-25 May analysis stall merges into 26 May, so the window
## grid is slightly coarser than the raw confirmed history.
_conf_windows = BVDOutbreakSize.confirmed_positivity_windows(
    obs.confirmed_history, obs.lab_history, obs.lab_daily_history
);
## Oldest-first: early (no denominator) → observed (analysed Binomial) →
## late (post-28 May; trusted 24h-analysed days are Binomial windows, the
## rest unanchored windows scored against the modelled volume).
_conf_window_days = vcat(
    _conf_windows.early_days, _conf_windows.obs_days,
    _conf_windows.late_days
);
function _confirmed_at(day)
    i = searchsortedlast(obs.confirmed_history.days, day)
    return i == 0 ? 0 : Int(obs.confirmed_history.counts[i])
end;
_conf_early = _vintage_replicates(
    pp_joint, @varname(early_increments.increments)
);
_conf_obs = collect(
    first(
        pp_joint[k]
            for k in keys(pp_joint)
            if occursin("confirmed_state.confirmed_positives.positives", string(k))
    )
);
_conf_late = _vintage_replicates(
    pp_joint, @varname(late_increments.increments)
);
confirmed_panel = (;
    id = :confirmed_cases,
    title = "Confirmed cases",
    dates = _vintage_dates(_conf_window_days),
    replicates = [
        vcat(collect(e), collect(p), collect(l))
            for (e, p, l) in zip(vec(_conf_early), vec(_conf_obs), vec(_conf_late))
    ],
    observed = [_confirmed_at(d) for d in _conf_window_days],
    colour = :goldenrod,
);

## Confirmed deaths are a per-vintage stream, scored as increments of the
## modelled confirmed-death trajectory up to the cut-off, so they get the
## same cumulative conditional check.
confirmed_deaths_panel = (;
    id = :confirmed_deaths,
    title = "Confirmed deaths",
    dates = _vintage_dates(obs.confirmed_deaths_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(cdeath_increments.increments)
    ),
    observed = obs.confirmed_deaths_history.counts, colour = :purple,
);

## Recovered among confirmed ("cumul guéris") is a cumulative per-vintage
## stream fitted through the increments of the modelled recovered trajectory
## (the confirmation-to-recovery convolution of the daily confirmed cases) up
## to the cut-off, so it gets the same cumulative conditional check.
recovered_panel = (;
    id = :recovered,
    title = "Recovered (confirmed)",
    dates = _vintage_dates(obs.recovered_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(recovered_increments.increments)
    ),
    observed = obs.recovered_history.counts, colour = :mediumseagreen,
);

## Tableau 6 treatment-centre daily flows (the new patient-movement data
## sources): admissions and the discharge reasons (in-care deaths, rule-outs,
## absconded). Per-day counts, so drawn with `cumulative = false` — each
## replicate is the modelled daily flow on a report day against the observed
## Tableau 6 count.
admissions_panel = (;
    id = :treatment_admissions,
    title = "Admissions/day",
    dates = _vintage_dates(obs.treatment_admissions_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(admissions.obs)
    ),
    observed = obs.treatment_admissions_history.counts,
    colour = :teal, cumulative = false,
);
incare_deaths_panel = (;
    id = :treatment_deaths,
    title = "In-care deaths/day",
    dates = _vintage_dates(obs.treatment_deaths_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(incare_deaths.increments)
    ),
    observed = obs.treatment_deaths_history.counts,
    colour = :darkred, cumulative = false,
);
ruleouts_panel = (;
    id = :treatment_ruleouts,
    title = "Rule-outs/day",
    dates = _vintage_dates(obs.treatment_ruleout_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(ruleouts.increments)
    ),
    observed = obs.treatment_ruleout_history.counts,
    colour = :goldenrod, cumulative = false,
);
absconded_panel = (;
    id = :treatment_absconded,
    title = "Absconded/day",
    dates = _vintage_dates(obs.treatment_absconded_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(absconded.increments)
    ),
    observed = obs.treatment_absconded_history.counts,
    colour = :slategray, cumulative = false,
);

## Tableau 6 occupancy split (`dont confirmes` / `dont suspects`): the two
## in-care prevalence sub-stocks. Per-day census counts, so drawn with
## `cumulative = false` — each replicate is the modelled confirmed-in-care or
## suspect-in-care bed count on a report day against the observed sub-stock.
## On these split days the total-occupancy panel is not scored, so the two
## sub-stock panels carry the window instead.
confirmed_incare_panel = (;
    id = :treatment_beds,
    title = "Confirmed in care",
    dates = _vintage_dates(obs.treatment_confirmed_incare_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(confirmed_incare_obs.increments)
    ),
    observed = obs.treatment_confirmed_incare_history.counts,
    colour = :darkgoldenrod, cumulative = false,
    ylabel = "Beds occupied",
);
suspect_incare_panel = (;
    id = :suspect_beds,
    title = "Suspects in care",
    dates = _vintage_dates(obs.treatment_suspect_incare_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(suspect_incare_obs.increments)
    ),
    observed = obs.treatment_suspect_incare_history.counts,
    colour = :chocolate, cumulative = false,
    ylabel = "Beds occupied",
);

## Symptom-onset reporting triangle: one cell per (onset day, report day)
## pair, several onset dates per snapshot, unlike every panel above (one
## value per report day). To fit the same per-vintage panel shape, sum
## the cells sharing a report day into one net correction per snapshot —
## "one panel showing each snapshot as of its own report date" — rather
## than a full onset-by-report grid (shown separately as the fitted-vs-
## digitised snapshot figure in the Results section above). Per-day net
## correction, not a running total, so `cumulative = false`.
_onset_ppc_report_days = sort(unique(obs.onset_curve_history.report_days))
_onset_ppc_groups = [
    findall(==(r), obs.onset_curve_history.report_days)
        for r in _onset_ppc_report_days
]
_onset_ppc_replicates_raw = _vintage_replicates(
    pp_joint, @varname(onset_report_state.increments)
)
onset_panel = (;
    id = :onset_reports,
    title = "Onset reports (net correction/snapshot)",
    dates = _vintage_dates(_onset_ppc_report_days),
    replicates = [
        [sum(collect(rep)[g]) for g in _onset_ppc_groups]
            for rep in vec(_onset_ppc_replicates_raw)
    ],
    observed = [
        sum(obs.onset_curve_history.increments[g])
            for g in _onset_ppc_groups
    ],
    colour = :mediumpurple, cumulative = false,
);

## Each panel runs to its own last vintage, so a stream that keeps
## reporting shows the full series the model is fitting rather than the
## window the streams that stopped earlier cover. The per-stream
## calibration table reads this ordered list too, so it stays whole and
## the two stream groups are filtered out of it.
vintage_panels = [
    reported_panel, suspected_daily_panel, isolation_panel, confirmed_panel,
    deaths_panel, suspected_daily_deaths_panel, confirmed_deaths_panel,
    recovered_panel, tests_analysed_panel, tests_analysed_daily_panel,
    admissions_panel, incare_deaths_panel, ruleouts_panel, absconded_panel,
    confirmed_incare_panel, suspect_incare_panel, onset_panel,
];
## The incidence view drops the treatment-centre flow and occupancy-split
## panels, whose per-day counts are already their own incidence.
vintage_incidence_panels = [
    reported_panel, suspected_daily_panel, isolation_panel, confirmed_panel,
    deaths_panel, suspected_daily_deaths_panel, confirmed_deaths_panel,
    recovered_panel, tests_analysed_panel, tests_analysed_daily_panel,
    onset_panel,
];
## Whether a panel's stream was still being reported at the cut-off, from
## the shared registry rule (last vintage within a week of the cut-off)
## rather than a per-page list of dates that goes stale.
_still_reporting(p) = stream_reporting(obs, p.id);
reporting_panels = filter(_still_reporting, vintage_panels);
stopped_panels = filter(!_still_reporting, vintage_panels);
reporting_incidence_panels = filter(
    _still_reporting, vintage_incidence_panels
);
stopped_incidence_panels = filter(
    !_still_reporting, vintage_incidence_panels
);
joint_vintage_ppc_fig = plot_vintage_conditional_ppc(reporting_panels);

#md # ```@raw html
#md # </details>
#md # ```

joint_vintage_ppc_fig #hide

# ##### Per-vintage incidence
#
# This is the same check applied to per-vintage incidence: the count reported between consecutive situation reports, rather than the running cumulative.
# Plotting the increment lets a rise or a slowdown in each stream read directly off the height of each step, where the near-straight cumulative line would hide it.
# The replicates are the modelled per-vintage increments, shown as 30/60/90% credible ribbons with the observed increment overlaid.

#md # ```@raw html
#md # <details><summary>Per-vintage incidence posterior predictive plot</summary>
#md # ```

joint_vintage_incidence_fig = plot_vintage_incidence_ppc(
    reporting_incidence_panels
);

#md # ```@raw html
#md # </details>
#md # ```

joint_vintage_incidence_fig #hide

# #### Streams no longer reporting
#
# These streams stopped reporting before the cut-off, so their panels end earlier than the ones above.
#
# ##### Cumulative

#md # ```@raw html
#md # <details><summary>Joint posterior predictive plot</summary>
#md # ```

joint_vintage_ppc_stopped_fig = plot_vintage_conditional_ppc(stopped_panels);

#md # ```@raw html
#md # </details>
#md # ```

joint_vintage_ppc_stopped_fig #hide

# ##### Per-vintage incidence

#md # ```@raw html
#md # <details><summary>Per-vintage incidence posterior predictive plot</summary>
#md # ```

joint_vintage_incidence_stopped_fig = plot_vintage_incidence_ppc(
    stopped_incidence_panels
);

#md # ```@raw html
#md # </details>
#md # ```

joint_vintage_incidence_stopped_fig #hide

# #### Stream calibration
#
# We score each stream's per-vintage conditional predictions against the observed counts.
# `bias` is the mean forecast bias over the vintages (negative under-predicted, positive over-predicted, zero when the observed counts sit at the predictive median).
# `50%/90% coverage` are the fractions of vintages whose observed count falls inside the central 50% and 90% predictive intervals; a well-calibrated stream keeps these near the nominal levels.
# Streams with a large bias or coverage far from nominal are the ones the joint fit reproduces less well.

stream_calibration_table = stream_calibration(vintage_panels);

# The calibration plot's left panel marks each stream's empirical 50% and 90% coverage against dashed reference lines at the nominal levels.
# The right panel marks the mean forecast bias against a dashed line at zero.

#md # ```@raw html
#md # <details><summary>Per-stream calibration plot</summary>
#md # ```

stream_calibration_fig = plot_stream_calibration(stream_calibration_table);

#md # ```@raw html
#md # </details>
#md # ```

stream_calibration_fig #hide

#md # ```@raw html
#md # <details><summary>Per-stream calibration table</summary>
#md # ```

stream_calibration_table #hide

#md # ```@raw html
#md # </details>
#md # ```

# #### Exports
#
# The Uganda export and export-death streams are dated per-day series, each import or death scored as a Poisson at its detection day.
# The scalar posterior predictive sums each replicate's per-day count vector across the dated days, giving the cumulative export and death total to compare with the observed count.

#md # ```@raw html
#md # <details><summary>Scalar posterior predictive plot</summary>
#md # ```

## The dated counts are nested under their submodel prefix as a single
## per-day count vector `<prefix>.counts`; look it up by its VarName with
## `Prefixed` (matching the key by its `<obs>.counts` tail) so the
## deterministic `expected_*_T` quantities cannot be picked up by a loose
## substring, then sum each replicate's per-day vector into the total.
function _dated_total(pp, vn)
    return [sum(v) for v in vec(Array(pp[_Prefixed(vn)]))]
end;

pp_exports = _dated_total(pp_joint, @varname(export_obs.counts));
pp_exports_deaths = _dated_total(
    pp_joint, @varname(death_obs.counts)
);

joint_ppc_fig = plot_posterior_predictive(
    pp_exports, nothing,
    obs.exported_cases, nothing;
    pp_exports_deaths = pp_exports_deaths,
    obs_exports_deaths = obs.exports_deaths
);

#md # ```@raw html
#md # </details>
#md # ```

joint_ppc_fig #hide

# ### Posterior correlations and stream totals
#
# The heatmap is the posterior correlation between each pair of headline quantities: the outbreak size ($C_T$), the reproduction number ($R_T$), the outbreak age ($T$), the case-fatality ratio (CFR), the DRC and Uganda ascertainment fractions ($p_\text{drc}$, $p_\text{ug}$), the non-BVD background rate ($\lambda_\text{bg}$), the fraction tested ($\tau_\text{test}$), and the cut-off total expected for each stream.
# Blue is positive, red negative.

#md # ```@raw html
#md # <details><summary>Posterior correlation heatmap</summary>
#md # ```

correlation_fig = plot_correlation_heatmap(
    chn_joint,
    [
        :C_T, :R_T, :T, :CFR, :p_drc, :p_uganda, :lambda_bg, :tau_test,
        :expected_reports_T, :expected_deaths_T, :expected_confirmed_T,
    ];
    labels = Dict(
        :C_T => raw"C_T", :R_T => raw"R_T", :T => raw"T",
        :CFR => raw"\mathrm{CFR}", :p_drc => raw"p_\mathrm{drc}",
        :p_uganda => raw"p_\mathrm{ug}", :lambda_bg => raw"\lambda_\mathrm{bg}",
        :tau_test => raw"\tau_\mathrm{test}",
        :expected_reports_T => raw"\mathrm{susp.\ cases}",
        :expected_deaths_T => raw"\mathrm{susp.\ deaths}",
        :expected_confirmed_T => raw"\mathrm{conf.\ cases}"
    )
);

#md # ```@raw html
#md # </details>
#md # ```

correlation_fig #hide

# The stream-total plot takes each posterior draw, sums every stream over its own reporting dates, and marks the observed total with a crosshair.
# The diagonal panels are the predictive spread of each total against the observed value.
# The off-diagonal panels show whether the totals move together from draw to draw.

#md # ```@raw html
#md # <details><summary>Stream totals against observed</summary>
#md # ```

## Per-draw modelled total of each stream, summed over its own reporting
## vintages (the confirmed total adds the unscored first-vintage baseline),
## reusing the posterior-predictive replicates built for the vintage panels.
_stream_total(reps) = [sum(Float64.(collect(r))) for r in vec(reps)]
_conf_baseline = isempty(obs.confirmed_history.counts) ? 0 :
    Int(obs.confirmed_history.counts[1])
stream_totals = (;
    suspected_cases = _stream_total(reported_panel.replicates),
    suspected_deaths = _stream_total(deaths_panel.replicates),
    confirmed_cases = _stream_total(confirmed_panel.replicates) .+ _conf_baseline,
    confirmed_deaths = _stream_total(confirmed_deaths_panel.replicates),
    analysed = _stream_total(tests_analysed_panel.replicates),
);
stream_observed = (;
    suspected_cases = Float64(obs.reported_history.counts[end]),
    suspected_deaths = Float64(obs.deaths_history.counts[end]),
    confirmed_cases = Float64(obs.confirmed_cases),
    confirmed_deaths = Float64(obs.confirmed_deaths_history.counts[end]),
    analysed = Float64(obs.lab_history.counts[end]),
);
stream_pairs_fig = plot_stream_pairs(stream_totals, stream_observed);

#md # ```@raw html
#md # </details>
#md # ```

stream_pairs_fig #hide

# ## Forecast evaluation
#
# How the forecasts on the [forecasts](@ref "Forecasts") page have scored against the data that arrived afterwards.
# Scoring is the continuous ranked probability score against a persistence baseline, defined in the [forecast scoring](@ref "Forecast scoring against a persistence baseline") Methods section.

# ### Forecast validation
#
# How last week's forecast held up against the data since observed, using the frozen re-fit and one-week projection defined in [forecast-versus-frozen evaluation](@ref "Forecast-versus-frozen evaluation").
# Only the streams the situation reports are still updating are validated here.
# A stream that has stopped being reported carries a cumulative total that repeats its last reported value, so there is no observation for the past week to score against.
# The frozen fit also conditions on the isolation beds, so the projected bed occupancy is scored against the beds held a week later.
# The bed validation is weak at a one-week-back freeze.
# The reported occupancy rate starts only on 9 June, so the capacity has no implied-capacity anchor and rides its random walk back to the freeze date.
# Like the scores further down, the confirmed new-count rows here take out any retrospective harmonisation step the week contained.
# Such a step reattaches records notified earlier, so it is not something the forecast was predicting.
# The cumulative rows are scored against the published total, harmonisation included.

#md # ```@raw html
#md # <details><summary>Fit one week back and validate the one-week-ahead forecast</summary>
#md # ```

## frozen_lastweek and frozen_lastweek_streams are computed in the setup
## block above, and `validation_forecast_from` is defined there.
validation_forecast = validation_forecast_from(frozen_lastweek);

## Each frozen individual (single-stream) fit's own one-week-ahead new-count
## forecast at the same cut-off as `frozen_lastweek`, from
## [`forecast_stream`](@ref) (the same per-stream forecaster
## `stream_forecasts.csv` uses), so the validation plots below can show the
## individual fit alongside the joint rather than the joint alone. Recovered
## has no individual fit and is absent here, as it is throughout this report.
## Only the still-reported streams are fitted at the validation cut-off, so
## a stream the situation reports have stopped updating is absent from
## `frozen_lastweek_streams` and carries no individual series here.
function _validation_individual_new(sid, stream::Symbol, obs_field)
    haskey(frozen_lastweek_streams, sid) || return nothing
    f = frozen_lastweek_streams[sid]
    bp = f.o.n - f.o.who_first_sitrep_days
    return Float64.(
        forecast_stream(
            f.chn, stream; horizon = 7,
            obs_value = getproperty(f.o, obs_field), n = f.o.n, breakpoint = bp,
            rt_start = 1, rt_walk_start = 1
        )
    )
end
validation_individual = NamedTuple(
    k => v
        for (k, v) in pairs(
            (;
                cases_new = _validation_individual_new(
                    "cases", :reported_cases, :reported_cases
                ),
                deaths_new = _validation_individual_new(
                    "deaths", :suspected_deaths, :total_deaths
                ),
                confirmed_new = _validation_individual_new(
                    "confirmed", :confirmed_cases, :confirmed_cases
                ),
                confirmed_deaths_new = _validation_individual_new(
                    "confirmed_deaths", :confirmed_deaths, :confirmed_deaths
                ),
            )
        )
        if !isnothing(v)
)
## The frozen individual (treatment-only) fit's own bed-occupancy forecast,
## anchored on the beds occupied at ITS OWN cut-off (the frozen fit's own
## `o`, not the current `obs`), matching how the joint frozen forecast is
## itself anchored.
## `nothing` when the beds have stopped being reported, so the treatment fit
## is absent; the bed panel then draws the joint alone.
## A `let` block, not a bare `if`: a top-level `if` shares the script's
## global scope, so its working names would leak into the rest of the page.
validation_individual_isolation = let
    if haskey(frozen_lastweek_streams, "treatment")
        tf = frozen_lastweek_streams["treatment"]
        beds = isempty(tf.o.isolation_history.counts) ? 0.0 :
            Float64(tf.o.isolation_history.counts[end])
        Float64.(
            forecast_stream(
                tf.chn, :isolation_beds; horizon = 7,
                obs_value = beds, n = tf.o.n,
                breakpoint = tf.o.n - tf.o.who_first_sitrep_days,
                rt_start = 1, rt_walk_start = 1
            )
        )
    else
        nothing
    end
end

## The observed beds at the current cut-off (the forecast target), so the
## frozen-fit bed forecast is scored against what the beds actually held.
## Held back once the beds stop being reported, since the last count would
## then be carried forward rather than observed at the target date.
_obs_beds = stream_reporting(obs, :isolation_beds) ?
    obs.isolation_history.counts[end] : missing
## Same observed/baseline keying as the plot below, so the table covers every
## fitted count stream (cumulative and new-count rows) plus the bed level.
## A harmonisation-break day between the frozen cut-off and the current one
## puts records into the confirmed cumulative that were never notified in that
## week, so the new-count truth carries a step the forecast was never
## predicting. Take it out, the same correction `score_releases.jl` applies.
## Grid days are relative to a seeding date fixed by the genetic tmrca, so the
## frozen fit's own `n` and the current `obs.n` index the same grid.
validation_breaks = (
    confirmed_cum = confirmed_break_correction(
        obs, frozen_lastweek.o.n, obs.n
    ),
    confirmed_deaths_cum = confirmed_break_correction(
        obs, frozen_lastweek.o.n, obs.n; deaths = true
    ),
)

## Observed cumulative at the target date per stream, keyed by the forecast's
## cumulative column; `baseline` is each stream's origin cumulative (the
## frozen cut-off), so the new count is scored against observed minus origin,
## less any harmonisation the window carries (see `validation_breaks`). Both
## the table and the plot below take the still-reported streams
## (`reporting_cum_cols`, from the setup block): a stream the situation
## reports have stopped updating has an origin and a target reading the same
## repeated total, so its cumulative truth is stale and its new-count truth is
## a guaranteed zero.
validation_observed = (
    cases_cum = obs.reported_cases,
    deaths_cum = obs.total_deaths,
    confirmed_cum = obs.confirmed_cases,
    confirmed_deaths_cum = obs.confirmed_deaths,
    recovered_cum = obs.recovered_cases,
)
validation_baseline = (
    cases_cum = frozen_lastweek.o.reported_cases,
    deaths_cum = frozen_lastweek.o.total_deaths,
    confirmed_cum = frozen_lastweek.o.confirmed_cases,
    confirmed_deaths_cum = frozen_lastweek.o.confirmed_deaths,
    recovered_cum = frozen_lastweek.o.recovered_cases,
)

validation_table = forecast_vs_truth(
    validation_forecast;
    observed = keep_streams(validation_observed, reporting_cum_cols),
    baseline = keep_streams(validation_baseline, reporting_cum_cols),
    breaks = validation_breaks,
    isolation = _obs_beds
);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Forecast-versus-observed validation table</summary>
#md # ```

## `MarkdownTable` rather than a bare table expression: a DataFrame is #src
## `text/html`-showable, Literate prefers that mime, and the `@raw html` #src
## block it writes crosses Documenter's raw-block regex limit once the #src
## table grows. `MarkdownTable` is markdown-showable and not #src
## html-showable, so the table goes out as an ordinary markdown table #src
## rather than a fixed-width block of printed output. See its docstring #src
## for the mechanism. The same treatment is applied to every DataFrame #src
## display in this file and in `analysis.jl`. #src
MarkdownTable(validation_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The observation panels histogram the one-week-ahead forecast made from the frozen fit: a cumulative and a new-count panel for each still-reported count stream the forecast carries.
# The 90% predictive interval is shaded, and the count observed by the current cut-off is a dashed black rule.
# Where a stream has its own individual (single-stream) fit, that fit's forecast from the same frozen cut-off is overlaid as a dotted step outline on the joint's own histogram bins.

#md # ```@raw html
#md # <details><summary>Forecast-versus-observed plot</summary>
#md # ```

validation_fig = plot_forecast_vs_truth(
    validation_forecast;
    observed = keep_streams(validation_observed, reporting_cum_cols),
    baseline = keep_streams(validation_baseline, reporting_cum_cols),
    breaks = validation_breaks,
    individual = keep_streams(validation_individual, reporting_cum_cols)
);

#md # ```@raw html
#md # </details>
#md # ```

validation_fig #hide

# The bed panel scores last week's projected occupancy against the beds occupied now (the dashed rule), with the individual (treatment-only) fit's own projection overlaid as a dotted step outline on the joint's own histogram bins.

#md # ```@raw html
#md # <details><summary>Bed forecast-versus-observed plot</summary>
#md # ```

validation_beds_fig = plot_forecast_beds_vs_truth(
    validation_forecast;
    isolation = _obs_beds, individual = validation_individual_isolation
);

#md # ```@raw html
#md # </details>
#md # ```

validation_beds_fig #hide

# The latent quantities are not observed, so they are scored distribution against distribution: what the frozen fit forecast for the past week's new infections, onsets and deaths against what the current fit now estimates for the same window.

#md # ```@raw html
#md # <details><summary>Forecast-versus-now latent plot</summary>
#md # ```

## Current fit's draws of the new latent counts over the past week, the last
## seven days of each cumulative-trajectory deterministic.
function _now_new(chn, key)
    mat = chn[key]
    trajs = [collect(v) for v in vec(collect(mat))]
    return Float64[t[end] - t[max(1, length(t) - 7)] for t in trajs]
end
now_latent = (;
    infections_new = _now_new(chn_joint, :cumulative_infections),
    onsets_new = _now_new(chn_joint, :cumulative_onsets),
    deaths_latent_new = _now_new(chn_joint, :cumulative_expected_deaths),
)

validation_latent_fig = plot_forecast_vs_truth_latent(
    validation_forecast; now = now_latent
);

#md # ```@raw html
#md # </details>
#md # ```

validation_latent_fig #hide

# #### Streams no longer reported
#
# The situation reports have stopped updating some of the streams the model fits, listed with the date each was last reported below.
# The panels show what the frozen fit projected for those streams over the same week, without an observed rule, since the count they would be scored against has not moved since the stream stopped.

#md # ```@raw html
#md # <details><summary>Forecast for the streams no longer reported</summary>
#md # ```

## The last-reported date per stopped stream, and the frozen fit's own
## projection for them. `plot_forecast` draws a panel per new-count column
## the frame carries, so passing the stopped streams' columns alone gives the
## projection without the fabricated truth rule the validation figure would
## otherwise draw against a repeated total.
validation_stopped_streams = let s = stream_report_status(obs),
        ids = [stream_id(c) for c in stopped_cum_cols]

    keep = [r.stream in ids for r in eachrow(s)]
    DataFrame(
        "Stream" => s[keep, :label],
        "Last reported" => s[keep, :last_date]
    )
end
_stopped_new_cols = [
    c
        for c in new_cols(stopped_cum_cols)
        if c in propertynames(validation_forecast)
]
validation_stopped_fig = plot_forecast(
    validation_forecast[!, _stopped_new_cols]
);

#md # ```@raw html
#md # </details>
#md # ```

## See the comment above `validation_table`'s display for why this wraps #src
## the table in `MarkdownTable` instead of showing it directly. #src
## No blank line between the two. A blank line counts as visible, so #src
## Literate would write an empty code fence where the comment was. #src
MarkdownTable(validation_stopped_streams) #hide
validation_stopped_fig #hide

# ### Forecast scoring across releases
#
# Every release's saved one- to four-week-ahead forecast is scored against the data observed since, against a persistence baseline and, where one exists, the stream's own individual fit as well as the joint.
# The tables in this section are the joint model's, one row per stream.
# See [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline") for how the scores, the relative skill and the baseline are built.
# Recovered has no individual fit of its own, so its comparison is the baseline against the joint only.
# Reported cases and suspected deaths stopped being updated by the situation reports partway through the outbreak, and exports' confirmed-detection series is anchored to an earlier cut-off.
# Exports therefore contributes no scored forecast, and reported cases and suspected deaths each rest on exactly one matched forecast, a single window rather than a settled sample.
#
# Only a minority of the daily releases examined contribute a row to the table below, each a reconstruction of an earlier model version rather than the current fit.
# Only the newest few releases carry the current model's own individual-stream forecasts, and the backfilled reconstructions carry none at all.
# Every row also rests on one to a handful of matched forecasts, shown as its own count rather than rounded away.
#
# Four things are excluded from the scores here and in the frozen section below, each for a stated reason rather than for scoring badly, and the second of them applies to the frozen section alone.
#
# - One whole reconstruction (`results-v1.6.0`): its chain forecasts a near-zero median at every horizon and stream, with the upper predictive tail occasionally reaching five- and six-digit values, which is the signature of a chain that failed to sample rather than a forecast.
# - Frozen section only: the confirmed-death rows of the fourteen frozen reconstructions cut between 16 July and 15 August 2026, whose forecaster could not project that stream from its own trajectory and floored it at zero, so each carries a one-week median of exactly zero against an observed 250 to 370. Reconstructions cut after that window project the stream normally.
# - An onset window containing a vintage whose reread total falls, since its increment is not what the situation reports added.
# - A stream that carries no persistence baseline, which is what makes a window scoreable at all.
#
# Nothing is dropped from the archive itself; `data/forecast_scores*.csv` and `data/forecast_overlay*.csv` record everything that was scored.
#
# The symptom-onset stream is scored on the new reported count each vintage adds rather than on its level, because every vintage rereads the whole figure.
# Its printed total therefore moves with the scan error as well as with late reporting.
# On fourteen vintages the reread total falls, which a cumulative onset curve cannot do, and on many others it repeats unchanged.
# The fit absorbs that with a per-vintage scan level; the scored truth cannot, since it is the increment between the vintages at the two ends of a window.
# A window containing a falling vintage is therefore left unscored, the rule the [province scores](@ref "Forecast by province across releases") already apply to a window holding a harmonisation-break day.
# It bites hardest at the longer horizons, a four-week window being more likely to contain a reread than a one-week one: the frozen onset row keeps three of its twenty-nine windows, all at one week, and the cross-release row six of thirty-eight.
# Read the onset row's skill against the baseline rather than its coverage, and read it as resting on a handful of windows.

#md # ```@raw html
#md # <details><summary>Load and summarise the cross-release forecast scores</summary>
#md # ```

forecast_scores_df = _release_data(
    "forecast_scores.csv",
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, crps = Float64,
        log_crps = Float64, dispersion = Float64, overprediction = Float64,
        underprediction = Float64, coverage_50 = Float64,
        coverage_90 = Float64,
        bias = Float64, n_samples = Int,
        log_rel_to_baseline = Float64,
    )
)
forecast_overlay_df = _release_data(
    "forecast_overlay.csv",
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, observed = Float64,
        median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
        hi60 = Float64, lo90 = Float64, hi90 = Float64,
    )
)
## The digitised onset triangle's own per-vintage total, as calendar dates,
## and the windows it makes unscoreable. A vintage that rereads the figure
## lower gives a window an increment the reports did not add, so the window
## is dropped from the scores and the overlay alike (see
## `drop_rescanned_onset_windows`). The current triangle is used to judge
## every release, since the scored truth is read off this one series.
_onset_vintage_dates = grid_date.(obs.onset_report_history.days)
_onset_vintage_totals = obs.onset_report_history.counts
_drop_rescanned(tbl) = drop_rescanned_onset_windows(
    tbl; vintage_dates = _onset_vintage_dates,
    vintage_totals = _onset_vintage_totals
)
forecast_scores_df = _drop_rescanned(forecast_scores_df)
forecast_overlay_df = _drop_rescanned(forecast_overlay_df)

## One row per (stream, fit) pooled over every horizon and release. The
## by-horizon and by-release detail tables carry the same columns at a finer
## grain (see src/scoring.jl). Every fit is kept here, since the
## relative-skill figure below compares the roles against each other. The
## tables rendered in this section select the joint role, and the individual
## fits are tabulated in their own section.
forecast_score_overview_table = forecast_score_overview(forecast_scores_df)
forecast_score_by_horizon_table = forecast_score_by_horizon(forecast_scores_df)
forecast_score_by_release_table = forecast_score_by_release(forecast_scores_df)

joint_score_overview_table = select_fit_role(
    forecast_score_overview_table, "joint"
)
joint_score_by_horizon_table = select_fit_role(
    forecast_score_by_horizon_table, "joint"
)
## The trailing `;` on this last assignment matters: without it, this whole
## setup chunk's last statement (the DataFrame it assigns) is Literate's
## implicitly displayed "result" for the chunk, on top of the deliberate
## display further down -- and a bare DataFrame is html-showable, so it
## goes out as a second, undisplayed-in-source `@raw html` block that (for
## a table this size) can itself hit the PCRE limit described above.
joint_score_by_release_table = select_fit_role(
    forecast_score_by_release_table, "joint"
);

#md # ```@raw html
#md # </details>
#md # ```

# The headline pools every horizon and release into one row per stream for the joint model: the mean CRPS and its decomposition, coverage, bias, and the relative skill against the persistence baseline, on both the natural and the log scale.
# Each row also carries relative skill against the stream's own individual fit where one exists.

MarkdownTable(joint_score_overview_table) #hide

# The same relative skill against the baseline, by horizon: one panel per stream, one series per fit role, on a log-scaled skill axis with the reference line at one.

forecast_relative_skill_fig = plot_forecast_relative_skill(
    forecast_score_by_horizon_table
);

forecast_relative_skill_fig #hide

# What that error is made of, by horizon: the mean CRPS split into its width, its overprediction and its underprediction, one stacked bar per horizon and fit role.

forecast_crps_by_horizon_fig = plot_forecast_crps_by_horizon(
    forecast_score_by_horizon_table
);

forecast_crps_by_horizon_fig #hide

#md # ```@raw html
#md # <details><summary>Scores by horizon</summary>
#md # ```

MarkdownTable(joint_score_by_horizon_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The same relative skill against the baseline, release by release, so a run of releases that lost to the baseline reads as a run rather than as an average.

forecast_skill_by_cutoff_fig = plot_forecast_skill_by_cutoff(
    forecast_score_by_release_table;
    title = "Relative skill against the baseline, by release"
);

forecast_skill_by_cutoff_fig #hide

#md # ```@raw html
#md # <details><summary>Scores by release</summary>
#md # ```

MarkdownTable(joint_score_by_release_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# Forecasts made at each release against the value observed since, one panel per stream and horizon, the observed value in black.
# The median and 90% interval are coloured by fit role: the persistence baseline, the stream's individual fit and the joint.
# The x-axis is the date each forecast was made, so an incident stream's observed window pairs unambiguously with the forecast that made it.
# Each panel's axis is cropped to a small multiple of what that stream actually reached, so one very wide interval cannot squash every other series flat.
# An interval or median too wide for the panel is clamped at the top and marked with an open triangle rather than silently cut off.

#md # ```@raw html
#md # <details><summary>Forecasts-versus-now overlay</summary>
#md # ```

forecast_overlay_fig = plot_forecast_overlay(
    scored_overlay(forecast_overlay_df)
);

#md # ```@raw html
#md # </details>
#md # ```

forecast_overlay_fig #hide

# ### Frozen-fit forecast evaluation
#
# The current model, frozen at earlier data cut-offs (see [Forecast-versus-frozen evaluation](@ref "Forecast-versus-frozen evaluation")), is scored the same way as the cross-release forecasts above, against the same persistence baseline.
# The tables in this section are the frozen joint model's, one row per stream.
# Each stream's own frozen fit is also scored where one exists, at the one-week-back cut-off and for still-reported streams only, and is carried by the skill figures rather than by the tables.
# Every other cut-off carries the joint alone.
# The May cut-offs predate the first reported bed occupancy and the first reported recoveries, so those windows are left unscored rather than scored against a series that had not started.
# The baseline carries a weaker data-vintage guarantee than the cross-release one, since its snapshot was taken weeks after the frozen cut-off and can hold later revisions to earlier days (see [forecast scoring against a persistence baseline](@ref "Forecast scoring against a persistence baseline")).

#md # ```@raw html
#md # <details><summary>Load and summarise the frozen-fit forecast scores</summary>
#md # ```

frozen_scores_df = _release_data(
    "forecast_scores_frozen.csv",
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, crps = Float64,
        log_crps = Float64, dispersion = Float64, overprediction = Float64,
        underprediction = Float64, coverage_50 = Float64,
        coverage_90 = Float64,
        bias = Float64, n_samples = Int,
        log_rel_to_baseline = Float64,
    )
)
frozen_overlay_df = _release_data(
    "forecast_overlay_frozen.csv",
    (;
        release = String, made_date = Date, stream = String, horizon = Int,
        target_date = Date, fit = String, observed = Float64,
        median = Float64, lo30 = Float64, hi30 = Float64, lo60 = Float64,
        hi60 = Float64, lo90 = Float64, hi90 = Float64,
    )
)
## Rows a superseded forecaster produced are dropped before anything is
## summarised or drawn, from the scores and the overlay alike, so the
## tables and the figures rest on one set of rows. See
## `drop_superseded_forecasts` for the one exclusion in force and why.
frozen_scores_df = _drop_rescanned(
    drop_superseded_forecasts(frozen_scores_df)
)
frozen_overlay_df = _drop_rescanned(
    drop_superseded_forecasts(frozen_overlay_df)
)

## The frozen joint carries `FROZEN_FIT`, so it is named as the joint role
## here and compared against each stream's own frozen fit. A release
## published before the archive had a `fit` column carries joint rows only.
frozen_score_overview_table = forecast_score_overview(
    frozen_scores_df; joint_fit = FROZEN_FIT
)
frozen_score_by_horizon_table = forecast_score_by_horizon(
    frozen_scores_df; joint_fit = FROZEN_FIT
)
frozen_score_by_release_table = forecast_score_by_release(
    frozen_scores_df; joint_fit = FROZEN_FIT
)

## The tables show the frozen joint alone, as the cross-release tables show
## the joint alone, so a row reads as one model at one cut-off rather than a
## stream interleaving two fits. The archive's single-stream frozen fits stay
## in the scored data and in the figures, which compare the roles against
## each other. Selecting the joint role leaves `fit` single-valued, so it is
## dropped and the model named in the prose instead.
_frozen_joint_only(tbl) = drop_degenerate_fit_column(
    select_fit_role(tbl, "joint")
)
frozen_score_overview_display = _frozen_joint_only(
    frozen_score_overview_table
)
frozen_score_by_horizon_display = _frozen_joint_only(
    frozen_score_by_horizon_table
)
frozen_score_by_release_display = _frozen_joint_only(
    frozen_score_by_release_table
)

## One row per release for the cut-offs more than one release forecast.
## See the comment above `joint_score_by_release_table`'s assignment for why
## this setup chunk's last statement needs a trailing `;`.
frozen_score_by_vintage_table = forecast_score_by_vintage(
    frozen_scores_df; joint_fit = FROZEN_FIT
)
frozen_score_by_vintage_display = _frozen_joint_only(
    frozen_score_by_vintage_table
);

#md # ```@raw html
#md # </details>
#md # ```

MarkdownTable(frozen_score_overview_display) #hide

# The same relative skill against the baseline, by horizon, for the frozen cut-offs.

frozen_relative_skill_fig = plot_forecast_relative_skill(
    frozen_score_by_horizon_table
);

frozen_relative_skill_fig #hide

# What that error is made of, by horizon, as in the cross-release section above.

frozen_crps_by_horizon_fig = plot_forecast_crps_by_horizon(
    frozen_score_by_horizon_table;
    title = "CRPS decomposition by horizon, frozen cut-offs"
);

frozen_crps_by_horizon_fig #hide

#md # ```@raw html
#md # <details><summary>Scores by horizon</summary>
#md # ```

MarkdownTable(frozen_score_by_horizon_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The same relative skill against the baseline, cut-off by cut-off, pooled over the horizons each cut-off forecast.

frozen_skill_by_cutoff_fig = plot_forecast_skill_by_cutoff(
    frozen_score_by_release_table;
    xlabel = "Frozen cut-off",
    title = "Relative skill against the baseline, by frozen cut-off"
);

frozen_skill_by_cutoff_fig #hide

#md # ```@raw html
#md # <details><summary>Scores by frozen cut-off</summary>
#md # ```

MarkdownTable(frozen_score_by_release_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# #### Frozen skill by release
#
# Skill at each cut-off more than one release forecast, one point per release rather than pooled across releases.
# Releases run in the order they were cut, evenly spaced rather than to calendar scale.

#md # ```@raw html
#md # <details><summary>Frozen skill per release</summary>
#md # ```

frozen_skill_by_vintage_fig = plot_forecast_skill_by_vintage(
    frozen_score_by_vintage_table
);

#md # ```@raw html
#md # </details>
#md # ```

frozen_skill_by_vintage_fig #hide

#md # ```@raw html
#md # <details><summary>Frozen scores by release</summary>
#md # ```

MarkdownTable(frozen_score_by_vintage_display) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The frozen forecasts made at each cut-off against the value observed since, one panel per stream and horizon, the observed value in black.
# Each panel carries the frozen forecast and the persistence baseline, coloured as in the cross-release overlay above, and the x-axis is the cut-off each forecast was made from.

#md # ```@raw html
#md # <details><summary>Frozen-fit forecasts-versus-now overlay</summary>
#md # ```

## The frozen joint and the persistence baseline only. A single-stream
## frozen fit exists at the one-week-back cut-off alone, so its series lands
## on one made date of a panel spanning every cut-off, overplotting the
## joint point it sits beside rather than reading as a second series. The
## single-stream frozen fits are compared against the joint in the skill
## figures and in the validation plot at that cut-off.
frozen_overlay_fig = plot_forecast_overlay(
    scored_overlay(
        vcat(
            select_fit_role(frozen_overlay_df, "joint"),
            select_fit_role(frozen_overlay_df, "baseline")
        )
    )
);

#md # ```@raw html
#md # </details>
#md # ```

frozen_overlay_fig #hide

# The frozen re-fits below freeze the renewal data to an earlier cut-off and re-fit, so that a change driven by newer data can be distinguished from one driven by a change of method.
# Each uses the full headline settings: 1000 draws across two chains.

#md # ```@raw html
#md # <details><summary>Freeze the renewal data to a cut-off and re-fit</summary>
#md # ```

## Frozen re-fits and released_df are prepared in the setup block above.

#md # ```@raw html
#md # </details>
#md # ```

# ### Individual fits against the baseline
#
# This section carries the same cross-release forecast scoring as [Forecast scoring across releases](@ref "Forecast scoring across releases") above, for each stream's own individual fit rather than the joint, against the same persistence baseline.

#md # ```@raw html
#md # <details><summary>Individual-fit rows of the cross-release scores</summary>
#md # ```

## The relative skill against a stream's individual fit is only ever
## computed on the joint model's row, so on these rows it is missing by
## construction and the column is dropped rather than shown empty.
individual_score_overview_table = drop_individual_fit_columns(
    select_fit_role(forecast_score_overview_table, "individual")
)
individual_score_by_horizon_table = drop_individual_fit_columns(
    select_fit_role(forecast_score_by_horizon_table, "individual")
)
## See the comment above `joint_score_by_release_table`'s assignment for why
## this setup chunk's last statement needs a trailing `;`.
individual_score_by_release_table = drop_individual_fit_columns(
    select_fit_role(forecast_score_by_release_table, "individual")
);

#md # ```@raw html
#md # </details>
#md # ```

MarkdownTable(individual_score_overview_table) #hide

# The same relative skill against the baseline, by horizon, one panel per stream (dataset), for each stream's own individual fit.

individual_relative_skill_fig = plot_forecast_relative_skill(
    individual_score_by_horizon_table;
    empty_message = "Empty: no release old enough for its targets to " *
        "have been observed carries an individual-stream " *
        "forecast. Not a missing forecast."
);

individual_relative_skill_fig #hide

#md # ```@raw html
#md # <details><summary>Scores by horizon</summary>
#md # ```

MarkdownTable(individual_score_by_horizon_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Scores by release</summary>
#md # ```

MarkdownTable(individual_score_by_release_table) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ## Saving evaluation outputs
#
# The in-sample dashboard asset, and the one-week-back validation forecast in the same archive format as the release forecast, so the frozen "last week versus now" forecast is recorded as a release asset alongside the forecast it is scored against.

#md # ```@raw html
#md # <details><summary>Write the in-sample dashboard asset</summary>
#md # ```

dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
## The report splits the surveillance panels by whether the stream was still
## reporting at the cut-off. The dashboard shows one grid, so it is drawn here
## over the full ordered panel set.
CairoMakie.save(
    joinpath(dashboard_dir, "reported_cases.png"),
    plot_vintage_conditional_ppc(vintage_panels)
)

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Write forecast outputs</summary>
#md # ```

output_dir = get(
    ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output")
)
mkpath(output_dir)
CSV.write(
    joinpath(output_dir, "forecast_validation.csv"),
    forecast_archive(
        [(7, validation_forecast)];
        made_date = frozen_lastweek.o.cutoff, thin = 5
    )
)

#md # ```@raw html
#md # </details>
#md # ```
