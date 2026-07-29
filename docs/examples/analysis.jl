#md # ```@eval
#md # using BVDOutbreakSize, Markdown, Dates
#md # readme = read(joinpath(pkgdir(BVDOutbreakSize), "README.md"), String)
#md # body = strip(match(r"^(.*?)<!-- SHARED:END -->"s, readme).captures[1])
#md # # Keep the dates current automatically: "Last updated" is the build
#md # # date and "Data as of" is the loaded cut-off, so a rebuild always
#md # # refreshes them without a manual edit to README.md.
#md # built = Dates.format(Dates.today(), "d U yyyy")
#md # asof = Dates.format(load_observations().cutoff, "d U yyyy")
#md # body = replace(body,
#md #     r"\*\*Last updated:\*\* [^.]*\." => "**Last updated:** $built.",
#md #     r"\*\*Data as of:\*\* [^.]*\." => "**Data as of:** $asof.",
#md #     r"https://epiforecasts\.io/BVDOutbreakSize/stable/analysis" => "",
#md #     "https://epiforecasts.io/BVDOutbreakSize/stable/contributing" => "contributing.md")
#md # Markdown.parse(body)
#md # ```
#
# This page is generated from
# [`docs/examples/analysis.jl`](https://github.com/epiforecasts/BVDOutbreakSize/blob/main/docs/examples/analysis.jl);
# the model code it calls is in
# [`src/`](https://github.com/epiforecasts/BVDOutbreakSize/tree/main/src).
# See the *LLM-driven reimplementation* limitation below for the
# oversight context behind the Use of AI note.
#
# **Offline copy.** A self-contained single-file HTML version of this
# report, built from the same run, is attached to each results release:
# [download the latest](https://github.com/epiforecasts/BVDOutbreakSize/releases/latest/download/analysis.html).
#
# ## Origins of this work
#
# This work began as a replication of the McCabe et al. [mccabe2026](@cite)
# report.
# It has since evolved into a real-time joint Bayesian estimate of the
# current outbreak size, a discrete-time renewal process with a time-varying
# reproduction number fitted to more of the available data streams than the
# original.
# The points below summarise how it now differs from the report; the Methods
# section carries the full treatment, and the later
# [comparison with McCabe et al.](@ref "Comparison with McCabe et al.") sets
# the current estimates against theirs.
#
#md # ```@raw html
#md # <details><summary>Expand: differences from the report</summary>
#md # ```
#md #
# **Latent process and parameters**
#
# - *Discrete-time renewal model.* The whole model runs on a daily grid.
#   Infections follow the discrete renewal equation $I_t = R_t \sum_{s
#   \ge 1} I_{t-s} g_s$, where $g$ is the discretised generation-interval
#   PMF, and every delay is applied as a discrete convolution. McCabe et
#   al. [mccabe2026](@cite) use continuous-time closed forms.
# - *Time-varying reproduction number.* $R_t$ is held flat at the
#   established $R_0$ until the first WHO situation report (18 May 2026),
#   then follows a weekly Gaussian random walk on the log scale,
#   interpolated within weeks, with a logistic outbreak-response ramp of
#   about three weeks from that report. McCabe et al. use one constant
#   exponential growth rate.
# - *Joint posterior rather than scenario estimates.* The reproduction
#   number, case-fatality ratio, all delays, traveller volume and
#   surveillance dispersion have priors and are sampled together. McCabe
#   et al. [mccabe2026](@cite) fix each and report a set of scenarios.
# - *Two-phase seeding with a wide, genetically-floored outbreak age.*
#   A single import grows through an unobserved cryptic exponential phase
#   to a magnitude set by a wide prior on the doubling count, at the growth
#   rate the genetic estimate informs, before the renewal process takes
#   over. The established reproduction number is derived forward from that
#   growth rate. The genetic time to the most recent common ancestor floors
#   the cryptic duration from below. McCabe et al. fix the start from a
#   single seed.
#
# **Delays and convolutions**
#
# - *Delays re-estimated with uncertainty.* McCabe et al.
#   [mccabe2026](@cite) take the onset-to-death delay from the Isiro 2012
#   point estimate of Rosello et al. [rosello2015](@cite). We instead use
#   a Bayesian reanalysis of the same line list
#   [bdbv_linelist_analysis_2026](@cite) that re-estimates the delay with
#   uncertainty, and we sample every other delay (generation interval,
#   incubation period, onset-to-report, onset-to-confirmation and
#   onset-to-detection abroad) from a prior centred on published Ebola
#   estimates, discretised with double interval censoring
#   [charniga2024](@cite), so the delay uncertainty propagates.
#
# **Likelihoods and data streams**
#
# - *More streams fitted.* McCabe et al. [mccabe2026](@cite) fit the
#   Uganda export cases and deaths. We add the DRC suspected cases, the
#   laboratory-confirmed cases, the confirmed deaths and the deaths among
#   the Uganda exports.
# - *Per-vintage time-series fitting.* The DRC streams are fitted on the
#   incidence scale, as the between-vintage increments across successive
#   sitreps (the first vintage being the cumulative count to that date),
#   which sharpens $R_t$. McCabe et al. condition on a single cumulative
#   total.
# - *Ascertainment estimated.* We jointly estimate the outbreak size and the
#   fraction of cases each surveillance system reports. McCabe et
#   al. have no ascertainment component.
# - *Comparison against published scenarios.* The model is set beside the
#   McCabe et al. [mccabe2026](@cite) scenario estimates as an external
#   sense-check, matched in time at the cut-off each scenario was computed,
#   while the cumulative infection count, the running sum of the daily
#   infections, is the headline quantity reported separately. A forward
#   projection from a frozen fit is also set against the Chamla et al.
#   [chamla2026](@cite) confirmed-case projection and the data observed since.
#
# **Extensions**
#
# - *No-onward-transmission counterfactual and one-week-ahead
#   forecasts.* Future expected deaths from infections already seeded,
#   and a posterior-predictive projection of each stream.
#md #
#md # ```@raw html
#md # </details>
#md # ```
#
# ## Limitations
#
# The limitations are grouped by the data, the model assumptions and
# design, and the implementation, with the most consequential first in
# each group.
#
#md # ```@raw html
#md # <details><summary>Expand: limitations</summary>
#md # ```
#md #
# **Data and what it can support**
#
# - *Most quantities rest on weakly-informed priors.* Nearly all of the
#   delays, the case-fatality ratio and the laboratory assumptions are
#   set by priors informed at best by a handful of literature sources,
#   often from other outbreaks, and in places by our own prior judgement
#   rather than anything from this outbreak. The data do little to move
#   them, so these posteriors largely track their priors. We fit the
#   between-report increments, so the trajectory informs the change in
#   the reproduction number over the window but is uninformative about the
#   delays, the surveillance dispersion or the reporting fractions on
#   their own.
# - *Only report-date totals, no epidemiological dating.* We have no
#   counts by symptom onset or any other epidemiologically relevant date,
#   only cumulative totals at the report date. The timing of the
#   underlying epidemic is therefore weakly identified, and we recover it
#   only through the assumed delays.
# - *Fitted to aggregate counts.* The DRC data are situation-report
#   totals of suspected cases and deaths, laboratory-confirmed cases and
#   deaths, and specimens received and analysed; the Uganda data are
#   three export cases with one death. We do not have a line list or
#   information on case definitions or reporting completeness. The
#   laboratory testing series gives partial information on testing
#   capacity, but it is incomplete and stops at the cut-off. Every
#   estimate is a model-based extrapolation under strong assumptions, not
#   a measurement.
# - *Later sitreps revise earlier figures.* A later situation report can
#   revise an earlier total up or down as suspects are reclassified and
#   newly-reporting health zones are added, and ascertainment probably
#   rose over the window. We do not model this revision process.
# - *Streams share one case pool.* They are fitted as conditionally
#   independent given latent incidence but observe overlapping people,
#   which can understate uncertainty. Whether the streams imply mutually
#   consistent outbreak sizes is not assessed here.
#
# **Model assumptions and design**
#
# - *Inherits McCabe et al.'s epidemiological assumptions.* A single
#   zoonotic seed, an assumed generation interval, no spatial structure
#   beyond the Ituri / Nord Kivu split, and no depletion of
#   susceptibles. The onset-to-death delay is grounded on Isiro 2012 and
#   the genetic seeding bound on an external clock rate, neither
#   propagating cross-outbreak or clock uncertainty.
# - *Intervention ramp is weakly identified.* With only a few sitreps
#   straddling it, the ramp effect and the pre-ramp reproduction number
#   are not well separated.
# - *Single national bed capacity.* The treatment-centre model carries one
#   national bed capacity and one national demand, so it cannot represent
#   local saturation — on 13 June Ituri was at 93.9% occupancy while
#   Sud-Kivu was at 21.9% — and the national bed shortfall understates
#   local unmet need.
#
# **Implementation**
#
# - *LLM-driven reimplementation.* The model code, priors and analysis
#   were drafted by a language model from the McCabe et al.
#   [mccabe2026](@cite) report and the companion delay reanalysis, then
#   reviewed and revised. Not independently replicated against the
#   authors' code.
#md #
#md # ```@raw html
#md # </details>
#md # ```
#
#md # ```@raw html
#md # <details><summary>Load packages, data and fitted chains</summary>
#md # ```

## Shared setup: packages, observations, the fit registry and every model fit
## (loaded from the content-addressed cache). See `docs/examples/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "examples", "_setup.jl"))

#md # ```@raw html
#md # </details>
#md # ```

# ## Methods
#
# ### [Data](@id methods-data)
#
# The DRC data come from the situation reports of the Institut National
# de Santé Publique [insp_sitrep_2026](@cite). Each report gives the
# national cumulative suspected cases and deaths, laboratory-confirmed
# cases and deaths, and the specimens received and analysed by the
# laboratory, at the report date. From SitRep 013 (27 May) INSP began
# reclassifying suspects, so the cumulative suspected count falls. We freeze
# it at its last stable vintage (26 May) and instead read the daily
# new-suspect count ("nouveaux cas suspects du jour") that the
# confirmed-based reports publish from 4 June, fitting it as a daily
# incidence where the cumulative series stops. The same reports print a daily
# new suspected-death count alongside it ("cas suspects du jour N (M deces)",
# from 7 June), fitted the same way where the cumulative suspected-death
# series stops. The confirmed-based reports
# also publish a daily "Patients en isolement" count, the number of patients
# (confirmed plus suspected) in an isolation/treatment bed at the end of the
# day; we fit it as the suspect inflow carried through a length-of-stay
# survival into a daily bed count. The fitted
# series runs from 1 June (SitRep 018), where the column is relabelled to the
# all-patients "Patients en isolement - hospitalisation"; the narrower
# suspects-only count in SitReps 016-017 is a different quantity and is left
# out. The reports also print a cumulative "cumul guéris" total of confirmed
# cases recorded as recovered, from 6 June; we fit it as survivors among the
# modelled confirmed cases (a scaled confirmation-to-recovery convolution,
# the incidence analogue of the isolation prevalence stream). From 13 June the
# reports add a Tableau 6 patient-movement table for the treatment centres, and
# we read its daily admissions, in-care deaths, rule-outs and absconded flows as
# four count streams feeding the same treatment-centre model. The same table
# breaks the occupancy into `dont confirmés (NC+AC)` and `dont suspects`
# sub-rows, two prevalence sub-stocks that sum to the total each day, which we
# read as two further census streams splitting the occupancy. We extracted
# these figures from the
# written situation-report PDFs (archived by INRB-UMIE
# [inrb_umie_2026](@cite)) using a language model, with a second pass to
# re-read them, rather than the published per-zone CSVs. The zone sums in
# the CSVs are inconsistent with the national headline totals because they
# drop counts not yet attributed to a zone, so they understate the national
# totals.
# The Uganda data are the cases and the one death exported across the
# border, taken from the WHO situation reports and Disease Outbreak News
# [who_don_2026_602](@cite). The cross-border traveller volume and source
# population come from McCabe et al. [mccabe2026](@cite); the source
# population is fixed and the traveller volume is given a Normal prior
# around the McCabe et al. figure.
#
# From SitRep 059 (12 July) the analytique-format situation reports also
# carry a raster figure of confirmed cases by symptom-onset date, split
# alive/deceased ("courbe épidémique par date de début des symptômes"),
# with no accompanying data table. We digitise it with a dependency-free
# pixel-recovery script (later ported to Python for the automated
# updater), self-calibrated from each figure's
# own axis ticks, giving one block of onset-date counts per situation-report
# vintage; this is the only published source for the onset-date
# distribution, so we fit it directly rather than leaving it as a data
# table (see the [symptom-onset reporting delay](@ref "Symptom-onset
# reporting delay") submodel below).
# Digitisation carries measured error: per-vintage totals run -3.0% to
# +1.6% against the figure's own printed total (the sign is not one-sided),
# a per-bar pixel-noise SD of roughly $\pm 2.1$ cases, and an independent
# per-scan level error with SD roughly $4.0\%$ shared by every bar in a
# single vintage's figure (a shared dilation of that one read of the whole
# curve). Because a case's onset date is sometimes corrected between
# vintages (a case moves from one bar to a neighbouring one), an
# individual bar can fall between two vintages even though the running
# total cannot; this is expected and is not scan error. Some situation
# reports reprint an earlier vintage's figure unchanged: SitReps 061-062,
# 069-071 and 073-074 are each byte-identical to the vintage before them,
# so we collapse each such group to its earliest report date before
# fitting, keeping eleven distinct digitised snapshots out of fifteen
# digitised vintages, with report dates 12, 13, 14, 17, 18, 19, 20, 21,
# 22, 25 and 26 July. The fit only ever sees snapshots at or before the
# current cut-off, the same rule as every other stream, so a vintage
# digitised ahead of the manifest enters automatically, with no code
# change, once the cut-off advances past its report date.
#
# Each figure also stops its horizontal axis short of its own report date,
# by anything from zero to eight days depending on the vintage, and the
# last bar it does print is often a substantial count rather than a tail
# fading to zero. An onset date past that axis is not a bar of height zero,
# it is a date the figure says nothing about, so we drop those cells rather
# than scoring them as zeros: the axis gap is about five days for most
# vintages, close to the reporting delay itself, and reading it as "nothing
# reported yet" would force the fitted hazard to near zero over the first
# five days and pile the missing mass onto the delay at which the axis
# first covers the date. The publisher's choice of axis limit does not
# depend on the counts it hides, so treating those cells as missing rather
# than zero is the conservative reading. The price is that a correction
# needs both vintages of a pair to print the date, which leaves the very
# shortest delays thinly observed: no cell in the current triangle reaches
# delay zero, and the fitted hazard there rests on pooling across delays
# rather than on data.
#
# The first table lists each figure at the cut-off, or at the date
# reporting stopped for that stream. The second table gives the per-date
# history of each situation-report stream; the model fits the
# between-report increments of these series, so a single date reduces to
# the cut-off total.

#md # ```@raw html
#md # <details><summary>Loading observations and building the data table</summary>
#md # ```

observations_table = DataFrame(
    field = [
        "exported_cases",
        "exports_deaths",
        "suspected_deaths",
        "suspected_cases",
        "confirmed_cases",
        "confirmed_deaths",
        "onset_curve_reported",
        "specimens_analysed",
        "treatment_admissions",
        "treatment_deaths",
        "treatment_ruleouts",
        "treatment_absconded",
        "genetic_tmrca_bound",
        "daily_outbound_travellers (prior mean)",
        "daily_outbound_travellers_sd (prior SD)",
        "source_population"
    ],
    date = [
        isempty(obs.export_case_days) ? missing :
        grid_date(maximum(obs.export_case_days)),
        isempty(obs.export_death_days) ? missing :
        grid_date(maximum(obs.export_death_days)),
        hist_last_date(obs.deaths_history),
        hist_last_date(obs.reported_history),
        hist_last_date(obs.confirmed_history),
        hist_last_date(obs.confirmed_deaths_history),
        isempty(obs.onset_curve_history.report_days) ? missing :
        grid_date(maximum(obs.onset_curve_history.report_days)),
        hist_last_date(obs.lab_history),
        hist_last_date(obs.treatment_admissions_history),
        hist_last_date(obs.treatment_deaths_history),
        hist_last_date(obs.treatment_ruleout_history),
        hist_last_date(obs.treatment_absconded_history),
        grid_date(obs.n - obs.tmrca_days),
        missing,
        missing,
        missing
    ],
    value = [
        obs.exported_cases,
        obs.exports_deaths,
        obs.total_deaths,
        obs.reported_cases,
        obs.confirmed_cases,
        obs.confirmed_deaths,
        obs.onset_curve_history.last_total,
        obs.tests_analysed,
        isempty(obs.treatment_admissions_history.counts) ? missing :
        obs.treatment_admissions_history.counts[end],
        isempty(obs.treatment_deaths_history.counts) ? missing :
        obs.treatment_deaths_history.counts[end],
        isempty(obs.treatment_ruleout_history.counts) ? missing :
        obs.treatment_ruleout_history.counts[end],
        isempty(obs.treatment_absconded_history.counts) ? missing :
        obs.treatment_absconded_history.counts[end],
        obs.tmrca_days,
        ITURI_DAILY_TRAVEL,
        ITURI_DAILY_TRAVEL_SD,
        ITURI_POPULATION
    ]
);

#md # ```@raw html
#md # </details>
#md # ```

observations_table #hide

# The per-date cumulative history of the DRC situation-report streams,
# the national totals at each report date. The joint model fits the
# between-report increments of these series, so a single date reduces to
# the cut-off total. Two columns are the exception. `suspected_new_daily`
# is a per-day new-suspect count (not a cumulative total), fitted directly
# as a daily incidence, and it picks up where the cumulative
# `suspected_cases` column freezes on 26 May. `patients_isolated` is a
# daily count of patients in an isolation/treatment bed, fitted as the
# suspect inflow carried through a length-of-stay survival. See
# `data/observations.toml` for the per-stream sources.

#md # ```@raw html
#md # <details><summary>Building the per-date time-series table</summary>
#md # ```

vintage_table = let
    ## Each history carries grid day-indices and counts; key the counts
    ## by calendar date so every stream lines up in one table.
    bydate(h) = Dict(grid_date(d) => c for (d, c) in zip(h.days, h.counts))
    streams = (
        suspected_cases = bydate(obs.reported_history),
        suspected_new_daily = bydate(obs.suspected_daily_history),
        patients_isolated = bydate(obs.isolation_history),
        suspected_deaths = bydate(obs.deaths_history),
        suspected_new_daily_deaths = bydate(obs.suspected_daily_deaths_history),
        confirmed_cases = bydate(obs.confirmed_history),
        confirmed_deaths = bydate(obs.confirmed_deaths_history),
        recovered_confirmed = bydate(obs.recovered_history),
        specimens_received = bydate(obs.tests_received_history),
        specimens_analysed = bydate(obs.lab_history)
    )
    dates = sort(collect(union((keys(s) for s in streams)...)))
    at(s) = [haskey(s, d) ? s[d] : missing for d in dates]
    DataFrame(
        date = dates,
        suspected_cases = at(streams.suspected_cases),
        suspected_new_daily = at(streams.suspected_new_daily),
        patients_isolated = at(streams.patients_isolated),
        suspected_deaths = at(streams.suspected_deaths),
        confirmed_cases = at(streams.confirmed_cases),
        confirmed_deaths = at(streams.confirmed_deaths),
        recovered_confirmed = at(streams.recovered_confirmed),
        specimens_received = at(streams.specimens_received),
        specimens_analysed = at(streams.specimens_analysed)
    )
end;

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Per-date situation-report data table</summary>
#md # ```

vintage_table #hide

#md # ```@raw html
#md # </details>
#md # ```

# ### Model
#
# #### Model overview
#
# We model a single outbreak seeded by a zoonotic introduction on a
# daily grid from a seeding date to the cut-off (day $n$). The
# generating infection process produces daily infection incidence via
# the discrete renewal equation
#
# ```math
# I_t = R_t \sum_{s = 1}^{L} I_{t-s}\, g_s, \tag{1}
# ```
#
# where $g$ is the discretised probability mass function (PMF) of the
# generation interval, indexed from lag 1 so an infectee is always infected
# strictly after its infector, and $R_t$ is the per-day reproduction
# number. A PMF gives the probability assigned to each whole-day lag. We
# never observe infections directly. Each data stream observes a thinned,
# delayed or transformed view of the same latent incidence. This is the
# class of time-varying renewal model used in EpiNow2 [epinow2](@cite),
# with the streams fitted jointly here rather than in a pipeline.
#
# The model is assembled from modular Turing [ge2018turing](@cite)
# submodels, each holding the maths and priors for one part of the
# generative process. We describe them in generative order, from the
# infection process through the epidemiological delays to the observation
# streams. The implementation uses Mooncake [mooncake_jl](@cite)
# reverse-mode automatic differentiation, CensoredDistributions for delay
# discretisation, FlexiChains for chain handling, and PairPlots
# [pairplots_jl](@cite) with AlgebraOfGraphics [danisch2021makie](@cite)
# for the figures. Each submodel's source is shown in the collapsible
# block beneath its prose.
#
# The table below shows which parameters inform each observation submodel.
# The *analysed* column is the analysed-specimen volume, the single
# laboratory stream fitted as a count; the *confirmed* positives are scored
# as a Binomial of the observed analysed denominator with a positivity
# linked to the composition of the suspected pool, so the laboratory data
# help identify the non-BVD background. The *conf. deaths* column mirrors the
# laboratory pipeline on the death side, with a death testing fraction and a
# death-pool composition positivity built from the same assay:
#
# | Parameter | Exports | Deaths | Cases | Analysed | Confirmed | Conf. deaths | Export deaths |
# |---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
# | Reproduction number $R_t$ | ● | ● | ● | ● | ● | ● | ● |
# | Generation interval | ● | ● | ● | ● | ● | ● | ● |
# | Incubation period | ● | ● | ● | ● | ● | ● | ● |
# | Seed $I_0$ | ● | ● | ● | ● | ● | ● | ● |
# | Onset-to-death delay |  | ● |  |  |  | ● | ● |
# | Case-fatality ratio |  | ● |  |  |  | ● | ● |
# | Death ascertainment $p_{\text{death}}$ |  | ● |  |  |  | ● |  |
# | Background CFR $\mathrm{cfr}_{\text{bg}}$ |  | ● |  |  |  | ● |  |
# | Onset-to-report delay |  |  | ● | ● | ● |  |  |
# | Receipt delay |  |  |  | ● | ● | ● |  |
# | Onset-to-detection delay | ● |  |  |  |  |  |  |
# | Assay sensitivity / specificity |  |  |  |  | ● | ● |  |
# | Severity enrichment $\delta_0$ |  |  |  |  | ● |  |  |
# | Death testing fraction $\tau_{\text{death}}$ |  |  |  |  |  | ● |  |
# | Testing fraction $\tau_{\text{test}}$ |  |  |  | ● | ● |  |  |
# | Background rate $\lambda_{\text{bg}}$ |  | ● | ● | ● | ● | ● |  |
# | Surveillance dispersion |  | ● | ● | ● |  |  |  |
# | Ascertainment | ● |  | ● | ● | ● |  | ● |
# | Traveller volume | ● |  |  |  |  |  | ● |

#md # ```@setup main
#md # using BVDOutbreakSize, CodeTracking, Revise
#md # ```

# #### Infections
#
# The infection process combines several components. These are the
# reproduction number, the generation interval that drives the renewal, the
# seeding that sets the initial infection count, the genetic bound on the
# outbreak age, the growth rate that fills the unobserved cryptic phase, and
# the renewal construction that grows the seed forward to the cut-off. Each
# is described in a subsection below.
#
# ##### Reproduction number
#
# The reproduction number is held flat at the established reproduction
# number $R_0$ until a month before the first WHO situation report, then
# follows a non-centred Gaussian random walk on the log scale with weekly
# knots to the cut-off. The month-long lead lets $R_t$ start moving before the
# first report, since transmission may already have turned before the outbreak
# was formally reported; the walk start is floored at the renewal start. The
# walk starts from $R_0$ at its first knot:
#
# ```math
# \log R_k = \log R_0 + \sigma_{\text{rw}}
#            \sum_{j=1}^{k} z_j, \quad
# z_j \sim \mathrm{Normal}(0, 1), \qquad
# \sigma_{\text{rw}} \sim \mathrm{Normal}^{+}(0,\ 0.1). \tag{2}
# ```
#
# We do not place a prior on $R_0$ directly. We put the prior on the initial
# growth rate $r$ instead, given in the seeding and growth subsection below,
# and derive the established reproduction number forward from it through the
# Euler–Lotka relation under our generation interval $g$:
#
# ```math
# R_0 = \left( \sum_{s \ge 1} g_s\, e^{-r s} \right)^{-1}. \tag{3}
# ```
#
# The step-size prior keeps weekly changes in the reproduction number
# moderate. We set the half-normal on $\sigma_{\text{rw}}$ so that the
# reproduction number is unlikely to change by more than about 20% from one
# week to the next: two standard deviations of the weekly log-step is around
# $0.20$.
#
# Daily $\log R_t$ is the linear interpolation between the weekly knots, so
# the reproduction number varies piecewise linearly within each week; before
# the first knot it is held flat at $R_0$ (the interpolation clamps below the
# first knot rather than extrapolating):
#
# ```math
# \log R_t = \log R_k +
#     \frac{t - d_k}{d_{k+1} - d_k}\,(\log R_{k+1} - \log R_k),
# \qquad d_k \le t \le d_{k+1}, \tag{4}
# ```
#
# with $d_k$ the day of knot $k$. The outbreak response adds a sampled
# effect shaped by a logistic ramp at the first WHO situation report on 18
# May 2026. We assume the response takes about three weeks (21 days) to take
# effect, and that it can only reduce transmission, so the effect is
# constrained to be non-positive:
#
# ```math
# \log R_t \mathrel{+}= \delta \cdot
#     \mathrm{logistic}\!\left(\frac{t - t_{\text{bp}}}{21}\right),
# \qquad
# \delta \sim \mathrm{Normal}^{-}(0,\ 0.4). \tag{5}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: rt_walk_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.rt_walk_model(7, 0.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Generation interval
#
# We assume the generation interval $g$ is a Gamma distribution with a
# sampled shape $\alpha$ and scale $\theta$, taken from the Ebola virus
# disease serial interval used as a generation-time proxy (mean 15.3 d, SD
# 9.3 d; WHO Ebola Response Team 2014). That distribution maps once to a
# Gamma shape near $2.71$ and scale near $5.65$, and the priors are centred
# on those values, with spreads that carry the source's reported uncertainty
# on the mean rather than a spread we assign ourselves:
#
# ```math
# \alpha \sim \mathrm{Normal}^{+}(2.71,\ 0.70), \qquad
# \theta \sim \mathrm{Normal}^{+}(5.65,\ 1.50). \tag{6}
# ```
#
# The Gamma is discretised through the same double-interval-censoring route
# as every delay, described with the first epidemiological process model
# below, and the lag-0 bin is dropped and the remainder renormalised so the
# generation interval starts at one day and an infectee is infected strictly
# after its infector.

#md # ```@raw html
#md # <details><summary>Submodel: generation_interval_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.generation_interval_model(40)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Seeding and growth
#
# We assume the outbreak started from a single seed case introduced by a
# zoonotic spillover. The initial infection count $I_0$ on the last day of
# the seeding window has a prior centred on a single seed:
#
# ```math
# I_0 \sim \mathrm{Normal}^{+}(0.1,\ 0.1). \tag{7}
# ```
#
# From that seed we assume the outbreak grew deterministically through an
# unobserved cryptic exponential phase, doubling $m$ times before sustained
# transmission was established. The cryptic phase grows the seed to $2^m$
# infections at the renewal start, the day the renewal takes over, over a
# duration $m\,\tau$ with $\tau$ the doubling time. The doubling count has a
# wide prior centred on three cryptic doublings:
#
# ```math
# m \sim \mathrm{Normal}^{+}(3,\ 3), \qquad
# \tau = \frac{\log 2}{r}, \qquad
# T_{\text{cryptic}} = m\,\tau. \tag{8}
# ```
#
# The growth rate $r$ carries the prior the genetic source informs. The
# genetic reanalysis reports the epidemic doubling time as 15.2 to 24.5 d
# across substitution-rate assumptions, with a centre near 20 d. We put a
# log-normal prior on $r$ equivalent to a log-normal prior on the doubling
# time centred on 20 d, with its log spread read from that range and
# inflated a little, so the prior is slightly wider in spread than the
# source but unbiased relative to it:
#
# ```math
# r \sim \mathrm{LogNormal}\!\left(\log\tfrac{\log 2}{20},\ 0.15\right). \tag{9}
# ```
#
# This single growth rate fills the cryptic phase and, through the forward
# Euler–Lotka derivation above, sets the established reproduction number, so
# the cryptic phase and the renewal share one growth source. The genetic
# report's own established reproduction number of about $1.31$ to $1.55$ uses
# its own generation interval; deriving $R_0$ forward from the shared growth
# rate under our generation interval is the consistent choice.

#md # ```@raw html
#md # <details><summary>Submodel: seed_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.seed_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: exponential_growth_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exponential_growth_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Genetic bound on outbreak age
#
# A BEAST time tree of the first ten sequenced genomes
# [virological2026](@cite) places the TMRCA, the age of the oldest
# internal node of the tree, at a mean of 25 March 2026. The temporal
# sampling range is too short to estimate the molecular clock, so we fix it
# to the $1.2\times10^{-3}$ substitutions/site/year rate of the 2013-2016
# West African Ebola epidemic [holmes2016](@cite). The TMRCA is a lower
# bound on the outbreak age. Adding sequences, or more geographically
# representative ones, can only push it earlier, never later, as the
# sampled tree is almost entirely from Bunia. Using the genetic TMRCA as a
# one-sided seeding bound rather than a point estimate follows a suggestion
# of N. Ferguson [ferguson2026](@cite).
#
# We treat the TMRCA day as a right-censored, noisy reading of the total
# outbreak age $T$ (the cryptic duration plus the observed window, defined
# in the infection process below):
#
# ```math
# \text{tmrca}_{\text{days}} \sim
#   \mathrm{censored}\!\bigl(\mathrm{Normal}(T,\ \sigma);\
#   \text{upper} = \text{tmrca}_{\text{days}}\bigr),
# \qquad \sigma = 15\ \text{d}. \tag{10}
# ```
#
# The renewal starts on the grid day on which the renewal recursion begins
# and sustained transmission is treated as established. We place it 14 days
# after the genetic TMRCA day, past the molecular-clock uncertainty, so the
# observed window from the renewal start to the cut-off is shorter than the
# TMRCA age. The bound therefore stays informative on the cryptic duration,
# pulling the origin to sit at or before the most recent common ancestor and
# bounding the cryptic phase from below. It is one-sided, leaving the age
# free above the TMRCA. We fix the clock and do not propagate cross-outbreak
# or clock uncertainty.

#md # ```@raw html
#md # <details><summary>Submodel: genetic_seeding_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.genetic_seeding_model(100.0, 50.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Infection process
#
# The renewal start and observed window from the genetic bound above are
#
# ```math
# \text{renewal start} = n - \text{tmrca}_{\text{days}} + 14, \qquad
# \tau_{\text{obs}} = n - \text{renewal start}. \tag{11}
# ```
#
# The grid days before the renewal start are filled by the cryptic
# exponential curve at rate $r$ ending at $2^m$, giving the recursion a full
# generation interval of history. The renewal then grows the trajectory
# forward under the time-varying reproduction number. The total outbreak age
# is the cryptic duration plus the observed window:
#
# ```math
# T = m\,\tau + \tau_{\text{obs}}. \tag{12}
# ```
#
# Cumulative infections are the running sum of the daily infection series.
# The cumulative infection count at the cut-off is the headline outbreak
# size. The current growth rate is the exponential growth implied by the
# cut-off reproduction number and the generation interval through forward
# Euler–Lotka, so it is sign-consistent with that number by construction, and
# the current doubling time is $\log 2$ divided by that rate.

#md # ```@raw html
#md # <details><summary>Submodel: infection_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.infection_model(40)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# #### Epidemiological process models
#
# We model each observed stream as a delayed and thinned view of the daily
# onset incidence. This section gives the delays that map infections to
# onsets and onsets to each observed endpoint, and the case-fatality ratio
# that maps onsets to deaths. The incubation period comes first, then the
# onset-to-report delay (also used for export detection), the onset-to-death
# delay and the report-to-receipt delay, then the case-fatality ratio.
#
# ##### Incubation period
#
# Infections are convolved with the incubation-period PMF to give daily
# symptom-onset incidence, computed once and consumed by every downstream
# observation stream. We use the Bundibugyo virus incubation-period estimate
# from the 2007 Uganda outbreak (mean 6.3 d, 95% CI 5.2-7.3,
# $n = 24$; [macneil2010](@cite)). The mean prior reproduces that 95% CI;
# the source reports no interval on the spread, so the SD prior is our own
# choice:
#
# ```math
# \mu_{\text{inc}} \sim \mathrm{Normal}^{+}(6.3,\ 0.54), \qquad
# \sigma_{\text{inc}} \sim \mathrm{Normal}^{+}(3.5,\ 0.8). \tag{13}
# ```
#
# Every delay is discretised to a daily PMF over lags $0,\dots,n_{\max}$ by
# double interval censoring [charniga2024](@cite). The delays the companion
# line-list reanalysis reports, namely the onset-to-admission delay (used for
# both suspected-case reporting and export detection) and the two
# onset-to-death components, are carried through on their natural Gamma shape
# and scale,
# with the reanalysis's reported uncertainty, like the generation interval
# above. The incubation period and the laboratory receipt delay are not in the
# line list, so they keep a mean-and-SD prior moment-matched to a LogNormal.
# The LogNormal and Gamma CDFs both differentiate cleanly under the
# reverse-mode automatic differentiation. The maximum lag $n_{\max}$ is not
# hand-set: for each delay it is the 98th percentile of the prior-centre
# distribution, computed once outside the model, so every delay captures a
# consistent 98% of its mass before the truncated PMF is renormalised.
#
# Both the primary event (the onset, say) and the secondary event (the
# report) are observed only to the day, so the discretisation censors both.
# The primary event is taken uniform over its day and the secondary event is
# interval-censored to its day, giving the daily PMF
#
# ```math
# f_s = \int_0^1 \big[\, F(s + 1 - u) - F(s - u) \,\big]\, \mathrm{d}u,
# \qquad F = \text{the delay CDF}, \tag{14}
# ```
#
# which is then renormalised over lags $0,\dots,n_{\max}$.
#
# The incubation period also enters the infection-to-detection and
# infection-to-death delays for the export streams, where the survival clock
# runs from infection rather than onset.

#md # ```@raw html
#md # <details><summary>Submodel: onset_incidence_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.onset_incidence_model(Float64[])), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: censored_delay_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown, Distributions
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.censored_delay_model(30;
#md #         mean_prior = truncated(Normal(5, 1); lower = 1),
#md #         sd_prior = truncated(Normal(3, 1); lower = 1))), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Onset-to-report delay
#
# The delay from symptom onset to a suspected case being detected and
# reported into surveillance. We use a Bayesian reanalysis
# [bdbv_linelist_analysis_2026](@cite) of the 2012 Isiro Bundibugyo virus
# outbreak line list [rosello2015](@cite), taking its onset-to-admission delay
# as a Gamma sampled on its natural shape and scale, with priors centred on
# the reanalysis posterior (implied mean about 4 d) and carrying its reported
# uncertainty:
#
# ```math
# \alpha_{\text{rep}} \sim \mathrm{Normal}^{+}(1.18,\ 0.28), \qquad
# \theta_{\text{rep}} \sim \mathrm{Normal}^{+}(3.69,\ 1.20). \tag{15}
# ```
#
# We do not use the reanalysis onset-to-notification delay, a near-exponential
# Gamma with mean about 20 d.
# We assume that delay reflects a longer notification pathway, likely
# including laboratory confirmation and administrative processing, rather than
# the rapid surveillance report we model.
# This delay drives the suspected-case, laboratory and confirmed-death
# streams, and the export model uses the same onset-to-admission delay for
# detection abroad.
#
# ##### Onset-to-death delay
#
# McCabe et al. take the onset-to-death delay from the same line list as a
# point estimate [rosello2015](@cite), fitting a $t$-distributed delay. The
# reanalysis instead fits it as two atomic Gamma components, onset-to-admission
# and admission-to-death, and convolves them. We do the same: each component is
# a Gamma sampled on its natural shape and scale, with priors centred on the
# reanalysis posteriors:
#
# ```math
# \alpha_{\text{oa}} \sim \mathrm{Normal}^{+}(1.18,\ 0.28), \quad
# \theta_{\text{oa}} \sim \mathrm{Normal}^{+}(3.69,\ 1.20), \\
# \alpha_{\text{ad}} \sim \mathrm{Normal}^{+}(2.15,\ 0.60), \quad
# \theta_{\text{ad}} \sim \mathrm{Normal}^{+}(3.91,\ 1.38). \tag{16}
# ```
#
# and the onset-to-death PMF is the convolution of the two discretised
# components (implied mean about 13 d). The source is shown with the deaths
# submodel below, where the delay is injected.
#
# ##### Onset-to-detection delay (exports)
#
# An exported case is detected at a point of entry abroad when it first enters
# surveillance, the same event as a domestic suspected-case report, so the
# export model uses the same line-list onset-to-admission delay
# [bdbv_linelist_analysis_2026](@cite) as the onset-to-report delay above,
# with the same natural shape and scale priors:
#
# ```math
# \alpha_{\text{det}} \sim \mathrm{Normal}^{+}(1.18,\ 0.28), \qquad
# \theta_{\text{det}} \sim \mathrm{Normal}^{+}(3.69,\ 1.20). \tag{17}
# ```
#
# It drives the exports streams; its source is shown with the exports
# submodel below.
#
# ##### Report-to-analysed delay
#
# The delay from a suspected case being reported to its specimen being
# analysed by the laboratory, centred on a short turnaround with a heavy
# right tail allowing for specimen shipment to a confirmatory laboratory and
# the analysis queue. No per-sample outbreak data grounds this, so the prior
# is our own choice:
#
# ```math
# \mu_{\text{rec}} \sim \mathrm{Normal}^{+}(4.5,\ 1.0), \qquad
# \sigma_{\text{rec}} \sim \mathrm{Normal}^{+}(4.0,\ 0.75). \tag{18}
# ```
#
# It drives the laboratory analysed-specimen volume; its source is shown
# with the laboratory submodel below.
#
# ##### Case-fatality ratio
#
# The US Centers for Disease Control and Prevention (CDC) summary for
# the two previous BVD outbreaks is $55$ deaths in $169$ cases
# ($\approx 33\%$;
# [CDC outbreak history](https://www.cdc.gov/ebola/outbreaks/index.html)),
# with confidence bands spanning
# roughly $26$-$40\%$. The companion Bundibugyo virus (BDBV) reanalysis
# reports a baseline of $0.47$ ($95\%$ CrI $0.31$-$0.65$) for
# non-healthcare-worker (non-HCW) confirmed cases. Based on this we use a
# prior of
#
# ```math
# \mathrm{CFR} \sim \mathrm{Beta}(6.6,\ 13.4), \tag{19}
# ```
#
# with mean $0.33$ and $95\%$ interval roughly $0.15$-$0.54$. The mean
# matches the CDC $55/169 \approx 33\%$ figure and the corrected central
# CFR in the 20 May report [mccabe2026update](@cite).

#md # ```@raw html
#md # <details><summary>Submodel: cfr_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.cfr_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# The prior density, with the CDC $0.33$ figure marked.

cfr_prior_fig = plot_cfr_prior(Beta(6.6, 13.4)); #hide
cfr_prior_fig #hide

# #### Observation models
#
# Each observation submodel takes the shared daily onset incidence,
# convolves it with a sampled onset-to-event delay, scales by the relevant
# ascertainment, case-fatality ratio or positivity factor, and reads the
# modelled count off the daily series at each vintage day. Likelihoods score
# the between-vintage increments. The surveillance streams come first, then
# the geographic-spread exports.
#
# ##### Shared observation submodels
#
# Several parameters are assumed shared across the streams: the surveillance
# dispersion, the ascertainment fractions, the laboratory testing priors and
# the traveller volume. We assume the passive-surveillance count datasets are
# overdispersed and share a common dispersion, and that the laboratory
# testing priors are shared between the suspected-case and laboratory
# streams. More detail is given in the subsections below.
#
# ###### Surveillance dispersion
#
# Each passive-surveillance count stream has its own negative-binomial
# dispersion, partially pooled across the streams so the sparse ones borrow
# strength. Following Stan prior-choice recommendations
# [stan_prior_choice](@cite), the dispersion is sampled on the $1/\sqrt{k}$
# scale in non-centred log form:
#
# ```math
# \log\!\bigl(1/\sqrt{k_s}\bigr) = \mu + \tau\, z_s, \quad
# z_s \sim \mathrm{Normal}(0, 1), \qquad
# \mu \sim \mathrm{Normal}(\log 0.6,\ 0.33), \quad
# \tau \sim \mathrm{Normal}^{+}(0,\ 0.6), \tag{20}
# ```
#
# so $k_s = 1/\exp(\mu + \tau z_s)^2$ per stream, with $\tau$ setting the
# pooling ($\tau = 0$ collapses to one shared dispersion). The population
# value $k = 1/\exp(\mu)^2$ is the headline dispersion.

#md # ```@raw html
#md # <details><summary>Submodel: pooled_dispersion_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.pooled_dispersion_model(6)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ###### Ascertainment
#
# Two surveillance systems detect cases: DRC passive community
# surveillance (the reported suspected-case count) and Uganda's
# point-of-entry / hospital surveillance (the exported-case count). Each
# captures a fraction of the true cases passing through it. The two
# ascertainment fractions $p_{\text{DRC}}$ and $p_{\text{Uganda}}$ share a
# logit-scale hyperprior with mean $\mu$ and pooling strength $\tau$, centred
# on a reporting fraction of $75\%$, reflecting the active case-finding of a
# declared Ebola response rather than baseline passive surveillance:
#
# ```math
# \mu \sim \mathrm{Normal}(\mathrm{logit}(0.75),\ 1),
# \qquad
# \tau \sim \mathrm{Normal}^{+}(0,\ 0.5), \tag{21}
# ```
#
# ```math
# \mathrm{logit}(p_{\text{DRC}}) \sim \mathrm{Normal}(\mu,\ \tau),
# \qquad
# \mathrm{logit}(p_{\text{Uganda}}) \sim \mathrm{Normal}(\mu,\ \tau). \tag{22}
# ```
#
# The cases likelihood uses $p_{\text{DRC}}$; the two Uganda-side
# likelihoods use $p_{\text{Uganda}}$.

#md # ```@raw html
#md # <details><summary>Submodel: pooled_ascertainment_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.pooled_ascertainment_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ###### Laboratory priors
#
# We model the process of confirming cases via laboratory testing. The
# testing fraction $\tau_{\text{test}}$ is the share of suspected cases routed
# to the laboratory. A truly BVD specimen tests positive with the assay
# sensitivity $s$, and a non-BVD specimen tests positive with the
# false-positive rate $1 - \mathrm{spec}$ from the assay specificity. We assume
# that more severe cases, more likely to be Ebola, are preferentially tested,
# through an
# enrichment factor $\delta_0$ that raises the tested BVD share above the
# suspect-pool composition early on and relaxes towards it as testing
# broadens. The confirmed deaths mirror this laboratory pipeline rather than
# enriching the case composition: a fraction $\tau_{\text{death}}$ of suspected
# deaths reach the laboratory, and they confirm at the assay positivity
# $p = s\,q_{\text{death}} + (1-\mathrm{spec})(1-q_{\text{death}})$ built from
# the same sensitivity and specificity but the death-pool BVD share
# $q_{\text{death}}$. Confirmation runs
# on the altona RealStar Filovirus Screen RT-PCR rather than the
# Zaire-specific GeneXpert Ebola assay, which does not reliably detect
# Bundibugyo virus. Sensitivity for Bundibugyo virus is less well
# characterised than for Zaire ebolavirus, so we centre the sensitivity prior
# below the values reported for other strains and give it a wide
# spread. The specificity is high but imperfect; the
# severity enrichment is moderate and one-sided (triage upsamples BVD,
# never down); the death testing-intensity scaling is a tight log-normal
# centred on one, since no death-testing data grounds it:
#
# ```math
# \tau_{\text{test}} \sim \mathrm{Beta}(5,\ 2), \qquad
# s \sim \mathrm{Beta}(10,\ 1.76), \qquad
# \mathrm{spec} \sim \mathrm{Beta}(60,\ 2),
# ```
#
# ```math
# \delta_0 \sim \mathrm{Normal}^{+}(1.5,\ 0.75), \qquad
# \text{scaling} \sim \mathrm{LogNormal}(0,\ 0.25). \tag{23}
# ```
#
# The non-BVD background rate $\lambda_{\text{bg}}$ enters the suspected-case
# stream and is described with it below; the suspected deaths carry a death
# ascertainment $p_{\text{death}} \sim \mathrm{logit}^{-1}\mathrm{Normal}
# (\mathrm{logit}\,0.9,\ 0.5)$ and a non-BVD death background tied to the case
# background by a background CFR $\mathrm{cfr}_{\text{bg}} \sim
# \mathrm{Beta}(2,\ 6)$.

#md # ```@raw html
#md # <details><summary>Submodel: test_positivity_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.test_positivity_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: confirmed_positivity_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.confirmed_positivity_model(5)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: test_sensitivity_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.test_sensitivity_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: test_specificity_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.test_specificity_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: severity_enrichment_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.severity_enrichment_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: death_testing_fraction_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.death_testing_fraction_model()),
#md #     "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: death_ascertainment_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.death_ascertainment_model()),
#md #     "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: background_cfr_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.background_cfr_model()),
#md #     "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ###### Traveller volume
#
# The number of people crossing from the source area to Uganda each day
# sets the travel rate in the exports likelihood. We treat it as an
# estimated quantity rather than a fixed input. McCabe et al. Table 3
# records mean weekly passenger counts across seven points of entry; the
# Ituri-side daily total of $1871$ is a sample mean across roughly
# $15$-$21$ point-of-entry-weeks. We use a Normal prior centred on $1871$
# with SD $200$ ($\approx 10\%$ CV), truncated at zero, covering
# point-of-entry variation and the sitrep sampling uncertainty; the source
# population is kept fixed (census):
#
# ```math
# N_{\text{travel}} \sim \mathrm{Normal}^{+}(1871,\ 200). \tag{24}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: traveller_volume_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.traveller_volume_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Reported cases
#
# Reported suspected cases are the sum of two parts. The first is a
# BVD-driven component: the daily onsets convolved with the onset-to-report
# delay $f_{\text{rep}}$ and scaled by the DRC ascertainment $p_{\text{DRC}}$.
# The convolution of a daily series $x$ with a delay PMF $f$ is the lagged sum
#
# ```math
# (x * f)_t = \sum_{s \ge 0} x_{t-s}\, f_s,
# ```
#
# used for every delay below. We write the BVD onset-to-report series at unit
# ascertainment as
#
# ```math
# \text{bvd}_t = \sum_{s \ge 0} \text{onsets}_{t-s}\, f_{\text{rep},s}.
# ```
#
# The second part is an additive non-BVD background, so a suspected case need
# not be a true BVD infection. It is a per-day rate $\lambda_{\text{bg},t}$ that
# follows a lognormal random walk on weekly knots around a baseline
# $\lambda_\mu$, linearly interpolated to the daily grid,
#
# ```math
# \lambda_{\text{bg},t} = \lambda_\mu \exp\!\bigl(w_t\bigr), \qquad
# w_t = \mathrm{interp}\Bigl(\sigma_{\text{rw}} \sum_{s < k} z_s\Bigr),
# \qquad z_s \sim \mathcal N(0, 1),
# ```
#
# gated to zero before the surveillance onset (a report-to-receipt lead before
# the first suspected-case report — the background does not exist before
# surveillance began) and shared, with one tight innovation SD
# $\sigma_{\text{rw}}$, between the suspected-case and suspected-death streams.
# Weekly knots match the reproduction-number walk and keep the background a
# gentle drift over a small number of innovations. The baseline carries a
# half-normal $\mathrm{Normal}^{+}(0, 8)$ prior on the natural scale. A
# log-scale level would have a heavy right tail the background/outbreak-size
# degeneracy could exploit, whereas the natural-scale half-normal bounds it. It
# is wide enough that the laboratory positivity (only $210/755 \approx 0.28$ of
# analysed specimens are positive) identifies the background, which is inferred
# to be the majority of the suspect pool. The
# daily expected suspected case count is
#
# ```math
# c_t = p_{\text{DRC}}\, \text{bvd}_t + \lambda_{\text{bg},t}.
# ```
#
# The per-vintage increments are scored with a NegBinomial sharing the
# dispersion $k$:
#
# ```math
# Y_{\text{cases},i} - Y_{\text{cases},i-1} \sim \mathrm{NegBinomial}\!\Bigl(
#     \sum_{t = d_{i-1}+1}^{d_i} c_t,\ k\Bigr). \tag{25}
# ```
#
# From SitRep 013 (27 May) INSP reclassifies suspects, so the national
# cumulative suspected total falls. We freeze it at 26 May and instead fit the
# daily new-suspect count that the confirmed-based reports publish (the
# "nouveaux cas suspects du jour" $a_j$ on report day $t_j$, 4-7 June). This is
# a genuine daily incidence, not a cumulative total, so it is scored against
# the modelled daily suspected count $c_{t_j}$ on that day directly (a
# single-day mean, not a between-vintage sum) with a NegBinomial sharing $k$:
#
# ```math
# a_j \sim \mathrm{NegBinomial}(c_{t_j},\ k).
# ```
#
# The daily report days fall strictly after the frozen cumulative series ends,
# so the two suspected likelihoods cover disjoint days and do not double-count.
# The suspected-death stream is fitted the same way: the cumulative
# suspected-death total freezes at 26 May and the daily new suspected-death
# count ("cas suspects du jour N (M deces)", from 7 June) is scored against the
# modelled daily suspected-death count on each report day with a NegBinomial
# sharing $k$ (see `deaths_model`).

#md # ```@raw html
#md # <details><summary>Submodel: reported_cases_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.reported_cases_model(
#md #         (; days = Int[], counts = Int[]), missing,
#md #         Float64[], 1.0, 0.25)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Treatment-centre flow
#
# The treatment-centre stream models the daily patient flow through the
# isolation/treatment centres: the occupied-bed count ("Patients en isolement"),
# the daily admissions, and the daily discharges split by outcome — in-care
# deaths, rule-outs and absconded — read from the situation-report Tableau 6
# patient-movement table. Two parallel processes act on each patient. A clinical
# course governs how long a patient occupies a bed and how they leave it, and so
# sets the total occupancy and every discharge flow; a laboratory label runs
# alongside it and only relabels a patient from suspected to confirmed, carving
# the suspect/confirmed split of the census. Death is clinical and happens under
# either label, so a true case may die before its test confirms it. Separating
# the two keeps the operational churn in the suspected pool out of the part of
# the occupancy the infection estimate leans on.
#
# A proportion $p_{\text{iso}}$ of the reported suspects need a bed. These
# admissions split into a BVD true-case inflow, admitted at a severity-skewed
# rate $p_{\text{iso,bvd}} = \mathrm{logit}^{-1}(\mathrm{logit}\,p_{\text{iso}} +
# \delta_{\text{iso}})$ above the base rate, and a non-BVD background inflow at
# the base rate,
#
# ```math
# A_{\text{bvd},t} = p_{\text{iso,bvd}}\,p_{\text{DRC}}\,\text{bvd}_t,
# \qquad
# A_{\text{bg},t} = p_{\text{iso}}\,\lambda_{\text{bg},t}, \tag{26}
# ```
#
# each carried through a short suspected-to-admission delay that captures triage,
# transport and the wait for a bed. A patient then leaves by one of four routes.
# A BVD true case dies at the in-care case-fatality ratio $\text{CFR}_{\text{iso}}$
# over the admission-to-death stay or recovers over the longer
# admission-to-recovery stay; a non-case is ruled out by a negative test over the
# rule-out stay or absconds. The death-stay prior is the admission-to-death delay
# from the line-list reanalysis [bdbv_linelist_analysis_2026](@cite), and the
# non-BVD rule-out stay takes the report-to-receipt laboratory turnaround.
#
# Occupancy is the running balance of these latent events rather than a
# length-of-stay convolution: each day the bed stock is yesterday's stock plus
# the day's admissions less the day's deaths, recoveries, rule-outs and absconds.
# The abscond outflow drains the suspected pool at a small daily fraction
# $\kappa$ of the previous day's suspected occupancy,
#
# ```math
# \text{absconds}_t = \kappa\, O_{\text{susp},t-1}, \tag{27}
# ```
#
# so the total bed demand is the forward balance
#
# ```math
# D_t = D_{t-1} + A_t - \text{deaths}_t - \text{recoveries}_t
#       - \text{rule-outs}_t - \text{absconds}_t, \tag{28}
# ```
#
# with $A = A_{\text{bvd}} + A_{\text{bg}}$. The stay lives entirely in the
# discharge flows, each an admission stream convolved with its outcome density,
# so one inflow and one set of outcome timings generate the bed stock, the
# discharge flows and the demand together. The occupancy accumulates admissions
# over the stay, so it is a smooth integral of the infection signal: a high,
# sustained bed stock informs the reproduction number, and anything in the
# occupancy that is not infection, such as an overnight reclassification of who
# is counted, is modelled rather than left to bend the transmission estimate.
#
# In-care deaths combine the two labels: a true case who dies before its test
# returns is a suspected death and one who dies after is a confirmed death, and
# the report records the two together, so the death flow is the in-care fatality
# applied to the BVD inflow over the admission-to-death stay, scored against the
# combined deaths directly and never gated by confirmation. The in-care fatality
# is a sampled log-odds modifier $\beta_{\text{iso}}$ on the infection
# case-fatality ratio,
#
# ```math
# \text{CFR}_{\text{iso}} = \mathrm{logit}^{-1}\bigl(\mathrm{logit}\,\text{CFR}
#     + \beta_{\text{iso}}\bigr), \tag{29}
# ```
#
# identified by the in-care death flow relative to admissions and occupancy. It
# is a fatality conditional on admission rather than a causal treatment effect,
# sitting below the infection case-fatality ratio where care lowers mortality,
# and it is reported with $\beta_{\text{iso}}$ and the overall length of stay (the
# death/recovery mixture mean).
#
# The laboratory label carves the census into a confirmed and a suspected
# sub-stock. Confirmation relabels a true case already in a bed at the daily
# hazard $\rho\,\tau_{\text{test}}\,p_{\text{pos},t}$: the community confirmation
# hazard $\tau_{\text{test}}\,p_{\text{pos},t}$ — the share of suspects routed to
# the laboratory times the day's positivity — borrowed from the confirmed-case
# pipeline rather than re-estimated, so the in-care confirmed stock is a subset of
# the total confirmed by construction. The confirmed-in-care stock is tracked by
# admission cohort. Each true-case admission carries two clocks from the day it
# enters a bed, a confirmation clock and a clinical-stay clock, and it counts
# toward the confirmed census only once it has been confirmed and while it is
# still in a bed,
#
# ```math
# O_{\text{conf},t} = \sum_{u \le t} A_{\text{bvd},u}\,
#     F_{\text{conf}}(u, t)\, S_{\text{clin}}(t - u), \tag{30}
# ```
#
# with $F_{\text{conf}}$ the cumulative confirmation probability of a cohort
# admitted on day $u$,
#
# ```math
# F_{\text{conf}}(u, t) = 1 - \prod_{j = u+1}^{t}
#     \bigl(1 - \tau_{\text{test}}\,p_{\text{pos},j}\bigr),
# ```
#
# a cumulative product along the cohort's age rather than a fixed distribution
# because the hazard is time-varying, and $S_{\text{clin}}$ the clinical-stay
# survival,
#
# ```math
# S_{\text{clin}}(d) = 1 - \sum_{j=0}^{d}
#     \bigl(\text{CFR}_{\text{iso}}\, f^{\text{death}}_j
#     + (1 - \text{CFR}_{\text{iso}})\, f^{\text{rec}}_j\bigr),
# ```
#
# the probability an admitted case is still in a bed after $d$ days, the
# discharge-complement of the death/recovery mixture built from the same
# admission-to-death and admission-to-recovery stays the discharge flows use.
# Because $\sum_{u \le t} A_{\text{bvd},u}\, S_{\text{clin}}(t - u)$ reconstructs
# the occupied true-case stock, the confirmed-and-present cohort is a subset of
# it and $O_{\text{conf}} \le O_{\text{bvd}} \le D$ holds by construction.
#
# Cohort tracking is needed because deaths and recoveries are observed combined
# across the two labels, so the data do not say which departing patients had
# already been confirmed. Carrying the confirmation and stay clocks separately
# excludes cases that die before their test returns from the confirmed pool,
# which a single pool-average discharge rate would over-attribute. The suspected
# sub-stock is the remainder $O_{\text{susp},t} = D_t -
# O_{\text{conf},t}$, holding the not-yet-confirmed BVD occupancy together with
# the non-case occupancy awaiting rule-out, and the abscond outflow drains this
# suspected stock at the daily fraction $\kappa$. Recoveries among the confirmed
# (the published recovery total) are the confirmed subset of recoveries and are
# modelled as a separate confirmed-recovery stream (below).
#
# Capacity enters only as a censored observation; the latent demand is never
# capped, because the demand is the quantity of interest. The bed capacity is a
# non-decreasing random walk on weekly knots, since beds are added over the
# response and not taken away, pinned by the implied bed count — the reported
# occupancy divided by the reported occupancy rate (about $400$ rising to $452$
# beds over 9–13 June). The occupied beds are scored as the latent demand
# right-censored at the recorded implied capacity, so demand above a saturated
# capacity is left uncensored, and the daily admissions are right-censored at the
# recorded free-bed headroom, the implied capacity less the previous day's
# observed occupancy. Both censoring bounds are fixed recorded data, which keeps
# the admissions censor stable where a bound tied to the modelled, wandering
# capacity would drift. The occupancy and a flow stream are scored as
#
# ```math
# O_j \sim \mathrm{censored}\bigl(\mathrm{NegBinomial}(D_{t_j},\ k_{\text{iso}});\
#     \text{upper} = C^{\text{cap}}_j\bigr),
# \qquad
# F_j \sim \mathrm{NegBinomial}(\mu^{F}_{t_j},\ k_{\text{iso}}), \tag{31}
# ```
#
# with each flow mean $\mu^{F}_t$ the matching modelled event series — the
# admissions, the in-care deaths, the rule-outs and the absconds — all sharing
# the treatment dispersion $k_{\text{iso}}$, and the implied capacity carried by a
# NegBinomial of its own. Occupancy below capacity identifies the demand directly;
# demand above a saturated capacity is only partially identified, since the
# occupancy reveals that demand was at least the beds filled but not how much
# more, so the bed shortfall above capacity is informed by the demand model and
# its priors rather than measured. Bed demand is the uncapped diagnostic, and the
# model exposes the cut-off occupancy, the cut-off bed demand (the need under
# unconstrained supply), their difference (the bed shortfall) and the utilisation.
# Out of sample, where no occupancy is observed, the forecast caps admissions at
# the modelled free beds; this modelled bound is safe because the forecast is a
# forward simulation outside the likelihood.
#
# The fitted occupancy series is the all-patients column from 1 June (SitRep 018)
# onward, and from 13 June the report adds a two-row breakdown into confirmed and
# suspected beds that sums to the total each day. The total-occupancy term is the
# backbone present from 1 June; the daily flows and the confirmed/suspected census
# add likelihood on the days they exist, scored per day as either the total or the
# split so the total and its parts are never both counted on one day. The early
# window, with only the total occupancy reported, fits the backbone alone while
# the latent admissions still drive the stock, and the split is scored only where
# the borrowed confirmation hazard is non-zero, that is, where the laboratory
# pipeline of the full model supplies it.
#
# One reporting artefact is modelled, on identified days only: an overnight
# reclassification of the total, where the published start-of-day in-bed count is
# differenced against the previous report day's occupancy, and a day whose gap
# exceeds a threshold is flagged as a break day. One step is fitted per flagged
# day, with a prior centred on that day's observed gap but free to move, so the
# fit can attribute part of a gap to genuine change in demand; the steps
# accumulate into a persistent additive offset on the modelled total occupancy,
# carried forward to every later day, that absorbs the overnight gap without
# bending the reproduction number to chase it. The split does not change the
# occupancy before 13 June, since no breakdown is published there and the total
# backbone carries that window.

#md # ```@raw html
#md # <details><summary>Submodel: treatment_flow_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.treatment_flow_model(
#md #         (; days = Int[], counts = Int[]),
#md #         Float64[], Float64[], 0.25, 0.3)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Suspected deaths
#
# Suspected deaths are the ascertained, CFR-weighted convolution of the daily
# onsets with the onset-to-death PMF $f_d$, plus a non-BVD background, modelled
# on the incidence scale. The death history ends at the cut-off, so the cut-off
# total is the final increment and is not scored separately. A fatal BVD
# infection enters the suspected-death count only when ascertained, so the BVD
# deaths carry a death ascertainment $p_{\text{death}}$, the death analogue of
# the case ascertainment $p_{\text{DRC}}$, with an informative prior centred
# high (a death is more reliably reported than a living suspect). The non-BVD
# background suspected deaths are a background CFR $\mathrm{cfr}_{\text{bg}}$
# applied to the per-day non-BVD suspected-case background
# $\lambda_{\text{bg},t}$, lagged by the same onset-to-death delay so a
# background death follows its background case; the death background tracks the
# identified case background rather than a second free, outbreak-size-
# degenerate rate. The daily death series is
#
# ```math
# m_t = p_{\text{death}}\,\mathrm{CFR} \sum_{s \ge 0} \text{onsets}_{t-s}\,
#     f_{d,s} \; + \; \mathrm{cfr}_{\text{bg}} \sum_{s \ge 0}
#     \lambda_{\text{bg},t-s}\, f_{d,s}.
# ```
#
# The per-vintage increments are scored with a NegBinomial sharing the
# dispersion $k$:
#
# ```math
# Y_{\text{deaths},i} - Y_{\text{deaths},i-1} \sim \mathrm{NegBinomial}\!\Bigl(
#     \sum_{t = d_{i-1}+1}^{d_i} m_t,\ k\Bigr). \tag{32}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: deaths_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.deaths_model(
#md #         (; days = Int[], counts = Int[]), missing,
#md #         Float64[], 1.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Laboratory pipeline
#
# The laboratory pipeline fits a single analysed-specimen volume. There is no
# separately-modelled testing capacity: the analysed volume is a deterministic
# function of the suspected-case incidence. It is the suspected daily pipeline
# ($p_{\text{DRC}}\,\text{bvd}_t$ plus the non-BVD background
# $\lambda_{\text{bg}}$) carried through the report-to-analysed delay
# $f_{\text{rec}}$ and thinned by the testing fraction $\tau_{\text{test}}$ (the
# share of suspected cases routed to the laboratory),
#
# ```math
# v_t = \tau_{\text{test}} \sum_{s \ge 0}
#     \bigl(p_{\text{DRC}}\, \text{bvd}_{t-s} + \lambda_{\text{bg},t-s}\bigr)\,
#     f_{\text{rec},s}.
# ```
#
# This analysed volume is gated to zero before the testing onset: no
# specimens are analysed before the laboratory existed, so $v_t$ does not accrue
# over the pre-surveillance cryptic phase (modelling a pre-testing volume would
# both invent capacity and roll it into the first laboratory and early-confirmed
# bins, over-predicting the early confirmed counts). The first confirmed
# vintage is treated as the baseline and the early confirmed increments are
# scored from it. The suspected-case count itself is not gated, as those cases
# did accumulate over the cryptic phase.
#
# This construction, a testing fraction times the suspected pipeline carried
# to laboratory receipt, gives the modelled case analysed volume that the
# confirmed deaths reuse: the death volume scales it at the per-day suspected
# death-to-case ratio (see the confirmed-deaths section below), so the two
# share the laboratory capacity onset.
#
# The per-vintage increments are scored against the cumulative analysed
# series with a NegBinomial sharing the dispersion $k$:
#
# ```math
# Y_{\text{ana},i} - Y_{\text{ana},i-1} \sim \mathrm{NegBinomial}\!\Bigl(
#     \sum_{t = d_{i-1}+1}^{d_i} v_t,\ k\Bigr). \tag{33}
# ```
#
# The confirmed positives in each laboratory window $v$ are scored as a
# Binomial of the observed specimens-analysed denominator $A_v$ with a
# per-window tested-positive probability $p_{\text{pos},v}$, and where no
# analysed count is observed (the early and unanchored windows) the modelled
# volume $v_t$ is the denominator instead, so the fitted volume and the
# proxy denominator are the same quantity. We tie that probability to the
# composition of the tested pool, so the confirmed data help identify the
# non-BVD background. The suspect-pool composition $\varphi_v$ is the BVD
# share among the specimens analysed in the window, carried through the same
# delay as the volume so composition and volume share one clock:
#
# ```math
# \varphi_v = \frac{(p_{\text{DRC}}\,\text{bvd} * f_{\text{rec}})_v}
#     {(p_{\text{DRC}}\,\text{bvd} * f_{\text{rec}})_v +
#      (\lambda_{\text{bg}} * f_{\text{rec}})_v}.
# ```
#
# The tested BVD share $q_v$ raises $\varphi_v$ by the decaying severity
# enrichment $\delta_0$. A truly BVD specimen then tests positive with the
# sensitivity $s$, and a non-BVD specimen with the false-positive rate
# $1 - \mathrm{spec}$, so the false-positive term carries the non-BVD share
# and the laboratory data identify the background:
#
# ```math
# q_v = \mathrm{logistic}\!\bigl(\mathrm{logit}(\varphi_v) +
#     \delta_0\, e^{-c_v / \text{decay}}\bigr), \qquad
# p_{\text{pos},v} = s\, q_v + (1 - \mathrm{spec})(1 - q_v),
# ```
#
# ```math
# C_v \sim \mathrm{Binomial}(A_v,\ p_{\text{pos},v}), \tag{34}
# ```
#
# with $c_v$ the cumulative modelled laboratory volume at window $v$, the
# clock on which the enrichment decays. The confirmed vintages before the
# first and after the last laboratory date carry no observed analysed
# denominator. They are scored as NegBinomial counts against the modelled
# laboratory volume $V_v$, the daily modelled volume $v_t$ summed over the
# window, with the same composition-linked positivity, so all the confirmed
# data are used:
#
# ```math
# C_v^{\text{no-denom}} \sim
#     \mathrm{NegBinomial}(p_{\text{pos},v}\, V_v,\ k). \tag{35}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: lab_delay_model (receipt delay)</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.lab_delay_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: confirmed_cases_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.confirmed_cases_model(
#md #         (; days = Int[], counts = Int[]), missing, Float64[],
#md #         1.0, 0.25, Float64[], 0.5, Float64[])), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Confirmed deaths
#
# The confirmed deaths mirror the confirmed-case laboratory pipeline. The
# confirmed cases fit a modelled analysed-specimen volume and score the
# positives as that volume times a composition-linked positivity; the death
# side has no published analysed denominator, so we build the death analogue
# of that volume and score the confirmed-death increments as NegBinomial
# counts of it.
#
# Deaths are tested out of the same laboratory as cases, so the death analysed
# volume tracks the modelled case analysed volume $v^{\text{c}}_t$ at the
# per-day suspected death-to-case ratio, times a testing-intensity scaling,
#
# ```math
# v^{\text{d}}_t = \text{scaling}\; v^{\text{c}}_t\,
#     \frac{\sum_{s\ge 0} m^{\text{d}}_{t-s}\, f_{\text{rec},s}}
#          {\sum_{s\ge 0} m^{\text{c}}_{t-s}\, f_{\text{rec},s}},
# ```
#
# with $m^{\text{d}}$ and $m^{\text{c}}$ the modelled suspected-death and
# suspected-case series and $f_{\text{rec}}$ the report-to-receipt delay the
# confirmed cases use. The death-to-case ratio carries the suspect-pool
# severity and the suspected-death level, so the scaling is the per-suspect
# testing-intensity difference between deaths and living suspects alone; with
# no death-testing data it is a tight log-normal centred on one. Those
# specimens confirm at the assay positivity built from the death-pool BVD share
#
# ```math
# q_{\text{death},t} = \frac{\text{bvd}^{\text{d}}_t}
#     {\text{bvd}^{\text{d}}_t + \text{bg}^{\text{d}}_t}, \qquad
# p_t = s\,q_{\text{death},t} + (1-\mathrm{spec})(1-q_{\text{death},t}),
# ```
#
# with $\text{bvd}^{\text{d}}$ and $\text{bg}^{\text{d}}$ the BVD and non-BVD
# components of the suspected deaths (both at receipt) and $s$, $\mathrm{spec}$
# the same assay sensitivity and specificity as the confirmed cases. The
# false-positive term $(1-\mathrm{spec})(1-q_{\text{death}})$ makes the
# confirmed deaths respond to the non-BVD death share, the same structural link
# the confirmed cases use; the death background (the background CFR applied to
# the case background, lagged by the onset-to-death delay) keeps the
# composition below one. The daily confirmed deaths are the positivity times
# the death analysed volume,
#
# ```math
# \text{cd}_t = p_t\, v^{\text{d}}_t,
# ```
#
# and the per-vintage increments are scored with a NegBinomial sharing the
# dispersion $k$:
#
# ```math
# Y_{\text{cd},i} - Y_{\text{cd},i-1} \sim \mathrm{NegBinomial}\!\Bigl(
#     \sum_{t = d_{i-1}+1}^{d_i} \text{cd}_t,\ k\Bigr). \tag{36}
# ```
#
# The death analysed volume inherits the laboratory capacity onset from the
# case volume $v^{\text{c}}_t$, so $\text{cd}_t$ is zero before the first
# confirmed-case vintage: no deaths are confirmed before the laboratory
# existed.

#md # ```@raw html
#md # <details><summary>Submodel: confirmed_deaths_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.confirmed_deaths_model(
#md #         missing, missing, Float64[], Float64[], Float64[], 1.0)),
#md #     "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Recovered among confirmed
#
# Recoveries ("cumul guéris") are the survivors among laboratory-confirmed
# cases, the incidence analogue of the convolution-and-scaling
# secondary-observation model of EpiNow2 [epinow2](@cite). The modelled daily
# confirmed incidence $\text{confirmed}_t$ (the per-window tested-positive
# probability on the modelled analysed volume, the same daily series the
# cumulative-confirmed trajectory uses) is scaled by the recovery proportion
# $p_{\text{rec}}$ and convolved with a sampled confirmation-to-recovery
# delay $f_{\text{rec}}$,
#
# ```math
# \text{recovered}_t = p_{\text{rec}} \sum_{s \ge 0}
#     \text{confirmed}_{t-s}\, f_{\text{rec},s}.
# ```
#
# A recovered case is one that did not die, so the recovery proportion is
# grounded on the case-fatality ratio rather than estimated independently. It
# is the complement $1 - \mathrm{CFR}$ adjusted on the log-odds scale by a
# sampled offset $\delta_{\text{rec}} \sim \mathrm{Normal}(0, 0.5)$, since the
# confirmed cases are a slightly different population from the one the CFR is
# defined over,
#
# ```math
# p_{\text{rec}} = \operatorname{logistic}\!\bigl(
#     \operatorname{logit}(1 - \mathrm{CFR}) + \delta_{\text{rec}}\bigr).
# ```
#
# A case is taken to be confirmed before it is recorded as recovered (the
# report counts recoveries among confirmed cases); a positive result could in
# principle return after a patient has already recovered, but we assume the
# reported total reflects confirmed cases recorded as recovered.
# The cumulative recovered series ends at the cut-off, so its per-vintage
# increments are fitted, like the confirmed and confirmed-death streams, with
# a NegBinomial of an independent dispersion $k_{\text{rec}}$:
#
# ```math
# Y_{\text{rec},i} - Y_{\text{rec},i-1} \sim \mathrm{NegBinomial}\!\Bigl(
#     \sum_{t = d_{i-1}+1}^{d_i} \text{recovered}_t,\ k_{\text{rec}}\Bigr).
# ```
#
# The convolution right-censors recoveries that have not yet resolved by the
# cut-off, so the small observed totals (12 to 40 over 6-13 June) are
# consistent with a high eventual survival fraction and a multi-week recovery
# delay.

#md # ```@raw html
#md # <details><summary>Submodel: recovered_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.recovered_model(
#md #         (; days = Int[], counts = Int[]), missing, Float64[], 0.3)),
#md #     "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Exported cases
#
# The exports stream is travel-gated, so the at-risk clock runs from
# infection. An infected person travels to Uganda at the daily per-capita
# travel rate $q = N_{\text{travel}} / N_{\text{source}}$ and stays at risk
# of being exported and detected only until the infection-to-detection delay
# has elapsed. The daily at-risk export prevalence is the infections still
# infected and not yet detected, scaled by the Uganda ascertainment and the
# travel rate. The infection-to-detection delay is the onset-to-detection
# delay convolved with the incubation period, so the survival clock runs from
# infection. Write the cumulative infections and the infections that have
# completed the detection delay as
#
# ```math
# C_t = \sum_{u \le t} I_u, \qquad
# \text{det}_t = \sum_{s \ge 0} I_{t-s}\,
#     (f_{\text{inc}} * f_{\text{det}})_s.
# ```
#
# Then the daily export intensity and its running sum are
#
# ```math
# \lambda_t = p_{\text{Uganda}}\, q\, (C_t - \text{det}_t), \qquad
# \Lambda(t) = \sum_{u \le t} \lambda_u. \tag{37}
# ```
#
# We model outbound travel only, not return, so this term would overestimate
# the infections on its own. Each observed Uganda import is fitted at its
# reported detection date. An import detected on a given day is scored as a
# Poisson of the rise in cumulative export intensity between consecutive
# detection dates, with a term before the earliest detection $d_1$ observed
# at zero, since no export is expected then. After the last detection date we
# stop modelling exports rather than scoring further zeros: travellers'
# reasons for crossing the border change over the outbreak, so the baseline
# travel rate no longer applies beyond it and the export clock is truncated
# there:
#
# ```math
# Y_{\text{exports},i} \sim
#     \mathrm{Poisson}\!\bigl(\Lambda(d_i) - \Lambda(d_{i-1})\bigr),
# \qquad
# 0 \sim \mathrm{Poisson}\!\bigl(\Lambda(d_1 - 1)\bigr). \tag{38}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: exports_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_model(
#md #         missing, Float64[], 0.25; incubation_pmf = Float64[])), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Deaths among exports
#
# The expected deaths among exports weight the travelled at-risk prevalence
# by the infection-to-death delay (the onset-to-death PMF convolved with the
# incubation period) and scale by the CFR.
# The travelled prevalence is the export prevalence before the ascertainment
# factor $p_{\text{Uganda}}$, because a death among an exported case would be
# reported whether or not the case itself was ascertained as an import.
# Writing it $\ell_t = q\,(C_t - \text{det}_t)$, the daily export-death
# intensity is
#
# ```math
# \mu_t = \mathrm{CFR} \sum_{s \ge 0} \ell_{t-s}\, (f_{\text{inc}} * f_d)_s.
# ```
#
# Its running sum is the cumulative export-death intensity:
#
# ```math
# \Lambda_d(t) = \sum_{u \le t} \mu_u. \tag{39}
# ```
#
# Each dated Uganda export death is scored at its reported date with a
# per-day Poisson, the same dated-event likelihood the exports use, with a
# zero term before the first death day $\delta_1$:
#
# ```math
# Y_{\text{exp-deaths},i} \sim
#     \mathrm{Poisson}\!\bigl(\Lambda_d(\delta_i)
#     - \Lambda_d(\delta_{i-1})\bigr),
# \qquad
# 0 \sim \mathrm{Poisson}\!\bigl(\Lambda_d(\delta_1 - 1)\bigr). \tag{40}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: exports_deaths_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_deaths_model(
#md #         missing, Float64[], 0.33, Float64[], Float64[])), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Symptom-onset reporting delay
#
# Every other observation model sees the shared onset series only after a
# further convolution: a suspected-case report, a death, or a laboratory
# confirmation. The digitised onset epidemic curve (the [Data](@ref
# methods-data) section) is the only direct observation of it, so it can
# identify things the other streams cannot on their own, plausibly
# including the split between reporting and laboratory receipt that the
# laboratory pipeline currently pins with an external constraint.
#
# We model the onset-to-report delay as a discrete-time hazard over delay
# $d = 0,\dots,D-1$ days, $D = 28$: by then the triangle's between-vintage
# increments have decayed into digitisation noise. The baseline hazard is
# a non-centred logit random effect over the delay, free to rise and fall
# rather than forced monotone or parametric:
#
# ```math
# \eta_0 \sim \mathrm{Normal}(\mathrm{logit}(0.13),\ 0.7), \qquad
# \sigma_{h0} \sim \mathrm{Normal}^{+}(0,\ 1), \qquad
# \mathrm{logit}\,h_0(d) = \eta_0 + \sigma_{h0}\,z_{h0,d}. \tag{41}
# ```
#
# The hazard at delay $d$ for an onset on day $u$ is modified by a
# calendar-time effect indexed on the report day $u + d$: a weekly-knot
# non-centred random walk on the logit scale, the same construction as the
# reproduction-number walk above and concentrated near zero
# ($\sigma_\gamma \sim \mathrm{Normal}^{+}(0,\ 0.3)$), so a flat reporting
# profile stays the default the data has to argue away from while the walk
# can still follow a real drift in reporting speed:
#
# ```math
# \gamma_t = \mathrm{interp}\Bigl(\sigma_\gamma \sum_{s < k} z_{\gamma,s}\Bigr),
# \qquad
# h(d, t) = \mathrm{logistic}\bigl(\mathrm{logit}\,h_0(d) + \gamma_t\bigr).
# \tag{42}
# ```
#
# The cumulative reported proportion of onset date $u$'s eventual cases,
# reported within $\delta$ days, is the survival product of the daily
# hazards along that onset date's diagonal:
#
# ```math
# F(u, \delta) = \begin{cases} 0 & \delta < 0 \\
#     1 - \prod_{j=0}^{\min(\delta,\, D-1)} \bigl(1 - h(j, u + j)\bigr)
#     & \delta \ge 0. \end{cases} \tag{43}
# ```
#
# $\delta < 0$ is right truncation: a recent onset date's count in a given
# snapshot is only what its own delay allows, and rises as later snapshots
# see more of it, the idea EpiNow2 uses for its nowcast. $F$ deliberately
# does not reach $1$ as $\delta \to D-1$. The hazard's own asymptote
# absorbs ascertainment instead, so $F(u, D-1)$ is onset date $u$'s
# modelled ascertainment and $F(u,\delta) / F(u, D-1)$ its normalised delay
# profile, one mechanism read two ways rather than a second parameter block
# that would trade against the hazard total.
#
# The expected reported count is the onset series convolved with $F$,
# $\mathbb E[N(u, R_s)] = \mathrm{onsets}_u \cdot F(u, R_s - u)$. The
# likelihood scores the difference between consecutive snapshots at each
# onset date, in a trailing $D$-day window of the newer snapshot's report
# day, which avoids double-counting a case already reported earlier and
# drops the older onset dates that carry only noise by then. A count
# likelihood cannot be used, since a re-dated case can move a bar down in a
# later scan even though the true running total cannot fall, so the
# increment is scored with a Student-$t$ at fixed degrees of freedom
# ($\nu = 4$, a standard robust-regression choice):
#
# ```math
# y_u \sim \mathrm{Student}\text{-}t\Bigl(
#     \mathrm{onsets}_u\bigl(F(u, R_s{-}u) - F(u, R_{s-1}{-}u)\bigr),\
#     \sigma_u,\ \nu{=}4\Bigr). \tag{44}
# ```
#
# The likelihood therefore admits a negative increment, but the mean above
# cannot produce one: $F$ is non-decreasing in $\delta$, so the modelled
# increment is bounded below at zero. Re-dating is absorbed as observation
# noise rather than modelled, which issue #518 records along with how much
# of the movement between vintages is downward.
#
# $\sigma_u$ collects three sources: counting variation around the cell's
# own modelled mean; a $\pm 2.1$-case pixel-noise SD on the digitised bar,
# doubled for a correction since it differences two reads; and a $4.0\%$
# level error on each scan's own cumulative reading. Every magnitude
# entering $\sigma_u$ is the modelled one and never the observed count, so
# the likelihood's noise cannot feed into its own variance, and a sampled
# slack multiplier sits on top. That multiplier can only inflate the scale,
# because each term is a lower bound on the truth: a bar cannot be read off
# a raster more precisely than its pixels allow, and a count of newly
# reported cases carries at least its own counting variation.
#
# The first scored snapshot is differenced against an implicit empty
# predecessor, so its cells score levels rather than corrections. That is
# what anchors the asymptote: corrections only ever pin differences of $F$,
# so scaling $F$ up while scaling the onset series down leaves every
# correction cell unchanged, and something has to score a level.
#
# Three time-varying objects act on the same latent series: the
# reproduction-number walk on the onset axis, the calendar walk on the
# report axis, and ascertainment, which is the hazard's asymptote rather
# than a free object. Each leaves a different footprint on the grid of
# scored cells, which runs over onset date and snapshot pair. The onset
# series moves a column, the calendar walk moves a row, and the baseline
# hazard moves a diagonal band, and the three are distinguishable once
# there is more than one snapshot. That is the structural reason this
# stream is worth fitting, and also the reason the calendar walk is
# indexed on the report date: on the onset axis its footprint would
# coincide with the reproduction-number walk's and neither would be
# identified.
#
# Three things stay weak. The asymptote is pinned by the first snapshot's
# level cells and by the onset series the other streams supply, so in the
# onsets-only fit below it is confounded with the outbreak size and that
# fit's $C_T$ is close to prior-driven. The hazard below about two days'
# delay is barely observed, for the axis reason given in the Data section,
# and rests on pooling across delays. And a falling asymptote and a
# slowing hazard both suppress recent bars, where truncation self-corrects
# for the delay explanation but not for an ascertainment fall.
#
# The triangle counts confirmed cases by onset date, so its asymptote
# estimates the same quantity the confirmed pipeline builds from
# ascertainment, the tested fraction and test positivity. The two are
# deliberately not tied together, so that a digitised figure cannot pull
# on the ascertainment every other stream shares. The alive and dead split
# the raw figure carries is not modelled separately, since the
# confirmed-death stream already carries it from other data. An earlier
# line-list-independent reanalysis of this triangle put the median
# onset-to-report delay at around 6 days and the 7-day reporting fraction
# at 54-62%, an interval that wide because the digitisation noise is close
# in size to the increments the estimate rests on.

#md # ```@raw html
#md # <details><summary>Submodel: onset_report_hazard_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.onset_report_hazard_model(1, 28)),
#md #     "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: onset_reporting_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.onset_reporting_model(
#md #         (; onset_days = Int[], report_days = Int[],
#md #             prev_report_days = Int[], increments = Int[]),
#md #         Float64[])), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# #### Joint model
#
# The joint model runs the generating infection process once, stages it to
# daily onset incidence, and routes the shared onsets into every observation
# stream. It samples a single dispersion $k$ and the pooled ascertainment
# fractions, threading $p_{\text{DRC}}$ to the suspected-case, laboratory and
# confirmed-death likelihoods and $p_{\text{Uganda}}$ to the two Uganda-side
# likelihoods, and adds the genetic seeding bound on the outbreak age. Each
# observation stream argument may be dropped, so the same model structure
# generates prior- and posterior-predictive draws.
#
# The symptom-onset reporting triangle is threaded in the same way, as a
# standard stream rather than an opt-in one: the onset-curve input defaults
# to an empty history, so a missing input file degrades to a no-op rather
# than an error, but the production path fits it every time alongside the
# other streams.
#
# Alongside the joint model we write single-stream models for each
# count-based stream (exported cases, suspected deaths, suspected cases,
# laboratory-confirmed cases, confirmed deaths, deaths among exports and
# the symptom-onset reporting triangle), so each stream's posterior over
# the outbreak size can be compared with the joint. Other model variants
# reuse these models with different amounts of data, cutting the data to
# an earlier date or dropping the counts.

#md # ```@raw html
#md # <details><summary>Composer: exports-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_only_model(1, 1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Composer: deaths-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.deaths_only_model(1, 1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Composer: cases-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.cases_only_model(1, 1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Composer: confirmed-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.confirmed_only_model(1, 1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Composer: onsets-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.onsets_only_model(1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Composer: exports-deaths-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_deaths_only_model(1, 1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Composer: joint fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.bvd_joint(1, 1, 1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ### Model fitting and evaluation
#
# #### Prior predictive check
#
# Before any observation is taken into account, what does the prior
# imply about replicated exports, deaths and reported cases? Draws from
# the prior over the unobserved data should bracket the observed counts.

#md # ```@raw html
#md # <details><summary>Sample the joint prior</summary>
#md # ```

prior_chn = let
    breakpoint = obs.n - obs.who_first_sitrep_days
    m = bvd_joint(obs.n, missing, missing, missing, missing, missing;
        deaths_history = (; days = Int[], counts = Int[]),
        reported_history = (; days = Int[], counts = Int[]),
        confirmed_history = (; days = Int[], counts = Int[]),
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = breakpoint,
        background_re = true,
        confirmed_positivity_link = :composition,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days)
    sample(m, Prior(), 2_000; progress = false)
end;

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

prior_pair_fig = plot_pair(prior_chn,
    [:C_T, :R_T, :r, :T, :CFR, :k,
        :p_drc, :p_uganda]);

#md # ```@raw html
#md # </details>
#md # ```

prior_pair_fig #hide

# #### Fitting the models
#
# We sample with NUTS [hoffman2014nuts](@cite) and Mooncake
# [mooncake_jl](@cite) reverse-mode automatic differentiation, running two
# chains of 1000 post-warmup draws each after 1000 warmup adaptation steps,
# at a target acceptance probability of 0.85. Chains initialise from the
# prior. We fit the joint model and each single-stream model so the
# per-stream posteriors over the outbreak size can be compared with the
# joint.

# #### Fit diagnostics
#
# Fit-quality diagnostics for the joint and per-stream fits: the worst
# R-hat, the smallest bulk effective sample size, and the number of
# divergent transitions.

#md # ```@raw html
#md # <details><summary>Fit diagnostics</summary>
#md # ```

diagnostics_table( #hide
    "joint" => chn_joint, #hide
    "exports" => chn_exports, #hide
    "deaths (DRC)" => chn_deaths, #hide
    "cases (DRC)" => chn_cases, #hide
    "confirmed (DRC)" => chn_confirmed, #hide
    "confirmed deaths (DRC)" => chn_confirmed_deaths, #hide
    "isolation (DRC)" => chn_treatment, #hide
    "onsets (DRC)" => chn_onsets, #hide
    "frozen (1wk back)" => frozen_lastweek.chn, #hide
    (RUN_SENSITIVITY ? #hide
     ["delay sensitivity" => chn_joint_community_delay, #hide
        "clock sensitivity (ExpGrowth)" => chn_joint_exp_growth_clock] : [])...) #hide

#md # ```@raw html
#md # </details>
#md # ```

# #### No-onward-transmission counterfactual
#
# To bound the deaths already committed at the cut-off, we project the
# deaths that would still occur if all transmission stopped on the report
# date. Every infection present by the cut-off still dies with probability
# CFR, so the committed future deaths are the CFR-weighted cumulative
# infection count net of the deaths already expected,
# $\Delta D = \mathrm{CFR}\cdot I_T - \mathbb{E}[D_T]$, where $I_T$ is the
# cumulative infection count to the cut-off. The figure is shown in the
# counterfactual results below.
#
# #### Delay-corrected confirmed case-fatality ratio
#
# The case-fatality ratio above is the onset-level CFR, the share of
# symptomatic infections that die.
# It is hard to read directly off the data because the case and death streams
# are ascertained differently, so a reader who wants a figure anchored in the
# observed counts is left with the naive confirmed ratio, the cumulative
# confirmed deaths over the cumulative confirmed cases.
# That naive ratio is biased low in real time.
# A case confirmed close to the cut-off has not yet had time to die, so it
# enters the denominator before it can enter the numerator.
#
# We report a delay-corrected confirmed CFR that debiases the real-time ratio
# following [nishiura2009](@cite).
# The denominator is shrunk from all confirmed cases to those expected to have
# had their death confirmed by the cut-off, weighting each day of
# confirmed-case incidence by the probability that a case confirmed that day,
# if it is going to die, has had its death confirmed by the cut-off:
#
# ```math
# \mathrm{cCFR}_{\text{corr}}(T) =
#   \frac{D_{\text{conf}}(T)}
#        {\sum_{t} c_{\text{conf}}(t)\,
#         \Pr(X_d - X_c \le T - t)}, \tag{45}
# ```
#
# with $D_{\text{conf}}(T)$ the cumulative confirmed deaths, $c_{\text{conf}}(t)$
# the modelled daily confirmed-case incidence, and $X_d - X_c$ the residual
# delay between a confirmed case and its confirmed death.
# $X_d$ is the onset-to-death-confirmation lag (onset-to-death convolved with
# the report-to-receipt laboratory delay) and $X_c$ the onset-to-confirmation
# lag (onset-to-report convolved with the same laboratory delay), so the common
# receipt delay cancels in the mean and the residual centres on onset-to-death
# minus onset-to-report.
# Both lags and the confirmed trajectories are taken per posterior draw from
# the joint fit, so the corrected ratio carries the joint uncertainty.
# As the outbreak matures and recent incidence resolves the correction shrinks
# and the corrected ratio approaches the eventual confirmed CFR.
# It is the confirmed-case counterpart of the structural CFR, anchored in the
# confirmed counts rather than the latent infections, and the gap between the
# two reflects the difference in case and death ascertainment the structural
# CFR has to absorb.
# The result is shown in the
# [confirmed case-fatality ratio results](@ref "Confirmed case-fatality ratio")
# below.
#
# #### One-week-ahead forecast
#
# We project each DRC stream seven days beyond the cut-off, letting the
# reproduction number keep evolving over the horizon by continuing the
# recent trend of its trajectory rather than holding it fixed, with no
# further interventions and no saturation imposed.
# The projection carries both parameter and observation uncertainty.
# We forecast the DRC observation streams as forecast targets: the reported
# cases and suspected deaths, the laboratory-confirmed cases and confirmed
# deaths, the isolation/treatment beds and the recovered total. For the beds
# we project both the bed demand (the need a week ahead, under unconstrained
# supply, the cut-off demand grown by the horizon factor like the case
# inflow) and the supply-limited occupancy that demand produces against the
# bed capacity. The gap between them is the projected bed shortfall, the
# quantity of interest if bed occupancy is supply-constrained. The reported
# case and suspected death streams are no longer published, so their
# forecasts extend the last published cumulative total rather than a
# still-growing series. Exports are not forecast, since cross-border travel is
# unlikely to continue at its baseline rate, so the forward travel rate the
# export model relies on no longer holds. The figure is shown in the
# [one-week-ahead forecast results](@ref "One-week-ahead forecast results")
# below.
#
# ##### Symptom-onset nowcast and forecast
#
# The symptom-onset stream is projected differently from every stream above,
# because the reporting triangle lets us separate two things the other
# streams cannot tell apart: cases whose symptoms have already begun but
# whose report has not yet arrived, and cases whose symptoms have not begun
# at all. The first is a nowcast, the second a forecast, and a count of
# "cases still to come" that mixes them is not interpretable.
#
# The separation comes from the same cumulative reported proportion
# $F(u, \delta)$ the likelihood is built on
# (see [symptom-onset reporting delay](@ref "Symptom-onset reporting
# delay")). The expected reported total as of day $a$ is
#
# ```math
# S(a) = \sum_{u \le a} \text{onsets}_u\, F(u,\, a - u),
# ```
#
# so at the cut-off $T$ the triangle should have printed $S(T)$ while
# $\sum_{u \le T} \text{onsets}_u$ onsets have actually happened. Their
# difference is the nowcast, onsets already in the population but not in
# the figure. Not all of it will ever be reported, since $F$ carries
# ascertainment as well as delay, so it is the reporting backlog and the
# never-ascertained cases together.
#
# Splitting $S(T + h) - S(T)$ at the cut-off splits the coming week the same
# way:
#
# ```math
# \underbrace{\sum_{u \le T} \text{onsets}_u \bigl(F(u,\, T + h - u) -
#     F(u,\, T - u)\bigr)}_{\text{already happened, reported this week}}
# \; + \;
# \underbrace{\sum_{T < u \le T+h} \text{onsets}_u\, F(u,\, T + h -
#     u)}_{\text{not yet happened}} .
# ```
#
# Onsets past the cut-off are projected under the same evolving growth-rate
# path the other streams use, and the calendar-time effect $\gamma$ is held
# flat at its last fitted value across the horizon, the assumption any
# nowcast makes about reporting behaviour continuing.
#
# We score the sum of the two terms, the increment the triangle should add
# over the horizon, rather than its cumulative level. Every vintage rereads
# the whole figure, so the printed total moves with the roughly 4% per-scan
# level error as well as with genuine late reporting, and it falls between
# consecutive vintages more than once in the current data. Scoring the level
# would charge the forecast for a rescan of cases it had already predicted
# and would count the same revision again at every later horizon.
#
# The interval on that increment is dominated by the scan error rather than
# by epidemic uncertainty, so the forecast is worth more as a check that the
# fitted delay and ascertainment reproduce the next vintage than as a
# case-count prediction.
#
# Each release now saves its forecast as an asset so it can later be scored
# against what is observed.
# Earlier releases showed a forecast but did not store it, so those forecasts
# are reconstructed.
# A backfill script checks out each past release at its own tag,
# re-runs that release's own model code on that release's data through its own
# fit, and writes the forecast in the same archive schema, so a reconstructed
# forecast is the release's own output rather than a current-code
# approximation.
# The reconstruction re-resolves dependencies at current versions, since the
# release manifests were not pinned, so it reproduces the release's model and
# data but not its exact solver build.
# Reconstruction reaches back to the first release that carried any forecast,
# so it covers the whole release history, but the streams available grow over
# time as the model and data did.
# The renewal releases from v1.4.0 reconstruct the incident case and death
# streams, extending to all four streams from v1.6.0 once the recovered and
# isolation series entered the data.
# The integral-era v1.3.0 reconstructs the confirmed case and death streams.
# The earliest releases, v1.0.0 to v1.2.0, forecast the reported case,
# suspected death, and export streams, and are reconstructed from each tag's
# own inline model code.
# Reconstructed forecasts are published as a separate backfill release and
# scored alongside the stored ones by the
# [forecast scoring across releases](@ref "Forecast scoring across releases")
# section.
#
# #### Forecast-versus-frozen evaluation
#
# We assess the forecast against data observed since by freezing the data to
# roughly one week before the current cut-off, re-fitting, and projecting one
# week ahead with the same forecast machinery, then comparing that projection
# against the counts observed by the current cut-off. The frozen re-fit cuts
# the data to an earlier cut-off and re-fits the joint model, so that a change
# driven by newer data can be distinguished from one driven by a change of
# method. Each frozen re-fit uses the
# full headline settings (1000 draws across two chains). The same frozen re-fit is
# reused to compare against McCabe et al. at the cut-offs they used. The helper
# below performs one frozen joint re-fit and is reused by the forecast
# validation and matched-in-time results.

#md # ```@raw html
#md # <details><summary>Frozen-fit helper (reused by the forecast validation and matched-in-time sections)</summary>
#md # ```

## The frozen re-fits are defined in the fit registry (`docs/fits/registry.jl`) and
## loaded through the cache in the setup block above.

#md # ```@raw html
#md # </details>
#md # ```

# #### Forecast scoring against a persistence baseline
#
# Every forecast above, and every past release's saved forecast, is scored
# with the continuous ranked probability score (CRPS), on both the raw and
# the log scale (which downweights the largest counts rather
# than treating an error of a hundred cases the same at every scale).
# The CRPS decomposes into the ensemble's own spread (dispersion) and the
# extra cost of over- or under-predicting the observation, so a fit
# reading consistently high carries more overprediction than
# underprediction and the reverse for a fit reading low.
# Coverage is the share of forecasts whose 50% or 90% predictive interval
# contained the observation, which should sit near 0.5 and 0.9 for a
# well-calibrated fit; bias is signed, running from -1 (every forecast too
# low) to 1 (every forecast too high).
# A relative skill is the ratio of two fits' mean CRPS over the matched
# set of forecasts both scored, not a mean of per-forecast ratios, so a
# comparator that happened to score zero on one forecast cannot make the
# ratio infinite; a ratio below one beats the comparator.
#
# The count streams are running cumulative totals, so each is scored on
# its incidence, the increment between one vintage and the next, rather
# than the cumulative level; the isolation-bed occupancy is a level and is
# scored as one.
# A stream is scored only over the period its own reporting source
# covers, from the day it was first reported to the day it was last
# updated: outside that period an unmoved cumulative total is the absence
# of a series rather than an observed zero, and scoring against it would
# reward persistence for predicting a number nobody measured.
#
# Every scored forecast is also compared against a persistence baseline:
# the last observed occupancy carried forward unchanged for the isolation
# bed level, matching the COVID-19 Forecast Hub baseline
# (`COVIDhub-baseline`) directly, and for a count stream's own increments
# the observed count over the whole horizon-length window ending at the
# forecast's cut-off, rather than a single day's change; the streams here
# are sparse, irregularly-timed cumulative counts, so a single most recent
# increment would be a noisier centre than that window pooled.
# Its predictive uncertainty comes from an iterated random walk, day by
# day out to the forecast horizon, each day's step drawn from how much
# that same quantity has moved historically in the stream's own reporting
# record, so the interval widens with the square root of the horizon.
# The baseline is built from the release's own data vintage, frozen no
# later than the day the forecast was made, so a correction or backfill
# that landed after a forecast was made cannot leak into its baseline and
# make the comparison unfair.
# That guarantee is exact for forecasts scored at a release's own cut-off,
# where the vintage snapshot and the forecast's made date are the same
# day or close to it.
# For the frozen re-fit evaluation below, the made date is a fixed historical
# cut-off reused across many later release tags, so it can sit weeks
# before the snapshot's own cut-off; a correction that landed between
# the made date and the snapshot's cut-off is already baked into the
# snapshot and can still reach the frozen baseline.
# Closing that residual gap would need a made-date-specific vintage
# manifest for each frozen cut-off, which is not archived.

# ## Results
#
# ### Summary
#
# The numbers below are our estimate of the underlying infections to date,
# reported and unreported, from the joint posterior.
# Each is given as equal-tailed 30%, 60% and 90% credible intervals.

#md # ```@raw html
#md # <details><summary>Compute the headline ranges</summary>
#md # ```

summary_ranges = let
    med(x) = quantile(x, 0.5)
    iqr(x) = quantile(x, 0.75) - quantile(x, 0.25)
    ## Posterior-minus-prior shift in units of the parameter's prior IQR,
    ## reusing the prior draws so nothing is respecified here.
    shift(post, prior) = round((med(post) - med(prior)) / iqr(prior);
        digits = 2)

    C = posterior_C_joint
    Td = vec(Array(chn_joint[:T]))
    r0d = vec(Array(chn_joint[:r0]))
    rd = vec(Array(chn_joint[:r]))
    dt0 = log(2) ./ r0d
    dt = vec(Array(chn_joint[:doubling_time]))
    R0d = vec(Array(chn_joint[:R0]))
    RTd = vec(Array(chn_joint[:R_T]))
    cfrd = vec(Array(chn_joint[:CFR]))
    sC = posterior_summary(C)
    sT = posterior_summary(Td)
    sr0 = posterior_summary(r0d)
    sr = posterior_summary(rd)
    sdt0 = posterior_summary(dt0)
    sdt = posterior_summary(dt)
    sR0 = posterior_summary(R0d)
    sRT = posterior_summary(RTd)
    scfr = posterior_summary(cfrd)

    ints_i(s) = string(
        "30% ", round(Int, s.lo30), "–", round(Int, s.hi30),
        ", 60% ", round(Int, s.lo60), "–", round(Int, s.hi60),
        ", 90% ", round(Int, s.lo90), "–", round(Int, s.hi90))
    ints_f(s,
        d) = string(
        "30% ", round(s.lo30; digits = d), "–", round(s.hi30; digits = d),
        ", 60% ", round(s.lo60; digits = d), "–", round(s.hi60; digits = d),
        ", 90% ", round(s.lo90; digits = d), "–", round(s.hi90; digits = d))
    start_from(t) = obs.cutoff - Day(round(Int, t))
    ints_d(s) = string(
        "30% ", start_from(s.hi30), "–", start_from(s.lo30),
        ", 60% ", start_from(s.hi60), "–", start_from(s.lo60),
        ", 90% ", start_from(s.hi90), "–", start_from(s.lo90))
    f_lo = round(sC.lo90 / obs.confirmed_cases; digits = 1)
    f_hi = round(sC.hi90 / obs.confirmed_cases; digits = 1)

    ## How far the data has moved each estimate from its prior, in prior
    ## interquartile ranges, reusing the prior draws.
    moves = [
        "cumulative infection count" => shift(C, vec(Array(prior_chn[:C_T]))),
        "outbreak age" => shift(Td, vec(Array(prior_chn[:T]))),
        "doubling time" => shift(dt, vec(Array(prior_chn[:doubling_time])))]
    biggest = argmax(p -> abs(p.second), moves)

    Markdown.parse("""
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
    """)
end;

#md # ```@raw html
#md # </details>
#md # ```

summary_ranges #hide

# ### Joint model estimates
#
# This section reports the joint posterior over the cumulative infection
# count to date, fitting every data stream together.

#md # ```@raw html
#md # <details><summary>Cumulative infection count summary table</summary>
#md # ```

cumulative_cases_summary = summary_table(
    chn_joint, [:C_T]; digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_cases_summary #hide

# The figure below shows the cumulative trajectories and current-cut-off
# densities for three latent quantities: infections, symptom onsets and deaths.
# The infection density is the headline outbreak size, a count of infections
# rather than reported cases.

#md # ```@raw html
#md # <details><summary>Cumulative infections, onsets and deaths figure</summary>
#md # ```

cumulative_traj_fig = plot_cumulative_trajectories(chn_joint;
    n = obs.n, seeding = obs.seeding);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_traj_fig #hide

# The cumulative infection count is set by the reproduction number trajectory
# and the outbreak age, the elapsed time from the import that started the
# outbreak to the cut-off.
# The left panel below shows the posterior for that start date; the right
# panel shows the joint posterior of the outbreak age and the early
# doubling time.

#md # ```@raw html
#md # <details><summary>Outbreak start date and seeding-time posterior</summary>
#md # ```

start_date_fig = plot_start_date_pair(chn_joint;
    as_of_date = string(obs.cutoff));

#md # ```@raw html
#md # </details>
#md # ```

start_date_fig #hide

# The summary table reports the credible intervals on the infection-process
# parameters: the growth rate and doubling time, the reproduction number,
# the outbreak age, the case-fatality ratio and the cumulative infection
# count.
# The pair plot beside it shows their joint distribution, with the prior
# overlaid so the data's contribution to each marginal is visible.

#md # ```@raw html
#md # <details><summary>Infection-parameter summary table</summary>
#md # ```

infection_summary = summary_table(chn_joint,
    [:r, :doubling_time, :T, :R_T, :CFR, :C_T]; digits = 2);

#md # ```@raw html
#md # </details>
#md # ```

infection_summary #hide

#md # ```@raw html
#md # <details><summary>Infection-parameter pair plot (prior overlaid)</summary>
#md # ```

infection_pair_fig = plot_pair(chn_joint,
    [:R_T, :r, :T, :CFR,
        Symbol("rt_state.sigma_rw"), Symbol("rt_state.intervention_effect")];
    prior = prior_chn, labels = display_names);

#md # ```@raw html
#md # </details>
#md # ```

infection_pair_fig #hide

# The infection model carries two delays: the generation interval, the time
# between an infector's and an infectee's onset that drives the renewal
# recursion, and the incubation period, the time from infection to symptom
# onset that turns infections into onsets.
# The table reports their posterior means and standard deviations; the pair
# plot beside it shows their joint posterior with the prior overlaid.

#md # ```@raw html
#md # <details><summary>Infection-delay summary table</summary>
#md # ```

infection_delay_summary = summary_table(chn_joint,
    [Symbol("gi_state.α"), Symbol("gi_state.θ"),
        Symbol("inc_state.delay_mean"), Symbol("inc_state.delay_sd")];
    digits = 2, labels = display_names);

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

infection_delay_pair_fig = plot_pair(chn_joint,
    [Symbol("gi_state.α"), Symbol("gi_state.θ"),
        Symbol("inc_state.delay_mean"), Symbol("inc_state.delay_sd")];
    prior = prior_chn, labels = display_names);

#md # ```@raw html
#md # </details>
#md # ```

infection_delay_pair_fig #hide

# ### Reproduction number over time
#
# The daily reproduction number over the period we estimate it for, the
# established outbreak from the genetic bound to the cut-off.
# The 30%, 60% and 90% credible ribbons are shown with about a hundred
# sampled trajectories, and the no-growth threshold at one as a grey dashed
# line.
# The first situation report on 18 May 2026 marks the start of the response
# scale-up (red dashed) and the end of the three-week scale-up is the red
# dotted line; the data cut-off is grey dashed.

#md # ```@raw html
#md # <details><summary>Reproduction-number trajectory</summary>
#md # ```

## `rt_start` is the renewal/established-window start the plot shows from;
## `rt_walk_start` is where the random walk's knots begin — `RT_WALK_LEAD`
## days (a month) before the first situation report, matching `bvd_joint`'s
## `rt_walk_lead` — so the chain reconstruction uses the same knot grid the
## model did, floored at the renewal start. R_t is flat at R0 between the two.
rt_fig = plot_rt(chn_joint;
    n = obs.n, breakpoint = _BREAKPOINT,
    rt_start = _rt_start_plot,
    rt_walk_start = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n),
    as_of_date = string(obs.cutoff), seeding = obs.seeding,
    ramp = RT_INTERVENTION_RAMP);

#md # ```@raw html
#md # </details>
#md # ```

rt_fig #hide

# The table reports the posterior of the response effect on the
# reproduction number as a multiplier, where a value below one is the factor
# by which the response lowers the reproduction number once the scale-up
# completes.

#md # ```@raw html
#md # <details><summary>Intervention-effect summary table</summary>
#md # ```

intervention_effect = vec(Array(
    chn_joint[Symbol("rt_state.intervention_effect")]));
intervention_table = streams_table(
    "Rt multiplier exp(effect)" => exp.(intervention_effect);
    digits = 2);

#md # ```@raw html
#md # </details>
#md # ```

intervention_table #hide

# ### Observation delays
#
# The delays carry latent infections through to each observed event: reporting,
# death, detection abroad and laboratory receipt.
# The onset-to-report and onset-to-detection delays are the same line-list
# onset-to-admission delay, sampled on its natural Gamma shape and scale, and
# onset-to-death is the convolution of two atomic Gamma delays, onset to
# admission and admission to death, each with its own shape and scale.
# The report-to-receipt delay is sampled by its mean and standard deviation.
# The length-of-stay delays are also shown: the isolation-bed BVD treatment
# length-of-stay — for how long an admitted BVD patient occupies a bed, with
# the line-list admission-to-death delay as its prior; the non-BVD rule-out
# stay — how long a ruled-out suspect occupies a bed before discharge, with the
# report-to-receipt turnaround as its prior; and the confirmation-to-recovery
# delay (how long after confirmation a case is recorded as recovered).
# The table reports their posteriors; the pair plot beside it shows their
# joint posterior with the prior overlaid, so the data's contribution to
# each marginal is visible.

#md # ```@raw html
#md # <details><summary>Observation-delay summary table</summary>
#md # ```

obs_delay_summary = summary_table(chn_joint,
    [Symbol("cases_state.report_state.α"),
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
        :recovery_delay_mean];
    digits = 2, labels = display_names);

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

obs_delay_pair_fig = plot_pair(chn_joint,
    [Symbol("cases_state.report_state.α"),
        Symbol("deaths_state.od_state.oa.α"),
        Symbol("exports_state.detect_state.α"),
        Symbol("confirmed_state.receipt_state.d.delay_mean"),
        :isolation_bvd_los_mean,
        :isolation_ruleout_los_mean,
        :recovery_delay_mean];
    prior = prior_chn, labels = display_names);

#md # ```@raw html
#md # </details>
#md # ```

obs_delay_pair_fig #hide

# ### Surveillance parameters
#
# The surveillance-data parameters: the reporting fractions for the DRC and
# Uganda, the surveillance dispersions, and the laboratory pipeline (the
# testing fraction and receipt delay, the per-suspected and per-test
# positivity, the non-BVD background rate, and the death-confirmation
# probability). The six passive-surveillance count streams (suspected
# cases, suspected deaths, confirmed cases, confirmed deaths, isolation
# occupancy and recovered) each have their own negative-binomial dispersion
# partially pooled from a shared population: $k$ is the population-level
# dispersion, $k_{\text{cases}}$, $k_{\text{deaths}}$, $k_{\text{confirmed}}$ and
# $k_{\text{confirmed deaths}}$ the per-stream values for the four DRC count
# streams, and a pooling spread. The isolation and recovered streams add the
# proportion of suspects admitted to a bed and the recovery probability among
# confirmed cases, with their dispersions ($k_{\text{iso}}$, $k_{\text{rec}}$)
# drawn from the same pooled population (see the length-of-stay delays in the
# observation-delay table above).
# The table reports their credible intervals; the pair plot beside it shows
# their joint posterior with the prior overlaid.

#md # ```@raw html
#md # <details><summary>Surveillance-parameter summary table</summary>
#md # ```

surveillance_summary = summary_table(chn_joint,
    [:p_drc, :p_uganda, :k, :k_cases, :k_deaths, :k_confirmed,
        :k_confirmed_deaths, :dispersion_sd, :tau_test, :lambda_bg,
        :suspected_positivity, :test_positivity, :expected_confirmed_T,
        :expected_analysed_T, :death_ascertainment, :background_cfr,
        :tau_death, :death_composition,
        :death_confirmation, :expected_confirmed_deaths_T,
        :isolation_admission, :isolation_dispersion, :expected_isolation_T,
        :expected_bed_demand_T, :bed_capacity, :bed_shortfall_T,
        :incare_cfr, :incare_cfr_modifier, :incare_confirm_modifier,
        :isolation_death_los_mean,
        :isolation_recovery_los_mean, :abscond_fraction,
        :recovery_probability, :recovered_dispersion, :expected_recovered_T];
    digits = 3);

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

surveillance_pair_fig = plot_pair(chn_joint,
    [:p_drc, :p_uganda, :k, :tau_test, :lambda_bg, :test_positivity,
        :death_confirmation];
    prior = prior_chn);

#md # ```@raw html
#md # </details>
#md # ```

surveillance_pair_fig #hide

# ### Symptom-onset reporting delay and ascertainment
#
# The table reports the onset-report hazard's hyperparameters together with
# two derived quantities: the share of a representative onset date's
# eventual reports that arrive within 7 days, and the median modelled
# ascertainment across the onset dates whose asymptote is observable within
# the fitted grid.
# Both come from the one fitted calendar effect, so they are read together
# rather than as independent estimates (see the [symptom-onset reporting
# delay](@ref "Symptom-onset reporting delay") Methods section).
# The pair plot shows the hyperparameters alone, with the prior overlaid.
#
# The scale slack row is a diagnostic rather than a quantity of interest.
# Its prior is bounded below at one, because each term in the observation
# scale is a lower bound on the truth, so a posterior sitting on that bound
# says the fit would like a tighter likelihood than the figures can support.
#
# The stream is not underdispersed, and the likelihood is left unchanged
# here: individual cells are already scored wider than they need to be,
# while the per-snapshot aggregate is mis-centred rather than mis-scaled.
# Issue #507 carries the numbers and why no variance term fixes it.

#md # ```@raw html
#md # <details><summary>Reconstruct the onset-report hazard and calendar walk</summary>
#md # ```

## The report-date grid the calendar walk spans is a fixed function of the
## digitised triangle, not sampled, so it is recomputed once directly from
## `obs.onset_curve_history` rather than pulled from the chain (mirroring
## `onset_reporting_model`'s own `grid_start`/`grid_end` construction).
_onset_grid_start = isempty(obs.onset_curve_history.onset_days) ? 1 :
                    minimum(obs.onset_curve_history.onset_days)
_onset_grid_end = isempty(obs.onset_curve_history.report_days) ?
                  _onset_grid_start :
                  max(maximum(obs.onset_curve_history.report_days),
    _onset_grid_start)

## Every posterior draw's `logit_h0` (the baseline delay hazard) and `γ`
## (the report-date calendar walk), rebuilt from the non-centred
## innovations the chain stores. `reconstruct_onset_hazard` is the package
## function the onset forecast also uses, so the hazard plotted here and
## the one projected forward are the same object rather than two copies of
## the same reconstruction that could drift apart.
_onset_hazard = reconstruct_onset_hazard(chn_joint;
    grid_start = _onset_grid_start, grid_end = _onset_grid_end)

## A representative onset day (the median scored onset date), so the 7-day
## fraction below reflects a typical, not an edge, calendar day.
_onset_u_ref = isempty(obs.onset_curve_history.onset_days) ?
               _onset_grid_start :
               round(Int, quantile(obs.onset_curve_history.onset_days, 0.5))
## The delay profile is the cumulative F normalised by its own asymptote,
## since F itself is deliberately not made to reach 1 and so carries
## ascertainment as well as delay. The two are reported separately below.
_onset_7d_fraction = [begin
                          D = length(_onset_hazard.logit_h0[i])
                          f7 = onset_report_cdf(6,
                              _onset_hazard.logit_h0[i],
                              _onset_hazard.γ[i], _onset_u_ref,
                              _onset_grid_start)
                          fa = onset_report_cdf(D - 1,
                              _onset_hazard.logit_h0[i],
                              _onset_hazard.γ[i], _onset_u_ref,
                              _onset_grid_start)
                          fa > 0 ? f7 / fa : NaN
                      end
                      for i in eachindex(_onset_hazard.logit_h0)]
filter!(isfinite, _onset_7d_fraction)

## Per-draw median ascertainment over every onset date whose asymptote is
## observable within the fitted grid (`onset_report_ascertainment`), then
## summarised across draws: the one fitted calendar effect read as an
## ascertainment level rather than a delay profile (see the identifiability
## discussion in the Methods section above).
_onset_ascertainment_draws = Float64[]
for i in eachindex(_onset_hazard.logit_h0)
    a = onset_report_ascertainment(_onset_hazard.logit_h0[i],
        _onset_hazard.γ[i], _onset_grid_start, _onset_grid_end)
    isempty(a) || push!(_onset_ascertainment_draws, quantile(a, 0.5))
end

_onset_labels = merge(display_names,
    Dict(
        Symbol("onset_report_state.η0") => "onset-report hazard baseline (logit)",
        Symbol("onset_report_state.σ_h0") => "onset-report hazard pooling SD",
        Symbol("onset_report_state.σ_γ") => "onset-report calendar-walk step size",
        Symbol("onset_report_state.σ_mult") => "onset-report scale slack"));

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
    upper_60 = Float64[], upper_90 = Float64[])
for (label, draws) in [
    ("share of eventual reports arriving within 7 days (median onset date)",
        _onset_7d_fraction),
    ("median modelled ascertainment", _onset_ascertainment_draws)]
    s = posterior_summary(draws)
    push!(onset_derived_raw,
        (label, round(s.lo90; digits = 3), round(s.lo60; digits = 3),
            round(s.lo30; digits = 3), round(s.hi30; digits = 3),
            round(s.hi60; digits = 3), round(s.hi90; digits = 3)))
end
onset_derived_table = BVDOutbreakSize._prettify(onset_derived_raw)

onset_summary = vcat(
    summary_table(chn_joint,
        [Symbol("onset_report_state.η0"), Symbol("onset_report_state.σ_h0"),
            Symbol("onset_report_state.σ_γ"),
            Symbol("onset_report_state.σ_mult")];
        digits = 3, labels = _onset_labels),
    onset_derived_table);

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

#md # ```@raw html
#md # <details><summary>Symptom-onset reporting-delay pair plot (prior overlaid)</summary>
#md # ```

onset_pair_fig = plot_pair(chn_joint,
    [Symbol("onset_report_state.η0"), Symbol("onset_report_state.σ_h0"),
        Symbol("onset_report_state.σ_γ"),
        Symbol("onset_report_state.σ_mult")];
    prior = prior_chn, labels = _onset_labels);

#md # ```@raw html
#md # </details>
#md # ```

onset_pair_fig #hide

# Each panel below is one digitised snapshot, and each point is one onset
# date's bar as that snapshot printed it, against the model's count for the
# same onset date and the same reporting delay.
# The dark band is the modelled count itself and the pale band adds the
# measurement error the likelihood gives a digitised bar, so the pale band
# is the one the points should fall inside: they do for 95% of cells,
# against 42% for the dark band alone.
# A recent onset date sits below its eventual value in its own snapshot's
# panel and catches up in a later panel, the right-truncation behaviour the
# model relies on (see the [symptom-onset reporting delay](@ref
# "Symptom-onset reporting delay") Methods section).
#
# The panels also show two misses the fitted delay does not absorb: it is
# too slow between 6 and 19 days, and the miss grows across the snapshots.
# Issue #517 carries both.

#md # ```@raw html
#md # <details><summary>Fits to the digitised reporting-triangle snapshots</summary>
#md # ```

## The digitised, deduplicated, cut-off-filtered snapshot blocks
## `load_onset_curve` scores, kept here for their raw cumulative
## onset-date counts: the fitted stream only ever sees between-vintage
## increments, so the observed cumulative levels this figure plots are
## read back from the source blocks directly rather than reconstructed
## from the fitted increments.
_onset_path = joinpath(pkgdir(BVDOutbreakSize), "data",
    "onset_curve_scanned.csv")
_onset_snaps = filter(b -> b.report_date <= obs.cutoff,
    BVDOutbreakSize._dedup_onset_blocks(
        BVDOutbreakSize._read_onset_curve_blocks(_onset_path)))
## Keyed by report day rather than kept in order: a snapshot whose printed
## extent misses the scored window contributes no cells, so the panels and
## the snapshot blocks are not guaranteed to line up positionally.
_onset_snap_by_day = Dict(
    obs.n - value(obs.cutoff - b.report_date) => b for b in _onset_snaps)

## Cell indices grouped by their snapshot's report day, and the daily
## onsets series per posterior draw (the diff of the chain's stored
## `cumulative_onsets` trajectory; the model does not store the daily
## series itself, only its running sum).
_onset_cells_by_report = Dict{Int, Vector{Int}}()
for (i, r) in enumerate(obs.onset_curve_history.report_days)
    push!(get!(_onset_cells_by_report, r, Int[]), i)
end
_onset_report_grid_days = sort(collect(keys(_onset_cells_by_report)))
_onset_daily_draws = [vcat(v[1], diff(v))
                      for v in vec(collect(chn_joint[:cumulative_onsets]))]

## Modelled count at onset day `u` as of report day `R`, per posterior
## draw: `onsets[u] * F(u, R - u)`, the same expected value
## `onset_report_expected_total` sums, evaluated at a single onset day.
function _onset_modelled_cumulative(u::Integer, R::Integer)
    return [_onset_daily_draws[i][u] *
            onset_report_cdf(R - u, _onset_hazard.logit_h0[i],
                _onset_hazard.γ[i], u, _onset_grid_start)
            for i in eachindex(_onset_daily_draws)]
end

## The same counts put through the stream's own observation model, so the
## band is a posterior predictive of a digitised bar rather than of the
## latent count behind it. A bar is one read off one scan with no previous
## level to difference against, which is `onset_report_scale`'s level case
## (`reads = 1`, `level_prev = 0`) — the case the first scored snapshot's
## own cells carry. Without this the band is the modelled count alone and
## covers 42% of the observed bars at a nominal 90%.
_onset_σ_mult = vec(collect(chn_joint[Symbol("onset_report_state.σ_mult")]))
_onset_ppc_rng = Random.MersenneTwister(20260729)
## Four replicates per draw rather than one: the band is a 90% interval of
## a heavy-tailed replicate, and at one per draw its edge is visibly ragged
## from Monte Carlo error alone.
function _onset_replicated(draws::AbstractVector)
    return [begin
                μ = draws[i]
                σ = _onset_σ_mult[i] * onset_report_scale(μ, μ, 0.0, 1)
                μ + σ * rand(_onset_ppc_rng, TDist(4.0))
            end
            for _ in 1:4 for i in eachindex(draws)]
end

onset_fit_fig = let
    ncol = 3
    nrow = cld(length(_onset_report_grid_days), ncol)
    fig = CairoMakie.Figure(; size = (330 * ncol, 260 * nrow + 60))
    for (k, R) in enumerate(_onset_report_grid_days)
        r, c = fldmod1(k, ncol)
        snap = _onset_snap_by_day[R]
        idx = sort(_onset_cells_by_report[R];
            by = i -> obs.onset_curve_history.onset_days[i])
        us = obs.onset_curve_history.onset_days[idx]
        observed = Float64[get(snap.onsets, grid_date(u), 0) for u in us]
        draws = [_onset_modelled_cumulative(u, R) for u in us]
        reps = [_onset_replicated(d) for d in draws]
        med = [quantile(d, 0.5) for d in draws]
        lo90 = [quantile(d, 0.05) for d in draws]
        hi90 = [quantile(d, 0.95) for d in draws]
        plo90 = [quantile(d, 0.05) for d in reps]
        phi90 = [quantile(d, 0.95) for d in reps]
        xs = Float64.(1:length(us))
        ax = CairoMakie.Axis(fig[r, c]; title = string(snap.report_date),
            xlabel = r == nrow ? "onset day (oldest to newest)" : "",
            ylabel = c == 1 ? "cases at this onset date" : "")
        CairoMakie.band!(ax, xs, plo90, phi90; color = (:steelblue, 0.12))
        CairoMakie.band!(ax, xs, lo90, hi90; color = (:steelblue, 0.30))
        CairoMakie.lines!(ax, xs, med; color = :steelblue, linewidth = 2,
            label = "modelled")
        CairoMakie.scatter!(ax, xs, observed; color = :black,
            markersize = 6, label = "digitised")
    end
    CairoMakie.Label(fig[0, 1:ncol],
        "Symptom-onset reporting triangle: fitted vs digitised";
        font = :bold, tellwidth = false)
    CairoMakie.Legend(fig[nrow + 1, 1:ncol],
        [CairoMakie.LineElement(color = :steelblue),
            CairoMakie.PolyElement(color = (:steelblue, 0.30)),
            CairoMakie.PolyElement(color = (:steelblue, 0.12)),
            CairoMakie.MarkerElement(color = :black, marker = :circle)],
        ["modelled median", "modelled count, 90%",
            "with measurement error, 90%", "digitised"];
        orientation = :horizontal, tellwidth = false)
    fig
end;

#md # ```@raw html
#md # </details>
#md # ```

onset_fit_fig #hide

# ### Symptom onsets by date of onset
#
# The panels above read one snapshot at a time.
# This figure reads the other way, along the onset date, and puts three
# quantities for the same day side by side.
#
# The open circles are the bar as the first figure to cover that onset date
# printed it, which is what a reader had at the time.
# The filled circles are the same bar as the most recent figure prints it,
# after every later revision.
# The gap between them is what late reporting and re-dating have added
# since, and it widens towards the right, where the first figure to print a
# date caught only a few days of its reporting.
#
# The two bands are the model's current view of the same days.
# The lower one is the count the triangle should eventually print for that
# onset date, the modelled onsets multiplied by the fitted asymptote, and
# is the band the filled circles are converging towards.
# The upper one is the modelled symptom onsets themselves, which sits above
# it because the asymptote is the ascertainment and does not reach one.
# The distance between the two bands is the part of the epidemic this
# stream will never print, and the distance between the upper band and the
# filled circles at the right-hand end is the nowcast.

#md # ```@raw html
#md # <details><summary>Onsets by onset date: as first printed, as printed now, and modelled</summary>
#md # ```

## First and most recent printed value for each onset date the digitised
## figures cover. Ordered by report date already, so the first block
## carrying a date is the earliest to print it and the last is the current
## reading. A date inside a block's own printed extent but with no row is a
## zero-height bar and does count; a date outside that extent is not
## covered by that figure at all and is skipped (the same rule the loader
## applies, see the [Data](@ref methods-data) section).
_onset_first_printed = Dict{Int, Float64}()
_onset_last_printed = Dict{Int, Float64}()
for snap in _onset_snaps
    lo, hi = extrema(keys(snap.onsets))
    for d in lo:Day(1):hi
        u = obs.n - value(obs.cutoff - d)
        (1 <= u <= obs.n) || continue
        v = Float64(get(snap.onsets, d, 0))
        haskey(_onset_first_printed, u) || (_onset_first_printed[u] = v)
        _onset_last_printed[u] = v
    end
end
_onset_by_date_days = sort(collect(keys(_onset_last_printed)))

## Modelled onsets on each of those days, and the count the triangle should
## eventually print for them: the same onsets times the hazard's asymptote,
## which is this stream's ascertainment. The extrapolated form of the
## cumulative reported proportion is used because the most recent onset
## dates sit past the calendar walk's fitted support.
_onset_by_date_onsets = [[_onset_daily_draws[i][u]
                          for i in eachindex(_onset_daily_draws)]
                         for u in _onset_by_date_days]
_onset_by_date_eventual = [[_onset_daily_draws[i][u] *
                            onset_report_cdf_extrapolated(
                                length(_onset_hazard.logit_h0[i]) - 1,
                                _onset_hazard.logit_h0[i], _onset_hazard.γ[i],
                                u, _onset_grid_start)
                            for i in eachindex(_onset_daily_draws)]
                           for u in _onset_by_date_days]

onset_by_date_fig = let
    fig = CairoMakie.Figure(; size = (900, 430))
    ax = CairoMakie.Axis(fig[1, 1];
        title = "Symptom onsets by date of onset: printed then, printed " *
                "now, and modelled",
        xlabel = "onset date", ylabel = "cases")
    xs = Float64.(_onset_by_date_days)
    q(ds, p) = [quantile(d, p) for d in ds]
    for (ds, col) in ((_onset_by_date_onsets, :seagreen),
        (_onset_by_date_eventual, :mediumpurple))
        CairoMakie.band!(ax, xs, q(ds, 0.05), q(ds, 0.95);
            color = (col, 0.20))
        CairoMakie.band!(ax, xs, q(ds, 0.25), q(ds, 0.75);
            color = (col, 0.30))
        CairoMakie.lines!(ax, xs, q(ds, 0.5); color = col, linewidth = 2)
    end
    CairoMakie.scatter!(ax, xs,
        [_onset_first_printed[u] for u in _onset_by_date_days];
        color = :transparent, marker = :circle, markersize = 8,
        strokewidth = 1, strokecolor = :grey40)
    CairoMakie.scatter!(ax, xs,
        [_onset_last_printed[u] for u in _onset_by_date_days];
        color = :black, marker = :cross, markersize = 9)
    ## Calendar labels on a grid-day axis, at weekly ticks so they do not
    ## collide at this width.
    _ticks = _onset_by_date_days[1:7:end]
    ax.xticks = (Float64.(_ticks), string.(grid_date.(_ticks)))
    ax.xticklabelrotation = pi / 4
    CairoMakie.Legend(fig[2, 1],
        [
            CairoMakie.MarkerElement(color = :transparent, marker = :circle,
                strokewidth = 1, strokecolor = :grey40),
            CairoMakie.MarkerElement(color = :black, marker = :cross),
            CairoMakie.PolyElement(color = (:mediumpurple, 0.30)),
            CairoMakie.PolyElement(color = (:seagreen, 0.30))],
        ["as first printed", "as printed now",
            "modelled eventual reports", "modelled onsets"];
        orientation = :horizontal, tellwidth = false, tellheight = true)
    fig
end;

#md # ```@raw html
#md # </details>
#md # ```

onset_by_date_fig #hide
# ### Posterior predictive checks
#
# A posterior predictive check draws replicated observations from the
# fitted joint model and compares them to the observed counts.
# The checks cover two groups: the dated DRC surveillance streams and the
# Uganda exports.
# The latent infection process is not checked here, as it carries no direct
# observation, and is shown instead as the estimated cumulative trajectories
# in the [joint model estimates](@ref "Joint model estimates") figure.
#
# The surveillance group is checked first.
# Each panel is shown over its own reporting dates with the observed series
# overlaid: the cumulative streams as replicated cumulative trajectories, and
# the daily new-suspect inflow and the daily isolation-bed occupancy on a
# daily scale (each day's replicated count against the observed count). The
# cumulative suspected case and death streams stop at their last stable
# vintage on 26 May; the daily new-suspect inflow then runs 4-11 June, where
# the cumulative suspected series freezes, and the isolation occupancy runs
# 1-11 June; the laboratory-confirmed streams keep reporting to the cut-off.

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
            increments = missing),
        breakpoint = _BREAKPOINT,
        background_re = true,
        confirmed_positivity_link = :composition,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days),
    chn_joint);

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
    title = "Suspected cases",
    dates = _vintage_dates(obs.reported_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(reported_increments.increments)),
    observed = obs.reported_history.counts, colour = :steelblue);
## Daily new-suspect inflow: a per-day count (not cumulative), so the panel
## is drawn with `cumulative = false` — each replicate is its own daily
## count against the observed daily count rather than a running total. Its
## days (4-7 June) pick up where the cumulative suspected panel freezes on
## 26 May.
suspected_daily_panel = (;
    title = "New suspects/day",
    dates = _vintage_dates(obs.suspected_daily_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(suspected_daily.increments)),
    observed = obs.suspected_daily_history.counts,
    colour = :slateblue, cumulative = false);
## Isolation/treatment-bed occupancy: a daily count, so the panel is drawn
## with `cumulative = false` — each replicate is the modelled bed count on a
## report day against the observed "Patients en isolement" count. The count
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
    title = "Patients in isolation",
    dates = _vintage_dates(obs.isolation_history.days[_iso_keep]),
    replicates = _vintage_replicates(
        pp_joint, @varname(isolation.obs)),
    observed = obs.isolation_history.counts[_iso_keep],
    colour = :darkorange, cumulative = false);
deaths_panel = (;
    title = "Suspected deaths",
    dates = _vintage_dates(obs.deaths_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(death_increments.increments)),
    observed = obs.deaths_history.counts, colour = :firebrick);
## Daily new suspected deaths: a per-day count (not cumulative), so the panel
## is drawn with `cumulative = false` — each replicate is its own daily count
## against the observed daily count rather than a running total. Its days
## (7-14 June) pick up where the cumulative suspected-death panel freezes on
## 26 May, the deaths analogue of the new-suspects-per-day panel.
suspected_daily_deaths_panel = (;
    title = "New suspected deaths/day",
    dates = _vintage_dates(obs.suspected_daily_deaths_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(suspected_daily_deaths.increments)),
    observed = obs.suspected_daily_deaths_history.counts,
    colour = :indianred, cumulative = false);
## Specimens analysed is the single modelled laboratory volume (the
## report-to-analysed delay and tested-fraction throughput), fit to the
## cumulative analysed series, so it gets the same cumulative conditional
## check as the suspected streams. This is the testing volume the
## confirmed-positivity denominator is built from.
tests_analysed_panel = (;
    title = "Specimens analysed (cumulative)",
    dates = _vintage_dates(obs.lab_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(analysed_increments.increments)),
    observed = obs.lab_history.counts, colour = :seagreen);
## Post-cutoff 24h analysed volume: once the cumulative series stops, INSP
## reports a 24h analysed count on some days. These are fitted as per-day
## volumes (not cumulative), so the panel is a standalone daily check
## (`cumulative = false`): the modelled daily analysed volume against the
## observed 24h count on each reported day.
tests_analysed_daily_panel = (;
    title = "Specimens analysed (24h)",
    dates = _vintage_dates(obs.lab_daily_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(analysed_daily_increments.increments)),
    observed = obs.lab_daily_history.counts, colour = :teal,
    cumulative = false);

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
    obs.confirmed_history, obs.lab_history, obs.lab_daily_history);
## Oldest-first: early (no denominator) → observed (analysed Binomial) →
## late (post-28 May; trusted 24h-analysed days are Binomial windows, the
## rest unanchored windows scored against the modelled volume).
_conf_window_days = vcat(_conf_windows.early_days, _conf_windows.obs_days,
    _conf_windows.late_days);
function _confirmed_at(day)
    i = searchsortedlast(obs.confirmed_history.days, day)
    return i == 0 ? 0 : Int(obs.confirmed_history.counts[i])
end;
_conf_early = _vintage_replicates(
    pp_joint, @varname(early_increments.increments));
_conf_obs = collect(first(pp_joint[k]
for k in keys(pp_joint)
if occursin("confirmed_state.confirmed_positives.positives", string(k))));
_conf_late = _vintage_replicates(
    pp_joint, @varname(late_increments.increments));
confirmed_panel = (;
    title = "Confirmed cases",
    dates = _vintage_dates(_conf_window_days),
    replicates = [vcat(collect(e), collect(p), collect(l))
                  for (e, p, l) in zip(vec(_conf_early), vec(_conf_obs), vec(_conf_late))],
    observed = [_confirmed_at(d) for d in _conf_window_days],
    colour = :goldenrod);

## Confirmed deaths are a per-vintage stream, scored as increments of the
## modelled confirmed-death trajectory up to the cut-off, so they get the
## same cumulative conditional check.
confirmed_deaths_panel = (;
    title = "Confirmed deaths",
    dates = _vintage_dates(obs.confirmed_deaths_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(cdeath_increments.increments)),
    observed = obs.confirmed_deaths_history.counts, colour = :purple);

## Recovered among confirmed ("cumul guéris") is a cumulative per-vintage
## stream fitted through the increments of the modelled recovered trajectory
## (the confirmation-to-recovery convolution of the daily confirmed cases) up
## to the cut-off, so it gets the same cumulative conditional check.
recovered_panel = (;
    title = "Recovered (confirmed)",
    dates = _vintage_dates(obs.recovered_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(recovered_increments.increments)),
    observed = obs.recovered_history.counts, colour = :mediumseagreen);

## Tableau 6 treatment-centre daily flows (the new patient-movement data
## sources): admissions and the discharge reasons (in-care deaths, rule-outs,
## absconded). Per-day counts, so drawn with `cumulative = false` — each
## replicate is the modelled daily flow on a report day against the observed
## Tableau 6 count.
admissions_panel = (;
    title = "Admissions/day",
    dates = _vintage_dates(obs.treatment_admissions_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(admissions.obs)),
    observed = obs.treatment_admissions_history.counts,
    colour = :teal, cumulative = false);
incare_deaths_panel = (;
    title = "In-care deaths/day",
    dates = _vintage_dates(obs.treatment_deaths_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(incare_deaths.increments)),
    observed = obs.treatment_deaths_history.counts,
    colour = :darkred, cumulative = false);
ruleouts_panel = (;
    title = "Rule-outs/day",
    dates = _vintage_dates(obs.treatment_ruleout_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(ruleouts.increments)),
    observed = obs.treatment_ruleout_history.counts,
    colour = :goldenrod, cumulative = false);
absconded_panel = (;
    title = "Absconded/day",
    dates = _vintage_dates(obs.treatment_absconded_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(absconded.increments)),
    observed = obs.treatment_absconded_history.counts,
    colour = :slategray, cumulative = false);

## Tableau 6 occupancy split (`dont confirmes` / `dont suspects`): the two
## in-care prevalence sub-stocks. Per-day census counts, so drawn with
## `cumulative = false` — each replicate is the modelled confirmed-in-care or
## suspect-in-care bed count on a report day against the observed sub-stock.
## On these split days the total-occupancy panel is not scored, so the two
## sub-stock panels carry the 13-23 June window.
confirmed_incare_panel = (;
    title = "Confirmed in care",
    dates = _vintage_dates(obs.treatment_confirmed_incare_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(confirmed_incare_obs.increments)),
    observed = obs.treatment_confirmed_incare_history.counts,
    colour = :darkgoldenrod, cumulative = false);
suspect_incare_panel = (;
    title = "Suspects in care",
    dates = _vintage_dates(obs.treatment_suspect_incare_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(suspect_incare_obs.increments)),
    observed = obs.treatment_suspect_incare_history.counts,
    colour = :chocolate, cumulative = false);

## Symptom-onset reporting triangle: one cell per (onset day, report day)
## pair, several onset dates per snapshot, unlike every panel above (one
## value per report day). To fit the same per-vintage panel shape, sum
## the cells sharing a report day into one net correction per snapshot —
## "one panel showing each snapshot as of its own report date" — rather
## than a full onset-by-report grid (shown separately as the fitted-vs-
## digitised snapshot figure in the Results section above). Per-day net
## correction, not a running total, so `cumulative = false`.
_onset_ppc_report_days = sort(unique(obs.onset_curve_history.report_days))
_onset_ppc_groups = [findall(==(r), obs.onset_curve_history.report_days)
                     for r in _onset_ppc_report_days]
_onset_ppc_replicates_raw = _vintage_replicates(
    pp_joint, @varname(onset_report_state.increments))
onset_panel = (;
    title = "Onset reports (net correction/snapshot)",
    dates = _vintage_dates(_onset_ppc_report_days),
    replicates = [[sum(collect(rep)[g]) for g in _onset_ppc_groups]
                  for rep in vec(_onset_ppc_replicates_raw)],
    observed = [sum(obs.onset_curve_history.increments[g])
                for g in _onset_ppc_groups],
    colour = :mediumpurple, cumulative = false);

## Each panel runs to its own last vintage: the suspected case and death
## streams freeze at 26 May (their last stable vintage) while the
## laboratory-confirmed streams keep reporting to the cut-off, so the
## confirmed panels show the full series the model is fitting, not just the
## window the suspected streams cover.
vintage_panels = [
    reported_panel, suspected_daily_panel, isolation_panel, confirmed_panel,
    deaths_panel, suspected_daily_deaths_panel, confirmed_deaths_panel,
    recovered_panel, tests_analysed_panel, tests_analysed_daily_panel,
    admissions_panel, incare_deaths_panel, ruleouts_panel, absconded_panel,
    confirmed_incare_panel, suspect_incare_panel, onset_panel];
joint_vintage_ppc_fig = plot_vintage_conditional_ppc(vintage_panels);

#md # ```@raw html
#md # </details>
#md # ```

joint_vintage_ppc_fig #hide

# The same check as per-vintage incidence: the count reported between
# consecutive situation reports rather than the running cumulative. Plotting
# the increment makes the trend in each stream read directly off the height
# of each step, so a rise or a slowdown is visible where the near-straight
# cumulative line hides it. The replicates are the modelled per-vintage
# increments, shown as 30/60/90% credible ribbons with the observed
# increment overlaid.

#md # ```@raw html
#md # <details><summary>Per-vintage incidence posterior predictive plot</summary>
#md # ```

joint_vintage_incidence_fig = plot_vintage_incidence_ppc(
    [reported_panel, suspected_daily_panel, isolation_panel, confirmed_panel,
    deaths_panel, suspected_daily_deaths_panel, confirmed_deaths_panel,
    recovered_panel, tests_analysed_panel, tests_analysed_daily_panel,
    onset_panel]);

#md # ```@raw html
#md # </details>
#md # ```

joint_vintage_incidence_fig #hide

# We score each stream's per-vintage conditional predictions against the
# observed counts. `bias` is the mean forecast bias over the vintages
# (negative = under-predicted, positive = over-predicted, zero = the observed
# counts sit at the predictive median); `50%/90% coverage` are the fractions
# of vintages whose observed count falls inside the central 50% and 90%
# predictive intervals, which a well-calibrated stream keeps near those nominal
# levels. Streams with a large bias or coverage far from nominal are the ones
# the joint fit reproduces less well.

stream_calibration_table = stream_calibration(vintage_panels);

# The calibration plot reads the table at a glance: the left panel marks each
# stream's empirical 50% and 90% coverage against dashed reference lines at the
# nominal levels, and the right panel marks the mean forecast bias against a
# dashed line at zero.

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

# The exports group is checked next.
# The Uganda export and export-death streams are dated per-day series, each
# import or death scored as a Poisson at its detection day. The scalar
# posterior predictive sums each replicate's per-day count vector across the
# dated days, giving the cumulative export and death total to compare with
# the observed count.

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
    pp_joint, @varname(death_obs.counts));

joint_ppc_fig = plot_posterior_predictive(
    pp_exports, nothing,
    obs.exported_cases, nothing;
    pp_exports_deaths = pp_exports_deaths,
    obs_exports_deaths = obs.exports_deaths);

#md # ```@raw html
#md # </details>
#md # ```

joint_ppc_fig #hide

# ### Posterior correlations and stream totals
#
# The heatmap is the posterior correlation between each pair of headline
# quantities: the outbreak size ($C_T$), the reproduction number ($R_T$), the
# outbreak age ($T$), the case-fatality ratio (CFR), the DRC and Uganda
# ascertainment fractions ($p_\text{drc}$, $p_\text{ug}$), the non-BVD
# background rate ($\lambda_\text{bg}$), the fraction tested ($\tau_\text{test}$),
# and the cut-off total expected for each stream. Blue is positive, red negative.

#md # ```@raw html
#md # <details><summary>Posterior correlation heatmap</summary>
#md # ```

correlation_fig = plot_correlation_heatmap(chn_joint,
    [:C_T, :R_T, :T, :CFR, :p_drc, :p_uganda, :lambda_bg, :tau_test,
        :expected_reports_T, :expected_deaths_T, :expected_confirmed_T];
    labels = Dict(:C_T => raw"C_T", :R_T => raw"R_T", :T => raw"T",
        :CFR => raw"\mathrm{CFR}", :p_drc => raw"p_\mathrm{drc}",
        :p_uganda => raw"p_\mathrm{ug}", :lambda_bg => raw"\lambda_\mathrm{bg}",
        :tau_test => raw"\tau_\mathrm{test}",
        :expected_reports_T => raw"\mathrm{susp.\ cases}",
        :expected_deaths_T => raw"\mathrm{susp.\ deaths}",
        :expected_confirmed_T => raw"\mathrm{conf.\ cases}"));

#md # ```@raw html
#md # </details>
#md # ```

correlation_fig #hide

# The stream-total plot takes each posterior draw, sums every stream over its
# own reporting dates, and marks the observed total with a crosshair. The
# diagonal panels are the predictive spread of each total against the observed
# value; the off-diagonal panels show whether the totals move together from
# draw to draw.

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
    analysed = _stream_total(tests_analysed_panel.replicates));
stream_observed = (;
    suspected_cases = Float64(obs.reported_history.counts[end]),
    suspected_deaths = Float64(obs.deaths_history.counts[end]),
    confirmed_cases = Float64(obs.confirmed_cases),
    confirmed_deaths = Float64(obs.confirmed_deaths_history.counts[end]),
    analysed = Float64(obs.lab_history.counts[end]));
stream_pairs_fig = plot_stream_pairs(stream_totals, stream_observed);

#md # ```@raw html
#md # </details>
#md # ```

stream_pairs_fig #hide

# ### Counterfactual: lower bound under no further transmission
#
# The committed future deaths $\Delta D$ if transmission stopped at the
# report date, defined in the methods
# [counterfactual](@ref "No-onward-transmission counterfactual").

#md # ```@raw html
#md # <details><summary>Project no-onward deaths and summarise</summary>
#md # ```

no_onward = predict_no_onward_deaths(
    chn_joint; obs_deaths = obs.total_deaths);

no_onward_table = streams_table(
    "no-onward total" => no_onward.total_projected;
    digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

no_onward_table #hide

# Two panels: the *still expected* deaths $\Delta D$ (future deaths in
# cases already infected by $T$, net of those already observed) on the
# left, and the *projected total* $D(T) + \Delta D$ on the right with
# the observed death count marked as a dashed black rule.

#md # ```@raw html
#md # <details><summary>No-onward projected-deaths plot</summary>
#md # ```

no_onward_fig = plot_no_onward_deaths(
    no_onward; obs_deaths = obs.total_deaths);

#md # ```@raw html
#md # </details>
#md # ```

no_onward_fig #hide

# ### Confirmed case-fatality ratio
#
# The delay-corrected confirmed CFR defined in the methods
# [delay-corrected confirmed CFR](@ref "Delay-corrected confirmed case-fatality ratio"),
# set against the structural (infection-based) CFR and the naive confirmed
# ratio.
# The corrected ratio debiases the naive confirmed ratio for the real-time
# delay between a case being confirmed and a death being confirmed; the
# structural CFR is the onset-level estimate the joint model fits.
# Reading the three together separates the real-time delay bias (naive versus
# corrected) from the case/death ascertainment difference (corrected versus
# structural).

#md # ```@raw html
#md # <details><summary>Compute the confirmed-CFR comparison</summary>
#md # ```

confirmed_cfr = delay_corrected_confirmed_cfr(chn_joint;
    obs_confirmed = obs.confirmed_cases,
    obs_confirmed_deaths = obs.confirmed_deaths);

confirmed_cfr_summary = confirmed_cfr_table(confirmed_cfr);

## Summary line carrying the data-anchored corrected estimate alongside the
## infection-based structural CFR, so both can be quoted together.
confirmed_cfr_line = let r = confirmed_cfr
    pct(x) = round(100 * x; digits = 1)
    corr = filter(isfinite, r.corrected)
    struc = filter(isfinite, r.structural)
    cs = posterior_summary(corr)
    ss = posterior_summary(struc)
    Markdown.parse(string(
        "**Delay-corrected confirmed CFR:** ",
        pct(quantile(corr, 0.5)), "% (90% CrI ",
        pct(cs.lo90), "–", pct(cs.hi90), "%), versus a naive confirmed ratio ",
        "of ", pct(r.naive_observed), "% and a structural (infection-based) ",
        "CFR of ", pct(quantile(struc, 0.5)), "% (90% CrI ",
        pct(ss.lo90), "–", pct(ss.hi90), "%)."))
end;

#md # ```@raw html
#md # </details>
#md # ```

confirmed_cfr_line #hide

confirmed_cfr_summary #hide

# The posterior densities of the delay-corrected confirmed CFR and the
# structural CFR, with the naive observed confirmed ratio drawn as a solid
# vertical rule and the median uncorrected modelled confirmed ratio as a
# dashed rule. The gap from the naive rule to the corrected density is the
# real-time delay debiasing; the gap to the structural density is the residual
# case/death ascertainment difference.

#md # ```@raw html
#md # <details><summary>Confirmed-CFR density plot</summary>
#md # ```

confirmed_cfr_fig = plot_confirmed_cfr(confirmed_cfr);

#md # ```@raw html
#md # </details>
#md # ```

confirmed_cfr_fig #hide

# ### One-week-ahead forecast results
#
# The cumulative and new expected counts by $T + 7$ from the no-change
# projection defined in the methods
# [one-week-ahead forecast](@ref "One-week-ahead forecast"). The summary
# table reports the confirmed case and death streams, the recovered total and
# the isolation-bed levels and daily flows; the observed-forecast plot below
# additionally shows the reported cases and suspected deaths, so every
# projected stream appears.

#md # ```@raw html
#md # <details><summary>Generate the one-week-ahead forecast</summary>
#md # ```

forecast = forecast_reported(chn_joint;
    horizon = 7,
    obs_cases = obs.reported_cases,
    obs_deaths = obs.total_deaths,
    obs_confirmed = obs.confirmed_cases,
    obs_confirmed_deaths = obs.confirmed_deaths,
    obs_recovered = obs.recovered_cases);
forecast_summary = forecast_table(forecast);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>One-week-ahead forecast summary table</summary>
#md # ```

forecast_summary #hide

#md # ```@raw html
#md # </details>
#md # ```

# The coming week at a glance, split into the latent quantities and the
# observations.
# The latent figure shows the new infections, symptom onsets and deaths over
# the horizon, with the reproduction number left to keep evolving across it.

#md # ```@raw html
#md # <details><summary>One-week-ahead latent forecast plot</summary>
#md # ```

forecast_latent_fig = plot_forecast_latent(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_latent_fig #hide

# The observed figure shows the new count each reported stream adds over the
# horizon: reported cases, suspected deaths, laboratory-confirmed cases,
# confirmed deaths and recovered, one panel per stream the forecast carries.

#md # ```@raw html
#md # <details><summary>One-week-ahead observed forecast plot</summary>
#md # ```

forecast_fig = plot_forecast(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_fig #hide

# The bed figure shows the projected isolation/treatment-bed demand (the need
# a week ahead, under unconstrained supply) against the supply-limited
# occupancy the beds can actually meet; the gap between the two is the
# projected bed shortfall, shown in the right panel. The reported "Patients en
# isolement" count is the occupied-bed count (the report computes the "Taux
# d'occupation" as that count over the bed capacity), so isolation is bed
# usage, gated by supply; the demand is its unobserved counterpart, the number
# who need a bed. Because the model carries a single national bed capacity it
# cannot represent local saturation, so the national shortfall understates
# local unmet need. On 13 June Ituri was at 93.9% occupancy while Sud-Kivu was
# at 21.9%, and beds free in one province cannot serve patients in another.

#md # ```@raw html
#md # <details><summary>One-week-ahead isolation-bed forecast plot</summary>
#md # ```

forecast_beds_fig = plot_forecast_beds(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_beds_fig #hide

# The flow figure projects the daily isolation/treatment flows a week ahead:
# new admissions, in-care deaths and rule-outs, each grown from its cut-off
# daily rate and replicated through the isolation dispersion. These are the
# daily-flow counterparts of the bed-stock forecast above.

#md # ```@raw html
#md # <details><summary>One-week-ahead treatment-flow forecast plot</summary>
#md # ```

forecast_flows_fig = plot_forecast_flows(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_flows_fig #hide

# ### Symptom-onset nowcast and forecast results
#
# The onset stream's projection, built as described in the methods
# [symptom-onset nowcast and forecast](@ref "Symptom-onset nowcast and
# forecast").
# The table reads top to bottom as the nowcast first and the forecast
# second, and the two halves must not be added together: the first three
# rows are the state of the outbreak at the cut-off, the next three the
# coming week.
#
# One row deserves care.
# "Onsets not yet reported at T" is not a backlog that will all arrive,
# because the fitted hazard's asymptote does not reach one, so the row
# holds the reporting backlog and the cases surveillance will never
# confirm, together.
# The "reports this week of onsets before T" row is the part of it the
# coming week should actually clear, and is the smaller number.

#md # ```@raw html
#md # <details><summary>Generate the symptom-onset nowcast and forecast</summary>
#md # ```

## The triangle's own grid, recomputed from the observations exactly as
## `onset_reporting_model` derives it, since it is data rather than chain
## contents. `_onset_grid_start`/`_onset_grid_end` are already built for the
## reporting-delay section above and reused here.
onset_forecast = forecast_onsets(chn_joint;
    grid_start = _onset_grid_start, grid_end = _onset_grid_end,
    n = obs.n, horizon = 7,
    obs_value = something(obs.onset_curve_history.last_total, 0));
onset_forecast_summary = onset_forecast_table(onset_forecast);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Symptom-onset nowcast and forecast summary table</summary>
#md # ```

onset_forecast_summary #hide

#md # ```@raw html
#md # </details>
#md # ```

# The left panel splits the coming week's new onset reports into reports of
# onsets that had already happened by the cut-off and reports of onsets
# still to come, and shows their sum.
# The fourth bar is the same sum after it has been through the observation
# model, which is what the next vintage will actually print.
# It is much wider than the sum it replicates, and the gap is the rescan:
# at a printed total of a couple of thousand cases the per-scan level error
# alone is worth tens of cases either way, well above the epidemic
# uncertainty on a week of new reports.
# Only the fourth bar is comparable to a digitised figure, and only it is
# scored.
#
# The right panel puts the nowcast itself on the same axes, the onsets that
# have happened against the share of them the triangle has printed.

#md # ```@raw html
#md # <details><summary>Symptom-onset nowcast and forecast plot</summary>
#md # ```

onset_forecast_fig = let
    fig = CairoMakie.Figure(; size = (960, 420))
    ## Two-line tick labels rather than rotated ones: the leftmost rotated
    ## label overhangs the axis and is clipped at the figure edge.
    ax1 = CairoMakie.Axis(fig[1, 1];
        title = "New onset reports over the coming week",
        ylabel = "cases", xticks = (1:4,
            ["already\nhappened", "not yet\nhappened", "sum of\nthe two",
                "as the next\nfigure reads it"]))
    ## The first three bars are latent, so the third is exactly the first
    ## two added. The fourth is that same sum replicated through the
    ## observation model, which is the scored quantity and the only one
    ## comparable to a digitised figure; it is wider by the per-scan level
    ## error, which is why the three latent bars are shown as well rather
    ## than a decomposition that appears not to add up.
    _latent_total = onset_forecast.onset_reports_backfill .+
                    onset_forecast.onset_reports_future
    for (i, d, col) in ((1, onset_forecast.onset_reports_backfill,
        :mediumpurple),
        (2, onset_forecast.onset_reports_future, :mediumpurple),
        (3, _latent_total, :mediumpurple),
        (4, Float64.(onset_forecast.onset_reports_new), :slategray))
        s = posterior_summary(d)
        CairoMakie.rangebars!(ax1, [Float64(i)], [s.lo90], [s.hi90];
            color = col, linewidth = 3)
        CairoMakie.rangebars!(ax1, [Float64(i)], [s.lo60], [s.hi60];
            color = col, linewidth = 8)
        CairoMakie.scatter!(ax1, [Float64(i)], [quantile(d, 0.5)];
            color = :black, markersize = 9)
    end
    ax2 = CairoMakie.Axis(fig[1, 2];
        title = "Symptom onsets by the cut-off",
        ylabel = "cases", xticks = (1:3,
            ["onsets\nto date", "reported\nby T", "not yet\nreported"]))
    for (i, d) in enumerate((onset_forecast.onsets_to_date,
        onset_forecast.onset_reports_to_date,
        onset_forecast.onsets_unreported))
        s = posterior_summary(d)
        CairoMakie.rangebars!(ax2, [Float64(i)], [s.lo90], [s.hi90];
            color = :seagreen, linewidth = 3)
        CairoMakie.rangebars!(ax2, [Float64(i)], [s.lo60], [s.hi60];
            color = :seagreen, linewidth = 8)
        CairoMakie.scatter!(ax2, [Float64(i)], [quantile(d, 0.5)];
            color = :black, markersize = 9)
    end
    ## The digitised total the "reported by T" bar is a model of, so the
    ## reader can see the fitted reported level against the figure itself.
    ismissing(obs.onset_curve_history.last_total) ||
        CairoMakie.hlines!(ax2,
            [Float64(obs.onset_curve_history.last_total)];
            color = :black, linestyle = :dash)
    fig
end;

#md # ```@raw html
#md # </details>
#md # ```

onset_forecast_fig #hide

# ## Saving results
#
# The tables above are written to an output directory at the repo
# root so they can be archived and shared. On every push to the main branch a
# GitHub Actions workflow regenerates these files and publishes them
# as a GitHub Release, downloadable from the repository's releases
# page (<https://github.com/epiforecasts/BVDOutbreakSize/releases>).
# The release bundles the summary tables, a thinned set of
# posterior draws, the latent symptom-onset ("symptomatic cases")
# trajectory over time, the one- to four-week-ahead forecasts of the
# observed streams (so each release records the forecast it made, for later
# scoring), and a copy of the input data manifest so
# the exact data that produced each result is recorded alongside it.

#md # ```@raw html
#md # <details><summary>Write outputs to output/</summary>
#md # ```

## Outputs default to `output/` in the package directory (where the
## docs build and Release workflow expect them). Set `BVD_OUTPUT_DIR`
## to redirect them, e.g. when running from a read-only package
## install.
output_dir = get(ENV, "BVD_OUTPUT_DIR",
    joinpath(pkgdir(BVDOutbreakSize), "output"))
mkpath(output_dir)

## Full parameter summary for the published CSV (infection, surveillance and
## export parameters together).
joint_summary = summary_table(chn_joint,
    [:r, :r0, :doubling_time, :T, :R_T, :CFR, :C_T,
        :p_drc, :p_uganda, :k, :tau_test, :lambda_bg,
        Symbol("exports_state.travel_state.daily_travellers")]; digits = 2)
CSV.write(joinpath(output_dir, "posterior_summary.csv"), joint_summary)
CSV.write(joinpath(output_dir, "confirmed_cfr_summary.csv"),
    confirmed_cfr_summary)

## Copy the input data so the release records what produced these
## results.
cp(joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml"),
    joinpath(output_dir, "observations.toml"); force = true)

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
forecast_runs = [(h,
                     forecast_reported(chn_joint; horizon = h,
                         obs_cases = obs.reported_cases,
                         obs_deaths = obs.total_deaths,
                         obs_confirmed = obs.confirmed_cases,
                         obs_confirmed_deaths = obs.confirmed_deaths,
                         obs_recovered = obs.recovered_cases,
                         grid_n = obs.n,
                         onset_grid_start = _onset_grid_start,
                         onset_grid_end = _onset_grid_end))
                 for h in forecast_horizons]
CSV.write(joinpath(output_dir, "forecast.csv"),
    forecast_archive(forecast_runs; made_date = obs.cutoff, thin = 5));

## The same one- to four-week-ahead forecast made from each FROZEN joint
## re-fit (the McCabe-matched cut-offs, the Chamla anchor and the one-week-back
## validation fit), each stamped with its OWN cut-off as the made date. This
## gives historical forecast evaluation using the current model at past data
## cut-offs, scored against what has since been observed, without
## reconstructing old release tags. Each frozen fit uses its own frozen
## observations for the cut-off counts. The May cut-offs predate the isolation
## and recovered streams, so those are simply absent for them; the per-stream
## guard in `forecast_archive` skips a stream a fit does not carry.
frozen_forecast_fits = unique(f -> f.o.cutoff,
    [frozen_results; frozen_by_cutoff[chamla_cutoff]; frozen_lastweek])
frozen_forecast_archive = DataFrame(made_date = Date[], horizon = Int[],
    target_date = Date[], stream = String[], draw = Int[], value = Float64[])
## The onset grid belongs to the triangle each frozen fit actually saw, not
## to the live one: the May cut-offs predate the digitised figure entirely,
## so their grid is empty and the onset block is simply absent for them.
function _frozen_onset_grid(o)
    isempty(o.onset_curve_history.onset_days) && return (nothing, nothing)
    gs = minimum(o.onset_curve_history.onset_days)
    return (gs, max(maximum(o.onset_curve_history.report_days), gs))
end
for f in frozen_forecast_fits
    _fgs, _fge = _frozen_onset_grid(f.o)
    runs = [(h,
                forecast_reported(f.chn; horizon = h,
                    obs_cases = f.o.reported_cases,
                    obs_deaths = f.o.total_deaths,
                    obs_confirmed = f.o.confirmed_cases,
                    obs_confirmed_deaths = f.o.confirmed_deaths,
                    obs_recovered = f.o.recovered_cases,
                    grid_n = f.o.n,
                    onset_grid_start = _fgs, onset_grid_end = _fge))
            for h in forecast_horizons]
    append!(frozen_forecast_archive,
        forecast_archive(runs; made_date = f.o.cutoff, thin = 5))
end
CSV.write(joinpath(output_dir, "forecast_frozen.csv"),
    frozen_forecast_archive);

## Per-fit release assets: the reproduction number, outbreak size and forecasts
## for every fit rather than the joint alone, so a release records what each
## dataset implies on its own and can later be scored against the joint. The
## single-stream fits walk Rt from day 1 while the joint walks from
## `RT_WALK_LEAD` days before the first situation report, so each fit carries
## the starts its own fit used. Each single-stream fit forecasts only the
## dataset it observes; the joint forecasts every shared stream. Recovered has
## no single-stream fit, so it stays a joint-only stream in `forecast.csv`.
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
stream_fits = [
    (; fit = "joint", chn = chn_joint, rt_start = _rt_start_plot,
        rt_walk_start = _rt_walk_start_joint,
        streams = [(:reported_cases, "reported cases", obs.reported_cases),
            (:suspected_deaths, "suspected deaths", obs.total_deaths),
            (:confirmed_cases, "confirmed cases", obs.confirmed_cases),
            (:confirmed_deaths, "confirmed deaths", obs.confirmed_deaths),
            (:isolation_beds, "isolation beds", _iso_at_cutoff),
            (:exports, "exports", obs.exported_cases),
            (:onset_reports, "onset reports", _onset_at_cutoff)]),
    (; fit = "cases", chn = chn_cases, rt_start = 1, rt_walk_start = 1,
        streams = [(:reported_cases, "reported cases", obs.reported_cases)]),
    (; fit = "deaths", chn = chn_deaths, rt_start = 1, rt_walk_start = 1,
        streams = [(:suspected_deaths, "suspected deaths", obs.total_deaths)]),
    (; fit = "confirmed", chn = chn_confirmed, rt_start = 1, rt_walk_start = 1,
        streams = [(:confirmed_cases, "confirmed cases", obs.confirmed_cases)]),
    (; fit = "confirmed_deaths", chn = chn_confirmed_deaths, rt_start = 1,
        rt_walk_start = 1,
        streams = [(:confirmed_deaths, "confirmed deaths",
            obs.confirmed_deaths)]),
    (; fit = "treatment", chn = chn_treatment, rt_start = 1,
        rt_walk_start = 1,
        streams = [(:isolation_beds, "isolation beds", _iso_at_cutoff)]),
    (; fit = "exports", chn = chn_exports, rt_start = 1, rt_walk_start = 1,
        streams = [(:exports, "exports", obs.exported_cases)]),
    (; fit = "onsets", chn = chn_onsets, rt_start = 1, rt_walk_start = 1,
        streams = [(:onset_reports, "onset reports", _onset_at_cutoff)])
]

## Cut-off reproduction number per fit. The joint exposes it as `R_T`; the
## single-stream composers do not (the alias lives in `bvd_joint`, not the
## shared latent submodel), so theirs is rebuilt from the walk parameters every
## chain carries and read at the cut-off, the last day of the reconstructed
## path. `ramp` matches the model's 21-day intervention scale-up.
function _fit_rt_draws(f)
    f.fit == "joint" && return vec(Array(f.chn[:R_T]))
    rt = reconstruct_rt(f.chn; n = obs.n, breakpoint = _BREAKPOINT,
        rt_start = f.rt_start, rt_walk_start = f.rt_walk_start, ramp = RT_INTERVENTION_RAMP)
    return Float64[rt[i, obs.n] for i in axes(rt, 1)]
end

_stream_quantities = [(f.fit, _fit_rt_draws(f), vec(Array(f.chn[:C_T])))
                      for f in stream_fits]

## One row per fit and quantity, with the median and the 30/60/90% credible
## bounds the report's tables use.
function _stream_estimate_row(fit, quantity, draws)
    s = posterior_summary(draws)
    return (fit = fit, quantity = quantity, median = quantile(draws, 0.5),
        lo30 = s.lo30, hi30 = s.hi30, lo60 = s.lo60, hi60 = s.hi60,
        lo90 = s.lo90, hi90 = s.hi90)
end
stream_estimates = DataFrame([_stream_estimate_row(fit, q, d)
                              for (fit, rt, ct) in _stream_quantities
                              for (q, d) in (("R_T", rt), ("C_T", ct))])
CSV.write(joinpath(output_dir, "stream_estimates.csv"), stream_estimates);

## Thinned reproduction-number and outbreak-size draws per fit, so downstream
## scoring can recompute its own summaries rather than reuse the intervals.
stream_draws = DataFrame([(fit = fit, quantity = q, draw = d, value = v)
                          for (fit, rt, ct) in _stream_quantities
                          for (q, vals) in (("R_T", rt), ("C_T", ct))
                          for (d, v) in enumerate(vals[1:stream_thin:end])])
CSV.write(joinpath(output_dir, "stream_draws.csv"), stream_draws);

## Per-fit forecasts of each fit's own observed stream, in the `forecast.csv`
## long schema plus the fit that made them. Rebuilding a single-stream fit's
## cut-off growth rate needs the grid length and the breakpoint, which are data
## rather than chain contents, so both are passed.
stream_forecasts = DataFrame(made_date = Date[], horizon = Int[],
    target_date = Date[], stream = String[], draw = Int[], value = Float64[],
    fit = String[])
for f in stream_fits, (stream, label, obs_value) in f.streams,
    h in forecast_horizons
    ## The onset grid is ignored by every other stream, so it is passed
    ## unconditionally rather than branching the loop on the stream name.
    _vals = forecast_stream(f.chn, stream; horizon = h,
        obs_value = obs_value, n = obs.n, breakpoint = _BREAKPOINT,
        rt_start = f.rt_start, rt_walk_start = f.rt_walk_start,
        onset_grid_start = _onset_grid_start,
        onset_grid_end = _onset_grid_end)
    for (d, i) in enumerate(1:stream_thin:length(_vals))
        push!(stream_forecasts, (obs.cutoff, h, obs.cutoff + Day(h), label,
            d, Float64(_vals[i]), f.fit))
    end
end

## Confirmed/suspect ward-bed occupancy for the joint, taken from the joint's
## own `forecast_reported` runs (partitioned by the cut-off confirmed share in
## `forecast_reported`, not re-derived here) so the ward beds are scored on the
## same footing as the total occupancy in the preferred `stream_forecasts.csv`
## asset. `forecast_stream` cannot project these — the split is not a growable
## stream but a partition of the total — so they are read from the archive
## runs. Dormant until the chain carries the confirmed in-care split:
## `forecast_reported` emits these columns only when the confirmed in-care
## prevalence `expected_confirmed_incare_T` is present, so the guard skips
## them otherwise.
for (h, fc) in forecast_runs,
    (col, label) in ((:suspect_occupancy, "isolation beds (suspected)"),
        (:confirmed_occupancy, "treatment beds"))

    col in propertynames(fc) || continue
    _wvals = fc[!, col]
    for (d, i) in enumerate(1:stream_thin:length(_wvals))
        push!(stream_forecasts, (obs.cutoff, h, obs.cutoff + Day(h), label,
            d, Float64(_wvals[i]), "joint"))
    end
end
CSV.write(joinpath(output_dir, "stream_forecasts.csv"), stream_forecasts);

## Latent symptom-onset trajectory over time, the "symptomatic cases" curve,
## showing outbreak growth: one row per grid day with the 30/60/90%
## credible intervals of both the daily new and cumulative onsets.
onsets_over_time_table = onsets_over_time(chn_joint;
    n = obs.n, seeding = obs.seeding)
CSV.write(joinpath(output_dir, "onsets_over_time.csv"),
    onsets_over_time_table);

#md # ```@raw html
#md # </details>
#md # ```

# ### Summary-page assets
#
# The one-page [Summary dashboard](@ref) reuses the results computed above
# rather than re-fitting. Here we save its headline text, headline tables and
# the four figures it shows (reproduction number, the reproduction number each
# data stream implies on its own, infections over time, and modelled versus
# observed reported cases) into `docs/src/summary_assets/`,
# so the static dashboard page can embed them after this build step has run.

#md # ```@raw html
#md # <details><summary>Write the dashboard assets</summary>
#md # ```

dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets")
mkpath(dashboard_dir)

## Figures: estimated R(t), the R(t) each data stream implies on its own,
## latent infections over time, and the modelled versus observed reported
## cases. All are produced in the Results sections above; here we just write
## them out at the dashboard size.
CairoMakie.save(joinpath(dashboard_dir, "rt.png"), rt_fig)
CairoMakie.save(joinpath(dashboard_dir, "infections.png"),
    cumulative_traj_fig)
CairoMakie.save(joinpath(dashboard_dir, "reported_cases.png"),
    joint_vintage_ppc_fig)

## Headline prose: the same bullet summary shown at the top of the Results
## section, serialised to markdown so the dashboard renders it verbatim.
open(joinpath(dashboard_dir, "headline.md"), "w") do io
    print(io, sprint(Markdown.plain, summary_ranges))
end

## Headline tables: outbreak size and timing as whole numbers, and the
## growth and severity parameters to two decimals, each with reader-friendly
## quantity names.
dashboard_counts = summary_table(chn_joint, [:C_T, :T]; digits = 0,
    labels = Dict(:C_T => "Cumulative infections",
        :T => "Outbreak age (days)"))
dashboard_rates = summary_table(chn_joint,
    [:R0, :R_T, :r, :doubling_time, :CFR]; digits = 2,
    labels = Dict(:R0 => "Initial reproduction number",
        :R_T => "Latest reproduction number",
        :r => "Latest growth rate (per day)",
        :doubling_time => "Latest doubling time (days)",
        :CFR => "Case-fatality ratio"))
open(joinpath(dashboard_dir, "headline_counts.md"), "w") do io
    print(io, markdown_table(dashboard_counts))
end
open(joinpath(dashboard_dir, "headline_rates.md"), "w") do io
    print(io, markdown_table(dashboard_rates))
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
# The full analysis code, data and model definitions are in the
# [epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
# repository. Issues, corrections and suggestions are welcome there.
# Maintained by Sam Abbott, Kath Sherratt, Samuel Brand and Sebastian
# Funk.
