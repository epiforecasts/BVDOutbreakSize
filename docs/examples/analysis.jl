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
#   infections, is the headline quantity reported separately.
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
#md # <details><summary>Load packages and seed the RNG</summary>
#md # ```

using Turing
using Distributions
using StatsFuns: logistic
using DataFrames: DataFrame, eachrow
import CSV
using Random
using Markdown
using Dates: Date, Day, value
using BVDOutbreakSize
import CairoMakie
## Loading TensorBoardLogger activates the `tensorboard_callback` extension
## so `fit_callback` can stream each fit to TensorBoard as well as a progress
## log. Set `BVD_FIT_LOG=none` to disable all fit logging (CI release builds).
using TensorBoardLogger

## Render figures at higher resolution so they stay crisp in the docs.
CairoMakie.activate!(type = "png", px_per_unit = 3)

Random.seed!(20260518)

#md # ```@raw html
#md # </details>
#md # ```

# ## Methods
#
# ### Data
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
# the incidence analogue of the isolation prevalence stream). We extracted
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
# The first table lists each figure at the cut-off, or at the date
# reporting stopped for that stream. The second table gives the per-date
# history of each situation-report stream; the model fits the
# between-report increments of these series, so a single date reduces to
# the cut-off total.

#md # ```@raw html
#md # <details><summary>Loading observations and building the data table</summary>
#md # ```

obs = load_observations()
## Grid day-index (day n is the cut-off) back to a calendar date.
grid_date(day) = obs.cutoff - Day(obs.n - day)
## Date a cumulative history last reports (the cut-off for streams that
## run to it, or the freeze date for streams that stop earlier).
hist_last_date(h) = isempty(h.days) ? missing : grid_date(maximum(h.days))
observations_table = DataFrame(
    field = [
        "exported_cases",
        "exports_deaths",
        "suspected_deaths",
        "suspected_cases",
        "confirmed_cases",
        "confirmed_deaths",
        "specimens_analysed",
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
        hist_last_date(obs.lab_history),
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
        obs.tests_analysed,
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
# \tau \sim \mathrm{Normal}^{+}(0,\ 0.3), \tag{20}
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

# ##### Isolation occupancy
#
# The "Patients en isolement" figure is the daily count of occupied
# isolation/treatment beds. Bed occupancy may be supply-driven, with demand
# for beds able to outstrip supply and occupancy catching up as capacity
# expands. To allow for that we model a latent bed demand and the
# supply-limited occupancy it produces rather than occupancy directly.
#
# The latent demand is the suspect inflow carried through a length-of-stay
# survival $S(\tau) = P(\text{LOS} \ge \tau)$ (the renewal analogue of the
# convolution secondary-observation model of EpiNow2 [epinow2](@cite)). A
# proportion $p_{\text{iso}}$ of the reported suspects need a bed. The
# suspects are a BVD/background mixture leaving on different clocks, so the
# demand is the sum of two survival convolutions. The BVD demand uses the
# treatment length-of-stay $S_{\text{BVD}}$, the time an admitted BVD case
# occupies a bed, with the admission-to-death delay from the line-list
# reanalysis [bdbv_linelist_analysis_2026](@cite) as its prior. The
# non-BVD demand uses the rule-out stay $S_{\text{ruleout}}$, how long a
# ruled-out suspect occupies a bed before discharge, with the report-to-receipt
# laboratory turnaround as its prior,
#
# ```math
# D_t = p_{\text{iso}}\left[ \sum_{s \ge 0} p_{\text{DRC}}\,
#       \text{bvd}_{t-s}\, S_{\text{BVD}}(s) + \sum_{s \ge 0}
#       \lambda_{\text{bg},\,t-s}\, S_{\text{ruleout}}(s) \right].
# ```
#
# The bed capacity $C(t)$ is a non-decreasing random walk on weekly knots,
# since beds are added over the response and not taken away, pinned by the
# implied bed count, the reported occupancy (the "Patients en isolement" count)
# divided by the reported "Taux d'occupation" rate ($\approx 400 \to 452$ beds
# over 9–13 June).
#
# The occupied-bed count is the latent demand right-censored at the bed
# capacity: while demand is below capacity the count tracks it, and once demand
# reaches capacity the count is censored there. The censoring bound is fixed at
# the recorded implied capacity $C^{\text{cap}}_j$ so the latent demand is left
# uncensored,
#
# ```math
# O_j \sim \mathrm{censored}\bigl(\mathrm{NegBinomial}(D_{t_j},\ k_{\text{iso}});\
#     \text{upper} = C^{\text{cap}}_j\bigr),
# \qquad
# C^{\text{obs}}_j \sim \mathrm{NegBinomial}(C_{t_j},\ k_{\text{iso}}),
# ```
#
# with a dispersion $k_{\text{iso}}$ of its own (not shared with the other
# streams). Occupancy below capacity identifies the demand directly; the part
# of demand above a saturated capacity is only partially identified, since
# occupancy says demand was at least the beds filled and not how much more, so
# the bed shortfall above capacity is informed by the demand model and its
# priors rather than measured by the occupancy. The model exposes the cut-off
# occupancy, the cut-off bed demand (the need under unconstrained supply),
# their difference (the bed shortfall) and the utilisation $O_T / C$.
#
# One limitation is that this is a single national model, with one national
# bed capacity and one national demand, so it cannot represent local
# saturation. On 13 June Ituri was at 93.9% occupancy while Sud-Kivu was at
# 21.9%, so beds free in one province cannot serve patients who need them in
# another, and the national capacity averages over a saturated epicentre and
# slack elsewhere. The national shortfall therefore understates the local
# unmet need. The renewal model does not carry per-province inflow, so it
# cannot be split into the per-province bed model at which the supply
# constraint actually operates. A second
# limitation is that capacity is taken as a single (slowly varying) national
# quantity even though beds are being added.
#
# The exposed BVD share is the true-BVD fraction of demand (BVD-confirmed plus
# BVD-suspect), not the report's confirmed/suspect split. The fitted occupancy
# series is the all-patients column from 1 June (SitRep 018) onward.

#md # ```@raw html
#md # <details><summary>Submodel: treatment_admission_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.treatment_admission_model(
#md #         (; days = Int[], counts = Int[]),
#md #         Float64[], Float64[], 0.25)), "\n```"))
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
#     \sum_{t = d_{i-1}+1}^{d_i} m_t,\ k\Bigr). \tag{26}
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
#     \sum_{t = d_{i-1}+1}^{d_i} v_t,\ k\Bigr). \tag{27}
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
# C_v \sim \mathrm{Binomial}(A_v,\ p_{\text{pos},v}), \tag{28}
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
#     \mathrm{NegBinomial}(p_{\text{pos},v}\, V_v,\ k). \tag{29}
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
#     \sum_{t = d_{i-1}+1}^{d_i} \text{cd}_t,\ k\Bigr). \tag{30}
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
# a NegBinomial whose dispersion $k_{\text{rec}}$ is its own rather than
# shared with the other streams:
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
# \Lambda(t) = \sum_{u \le t} \lambda_u. \tag{31}
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
# 0 \sim \mathrm{Poisson}\!\bigl(\Lambda(d_1 - 1)\bigr). \tag{32}
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
# \Lambda_d(t) = \sum_{u \le t} \mu_u. \tag{33}
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
# 0 \sim \mathrm{Poisson}\!\bigl(\Lambda_d(\delta_1 - 1)\bigr). \tag{34}
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
# Alongside the joint model we write single-stream models for each
# count-based stream (exported cases, suspected deaths, suspected cases,
# laboratory-confirmed cases, confirmed deaths and deaths among exports), so
# each stream's posterior over the outbreak size can be compared with the
# joint. Other model variants reuse these models with different amounts of
# data, cutting the data to an earlier date or dropping the counts.

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

#md # ```@raw html
#md # <details><summary>Run the joint and per-stream NUTS fits</summary>
#md # ```

const _BREAKPOINT = obs.n - obs.who_first_sitrep_days

## The joint and six single-stream fits are independent, so they run through
## `fit_parallel`: model-level parallelism bounded by the available threads
## (sequential on CI's two threads, several at once on a many-core machine).
## The exports-deaths composer keeps the deaths and exports submodels only
## for their CFR, onset-to-death PMF and export onsets, leaving their own
## counts missing, which leaves two redundant sampled discrete draws, so its
## model check is disabled (see `nuts_sample`).
## Setup for the single master fit pool: relocated cut-offs, frozen-fit and
## sensitivity-variant helpers so every independent fit can run in one pool.
validation_cutoff = string(obs.cutoff - Day(7))

## Published per-release estimates, pulled from the tagged results
## releases by `scripts/refresh_releases.jl` into
## `data/released_estimates.csv`. Columns: tag, date, model (integral or
## renewal), median and the 30/60/90% bounds.
released_df = CSV.read(
    joinpath(pkgdir(BVDOutbreakSize), "data", "released_estimates.csv"),
    DataFrame)

## Integral-era release cut-offs for the frozen renewal overlay: a
## release whose data cut-off already has a renewal release needs no
## re-fit, since that release already is the renewal estimate. 20, 23 and
## 27 May are additionally fit for the matched-cutoff comparison further
## down (27 May matches the Lancet publication's cut-off).
renewal_release_dates = Set(string(r.date)
for r in eachrow(released_df) if r.model == "renewal")
frozen_evolution_cutoffs = sort(unique(string(r.date)
for r in eachrow(released_df)
if r.model == "integral" && string(r.date) ∉ renewal_release_dates))
frozen_cutoffs = sort(union(frozen_evolution_cutoffs,
    ["2026-05-20", "2026-05-23", "2026-05-27"]))

## A joint fit at the full headline settings (1000 draws × 2 chains) to the
## data frozen at `cutoff_date`. The frozen named tuple has the same shape as
## the full `obs`, so the model call mirrors the headline joint fit.
function fit_frozen_joint(cutoff_date; samples = 1000, chains = 2)
    o = freeze_observations(cutoff_date)
    bp = o.n - o.who_first_sitrep_days
    chn = nuts_sample(
        bvd_joint(
            o.n, o.exported_cases, o.total_deaths,
            o.reported_cases, o.exports_deaths, o.confirmed_cases,
            o.tests_analysed;
            confirmed_deaths = o.confirmed_deaths,
            deaths_history = o.deaths_history,
            reported_history = o.reported_history,
            confirmed_history = o.confirmed_history,
            confirmed_deaths_history = o.confirmed_deaths_history,
            lab_history = o.lab_history,
            lab_daily_history = o.lab_daily_history,
            isolation_history = o.isolation_history,
            bed_capacity_history = o.bed_capacity_history,
            export_case_days = o.export_case_days,
            export_death_days = o.export_death_days,
            breakpoint = bp,
            background_re = true,
            confirmed_positivity_link = :composition,
            genetic = genetic_seeding_model,
            tmrca_days = o.tmrca_days);
        samples = samples, chains = chains,
        callback = fit_callback("frozen_$(cutoff_date)"))
    return (; cutoff = o.cutoff, o, chn)
end

## One joint re-fit on the live data at the full headline settings, with hooks
## to override the genetic-seeding bound and the deaths submodel. The deaths
## submodel is passed the same way the genetic-seeding override is, as a
## closure matching the joint's `deaths(history, total, onsets, k;
## background_re)` call, so an alternative onset-to-death delay can be injected
## without touching the package.
function refit_joint_variant(;
        deaths = deaths_model,
        tmrca_days = obs.tmrca_days,
        tmrca_days_sd = 15.0,
        samples = 1000, chains = 2)
    chn = nuts_sample(
        bvd_joint(
            obs.n, obs.exported_cases, obs.total_deaths,
            obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
            obs.tests_analysed;
            confirmed_deaths = obs.confirmed_deaths,
            recovered_cases = obs.recovered_cases,
            deaths_history = obs.deaths_history,
            reported_history = obs.reported_history,
            confirmed_history = obs.confirmed_history,
            confirmed_deaths_history = obs.confirmed_deaths_history,
            lab_history = obs.lab_history,
            lab_daily_history = obs.lab_daily_history,
            suspected_daily_history = obs.suspected_daily_history,
            suspected_daily_deaths_history =
            obs.suspected_daily_deaths_history,
            isolation_history = obs.isolation_history,
            bed_capacity_history = obs.bed_capacity_history,
            recovered_history = obs.recovered_history,
            export_case_days = obs.export_case_days,
            export_death_days = obs.export_death_days,
            breakpoint = _BREAKPOINT,
            background_re = true,
            confirmed_positivity_link = :composition,
            deaths = deaths,
            genetic = genetic_seeding_model,
            tmrca_days = tmrca_days,
            tmrca_days_sd = tmrca_days_sd);
        samples = samples, chains = chains,
        callback = fit_callback("variant"))
    return chn
end

## Community-pathway onset-to-death delay from the Isiro 2012 line-list
## reanalysis community-death model (a single Gamma, implied mean ≈ 8 d,
## shape ≈ 5.5), a closure that re-injects this delay on its natural Gamma
## shape and scale into the deaths submodel while keeping its other defaults.
deaths_community_delay = (history,
    total,
    onsets,
    k;
    kwargs...) -> deaths_model(history, total, onsets, k;
    onset_to_death = gamma_delay_model(40;
        alpha_prior = truncated(Normal(5.48, 2.0); lower = 0.01),
        theta_prior = truncated(Normal(1.49, 0.5); lower = 0.1)),
    kwargs...)

## The faster early-epidemic clock dates the common ancestor about 17 days
## more recently, so the bound on the outbreak age sits that many days
## closer to the cut-off, with a tighter spread from its narrower interval.
clock_alt_offset = value(Date("2026-04-11") - Date("2026-03-25"))
tmrca_days_alt = obs.tmrca_days - clock_alt_offset

## Sensitivity refits (onset-to-death delay, molecular clock) are slow extra
## joint fits. They are PAUSED due to compute constraints: the serial
## sensitivity re-fits pushed the main docs build toward the 6h CI cap, so they
## are not run on any build for now. The code below is retained unchanged;
## re-enable by restoring the `BVD_RUN_SENSITIVITY` env gate:
##   RUN_SENSITIVITY = lowercase(strip(get(ENV, "BVD_RUN_SENSITIVITY",
##       "false"))) in ("true", "1", "yes", "on")
RUN_SENSITIVITY = false

## Every independent fit runs as one work-stealing pool so the long joint fit
## overlaps the per-stream, frozen and (gated) sensitivity re-fits and keeps
## all cores busy, rather than the short fits idling cores while the joint
## finishes. The fits are data-only independent, so the schedule is free.
_headline_thunks = [
    () -> nuts_sample(
        bvd_joint(
            obs.n, obs.exported_cases, obs.total_deaths,
            obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
            obs.tests_analysed;
            confirmed_deaths = obs.confirmed_deaths,
            recovered_cases = obs.recovered_cases,
            deaths_history = obs.deaths_history,
            reported_history = obs.reported_history,
            confirmed_history = obs.confirmed_history,
            confirmed_deaths_history = obs.confirmed_deaths_history,
            lab_history = obs.lab_history,
            lab_daily_history = obs.lab_daily_history,
            suspected_daily_history = obs.suspected_daily_history,
            suspected_daily_deaths_history = obs.suspected_daily_deaths_history,
            isolation_history = obs.isolation_history,
            bed_capacity_history = obs.bed_capacity_history,
            recovered_history = obs.recovered_history,
            export_case_days = obs.export_case_days,
            export_death_days = obs.export_death_days,
            breakpoint = _BREAKPOINT,
            background_re = true,
            confirmed_positivity_link = :composition,
            genetic = genetic_seeding_model,
            tmrca_days = obs.tmrca_days);
        callback = fit_callback("joint")),
    ## Exports fit cases and deaths jointly as one "exports" stream, sharing
    ## the travel-gated at-risk prevalence, rather than as two separate fits.
    () -> nuts_sample(
        exports_joint_only_model(obs.n, obs.exported_cases,
            obs.exports_deaths;
            export_case_days = obs.export_case_days,
            export_death_days = obs.export_death_days,
            breakpoint = _BREAKPOINT);
        check_model = false, callback = fit_callback("exports")),
    () -> nuts_sample(
        deaths_only_model(obs.n, obs.total_deaths;
            deaths_history = obs.deaths_history,
            suspected_daily_deaths_history = obs.suspected_daily_deaths_history,
            breakpoint = _BREAKPOINT);
        callback = fit_callback("deaths")),
    () -> nuts_sample(
        cases_only_model(obs.n, obs.reported_cases;
            reported_history = obs.reported_history,
            suspected_daily_history = obs.suspected_daily_history,
            breakpoint = _BREAKPOINT);
        callback = fit_callback("cases")),
    () -> nuts_sample(
        confirmed_only_model(obs.n, obs.confirmed_cases;
            confirmed_history = obs.confirmed_history,
            lab_history = obs.lab_history,
            lab_daily_history = obs.lab_daily_history,
            breakpoint = _BREAKPOINT);
        callback = fit_callback("confirmed")),
    () -> nuts_sample(
        confirmed_deaths_only_model(obs.n, obs.confirmed_deaths,
            obs.total_deaths;
            deaths_history = obs.deaths_history,
            confirmed_deaths_history = obs.confirmed_deaths_history,
            breakpoint = _BREAKPOINT);
        callback = fit_callback("confirmed_deaths")),
    () -> nuts_sample(
        treatment_only_model(obs.n;
            isolation_history = obs.isolation_history,
            bed_capacity_history = obs.bed_capacity_history,
            breakpoint = _BREAKPOINT);
        callback = fit_callback("treatment"))
]
_frozen_thunks = [() -> fit_frozen_joint(c) for c in frozen_cutoffs]
_sens_thunks = RUN_SENSITIVITY ?
               [() -> refit_joint_variant(deaths = deaths_community_delay),
    () -> refit_joint_variant(
        tmrca_days = tmrca_days_alt, tmrca_days_sd = 9.0)] : []
_all_fits = fit_parallel(vcat(_headline_thunks,
    [() -> fit_frozen_joint(validation_cutoff)], _frozen_thunks, _sens_thunks))

_n_head = length(_headline_thunks)
_n_frozen = length(frozen_cutoffs)
(chn_joint, chn_exports, chn_deaths, chn_cases, chn_confirmed,
    chn_confirmed_deaths, chn_treatment) = _all_fits[1:_n_head]
frozen_lastweek = _all_fits[_n_head + 1]
frozen_results = _all_fits[(_n_head + 2):(_n_head + 1 + _n_frozen)]
frozen_by_cutoff = Dict(zip(frozen_cutoffs, frozen_results))
frozen_C(c) = vec(Array(frozen_by_cutoff[c].chn[:C_T]))
if RUN_SENSITIVITY
    chn_joint_community_delay = _all_fits[_n_head + 2 + _n_frozen]
    chn_joint_fast_clock = _all_fits[_n_head + 3 + _n_frozen]
end

posterior_C_joint = vec(Array(chn_joint[:C_T]));
posterior_C_exports = vec(Array(chn_exports[:C_T]));
posterior_C_deaths = vec(Array(chn_deaths[:C_T]));
posterior_C_cases = vec(Array(chn_cases[:C_T]));
posterior_C_confirmed = vec(Array(chn_confirmed[:C_T]));
posterior_C_confirmed_deaths = vec(Array(chn_confirmed_deaths[:C_T]));
posterior_C_treatment = vec(Array(chn_treatment[:C_T]));

## Clean display names for the summary tables and pair plots. The submodel
## prefixes (`rt_state.`, `gi_state.`, ...) are kept in the model so the
## nested submodels stay distinct; this map only relabels them for display.
display_names = Dict{Symbol, String}(
    Symbol("rt_state.sigma_rw") => "Rt step size",
    Symbol("rt_state.intervention_effect") => "intervention effect",
    Symbol("gi_state.α") => "generation interval shape",
    Symbol("gi_state.θ") => "generation interval scale",
    Symbol("inc_state.delay_mean") => "incubation period mean",
    Symbol("inc_state.delay_sd") => "incubation period SD",
    Symbol("cases_state.report_state.α") => "onset-to-report shape",
    Symbol("cases_state.report_state.θ") => "onset-to-report scale",
    Symbol("deaths_state.od_state.oa.α") => "onset-to-admission shape",
    Symbol("deaths_state.od_state.oa.θ") => "onset-to-admission scale",
    Symbol("deaths_state.od_state.ad.α") => "admission-to-death shape",
    Symbol("deaths_state.od_state.ad.θ") => "admission-to-death scale",
    Symbol("exports_state.detect_state.α") => "onset-to-detection shape",
    Symbol("exports_state.detect_state.θ") => "onset-to-detection scale",
    Symbol("confirmed_state.receipt_state.d.delay_mean") => "report-to-receipt mean",
    Symbol("confirmed_state.receipt_state.d.delay_sd") => "report-to-receipt SD",
    :isolation_bvd_los_mean => "isolation BVD length-of-stay mean",
    :isolation_ruleout_los_mean => "isolation non-BVD rule-out stay mean",
    :recovery_delay_mean => "confirmation-to-recovery mean",
    Symbol("exports_state.travel_state.daily_travellers") => "daily travellers");

#md # ```@raw html
#md # </details>
#md # ```

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
    "frozen (1wk back)" => frozen_lastweek.chn, #hide
    (RUN_SENSITIVITY ? #hide
     ["delay sensitivity" => chn_joint_community_delay, #hide
        "clock sensitivity" => chn_joint_fast_clock] : [])...) #hide

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
#         \Pr(X_d - X_c \le T - t)}, \tag{35}
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
# We forecast the two confirmed DRC streams (laboratory-confirmed cases and
# confirmed deaths) as the forecast targets, and also the isolation/treatment
# beds and the cumulative recovered total. For the beds we project both the
# bed demand (the need a week ahead, under unconstrained supply, the cut-off
# demand grown by the horizon factor like the case inflow) and the
# supply-limited occupancy that demand produces against the bed capacity. The
# gap between them is the projected bed shortfall, the quantity of interest
# if bed occupancy is supply-constrained. The suspected case and death
# streams are no longer published, so they are not shown as
# targets. Exports are not forecast either, since cross-border travel is
# unlikely to continue at its baseline rate, so the forward travel rate the
# export model relies on no longer holds. The figure is shown in the
# [one-week-ahead forecast results](@ref "One-week-ahead forecast results")
# below.
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
# reused to compare against McCabe et al. at the cut-offs they used, and to
# trace how the estimate moved as the situation reports accrued, re-fitting
# the renewal model frozen at each release date. The helper below performs
# one frozen joint re-fit and is reused by the forecast validation,
# estimate-evolution and matched-in-time results.

#md # ```@raw html
#md # <details><summary>Frozen-fit helper (reused by the forecast validation, evolution and matched-in-time sections)</summary>
#md # ```

## fit_frozen_joint and the frozen re-fits are defined and run in the setup
## block above.

#md # ```@raw html
#md # </details>
#md # ```

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

# The figure below shows three modelled cumulative quantities over time, one
# per row, all latent.
# The top row is infections, the middle row symptom onsets and the bottom
# row deaths.
# The left column is the modelled cumulative trajectory with its 50% and 90%
# intervals; the right column is the posterior density of the current cut-off
# cumulative.
# The infection density on the top right is the headline outbreak size, a
# count of infections rather than reported cases.
# No observed counts are overlaid, since each quantity sits upstream of the
# ascertainment, confirmation and reporting that produce the observations.

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
# Read as a calendar date, the age places the outbreak start at the
# cut-off date minus the age.
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
_rt_start_plot = clamp(
    obs.n - round(Int, obs.tmrca_days) + RENEWAL_START_LEAD, 1, obs.n);
rt_fig = plot_rt(chn_joint;
    n = obs.n, breakpoint = _BREAKPOINT,
    rt_start = _rt_start_plot,
    rt_walk_start = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n),
    as_of_date = string(obs.cutoff), seeding = obs.seeding,
    ramp = 21.0);

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
# The delays from symptom onset to each observed event: onset to a suspected
# case being reported, onset to death, onset to detection abroad (the export
# model uses the same onset-to-admission delay as the report), and the
# report-to-laboratory receipt delay.
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

# ### Posterior predictive checks
#
# A posterior predictive check draws replicated observations from the
# fitted joint model and compares them to the observed counts.
# The checks read in three groups, following the generative order of the
# model.
# The first is the infection process, the latent infections, symptom
# onsets and deaths that drive every stream; these are not observed
# directly, so they carry no genuine replicate and are shown as the
# estimated cumulative trajectories in the
# [joint model estimates](@ref "Joint model estimates") figure rather than
# checked against data here.
# The second is the surveillance data, the dated DRC streams that are real
# per-vintage observations: cumulative suspected cases, the daily new-suspect
# inflow, the daily isolation-bed occupancy, confirmed cases, suspected
# deaths, confirmed deaths, recovered-among-confirmed and specimens analysed.
# The third is the exports, the cross-border imported cases and deaths
# detected in Uganda.
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
        recovered_history = _days_only(obs.recovered_history),
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = _days_only(obs.confirmed_deaths_history),
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
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
isolation_panel = (;
    title = "Patients in isolation",
    dates = _vintage_dates(obs.isolation_history.days),
    replicates = _vintage_replicates(
        pp_joint, @varname(isolation.obs)),
    observed = obs.isolation_history.counts,
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

## Each panel runs to its own last vintage: the suspected case and death
## streams freeze at 26 May (their last stable vintage) while the
## laboratory-confirmed streams keep reporting to the cut-off, so the
## confirmed panels show the full series the model is fitting, not just the
## window the suspected streams cover.
joint_vintage_ppc_fig = plot_vintage_conditional_ppc(
    [reported_panel, suspected_daily_panel, isolation_panel, confirmed_panel,
    deaths_panel, suspected_daily_deaths_panel, confirmed_deaths_panel,
    recovered_panel, tests_analysed_panel, tests_analysed_daily_panel]);

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
    recovered_panel, tests_analysed_panel, tests_analysed_daily_panel]);

#md # ```@raw html
#md # </details>
#md # ```

joint_vintage_incidence_fig #hide

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

# The four quantities side by side: the delay-corrected confirmed CFR, the
# structural CFR, the uncorrected modelled confirmed ratio and the naive
# observed confirmed ratio.

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
# The cumulative and new expected counts by $T + 7$ for the two confirmed DRC
# streams (laboratory-confirmed cases and confirmed deaths), from the
# no-change projection defined in the methods
# [one-week-ahead forecast](@ref "One-week-ahead forecast").
# The suspected reported cases and deaths are no longer reported, so they are
# not shown as forecast targets.
# The forecast also projects the isolation/treatment beds — both the bed
# demand a week ahead (the need under unconstrained supply) and the
# supply-limited occupancy it produces, whose gap is the projected bed
# shortfall — and the cumulative recovered total.

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

forecast_summary #hide

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

# The observation figure shows the new confirmed cases and confirmed deaths
# over the horizon.

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

# ### Forecast validation (last week versus now)
#
# How last week's forecast held up against the data since observed, using the
# frozen re-fit and one-week projection defined in the methods
# [forecast-versus-frozen evaluation](@ref
# "Forecast-versus-frozen evaluation"). The frozen fit also conditions on
# the isolation beds, so the projected bed occupancy is scored against the
# beds held a week later. The bed validation is weak at a one-week-back freeze:
# the reported occupancy rate starts only on 9 June, so the capacity has no
# implied-capacity anchor and rides its random walk back to the freeze date,
# widening the projected bed interval.

#md # ```@raw html
#md # <details><summary>Fit one week back and validate the one-week-ahead forecast</summary>
#md # ```

## frozen_lastweek is computed in the setup block above.
validation_forecast = forecast_reported(frozen_lastweek.chn;
    horizon = 7,
    obs_cases = frozen_lastweek.o.reported_cases,
    obs_deaths = frozen_lastweek.o.total_deaths,
    obs_confirmed = frozen_lastweek.o.confirmed_cases,
    obs_confirmed_deaths = frozen_lastweek.o.confirmed_deaths);

## The observed beds at the current cut-off (the forecast target), so the
## frozen-fit bed forecast is scored against what the beds actually held.
_obs_beds = isempty(obs.isolation_history.counts) ? missing :
            obs.isolation_history.counts[end]
validation_table = forecast_vs_truth(validation_forecast;
    confirmed = obs.confirmed_cases,
    confirmed_deaths = obs.confirmed_deaths,
    isolation = _obs_beds);

#md # ```@raw html
#md # </details>
#md # ```

validation_table #hide

# The observation panels histogram the one-week-ahead cumulative forecast
# made from the frozen fit, with the 90% predictive interval shaded and the
# count actually observed by the current cut-off drawn as a dashed black
# rule.

#md # ```@raw html
#md # <details><summary>Forecast-versus-observed plot</summary>
#md # ```

validation_fig = plot_forecast_vs_truth(validation_forecast;
    confirmed = obs.confirmed_cases,
    confirmed_deaths = obs.confirmed_deaths,
    baseline_confirmed = frozen_lastweek.o.confirmed_cases,
    baseline_confirmed_deaths = frozen_lastweek.o.confirmed_deaths);

#md # ```@raw html
#md # </details>
#md # ```

validation_fig #hide

# The bed panel scores last week's projected isolation-bed occupancy against
# the beds actually occupied now (the dashed rule). At a one-week-back freeze
# the capacity has no implied-capacity anchor (the occupancy rate starts only
# on 9 June), so the projection rides the capacity random walk back to the
# freeze date and its interval is wide.

#md # ```@raw html
#md # <details><summary>Bed forecast-versus-observed plot</summary>
#md # ```

validation_beds_fig = plot_forecast_beds_vs_truth(validation_forecast;
    isolation = _obs_beds);

#md # ```@raw html
#md # </details>
#md # ```

validation_beds_fig #hide

# The latent quantities are not observed, so they are scored distribution
# against distribution: what the frozen fit forecast for the past week's new
# infections, onsets and deaths against what the current fit now estimates
# for the same window.

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
    deaths_latent_new = _now_new(chn_joint, :cumulative_expected_deaths))

validation_latent_fig = plot_forecast_vs_truth_latent(
    validation_forecast; now = now_latent);

#md # ```@raw html
#md # </details>
#md # ```

validation_latent_fig #hide

# ### Outbreak size estimated by each data stream
#
# Each data stream constrains the latent outbreak size differently.
# The table below puts the posteriors over the infection count side by side,
# the single-stream fits and the joint, to show what each stream implies on
# its own and what the joint adds.

#md # ```@raw html
#md # <details><summary>Per-stream infection-count table</summary>
#md # ```

streams_C_table = streams_table(
    "exports" => posterior_C_exports,
    "deaths (DRC)" => posterior_C_deaths,
    "cases (DRC)" => posterior_C_cases,
    "confirmed (DRC)" => posterior_C_confirmed,
    "isolation (DRC)" => posterior_C_treatment,
    "joint" => posterior_C_joint);

#md # ```@raw html
#md # </details>
#md # ```

streams_C_table #hide

# Each single-stream fit projects its outbreak size out to the cut-off, even
# for streams whose data stops earlier.
# The first figure shows each fit's cumulative-infection trajectory over the
# grid as 50% and 90% credible ribbons, with a dotted vertical rule in each
# stream's colour where that stream's data stops reporting.
# The suspected case and death streams freeze at 26 May while the exports and
# confirmed streams run on, so the part of each ribbon beyond its rule is the
# model projecting forward from the last data it saw.

#md # ```@raw html
#md # <details><summary>Per-stream projected-trajectory plot</summary>
#md # ```

## Per-draw cumulative-infection trajectory carried by each single-stream
## fit out to the cut-off on day `n`, so streams whose data ends earlier are
## still projected to today.
function _cuminf(chn)
    mat = chn[:cumulative_infections]
    return [collect(v) for v in vec(collect(mat))]
end
## Grid day a stream's data last reports, used for the dotted rule. The
## suspected case and death histories freeze at 26 May; exports and confirmed
## run to the cut-off.
_last_day(days) = isempty(days) ? nothing : maximum(days)

stream_traj_fig = plot_stream_trajectories(
    [
        (; label = "exports", trajs = _cuminf(chn_exports),
            last_day = _last_day(vcat(obs.export_case_days,
                obs.export_death_days)), colour = :seagreen),
        (; label = "deaths (DRC)", trajs = _cuminf(chn_deaths),
            last_day = _last_day(obs.deaths_history.days),
            colour = :firebrick),
        (; label = "cases (DRC)", trajs = _cuminf(chn_cases),
            last_day = _last_day(obs.reported_history.days),
            colour = :steelblue),
        (; label = "confirmed (DRC)", trajs = _cuminf(chn_confirmed),
            last_day = _last_day(obs.confirmed_history.days),
            colour = :goldenrod),
        (; label = "isolation (DRC)", trajs = _cuminf(chn_treatment),
            last_day = _last_day(obs.isolation_history.days),
            colour = :darkorange)];
    n = obs.n, seeding = obs.seeding);

#md # ```@raw html
#md # </details>
#md # ```

stream_traj_fig #hide

# The second figure is the posterior density of each fit's cumulative
# infection count at the cut-off.
# The confirmed-cases-only stream is ill-defined on its own, so its posterior
# runs far wider than the rest.
# The horizontal axis is scaled to a multiple of the joint-fit 90% upper
# bound, the estimate that constrains every stream together, so the bulk of
# the joint and the other streams stays visible rather than being flattened
# by the confirmed-only tail.

#md # ```@raw html
#md # <details><summary>Cut-off infection-count density plot</summary>
#md # ```

## Scale the x-axis to twice the joint-fit 90% upper bound, so the joint and
## the streams that track it read clearly while the confirmed-only tail runs
## off the axis rather than dominating it.
density_xmax = 2.0 * quantile(posterior_C_joint, 0.95)

cumulative_density_fig = plot_cumulative_cases(
    "exports" => posterior_C_exports,
    "deaths (DRC)" => posterior_C_deaths,
    "cases (DRC)" => posterior_C_cases,
    "confirmed (DRC)" => posterior_C_confirmed,
    "isolation (DRC)" => posterior_C_treatment,
    "joint" => posterior_C_joint;
    scenarios = [], xmax = density_xmax);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_density_fig #hide

# The third figure is the reproduction number each stream implies on its own,
# one panel per stream with the joint fit overlaid in grey as the reference.
# Each panel draws 30/60/90% credible ribbons with no median line, matching
# the band style used elsewhere; the window, the response markers and the
# cut-off match the joint Rt figure above.

#md # ```@raw html
#md # <details><summary>Per-stream implied-Rt plot</summary>
#md # ```

## The per-stream fits walk Rt from day 1 (the default `rt_start`), while the
## joint walks from `RT_WALK_LEAD` days before the first situation report; the
## shared `display_start` is the joint renewal start so every stream reads over
## the same established window. `ramp` matches the joint Rt figure.
_rt_walk_start_joint = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n);
stream_rt_fig = plot_rt_streams(
    [
        (; label = "exports", chn = chn_exports, rt_start = 1,
            rt_walk_start = 1, colour = :seagreen),
        (; label = "deaths (DRC)", chn = chn_deaths, rt_start = 1,
            rt_walk_start = 1, colour = :firebrick),
        (; label = "cases (DRC)", chn = chn_cases, rt_start = 1,
            rt_walk_start = 1, colour = :steelblue),
        (; label = "confirmed (DRC)", chn = chn_confirmed, rt_start = 1,
            rt_walk_start = 1, colour = :goldenrod),
        (; label = "isolation (DRC)", chn = chn_treatment, rt_start = 1,
            rt_walk_start = 1, colour = :darkorange)];
    joint = (; label = "joint", chn = chn_joint, rt_start = _rt_start_plot,
        rt_walk_start = _rt_walk_start_joint),
    n = obs.n, breakpoint = _BREAKPOINT,
    as_of_date = string(obs.cutoff), seeding = obs.seeding,
    display_start = _rt_start_plot, ramp = 21.0);

#md # ```@raw html
#md # </details>
#md # ```

stream_rt_fig #hide

# The frozen re-fits below freeze the renewal data to an earlier cut-off
# and re-fit, so that a change driven by newer data can be distinguished from
# one driven by a change of method.
# Each uses the full headline settings (1000 draws across two chains),
# reusing the frozen-fit helper defined above.

#md # ```@raw html
#md # <details><summary>Freeze the renewal data to a cut-off and re-fit</summary>
#md # ```

## Frozen re-fits and released_df are prepared in the setup block above.

#md # ```@raw html
#md # </details>
#md # ```

# ### Estimate evolution across releases
#
# How the outbreak-size estimate moved as the situation reports accrued.
#
# The project publishes a tagged results release at each data cut-off
# (<https://github.com/epiforecasts/BVDOutbreakSize/releases>), bundling
# the posterior draws and input data.
# The released series, in blue, is the project's published estimate at each
# release: the closed-form integral model up to v1.3.0, then the renewal
# model from v1.4.0 on.
# Each release is its own fit, so it is drawn as a discrete estimate, a
# median with nested 30/60/90% interval bars, rather than a ribbon.
# The renewal series, in red, is the renewal model re-fit frozen at each
# integral-era release cut-off.
# The renewal-era releases already are renewal fits, so they carry no frozen
# re-fit.
# The current-data, current-model estimate is drawn in green as the
# cumulative-infection trajectory over time, a single fit shown across the
# period so the latest estimate reads against the earlier ones.
# Each release date is marked with a dotted vertical rule.

#md # ```@raw html
#md # <details><summary>Released estimates and frozen renewal re-fits</summary>
#md # ```

## Released median and 30/60/90% intervals per release, from
## `data/released_estimates.csv`. Each tuple is
## `(date, median, lo30, hi30, lo60, hi60, lo90, hi90)`.
release_evolution = [(string(r.date), r.median, r.lo30, r.hi30, r.lo60, r.hi60,
                         r.lo90, r.hi90) for r in eachrow(released_df)]

## The current renewal model re-fit frozen at each integral-era release
## cut-off, each its own discrete estimate. Each tuple carries the median
## and 30/60/90% credible bounds.
function _ci369(xs)
    q(p) = round(Int, quantile(xs, p))
    (q(0.5), q(0.35), q(0.65), q(0.20), q(0.80), q(0.05), q(0.95))
end
## Reuse the one-week-back frozen fit (already run for forecast validation)
## as an additional recent renewal point, so the evolution plot shows the
## current-vintage frozen estimate without an extra re-fit.
frozen_by_cutoff[validation_cutoff] = frozen_lastweek
renewal_frozen = [(c, _ci369(frozen_C(c))...)
                  for c in sort(union(frozen_evolution_cutoffs,
    [validation_cutoff]))]

## The current-data, current-model estimate as the cumulative-infection
## trajectory over the day grid (one calendar date per grid day, day 1 is
## the seeding date), summarised by per-day 30/60/90% credible bounds. This
## is the same latent quantity the cumulative-trajectory figure shows, so
## the current estimate rises over time on the release-date axis instead of
## sitting flat. Drawn against calendar dates, it lines up with the
## release and frozen-renewal points.
infection_trajectory = let
    mat = chn_joint[:cumulative_infections]
    trajs = [collect(v) for v in vec(collect(mat))]
    ## Only over the comparison window — from the earliest release date to the
    ## cut-off — not back to the seeding date.
    start_day = obs.n - value(obs.cutoff - Date(release_evolution[1][1]))
    days = max(start_day, 1):obs.n
    dates = [obs.seeding + Day(d - 1) for d in days]
    q(d, p) = quantile(Float64[t[d] for t in trajs], p)
    (dates,
        [q(d, 0.35) for d in days], [q(d, 0.65) for d in days],
        [q(d, 0.20) for d in days], [q(d, 0.80) for d in days],
        [q(d, 0.05) for d in days], [q(d, 0.95) for d in days])
end

evolution_fig = plot_estimate_evolution(release_evolution;
    renewal = renewal_frozen,
    trajectory = infection_trajectory,
    title = "Outbreak-size estimate as data accrued");

#md # ```@raw html
#md # </details>
#md # ```

evolution_fig #hide

# ### Comparison with McCabe et al.
#
# Our model is a discrete-time renewal model with a time-varying
# reproduction number and every data stream fitted jointly.
# McCabe et al. published their estimates as scenarios at fixed
# situation-report cut-offs, each scenario carrying a 95% confidence
# interval.
# We show both their reports, the 18 May report and the 20 May update, with
# those intervals.
# The two reports share the same geographic-spread scenarios from exported
# cases and travel volume, so those appear once.
# Their back-calculation-from-deaths scenarios differ between the reports,
# since the 18 May report used 88 reported deaths and the 20 May update 131,
# with a corrected set of case-fatality ratios.
# McCabe's scenarios estimate cumulative cases at their report dates, though
# their report is not fully explicit about whether this is symptomatic cases
# or all infections.
# We take the like-for-like quantity to be our cumulative symptom onsets, the
# symptomatic cases, on the same dates, rather than the latent infections
# (which include the not-yet-symptomatic) or our current cut-off total.
# We read our value off the joint fit's cumulative-onset trajectory at the
# grid day for each report date and show it with its credible interval, the
# 18 May report against our 18 May value, the 20 May update against our
# 20 May value, and the 27 May Lancet publication against our 27 May value,
# so each scenario sits beside our estimate for the date it was made.

#md # ```@raw html
#md # <details><summary>McCabe scenarios with uncertainty against our estimates</summary>
#md # ```

function _ci90row(xs)
    (round(Int, quantile(xs, 0.5)),
        round(Int, quantile(xs, 0.05)),
        round(Int, quantile(xs, 0.95)))
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
    _ci90row(Float64[t[d] for t in _onset_trajs])
end

## McCabe scenarios with their reported 95% confidence intervals, grouped by
## report date, then our modelled cumulative onsets on the matching
## report date with its credible interval.
mccabe_rows = [(label, mean, lo, hi)
               for (_, label, mean, lo, hi) in REPORT_SCENARIOS_CI]
mccabe_groups = [date == "2026-05-18" ? "McCabe et al. (18 May)" :
                 date == "2026-05-20" ? "McCabe et al. (20 May update)" :
                 "McCabe et al. (27 May, Lancet)"
                 for (date, _, _, _, _) in REPORT_SCENARIOS_CI]

ours_rows = [("Renewal onsets on 18 May", _ours_on("2026-05-18")...),
    ("Renewal onsets on 20 May", _ours_on("2026-05-20")...),
    ("Renewal onsets on 27 May", _ours_on("2026-05-27")...)]
ours_groups = fill("Our renewal estimate", 3)

matched_rows = vcat(mccabe_rows, ours_rows)
matched_groups = vcat(mccabe_groups, ours_groups)

matched_comparison_fig = plot_estimate_comparison(matched_rows;
    xlabel = "Cumulative cases / infections",
    groups = matched_groups,
    group_colours = ["McCabe et al. (18 May)" => :grey,
        "McCabe et al. (20 May update)" => :black,
        "McCabe et al. (27 May, Lancet)" => :steelblue,
        "Our renewal estimate" => :firebrick]);

#md # ```@raw html
#md # </details>
#md # ```

matched_comparison_fig #hide

# The McCabe scenarios are outbreak-size estimates, the same quantity our
# renewal model and the released integral model report.
# Their 95% confidence intervals come from exact negative-binomial counts
# for the geographic-spread method and a Poisson likelihood profile for the
# back-calculation from deaths.
#
# Side-by-side outbreak-size intervals for the two frozen fits and the
# current-data fit, so the shift with the data cut-off reads off directly.

#md # ```@raw html
#md # <details><summary>Frozen-fit C_T table</summary>
#md # ```

frozen_streams_table = streams_table(
    "frozen 20 May" => frozen_C("2026-05-20"),
    "frozen 23 May" => frozen_C("2026-05-23"),
    "frozen 27 May" => frozen_C("2026-05-27"),
    "current data" => posterior_C_joint);

#md # ```@raw html
#md # </details>
#md # ```

frozen_streams_table #hide

# ### Delay sensitivity
#
# **Paused due to compute constraints.** This sensitivity re-fit is not run in
# the current build to keep the docs build within compute limits. The method
# and code are unchanged and it can be re-enabled (see the sensitivity gate in
# the setup block); the description below documents what it does.
#
# The death stream dates the outbreak from how far deaths lag symptom onset,
# so the assumed onset-to-death delay sets the implied infection count.
# The baseline uses the hospital-pathway delay from the Isiro 2012 line-list
# reanalysis (onset to admission then admission to death, implied mean about
# 12 d).
# We re-fit the joint model under the community-pathway delay from the same
# reanalysis, the delay for deaths that occur in the community without a
# recorded admission, which is shorter (implied mean about 8 d).
# Both pathways come from the line list, so this varies the actual delay
# assumption rather than an arbitrary scenario.
# The re-fit uses the full headline settings (1000 draws across two chains).
#
# The infection count to date shifts with the assumed delay, and the
# table and overlaid densities below show how far.

#md # ```@raw html
#md # <details><summary>Re-fit the joint under the community-pathway onset-to-death delay</summary>
#md # ```

## refit_joint_variant, deaths_community_delay and the sensitivity re-fits are
## defined and run (when enabled) in the setup block above.
posterior_C_community_delay = RUN_SENSITIVITY ?
                              vec(Array(chn_joint_community_delay[:C_T])) : nothing

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Delay-sensitivity infection-count table</summary>
#md # ```

delay_sensitivity_table = RUN_SENSITIVITY ?
                          streams_table("baseline (hospital pathway)" => posterior_C_joint,
    "community pathway" => posterior_C_community_delay) :
                          Markdown.md"_Delay sensitivity is paused due to compute constraints._"

#md # ```@raw html
#md # </details>
#md # ```

delay_sensitivity_table #hide

#md # ```@raw html
#md # <details><summary>Delay-sensitivity infection-count density plot</summary>
#md # ```

delay_sensitivity_fig = RUN_SENSITIVITY ?
                        plot_cumulative_cases(
    "baseline (hospital pathway)" => posterior_C_joint,
    "community pathway" => posterior_C_community_delay; scenarios = []) :
                        Markdown.md"_Delay sensitivity is paused due to compute constraints._"

#md # ```@raw html
#md # </details>
#md # ```

delay_sensitivity_fig #hide

# ### Clock-rate sensitivity
#
# **Paused due to compute constraints.** This sensitivity re-fit is not run in
# the current build to keep the docs build within compute limits. The method
# and code are unchanged and it can be re-enabled (see the sensitivity gate in
# the setup block); the description below documents what it does.
#
# The whole outbreak-age estimate rests on the genetic bound, the oldest
# date the common ancestor of the sequenced cases can sit, which is set by
# the assumed molecular clock rate.
# The baseline uses the slower clock rate of $1.2\times10^{-3}$
# substitutions per site per year, the rate of the 2013-2016 West African
# Ebola epidemic, which dates the common ancestor to 25 March 2026.
# The sequencing source also reports a faster early-epidemic rate of
# $1.9\times10^{-3}$ substitutions per site per year, which dates the common
# ancestor about two and a half weeks more recently, to 11 April 2026,
# without favouring either [virological2026](@cite).
# We re-fit the joint model under the faster clock and compare the
# infection count to date and the outbreak age.
# The re-fit uses the full headline settings (1000 draws across two chains).

#md # ```@raw html
#md # <details><summary>Re-fit the joint under the faster clock rate</summary>
#md # ```

## clock_alt_offset, tmrca_days_alt and the faster-clock re-fit are defined and
## run (when enabled) in the setup block above.
posterior_C_fast_clock = RUN_SENSITIVITY ?
                         vec(Array(chn_joint_fast_clock[:C_T])) : nothing
T_baseline_clock = vec(Array(chn_joint[:T]))
T_fast_clock = RUN_SENSITIVITY ? vec(Array(chn_joint_fast_clock[:T])) : nothing

#md # ```@raw html
#md # </details>
#md # ```

# The infection count to date under the two clock rates, side by side.

#md # ```@raw html
#md # <details><summary>Clock-rate infection-count table</summary>
#md # ```

clock_sensitivity_C_table = RUN_SENSITIVITY ?
                            streams_table("baseline clock" => posterior_C_joint,
    "faster clock" => posterior_C_fast_clock) :
                            Markdown.md"_Clock-rate sensitivity is paused due to compute constraints._"

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_C_table #hide

#md # ```@raw html
#md # <details><summary>Clock-rate infection-count density plot</summary>
#md # ```

clock_sensitivity_C_fig = RUN_SENSITIVITY ?
                          plot_cumulative_cases("baseline clock" => posterior_C_joint,
    "faster clock" => posterior_C_fast_clock; scenarios = []) :
                          Markdown.md"_Clock-rate sensitivity is paused due to compute constraints._"

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_C_fig #hide

# The outbreak age, the number of days from seeding to the cut-off, under
# the two clock rates.
# A more recent common ancestor permits a younger outbreak.

#md # ```@raw html
#md # <details><summary>Clock-rate outbreak-age table</summary>
#md # ```

clock_sensitivity_T_table = RUN_SENSITIVITY ?
                            streams_table("baseline clock" => T_baseline_clock,
    "faster clock" => T_fast_clock; digits = 0) :
                            Markdown.md"_Clock-rate sensitivity is paused due to compute constraints._"

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_T_table #hide

#md # ```@raw html
#md # <details><summary>Clock-rate outbreak-age density plot</summary>
#md # ```

clock_sensitivity_T_fig = RUN_SENSITIVITY ?
                          plot_density_overlay("baseline clock" => T_baseline_clock,
    "faster clock" => T_fast_clock;
    xlabel = "Outbreak age (days before cut-off)",
    title = "Posterior outbreak age by clock rate") :
                          Markdown.md"_Clock-rate sensitivity is paused due to compute constraints._"

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_T_fig #hide

# ## Saving results
#
# The tables above are written to an `output/` directory at the repo
# root so they can be archived and shared. On every push to `main` a
# GitHub Actions workflow regenerates these files and publishes them
# as a GitHub Release, downloadable from the repository's releases
# page (<https://github.com/epiforecasts/BVDOutbreakSize/releases>).
# The release bundles the summary tables, a thinned set of
# posterior draws, the latent symptom-onset ("symptomatic cases")
# trajectory over time, and a copy of the input `observations.toml` so
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
CSV.write(joinpath(output_dir, "cumulative_cases_by_stream.csv"),
    streams_C_table)
CSV.write(joinpath(output_dir, "frozen_matched_cutoffs.csv"),
    frozen_streams_table)

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
)[1:10:end, :]
CSV.write(joinpath(output_dir, "posterior_draws.csv"), posterior_draws);

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
CairoMakie.save(joinpath(dashboard_dir, "rt_streams.png"), stream_rt_fig)
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
