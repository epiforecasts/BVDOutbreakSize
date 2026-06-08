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
# - *Time-varying reproduction number.* $R_t$ follows a weekly Gaussian
#   random walk on the log scale, interpolated within weeks, with a
#   logistic outbreak-response ramp of about three weeks at the first WHO
#   situation report (18 May 2026). McCabe et al. use one constant
#   exponential growth rate.
# - *Joint posterior rather than scenario estimates.* The reproduction
#   number, case-fatality ratio, all delays, traveller volume and
#   surveillance dispersion have priors and are sampled together. McCabe
#   et al. [mccabe2026](@cite) fix each and report a set of scenarios.
# - *Two-phase seeding with a wide, genetically-floored outbreak age.*
#   A single import grows through an unobserved cryptic exponential phase
#   to a magnitude set by a wide prior on the doubling count, at a rate
#   derived from the established reproduction number, before the renewal
#   process takes over. The genetic time to the most recent common
#   ancestor floors the cryptic duration from below. McCabe et al. fix the
#   start from a single seed.
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
# - *Per-vintage time-series fitting.* The DRC streams are fitted as
#   cumulative series of between-vintage increments across successive
#   sitreps, which sharpens $R_t$. McCabe et al. condition on a single
#   cumulative total.
# - *Ascertainment estimated.* We estimate the fraction of cases each
#   surveillance system reports jointly with the outbreak size. McCabe et
#   al. have no ascertainment component.
# - *Comparison against published scenarios.* The model is set beside the
#   McCabe et al. [mccabe2026](@cite) scenario estimates as an external
#   sense-check, matched in time at the cut-off each scenario was computed,
#   while $C_T$ (the latent infection count, summed from the renewal
#   trajectory) is the headline quantity reported separately.
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
#   the reproduction number over the window, but says little about the
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
#   which can understate uncertainty. We have not checked whether the
#   streams imply conflicting outbreak sizes.
#
# **Model assumptions and design**
#
# - *Inherits McCabe et al.'s epidemiological assumptions.* A single
#   zoonotic seed, an assumed generation interval, no spatial structure
#   beyond the Ituri / Nord Kivu split, and no depletion of
#   susceptibles. The onset-to-death delay is anchored on Isiro 2012 and
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
using DataFrames: DataFrame
import CSV
using Random
using Markdown
using Dates: Date, Day, value
using BVDOutbreakSize
import CairoMakie

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
# laboratory, at the report date. We extracted these figures from the
# written situation-report PDFs (archived by INRB-UMIE
# [inrb_umie_2026](@cite)) using a language model, with a second pass to
# re-read them, rather than the published per-zone CSVs. The zone sums in
# the CSVs are inconsistent with the national headline totals because they
# drop counts not yet attributed to a zone, so they understate the count.
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
# the cut-off total. See `data/observations.toml` for the per-stream
# sources.

#md # ```@raw html
#md # <details><summary>Building the per-date time-series table</summary>
#md # ```

vintage_table = let
    ## Each history carries grid day-indices and counts; key the counts
    ## by calendar date so every stream lines up in one table.
    bydate(h) = Dict(grid_date(d) => c for (d, c) in zip(h.days, h.counts))
    streams = (
        suspected_cases = bydate(obs.reported_history),
        suspected_deaths = bydate(obs.deaths_history),
        confirmed_cases = bydate(obs.confirmed_history),
        confirmed_deaths = bydate(obs.confirmed_deaths_history),
        specimens_received = bydate(obs.tests_received_history),
        specimens_analysed = bydate(obs.lab_history)
    )
    dates = sort(collect(union((keys(s) for s in streams)...)))
    at(s) = [haskey(s, d) ? s[d] : missing for d in dates]
    DataFrame(
        date = dates,
        suspected_cases = at(streams.suspected_cases),
        suspected_deaths = at(streams.suspected_deaths),
        confirmed_cases = at(streams.confirmed_cases),
        confirmed_deaths = at(streams.confirmed_deaths),
        specimens_received = at(streams.specimens_received),
        specimens_analysed = at(streams.specimens_analysed)
    )
end;

#md # ```@raw html
#md # </details>
#md # ```

vintage_table #hide

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
# where $g$ is the discretised generation-interval PMF (indexed from
# lag 1, so an infectee is always infected strictly after its infector)
# and $R_t$ is the per-day reproduction number. We never observe
# infections directly. Each data stream observes a thinned, delayed or
# transformed view of the same latent incidence.
#
# The model is assembled from small reusable Turing [ge2018turing](@cite)
# submodels, each owning the maths and priors for one part of the
# generative process. We describe them in generative order, from the
# infection process through the epidemiological delays to the observation
# streams. The implementation uses Mooncake [mooncake_jl](@cite)
# reverse-mode automatic differentiation, CensoredDistributions for delay
# discretisation, FlexiChains for chain handling, and PairPlots
# [pairplots_jl](@cite) with AlgebraOfGraphics [danisch2021makie](@cite)
# for the figures. Each submodel's source is shown in the collapsible
# block beneath its prose.
#
# The table below shows which parameters feed each observation submodel.
# The *received* column is the received-specimen volume, the laboratory
# stream fitted as a count; the *confirmed* positives are scored as a
# Binomial of the observed analysed denominator with a positivity linked to
# the composition of the suspected pool, so the laboratory data help
# identify the non-BVD background. The *conf. deaths* column thins the
# suspected deaths by the case composition:
#
# | Parameter | Exports | Deaths | Cases | Received | Confirmed | Conf. deaths | Export deaths |
# |---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
# | Reproduction number $R_t$ | ● | ● | ● | ● | ● | ● | ● |
# | Generation interval | ● | ● | ● | ● | ● | ● | ● |
# | Incubation period | ● | ● | ● | ● | ● | ● | ● |
# | Seed $I_0$ | ● | ● | ● | ● | ● | ● | ● |
# | Onset-to-death delay |  | ● |  |  |  |  | ● |
# | Case-fatality ratio |  | ● |  |  |  |  | ● |
# | Onset-to-report delay |  |  | ● | ● | ● | ● |  |
# | Receipt delay |  |  |  | ● | ● |  |  |
# | Onset-to-detection delay | ● |  |  |  |  |  |  |
# | Assay sensitivity / specificity |  |  |  |  | ● |  |  |
# | Severity enrichment $\delta_0$ |  |  |  |  | ● |  |  |
# | Confirmed-death enrichment $m_{\text{death}}$ |  |  |  |  |  | ● |  |
# | Testing fraction $\tau_{\text{test}}$ |  |  |  | ● | ● |  |  |
# | Background rate $\lambda_{\text{bg}}$ |  |  | ● | ● | ● | ● |  |
# | Surveillance dispersion |  | ● | ● | ● |  |  |  |
# | Ascertainment | ● |  | ● | ● | ● | ● | ● |
# | Traveller volume | ● |  |  |  |  |  | ● |

#md # ```@setup main
#md # using BVDOutbreakSize, CodeTracking, Revise
#md # ```

# #### Infections
#
# The infection process is built in generative order: the reproduction
# number, the generation interval that drives the renewal, the seeding that
# sets the initial infection count, the genetic bound on the outbreak age,
# and the renewal construction that grows the seed forward to the cut-off.
#
# ##### Reproduction number
#
# The reproduction number follows a non-centred Gaussian random walk on the
# log scale, with knots at weekly intervals (day 1, 8, 15, …, $n$):
#
# ```math
# \log R_0 \sim \mathrm{Normal}(\log 1.6,\ 0.10), \qquad
# \sigma_{\text{rw}} \sim \mathrm{Normal}^{+}(0,\ 0.05), \tag{2}
# ```
#
# ```math
# \log R_k = \log R_0 + \sigma_{\text{rw}}
#            \sum_{j=1}^{k} z_j, \quad
# z_j \sim \mathrm{Normal}(0, 1). \tag{3}
# ```
#
# Here $R_0$ is the established reproduction number at the anchor, the value
# the walk starts from. Its prior is anchored on the molecular-clock growth
# estimate for this outbreak: the phylodynamic reanalysis of the first ten
# sequenced genomes [cuomodannenburg2026](@cite) puts the mean epidemic
# doubling time at 15.2-24.5 d (centre 20 d), which under the renewal
# generation interval implies $R_0$ near $1.6$. The single $R_0$ prior
# anchors both the established walk and, through the cryptic growth rate
# derived from it, the seeding, so the outbreak has one growth source. The
# step-size prior assumes a small weekly change in the log reproduction
# number; we use a tight half-normal so that, over the long unobserved
# stretch the sitrep window does not cover, the reproduction number stays
# near the seeding value rather than drifting.
#
# Daily $\log R_t$ is the piecewise-linear interpolation between knots. The
# outbreak response adds a sampled effect shaped by a logistic ramp at the
# first WHO situation report on 18 May 2026. We assume the response takes
# about three weeks (21 days) to take effect, rather than switching at a
# single date, and that it can only reduce transmission, so the effect is
# constrained to be non-positive:
#
# ```math
# \log R_t \mathrel{+}= \delta \cdot
#     \mathrm{logistic}\!\left(\frac{t - t_{\text{bp}}}{21}\right),
# \qquad
# \delta \sim \mathrm{Normal}^{-}(0,\ 0.4). \tag{4}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: rt_walk_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.rt_walk_model(7)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Generation interval
#
# The generation-interval PMF $g$ drives the renewal and, through the
# inverse Euler–Lotka relation, the cryptic growth rate used for seeding.
# It is sampled from a prior centred on the Ebola virus disease serial
# interval as a generation-time proxy (mean 15.3 d, SD 9.3 d; WHO Ebola
# Response Team 2014, NEJM). Both the mean and the SD are sampled, so the
# generation time is estimated around the published value. The lag-0 bin is
# dropped and the remainder renormalised, so an infectee is always strictly
# later than its infector:
#
# ```math
# \mu_g \sim \mathrm{Normal}^{+}(15.3,\ 3.0), \qquad
# \sigma_g \sim \mathrm{Normal}^{+}(9.3,\ 2.0). \tag{5}
# ```
#
# The generation interval uses the shared double-interval-censored LogNormal
# discretisation described with the epidemiological process models below,
# additionally dropping the lag-0 bin so it starts at one day. The source
# reports a generation-time distribution rather than a mean and SD; sampling
# the mean and SD with self-assigned spreads, and truncating the
# discretisation by hand rather than from the source distribution, are
# limitations flagged in issue
# [#224](https://github.com/epiforecasts/BVDOutbreakSize/issues/224).

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

# ##### Seeding
#
# We assume the outbreak started from a single seed case introduced by a
# zoonotic spillover. The initial infection count $I_0$ on the last day of
# the seeding window has a prior centred on a single seed:
#
# ```math
# I_0 \sim \mathrm{Normal}^{+}(0.1,\ 0.1). \tag{6}
# ```
#
# From that seed the outbreak grew deterministically through an unobserved
# cryptic exponential phase, doubling $m$ times before sustained
# transmission was established. The cryptic phase therefore grows the seed
# to $2^m$ infections at the anchor, the day the renewal takes over, over a
# duration $m\,\tau$ with $\tau$ the doubling time:
#
# ```math
# \tau = \frac{\log 2}{r}, \qquad
# T_{\text{cryptic}} = m\,\tau, \qquad
# \text{anchor seed} = 2^m. \tag{7}
# ```
#
# The seed magnitude $2^m$ does not depend on the growth rate $r$, so $r$
# enters only the renewal that follows. We infer $r$ from the established
# reproduction number through the inverse Euler–Lotka relation, the same
# derivation EpiNow2 uses, so the cryptic phase and the renewal share one
# growth source rather than carrying separate growth rates. We place the
# prior directly on the doubling count $m$:
#
# ```math
# m \sim \mathrm{Normal}^{+}(2,\ 3). \tag{8}
# ```
#
# We place the prior on the doubling count and derive the growth rate from
# the reproduction number; putting the prior on the growth rate itself is
# the cleaner formulation, flagged in issue
# [#223](https://github.com/epiforecasts/BVDOutbreakSize/issues/223).

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
#md #     (@code_string BVDOutbreakSize.exponential_growth_model(0.0)), "\n```"))
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
# \qquad \sigma = 15\ \text{d}. \tag{9}
# ```
#
# The anchor sits seven days after the TMRCA day, so the observed window is
# shorter than the TMRCA age. The bound therefore stays informative on the
# cryptic duration, pulling the origin to sit at or before the most recent
# common ancestor and bounding the cryptic phase from below. It is
# one-sided, leaving the age free above the TMRCA. We fix the clock and do
# not propagate cross-outbreak or clock uncertainty.

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
# The renewal recursion runs from the anchor to the cut-off. The anchor is
# the grid day seven days after the genetic TMRCA, past the clock
# uncertainty, where sustained transmission is treated as established; the
# observed window is the span from the anchor to the cut-off:
#
# ```math
# \text{anchor} = n - \text{tmrca}_{\text{days}} + 7, \qquad
# \tau_{\text{obs}} = n - \text{anchor}. \tag{10}
# ```
#
# The pre-anchor grid days are filled by the cryptic exponential curve at
# rate $r$ ending at $2^m$, giving the recursion a full generation interval
# of history. The renewal then grows the trajectory forward under the
# time-varying reproduction number, so the cut-off size is data-driven while
# the doubling count sets only the anchor scale. The total outbreak age is
# the cryptic duration plus the observed window:
#
# ```math
# T = m\,\tau + \tau_{\text{obs}}. \tag{11}
# ```
#
# Cumulative infections are the running sum of the daily series. The cut-off
# cumulative $C_T$ is the headline outbreak size, the current growth rate is
# the day-over-day log-ratio of infections at the cut-off, and the doubling
# time is $\log 2$ divided by that rate.

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

# #### Onset incidence
#
# Infections are convolved with the incubation-period PMF to give daily
# symptom-onset incidence, computed once and consumed by every downstream
# observation stream. The incubation period cannot be fitted from the BDBV
# line list, which has no exposure dates, so we use the Bundibugyo virus
# estimate from the 2007 Uganda outbreak (mean 6.3 d, 95% CI 5.2-7.3,
# $n = 24$; [macneil2010](@cite)). The mean prior reproduces that 95% CI;
# the source reports no interval on the spread, so the SD prior is our own
# weakly-informative choice:
#
# ```math
# \mu_{\text{inc}} \sim \mathrm{Normal}^{+}(6.3,\ 0.54), \qquad
# \sigma_{\text{inc}} \sim \mathrm{Normal}^{+}(3.5,\ 0.8). \tag{12}
# ```
#
# The incubation period uses the same censored discretisation as every
# delay. Onsets map to each observed endpoint through a stream-specific
# onset-to-event delay, defined with the epidemiological process models
# below.

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

# #### Epidemiological process models
#
# Onsets map to each observed endpoint through a stream-specific delay, and
# deaths additionally through the case-fatality ratio. Every delay is sampled
# by its mean and SD, moment-matched to a LogNormal, and discretised to a
# daily PMF over lags $0,\dots,n_{\max}$ by double interval censoring
# [charniga2024](@cite). The LogNormal CDF differentiates cleanly under the
# reverse-mode automatic differentiation, so this is the discretisation route
# for every delay in the convolutions:
#
# ```math
# f_s = F(s + 1) - F(s), \qquad
# F = \text{LogNormal CDF moment-matched to } (\mu, \sigma). \tag{13}
# ```
#
# Whether the censored discretisation is built with CensoredDistributions in
# every case should be confirmed; this is a verification item only and does
# not change behaviour here.

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

# ##### Incubation period
#
# The incubation period is the infection-to-onset delay defined with the
# onset incidence above (mean 6.3 d; [macneil2010](@cite)). It enters the
# infection-to-detection and infection-to-death delays for the export
# streams, where the survival clock runs from infection.
#
# ##### Onset-to-report delay
#
# The delay from symptom onset to a suspected case being reported, centred
# on Ebola surveillance reporting delays:
#
# ```math
# \mu_{\text{rep}} \sim \mathrm{Normal}^{+}(4.5,\ 1.5), \qquad
# \sigma_{\text{rep}} \sim \mathrm{Normal}^{+}(3.6,\ 1.2). \tag{14}
# ```
#
# This delay drives the suspected-case, laboratory and confirmed-death
# streams; its source is shown with the reported-cases submodel below.
#
# ##### Onset-to-death delay
#
# The McCabe et al. report uses the point estimate of
# [rosello2015](@cite). We instead anchor the onset-to-death delay on the
# companion Bayesian reanalysis of the same Isiro 2012 BDBV line list
# [bdbv_linelist_analysis_2026](@cite), which re-estimates the delay with
# uncertainty (posterior mean 11.2 d, SD 5.4 d). The mean and SD priors are
# centred on those values:
#
# ```math
# \mu_d \sim \mathrm{Normal}^{+}(11.2,\ 2.0), \qquad
# \sigma_d \sim \mathrm{Normal}^{+}(5.4,\ 1.5). \tag{15}
# ```
#
# The reanalysis reports a Gamma shape and scale with uncertainty; sampling
# the mean and SD with self-assigned spreads rather than carrying the source
# shape, scale and uncertainty is flagged in issue
# [#224](https://github.com/epiforecasts/BVDOutbreakSize/issues/224). The
# source is shown with the deaths submodel below, where the delay is
# injected.
#
# ##### Onset-to-detection delay
#
# The delay from symptom onset to detection at a point of entry abroad,
# centred on the Ebola onset-to-hospitalisation delay (mean 5.0 d, SD 4.7 d;
# WHO Ebola Response Team 2014, NEJM):
#
# ```math
# \mu_{\text{det}} \sim \mathrm{Normal}^{+}(5.0,\ 2.0), \qquad
# \sigma_{\text{det}} \sim \mathrm{Normal}^{+}(4.7,\ 1.5). \tag{16}
# ```
#
# It drives the exports streams; its source is shown with the exports
# submodel below.
#
# ##### Receipt delay
#
# The delay from a suspected case being reported to its specimen being
# received by the laboratory, centred on a short turnaround with a heavy
# right tail allowing for specimen shipment to a confirmatory laboratory.
# No per-sample outbreak data anchors this, so the prior is our own choice:
#
# ```math
# \mu_{\text{rec}} \sim \mathrm{Normal}^{+}(4.5,\ 2.0), \qquad
# \sigma_{\text{rec}} \sim \mathrm{Normal}^{+}(4.0,\ 1.5). \tag{17}
# ```
#
# It drives the laboratory received-specimen stream; its source is shown
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
# \mathrm{CFR} \sim \mathrm{Beta}(6.6,\ 13.4), \tag{18}
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

# The prior density, with the CDC $0.33$ figure marked, as a sense
# check.

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
# Several priors are shared across streams. The surveillance dispersion ties
# the count likelihoods together; the ascertainment fractions scale the DRC
# and Uganda streams; the laboratory priors govern the confirmation pipeline;
# and the traveller volume sets the export travel rate.
#
# ###### Surveillance dispersion
#
# Passive-surveillance counts (DRC suspected deaths, reported cases and
# received specimens) are modelled with negative-binomial observation error
# sharing one dispersion $k$. Following Stan prior-choice recommendations
# [stan_prior_choice](@cite), the dispersion is sampled on the $1/\sqrt{k}$
# scale:
#
# ```math
# 1/\sqrt{k} \sim \mathrm{Normal}^{+}(0.6,\ 0.2). \tag{19}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: surveillance_dispersion_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.surveillance_dispersion_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ###### Ascertainment
#
# Two surveillance systems detect cases: DRC passive community
# surveillance (the reported suspected-case count) and Uganda's
# point-of-entry / hospital surveillance (the exported-case count). Each
# captures a fraction of the true cases passing through it, and each
# fraction is informed by essentially a single aggregate data point. The
# two ascertainment fractions $p_{\text{DRC}}$ and $p_{\text{Uganda}}$
# share a logit-scale hyperprior with mean $\mu$ and pooling strength
# $\tau$, centred on a reporting fraction of $75\%$, reflecting the active
# case-finding of a declared Ebola response rather than baseline passive
# surveillance:
#
# ```math
# \mu \sim \mathrm{Normal}(\mathrm{logit}(0.75),\ 1),
# \qquad
# \tau \sim \mathrm{Normal}^{+}(0,\ 0.5), \tag{20}
# ```
#
# ```math
# \mathrm{logit}(p_{\text{DRC}}) \sim \mathrm{Normal}(\mu,\ \tau),
# \qquad
# \mathrm{logit}(p_{\text{Uganda}}) \sim \mathrm{Normal}(\mu,\ \tau). \tag{21}
# ```
#
# The cases likelihood uses $p_{\text{DRC}}$; the two Uganda-side
# likelihoods use $p_{\text{Uganda}}$. An independent alternative drops
# the shared hyperprior and gives each system its own fraction.

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
# We model the process of confirming cases explicitly. The testing fraction
# $\tau_{\text{test}}$ is the share of suspected cases routed to the
# laboratory. A truly BVD specimen tests positive with the assay sensitivity
# $s$, and a non-BVD specimen tests positive with the small false-positive
# rate $1 - \mathrm{spec}$ set by the assay specificity. The severity
# enrichment $\delta_0$ raises the tested BVD share above the suspect-pool
# composition early on (we assume the laboratory is biased towards testing
# cases thought more likely to be Ebola), and decays as testing broadens. The
# confirmed-death enrichment $m_{\text{death}}$ scales the death-confirmation
# odds relative to the case composition. The sensitivity prior sits just below
# the field whole-blood clinical sensitivity reported for the GeneXpert Ebola
# assay; the specificity is high but imperfect; the severity enrichment is
# moderate and one-sided (triage upsamples BVD, never down); the
# confirmed-death enrichment is centred on no enrichment:
#
# ```math
# \tau_{\text{test}} \sim \mathrm{Beta}(5,\ 2), \qquad
# s \sim \mathrm{Beta}(30,\ 2), \qquad
# \mathrm{spec} \sim \mathrm{Beta}(60,\ 2),
# ```
#
# ```math
# \delta_0 \sim \mathrm{Normal}^{+}(1.5,\ 0.75), \qquad
# m_{\text{death}} \sim \mathrm{LogNormal}(0,\ 1). \tag{22}
# ```
#
# The non-BVD background rate $\lambda_{\text{bg}}$ enters the suspected-case
# stream and is described with it below. How the per-window confirmation
# probability ties to the suspect-pool composition, rather than being a free
# random effect, is set out with the laboratory submodel below, after the
# composition is defined.

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
#md # <details><summary>Submodel: confirmed_death_enrichment_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.confirmed_death_enrichment_model()),
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
# N_{\text{travel}} \sim \mathrm{Normal}^{+}(1871,\ 200). \tag{23}
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
# delay $f_{\text{rep}}$ and scaled by the DRC ascertainment
# $p_{\text{DRC}}$. We write the BVD onset-to-report series at unit
# ascertainment as $\text{bvdrep} = \mathrm{conv}(\text{onsets},
# f_{\text{rep}})$, where $\mathrm{conv}$ denotes discrete convolution of a
# daily series with a delay PMF. The second is an additive non-BVD
# background accruing at $\lambda_{\text{bg}}$ per day, so a suspected case
# need not be a true BVD infection. We assume this background is a small
# minority of suspected reports, so the prior on $\lambda_{\text{bg}}$ is an
# informative half-normal $\mathrm{Normal}^{+}(0, 1)$ per day. The
# between-vintage increments are scored with a NegBinomial sharing $k$:
#
# ```math
# Y_{\text{cases}} \sim \mathrm{NegBinomial}\!\bigl(p_{\text{DRC}}\,
#     \text{bvdrep} + \lambda_{\text{bg}}\, n,\ k\bigr). \tag{24}
# ```

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

# ##### Suspected deaths
#
# Suspected deaths are the CFR-weighted convolution of the daily onsets with
# the onset-to-death PMF $f_d$, modelled on the incidence scale. The
# between-vintage increments are scored at each sitrep date with a NegBinomial
# sharing $k$. The death history ends at the cut-off, so the cut-off total is
# the final increment and is not scored separately. Writing the daily death
# series as $\mathrm{CFR}\,\mathrm{conv}(\text{onsets}, f_d)$, the increment
# at vintage $i$ is
#
# ```math
# Y_{\text{deaths},i} - Y_{\text{deaths},i-1} \sim
#     \mathrm{NegBinomial}(\mu_i - \mu_{i-1},\ k),
#     \qquad \mu_0 = 0, \tag{25}
# ```
#
# with $\mu_i$ the modelled daily death series cumulated to sitrep day $d_i$.

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
# The laboratory pipeline fits two streams. The received-specimen volume is
# the suspected daily pipeline ($p_{\text{DRC}}\,\text{bvdrep}$ plus the
# non-BVD background $\lambda_{\text{bg}}$) carried through the receipt delay
# $f_{\text{rec}}$ and thinned by the testing fraction $\tau_{\text{test}}$,
# scored per vintage with the shared $k$:
#
# ```math
# Y_{\text{received}} \sim \mathrm{NegBinomial}\!\bigl(\tau_{\text{test}}\,
#     \mathrm{conv}(p_{\text{DRC}}\,\text{bvdrep}
#     + \lambda_{\text{bg}},\, f_{\text{rec}}),\ k\bigr). \tag{26}
# ```
#
# The confirmed positives in each laboratory window $v$ are scored as a
# Binomial of the observed specimens-analysed denominator $A_v$ with a
# per-window tested-positive probability $p_{\text{pos},v}$. We tie that
# probability to the composition of the tested pool rather than leave it as a
# free random effect, so the confirmed data help identify the non-BVD
# background. The suspect-pool composition $\varphi_v$ is the BVD share among
# suspected reports in the window, $\varphi_v = (p_{\text{DRC}}\,\text{bvd})_v
# / ((p_{\text{DRC}}\,\text{bvd})_v + \lambda_{\text{bg},v})$, which falls as
# the background grows. The tested BVD share $q_v$ raises $\varphi_v$ by the
# decaying severity enrichment $\delta_0$. A truly BVD specimen then tests
# positive with the sensitivity $s$, and a non-BVD specimen with the
# false-positive rate $1 - \mathrm{spec}$, so the false-positive term carries
# the non-BVD share and the laboratory data identify the background:
#
# ```math
# p_{\text{pos},v} = s\, q_v + (1 - \mathrm{spec})(1 - q_v),
# \qquad
# C_v \sim \mathrm{Binomial}(A_v,\ p_{\text{pos},v}). \tag{27}
# ```
#
# The early confirmed vintages (18-23 May) have no per-vintage analysed
# denominator, so they are scored as NegBinomial counts against the modelled
# laboratory volume $V_v$ (the modelled received-specimen volume binned to
# the window) with the same composition-linked positivity, extending the use
# of the confirmed data to where no laboratory denominator is observed:
#
# ```math
# C_v^{\text{early}} \sim \mathrm{NegBinomial}(p_{\text{pos},v}\, V_v,\ k).
# \tag{28}
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
# We currently model confirmed deaths as a thinning of the modelled
# suspected deaths, not as a cases-like laboratory process. The confirmation
# probability is the suspected-case BVD composition $q_{\text{susp}} =
# p_{\text{DRC}}\,\text{bvd} / (p_{\text{DRC}}\,\text{bvd} +
# \lambda_{\text{bg}})$ enriched on the odds scale by $m_{\text{death}}$, so
# it stays in $(0, 1)$ without a hard clamp and a confirmed-death observation
# informs the background and ascertainment. The between-vintage increments
# are scored with a NegBinomial sharing $k$:
#
# ```math
# p_{\text{conf-death}} =
#     \operatorname{logistic}(\operatorname{logit}(q_{\text{susp}})
#     + \log m_{\text{death}}). \tag{29}
# ```
#
# A cases-like laboratory process for confirmed deaths (with windows that
# carry no testing data) would be the cleaner treatment; that is flagged in
# issue
# [#225](https://github.com/epiforecasts/BVDOutbreakSize/issues/225). The
# thinning ties confirmed deaths to the modelled suspected-death trajectory,
# so it may only cover deaths up to where suspected deaths exist.

#md # ```@raw html
#md # <details><summary>Submodel: confirmed_deaths_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.confirmed_deaths_model(
#md #         missing, missing, Float64[], Float64[], 1.0, Float64[], 1.0)),
#md #     "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Exported cases
#
# The exports stream is travel-gated, so the at-risk clock runs from
# infection. An infected person travels to Uganda at the daily per-capita
# travel rate and stays at risk of being exported and detected only until the
# infection-to-detection delay has elapsed. The daily at-risk export
# person-time scales the Uganda ascertainment and the travel rate by the
# infections that have not yet completed that delay, and accumulates into the
# cumulative export intensity $\Lambda(t)$. The infection-to-detection delay
# is the onset-to-detection delay convolved with the incubation period, so
# the survival clock runs from infection:
#
# ```math
# \text{detected} = \mathrm{conv}(\text{infections},\,
#     \mathrm{conv}(f_{\text{inc}}, f_{\text{det}})). \tag{30}
# ```
#
# Each observed Uganda import is fitted at its reported detection date. An
# import detected on a given day is scored as a Poisson of the rise in
# cumulative export intensity between consecutive detection dates, with a term
# before the earliest detection observed at zero, since no export is expected
# then. We model zero reports after the last detection date, assuming travel
# across the border beyond this was inconsistent with our assumed travel
# rates:
#
# ```math
# Y_{\text{exports},i} \sim
#     \mathrm{Poisson}\!\bigl(\Lambda(d_i) - \Lambda(d_{i-1})\bigr),
# \qquad \Lambda(d_0) = \Lambda(d_1 - 1). \tag{31}
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
# The expected deaths among detected exports weight the same daily at-risk
# export person-time by the infection-to-death delay (the onset-to-death PMF
# convolved with the incubation period) and scale by the CFR. Each dated
# Uganda export death is scored at its reported date with a per-day Poisson,
# the same dated-event likelihood the exports use:
#
# ```math
# Y_{\text{exp-deaths},i} \sim
#     \mathrm{Poisson}\!\bigl(\Lambda_d(\delta_i)
#     - \Lambda_d(\delta_{i-1})\bigr), \tag{32}
# ```
#
# with $\Lambda_d$ the cumulative export-death intensity. We assume the
# observation follows a Poisson.

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
# joint. The frozen and predictive variants reuse the same structure with the
# data cut to an earlier date or the counts dropped.

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

prior_C_table #hide

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
# chains of 1500 post-warmup draws each after 1000 warmup adaptation steps,
# at a target acceptance probability of 0.9. Chains initialise from the
# prior. We fit the joint model and each single-stream model so the
# per-stream posteriors over the outbreak size can be compared with the
# joint.

#md # ```@raw html
#md # <details><summary>Run the joint and per-stream NUTS fits</summary>
#md # ```

const _BREAKPOINT = obs.n - obs.who_first_sitrep_days

chn_joint = nuts_sample(
    bvd_joint(
    obs.n, obs.exported_cases, obs.total_deaths,
    obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
    obs.tests_analysed;
    confirmed_deaths = obs.confirmed_deaths,
    deaths_history = obs.deaths_history,
    reported_history = obs.reported_history,
    confirmed_history = obs.confirmed_history,
    confirmed_deaths_history = obs.confirmed_deaths_history,
    lab_history = obs.lab_history,
    tests_received_history = obs.tests_received_history,
    export_case_days = obs.export_case_days,
    export_death_days = obs.export_death_days,
    breakpoint = _BREAKPOINT,
    background_re = true,
    confirmed_positivity_link = :composition,
    genetic = genetic_seeding_model,
    tmrca_days = obs.tmrca_days));

chn_exports = nuts_sample(
    exports_only_model(obs.n, obs.exported_cases;
    export_case_days = obs.export_case_days,
    breakpoint = _BREAKPOINT));

chn_deaths = nuts_sample(
    deaths_only_model(obs.n, obs.total_deaths;
    deaths_history = obs.deaths_history,
    breakpoint = _BREAKPOINT));

chn_cases = nuts_sample(
    cases_only_model(obs.n, obs.reported_cases;
    reported_history = obs.reported_history,
    breakpoint = _BREAKPOINT));

chn_confirmed = nuts_sample(
    confirmed_only_model(obs.n, obs.confirmed_cases;
    confirmed_history = obs.confirmed_history,
    lab_history = obs.lab_history,
    tests_received_history = obs.tests_received_history,
    breakpoint = _BREAKPOINT));

chn_confirmed_deaths = nuts_sample(
    confirmed_deaths_only_model(obs.n, obs.confirmed_deaths,
    obs.total_deaths;
    deaths_history = obs.deaths_history,
    confirmed_deaths_history = obs.confirmed_deaths_history,
    breakpoint = _BREAKPOINT));

## This composer keeps the deaths and exports submodels only for their
## CFR, onset-to-death PMF and export onsets, leaving their own counts
## missing, which leaves two redundant sampled discrete draws; the model
## check is disabled so NUTS will run (see `nuts_sample`).
chn_exports_deaths = nuts_sample(
    exports_deaths_only_model(obs.n, obs.exports_deaths;
        breakpoint = _BREAKPOINT); check_model = false);

posterior_C_joint = vec(Array(chn_joint[:C_T]));
posterior_C_exports = vec(Array(chn_exports[:C_T]));
posterior_C_deaths = vec(Array(chn_deaths[:C_T]));
posterior_C_cases = vec(Array(chn_cases[:C_T]));
posterior_C_confirmed = vec(Array(chn_confirmed[:C_T]));
posterior_C_confirmed_deaths = vec(Array(chn_confirmed_deaths[:C_T]));
posterior_C_exports_deaths = vec(Array(chn_exports_deaths[:C_T]));

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
    "exports (cases)" => chn_exports, #hide
    "exports (deaths)" => chn_exports_deaths, #hide
    "deaths (DRC)" => chn_deaths, #hide
    "cases (DRC)" => chn_cases, #hide
    "confirmed (DRC)" => chn_confirmed, #hide
    "confirmed deaths (DRC)" => chn_confirmed_deaths) #hide

#md # ```@raw html
#md # </details>
#md # ```

# #### No-onward-transmission counterfactual
#
# To bound the deaths already committed at the cut-off, we project the
# deaths that would still occur if all transmission stopped on the report
# date. Every infection present by the cut-off still dies with probability
# CFR, so the committed future deaths are $\Delta D = \mathrm{CFR}\cdot C_T -
# \mathbb{E}[D_T]$, the CFR-weighted cumulative infections net of the deaths
# already expected. The figure is shown in the counterfactual results
# below.
#
# #### One-week-ahead forecast
#
# We project each DRC stream seven days beyond the cut-off as a no-change
# forecast: the current growth rate is carried forward with no further
# interventions and no saturation, and the projection carries both parameter
# and observation uncertainty. We forecast the four DRC streams (suspected
# reported cases, suspected deaths, laboratory-confirmed cases and confirmed
# deaths) but not exports, since cross-border travel is unlikely to continue
# at its baseline rate, so the forward travel rate the export model relies on
# no longer holds. The figure is shown in the
# [one-week-ahead forecast results](@ref "One-week-ahead forecast results")
# below.
#
# #### Forecast-versus-frozen evaluation
#
# We assess the forecast against data observed since by freezing the data to
# roughly one week before the current cut-off, re-fitting, and projecting one
# week ahead with the same forecast machinery, then comparing that projection
# against the counts observed by the current cut-off. The frozen re-fit cuts
# the data to an earlier cut-off and re-fits the joint model, so a later
# result can be told apart from a later method. Each frozen re-fit is a
# reduced fit of 500 draws across two chains, illustrative rather than a
# production result. The same frozen-refit approach underlies the
# matched-in-time comparison with McCabe et al., who computed their scenarios
# at fixed situation-report cut-offs: we freeze the renewal data to the same
# cut-off and re-fit, with any remaining gap the method read at the date the
# scenario was computed. It also underlies the estimate-evolution series,
# where the renewal is re-fit frozen at each release date. The helper below
# performs one frozen joint re-fit and is reused by the forecast validation,
# estimate-evolution and matched-in-time results.

#md # ```@raw html
#md # <details><summary>Frozen-fit helper (reused by the forecast validation, evolution and matched-in-time sections)</summary>
#md # ```

## A reduced joint fit (500 draws × 2 chains) to the data frozen at
## `cutoff_date`. The frozen named tuple has the same shape as the full
## `obs`, so the model call mirrors the headline joint fit.
function fit_frozen_joint(cutoff_date; samples = 500, chains = 2)
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
            tests_received_history = o.tests_received_history,
            export_case_days = o.export_case_days,
            export_death_days = o.export_death_days,
            breakpoint = bp,
            background_re = true,
            confirmed_positivity_link = :composition,
            genetic = genetic_seeding_model,
            tmrca_days = o.tmrca_days);
        samples = samples, chains = chains)
    return (; cutoff = o.cutoff, o, chn)
end

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
    ec = posterior_summary(vec(Array(chn_joint[:expected_confirmed_T])))

    Markdown.parse("""
    - **Cumulative infections:** the posterior is $(ints_i(sC))
      infections.
    - Against the $(obs.confirmed_cases) laboratory-confirmed cases by the
      cut-off that is roughly $(f_lo)–$(f_hi)× as many infections, so
      confirmed cases capture only a small share of the estimated outbreak.
    - **Confirmed-case fit:** the model expects $(ints_i(ec)) confirmed
      cases by the cut-off, against $(obs.confirmed_cases) observed.
    - **Outbreak start and age:** the outbreak began on a start date of
      $(ints_d(sT)), an elapsed age to the cut-off of $(ints_i(sT)) days.
    - **Growth rate and doubling time:** the initial growth rate was
      $(ints_f(sr0, 3)) per day, an initial doubling time of
      $(ints_f(sdt0, 1)) days.
      The latest growth rate is $(ints_f(sr, 3)) per day, a latest doubling
      time of $(ints_f(sdt, 1)) days.
    - **Reproduction number:** the initial reproduction number was
      $(ints_f(sR0, 2)) and the latest is $(ints_f(sRT, 2)).
    - **Case-fatality ratio:** the posterior is $(ints_f(scfr, 2)).
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
    [:r, :r0, :doubling_time, :T, :R_T, :CFR, :C_T]; digits = 2);

#md # ```@raw html
#md # </details>
#md # ```

infection_summary #hide

#md # ```@raw html
#md # <details><summary>Infection-parameter pair plot (prior overlaid)</summary>
#md # ```

infection_pair_fig = plot_pair(chn_joint,
    [:R_T, :r, :T, :CFR];
    prior = prior_chn);

#md # ```@raw html
#md # </details>
#md # ```

infection_pair_fig #hide

# ### Reproduction number over time
#
# The daily reproduction number over the period we estimate it for, the
# established outbreak from the genetic bound to the cut-off.
# The median is shown with 50% and 90% ribbons and about a hundred sampled
# trajectories.
# The first situation report on 18 May 2026 marks the start of the response
# scale-up (red dashed) and the end of the three-week scale-up is the red
# dotted line; the data cut-off is grey dashed.

#md # ```@raw html
#md # <details><summary>Reproduction-number trajectory</summary>
#md # ```

rt_fig = plot_rt(chn_joint;
    n = obs.n, breakpoint = _BREAKPOINT,
    rt_start = clamp(
        obs.n - round(Int, obs.tmrca_days) + SEEDING_ANCHOR_LEAD, 1, obs.n),
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

# ### Surveillance parameters
#
# The surveillance-data parameters: the reporting fractions for the DRC and
# Uganda, the surveillance dispersion, and the laboratory pipeline (the
# testing fraction and receipt delay, the per-suspected and per-test
# positivity, the non-BVD background rate, and the death-confirmation
# probability).
# The table reports their credible intervals; the pair plot beside it shows
# their joint posterior with the prior overlaid.

#md # ```@raw html
#md # <details><summary>Surveillance-parameter summary table</summary>
#md # ```

surveillance_summary = summary_table(chn_joint,
    [:p_drc, :p_uganda, :k, :tau_test, :lambda_bg,
        :suspected_positivity, :test_positivity, :expected_confirmed_T,
        :expected_received_T, :m_death, :death_composition,
        :death_confirmation, :expected_confirmed_deaths_T];
    digits = 3);

#md # ```@raw html
#md # </details>
#md # ```

surveillance_summary #hide

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

# ### Export parameters
#
# The export-stream parameters: the daily outbound traveller volume that
# sets the cross-border travel rate, and the implied expected exported cases
# by the cut-off.
# The table reports their credible intervals; the pair plot beside it shows
# their joint posterior with the prior overlaid.

#md # ```@raw html
#md # <details><summary>Export-parameter summary table</summary>
#md # ```

export_summary = summary_table(chn_joint,
    [Symbol("exports_state.travel_state.daily_travellers"), :expected_exports_T];
    digits = 2);

#md # ```@raw html
#md # </details>
#md # ```

export_summary #hide

#md # ```@raw html
#md # <details><summary>Export-parameter pair plot (prior overlaid)</summary>
#md # ```

export_pair_fig = plot_pair(chn_joint,
    [Symbol("exports_state.travel_state.daily_travellers"), :expected_exports_T];
    prior = prior_chn);

#md # ```@raw html
#md # </details>
#md # ```

export_pair_fig #hide

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
# The second is the surveillance data, the five dated DRC streams that are
# real per-vintage observations: suspected cases, confirmed cases,
# suspected deaths, confirmed deaths and specimens received.
# The third is the exports, the cross-border imported cases and deaths
# detected in Uganda.
#
# The surveillance group is checked first.
# Each replicated cumulative trajectory is shown across the
# situation-report dates with the observed series overlaid.

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
        deaths_history = _days_only(obs.deaths_history),
        reported_history = _days_only(obs.reported_history),
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = _days_only(obs.confirmed_deaths_history),
        lab_history = obs.lab_history,
        tests_received_history = _days_only(obs.tests_received_history),
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
## `replicates` shape `plot_vintage_conditional_ppc` anchors on each
## vintage's observed previous cumulative for the one-step-ahead
## predictive. Find it by its prefix among the predict keys.
function _vintage_replicates(pp, prefix)
    key = first(k for k in keys(pp)
    if occursin("$prefix.increments", string(k)))
    return collect(pp[key])
end;

## Grid day-index → INSP situation-report date label.
_vintage_dates(days) = string.(obs.seeding .+ Day.(days .- 1));

reported_panel = (;
    title = "Suspected cases",
    dates = _vintage_dates(obs.reported_history.days),
    replicates = _vintage_replicates(
        pp_joint, "cases_state.reported_increments"),
    observed = obs.reported_history.counts, colour = :steelblue);
deaths_panel = (;
    title = "Suspected deaths",
    dates = _vintage_dates(obs.deaths_history.days),
    replicates = _vintage_replicates(
        pp_joint, "deaths_state.death_increments"),
    observed = obs.deaths_history.counts, colour = :firebrick);
## Specimens received is also a per-vintage time series (the receipt-delay
## and tested-fraction throughput), so it gets the same cumulative
## conditional check as the suspected streams.
tests_received_panel = (;
    title = "Specimens received",
    dates = _vintage_dates(obs.tests_received_history.days),
    replicates = _vintage_replicates(
        pp_joint, "confirmed_state.received_increments"),
    observed = obs.tests_received_history.counts, colour = :seagreen);

## Confirmed cases are scored over two groups of laboratory windows: the
## early confirmed vintages (no per-vintage analysed denominator, scored
## as counts against the modelled laboratory volume) and the observed
## windows (a Binomial of the observed analysed denominator). Both groups
## produce per-window replicate increments in `predict`, so concatenating
## them oldest-first gives the per-vintage cumulative confirmed-case
## trajectory, anchored on the observed cumulative confirmed at each window
## end-day. The 24-25 May analysis stall merges into 26 May, so the window
## grid is slightly coarser than the raw confirmed history.
_conf_windows = BVDOutbreakSize.confirmed_positivity_windows(
    obs.confirmed_history, obs.lab_history);
## Oldest-first: early (no denominator) → observed (analysed Binomial) →
## late (post-28 May, no denominator, modelled volume).
_conf_window_days = vcat(_conf_windows.early_days, _conf_windows.obs_days,
    _conf_windows.late_days);
function _confirmed_at(day)
    i = searchsortedlast(obs.confirmed_history.days, day)
    return i == 0 ? 0 : Int(obs.confirmed_history.counts[i])
end;
_conf_early = _vintage_replicates(
    pp_joint, "confirmed_state.early_increments");
_conf_obs = collect(first(pp_joint[k]
for k in keys(pp_joint)
if occursin("confirmed_state.confirmed_positives.positives", string(k))));
_conf_late = _vintage_replicates(
    pp_joint, "confirmed_state.late_increments");
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
        pp_joint, "confirmed_deaths_state.cdeath_increments"),
    observed = obs.confirmed_deaths_history.counts, colour = :purple);

joint_vintage_ppc_fig = plot_vintage_conditional_ppc(
    [reported_panel, confirmed_panel, deaths_panel,
    confirmed_deaths_panel, tests_received_panel]);

#md # ```@raw html
#md # </details>
#md # ```

joint_vintage_ppc_fig #hide

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
## per-day count vector `<prefix>.counts`; match the full prefixed varname so
## the deterministic `expected_*_T` quantities are not picked up by a loose
## substring, then sum each replicate's per-day vector into the total.
function _dated_total(pp, name)
    key = first(k for k in keys(pp) if occursin("$name.counts", string(k)))
    return [sum(v) for v in vec(Array(pp[key]))]
end;

pp_exports = _dated_total(pp_joint, "exports_state.export_obs");
pp_exports_deaths = _dated_total(
    pp_joint, "exports_deaths_state.death_obs");

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

# ### One-week-ahead forecast results
#
# The cumulative and new expected counts by $T + 7$ for the four DRC streams
# (suspected reported cases, suspected deaths, laboratory-confirmed cases and
# confirmed deaths), from the no-change projection defined in the methods
# [one-week-ahead forecast](@ref "One-week-ahead forecast").

#md # ```@raw html
#md # <details><summary>Generate the one-week-ahead forecast</summary>
#md # ```

forecast = forecast_reported(chn_joint;
    horizon = 7,
    obs_cases = obs.reported_cases,
    obs_deaths = obs.total_deaths,
    obs_confirmed = obs.confirmed_cases,
    obs_confirmed_deaths = obs.confirmed_deaths);
forecast_summary = forecast_table(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_summary #hide

# The coming week at a glance, split into the latent quantities and the
# observations.
# The latent figure shows the new infections, symptom onsets and deaths over
# the horizon, with the reproduction number carried forward.

#md # ```@raw html
#md # <details><summary>One-week-ahead latent forecast plot</summary>
#md # ```

forecast_latent_fig = plot_forecast_latent(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_latent_fig #hide

# The observation figure shows the new reported cases, confirmed cases and
# confirmed deaths over the horizon.

#md # ```@raw html
#md # <details><summary>One-week-ahead observed forecast plot</summary>
#md # ```

forecast_fig = plot_forecast(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_fig #hide

# ### Forecast validation (last week versus now)
#
# How last week's forecast held up against the data since observed, using the
# frozen re-fit and one-week projection defined in the methods
# [forecast-versus-frozen evaluation](@ref
# "Forecast-versus-frozen evaluation").

#md # ```@raw html
#md # <details><summary>Fit one week back and validate the one-week-ahead forecast</summary>
#md # ```

validation_cutoff = string(obs.cutoff - Day(7))
frozen_lastweek = fit_frozen_joint(validation_cutoff)

validation_forecast = forecast_reported(frozen_lastweek.chn;
    horizon = 7,
    obs_cases = frozen_lastweek.o.reported_cases,
    obs_deaths = frozen_lastweek.o.total_deaths,
    obs_confirmed = frozen_lastweek.o.confirmed_cases,
    obs_confirmed_deaths = frozen_lastweek.o.confirmed_deaths);

validation_table = forecast_vs_truth(validation_forecast;
    cases = obs.reported_cases, deaths = obs.total_deaths,
    confirmed = obs.confirmed_cases,
    confirmed_deaths = obs.confirmed_deaths);

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
    cases = obs.reported_cases, deaths = obs.total_deaths,
    confirmed = obs.confirmed_cases,
    baseline_cases = frozen_lastweek.o.reported_cases,
    baseline_deaths = frozen_lastweek.o.total_deaths,
    baseline_confirmed = frozen_lastweek.o.confirmed_cases);

#md # ```@raw html
#md # </details>
#md # ```

validation_fig #hide

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
    "exports (cases)" => posterior_C_exports,
    "deaths (DRC)" => posterior_C_deaths,
    "cases (DRC)" => posterior_C_cases,
    "confirmed (DRC)" => posterior_C_confirmed,
    "joint" => posterior_C_joint);

#md # ```@raw html
#md # </details>
#md # ```

streams_C_table #hide

# Overlaid posterior densities of the infection count from each fit.
# The horizontal axis is set to the upper tail of the widest fit so the
# bulk of every stream's posterior mass is visible.

#md # ```@raw html
#md # <details><summary>Overlaid infection-count density plot</summary>
#md # ```

## Set the x-axis to the 95th percentile across the streams so the bulk of
## each posterior is visible rather than clipping out the upper mass.
density_xmax = 1.05 * maximum(
    quantile(s, 0.95)
for s in (posterior_C_exports, posterior_C_deaths,
    posterior_C_cases, posterior_C_confirmed, posterior_C_joint))

cumulative_density_fig = plot_cumulative_cases(
    "exports (cases)" => posterior_C_exports,
    "deaths (DRC)" => posterior_C_deaths,
    "cases (DRC)" => posterior_C_cases,
    "confirmed (DRC)" => posterior_C_confirmed,
    "joint" => posterior_C_joint;
    scenarios = [], xmax = density_xmax);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_density_fig #hide

# The frozen re-fits below freeze the renewal data to an earlier cut-off
# and re-fit, so a later result can be told apart from a later method.
# Each is a reduced fit of 500 draws across two chains, illustrative rather
# than a production result, reusing the frozen-fit helper defined above.

#md # ```@raw html
#md # <details><summary>Freeze the renewal data to a cut-off and re-fit (reduced)</summary>
#md # ```

frozen_20may = fit_frozen_joint("2026-05-20")
frozen_23may = fit_frozen_joint("2026-05-23")
frozen_26may = fit_frozen_joint("2026-05-26")
frozen_28may = fit_frozen_joint("2026-05-28")

frozen_C_20may = vec(Array(frozen_20may.chn[:C_T]))
frozen_C_23may = vec(Array(frozen_23may.chn[:C_T]))

frozen_reports_20may = vec(Array(frozen_20may.chn[:expected_reports_T]))
frozen_reports_23may = vec(Array(frozen_23may.chn[:expected_reports_T]))

#md # ```@raw html
#md # </details>
#md # ```

# ### Estimate evolution across releases
#
# How the outbreak-size estimate moved as the situation reports accrued.
#
# The project has published a tagged results release at each data cut-off
# (<https://github.com/epiforecasts/BVDOutbreakSize/releases>), bundling
# the posterior draws and input data.
# The releases to date are from the earlier closed-form integral model, so
# the released series is the project's published estimate over time.
# The current-model series is this renewal model re-fit frozen at each date
# up to the present cut-off, so it is the current method evaluated at each
# past date and rises as the outbreak grows rather than sitting flat.
# Both climb as the suspected-case count grows.
# The released values are read off the archived posterior draws so the
# figure builds without refetching the releases.
# Each estimate is a median with an alpha-shaded 90% credible band.

#md # ```@raw html
#md # <details><summary>Released estimates and frozen renewal re-fits</summary>
#md # ```

## Released median C_T and 90% interval per data cut-off, read from each
## release's archived `posterior_draws.csv` (one canonical build per
## cut-off date: builds 241, 350, 470 and 586). These are the integral
## model's published estimates; see the release page for provenance.
release_evolution = [
    ("2026-05-18", 925, 419, 2075),
    ("2026-05-23", 1364, 680, 3137),
    ("2026-05-26", 3041, 1961, 5172),
    ("2026-05-28", 3510, 2196, 6325)
]

## The current model evaluated at each past date, the renewal re-fit frozen
## at each cut-off plus the live current-data fit, as a rising series.
function _ci90(xs)
    (round(Int, quantile(xs, 0.5)),
        round(Int, quantile(xs, 0.05)), round(Int, quantile(xs, 0.95)))
end
renewal_frozen = [
    (string(frozen_20may.cutoff), _ci90(frozen_C_20may)...),
    (string(frozen_23may.cutoff), _ci90(frozen_C_23may)...),
    (string(frozen_26may.cutoff),
        _ci90(vec(Array(frozen_26may.chn[:C_T])))...),
    (string(frozen_28may.cutoff),
        _ci90(vec(Array(frozen_28may.chn[:C_T])))...),
    (string(obs.cutoff), _ci90(posterior_C_joint)...)
]

evolution_fig = plot_estimate_evolution(release_evolution;
    renewal = renewal_frozen,
    title = "Outbreak-size estimate as data accrued");

#md # ```@raw html
#md # </details>
#md # ```

evolution_fig #hide

# ### Comparison with McCabe et al.
#
# Earlier versions of this work reimplemented McCabe et al. closely as an
# exponential-growth model.
# This version is a discrete-time renewal model with a time-varying
# reproduction number and every data stream fitted jointly.
# McCabe et al. computed their scenarios at fixed situation-report cut-offs,
# so we match in time, freezing our data to the same cut-off and re-fitting.
# The McCabe scenarios shown are the 20 May update [mccabe2026update](@cite).
# They are deterministic point estimates carrying no uncertainty, so they
# appear as bare points rather than intervals.
# Alongside them we show our own estimates at the matched cut-offs, the
# renewal fit frozen at 20 May and at 23 May with their 95% intervals, so
# the comparison carries our trajectory and not only McCabe's.
# The 20 May freeze is the earliest matched cut-off with a coherent
# suspected-case series, since the earliest situation report is 18 May.
#
# The expected reported-case count moves sharply with the data.
# At the 20 May cut-off it sits close to McCabe et al.'s scenarios and well
# below the current-data fit.

#md # ```@raw html
#md # <details><summary>Matched-in-time reported-case comparison</summary>
#md # ```

function _ci95(xs)
    (round(Int, quantile(xs, 0.5)),
        round(Int, quantile(xs, 0.025)),
        round(Int, quantile(xs, 0.975)))
end

matched_rows = vcat(
    [(label, val, val, val) for (label, val) in REPORT_SCENARIOS],
    [("Renewal frozen 20 May", _ci95(frozen_reports_20may)...),
        ("Renewal frozen 23 May", _ci95(frozen_reports_23may)...)])

## Colour the rows by source so McCabe's deterministic scenarios and our own
## renewal estimates read apart.
matched_groups = vcat(
    fill("McCabe et al. (20 May update)", length(REPORT_SCENARIOS)),
    ["Our renewal estimate", "Our renewal estimate"])

matched_comparison_fig = plot_estimate_comparison(matched_rows;
    xlabel = "Cumulative reported cases",
    groups = matched_groups,
    group_colours = ["McCabe et al. (20 May update)" => :grey,
        "Our renewal estimate" => :firebrick]);

#md # ```@raw html
#md # </details>
#md # ```

matched_comparison_fig #hide

# The figure is on the reported-case scale, the ascertained quantity McCabe
# et al. report.
# The project's released outbreak-size estimates are on the infection scale
# and are shown in the
# [estimate evolution](@ref "Estimate evolution across releases") figure
# above.
#
# Side-by-side outbreak-size intervals for the two frozen fits and the
# current-data fit, so the shift with the data cut-off reads off directly.

#md # ```@raw html
#md # <details><summary>Frozen-fit C_T table</summary>
#md # ```

frozen_streams_table = streams_table(
    "frozen 20 May" => frozen_C_20may,
    "frozen 23 May" => frozen_C_23may,
    "current data" => posterior_C_joint);

#md # ```@raw html
#md # </details>
#md # ```

frozen_streams_table #hide

# ### Delay sensitivity
#
# The second method dates the outbreak from how far deaths lag symptom
# onset, so the assumed onset-to-death delay sets the implied infection
# count.
# We probe it by re-fitting the joint model under a shorter and a longer
# onset-to-death delay centre either side of the baseline, holding
# everything else fixed.
# Each variant is a reduced fit of 500 draws across two chains,
# illustrative rather than a production result.
#
# The infection count to date shifts with the assumed delay, and the
# table and overlaid densities below show how far.

#md # ```@raw html
#md # <details><summary>Re-fit the joint under shorter and longer onset-to-death delays (reduced)</summary>
#md # ```

## One reduced joint re-fit on the live data, with hooks to override the
## genetic-seeding bound and the deaths submodel. The deaths submodel is
## passed the same way the genetic-seeding override is, as a closure
## matching the joint's `deaths(history, total, onsets, k; background_re)`
## call, so an alternative onset-to-death delay can be injected without
## touching the package.
function refit_joint_variant(;
        deaths = deaths_model,
        tmrca_days = obs.tmrca_days,
        tmrca_days_sd = 15.0,
        samples = 500, chains = 2)
    chn = nuts_sample(
        bvd_joint(
            obs.n, obs.exported_cases, obs.total_deaths,
            obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
            obs.tests_analysed;
            confirmed_deaths = obs.confirmed_deaths,
            deaths_history = obs.deaths_history,
            reported_history = obs.reported_history,
            confirmed_history = obs.confirmed_history,
            confirmed_deaths_history = obs.confirmed_deaths_history,
            lab_history = obs.lab_history,
            tests_received_history = obs.tests_received_history,
            export_case_days = obs.export_case_days,
            export_death_days = obs.export_death_days,
            breakpoint = _BREAKPOINT,
            background_re = true,
            confirmed_positivity_link = :composition,
            deaths = deaths,
            genetic = genetic_seeding_model,
            tmrca_days = tmrca_days,
            tmrca_days_sd = tmrca_days_sd);
        samples = samples, chains = chains)
    return chn
end

## Shorter and longer onset-to-death delay centres, each a closure that
## re-injects the onset-to-death delay prior into the deaths submodel while
## keeping its other defaults.
deaths_short_delay = (history,
    total,
    onsets,
    k;
    kwargs...) -> deaths_model(history, total, onsets, k;
    onset_to_death = censored_delay_model(40;
        mean_prior = truncated(Normal(7.0, 2.0); lower = 1),
        sd_prior = truncated(Normal(5.4, 1.5); lower = 1)),
    kwargs...)
deaths_long_delay = (history,
    total,
    onsets,
    k;
    kwargs...) -> deaths_model(history, total, onsets, k;
    onset_to_death = censored_delay_model(40;
        mean_prior = truncated(Normal(16.0, 2.0); lower = 1),
        sd_prior = truncated(Normal(5.4, 1.5); lower = 1)),
    kwargs...)

chn_joint_short_delay = refit_joint_variant(deaths = deaths_short_delay)
chn_joint_long_delay = refit_joint_variant(deaths = deaths_long_delay)

posterior_C_short_delay = vec(Array(chn_joint_short_delay[:C_T]))
posterior_C_long_delay = vec(Array(chn_joint_long_delay[:C_T]))

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Delay-sensitivity infection-count table</summary>
#md # ```

delay_sensitivity_table = streams_table(
    "shorter delay" => posterior_C_short_delay,
    "baseline delay" => posterior_C_joint,
    "longer delay" => posterior_C_long_delay);

#md # ```@raw html
#md # </details>
#md # ```

delay_sensitivity_table #hide

#md # ```@raw html
#md # <details><summary>Delay-sensitivity infection-count density plot</summary>
#md # ```

delay_sensitivity_fig = plot_cumulative_cases(
    "shorter delay" => posterior_C_short_delay,
    "baseline delay" => posterior_C_joint,
    "longer delay" => posterior_C_long_delay;
    scenarios = []);

#md # ```@raw html
#md # </details>
#md # ```

delay_sensitivity_fig #hide

# ### Clock-rate sensitivity
#
# The whole outbreak-age estimate rests on the genetic bound, the oldest
# date the common ancestor of the sequenced cases can sit, which is set by
# the assumed molecular clock rate.
# The baseline uses the slower clock rate matching this analysis; the
# sequencing source also reports a faster early-epidemic rate that dates
# the common ancestor about two and a half weeks more recently, without
# favouring either [virological2026](@cite).
# We re-fit the joint model under the faster clock and compare the
# infection count to date and the outbreak age.
# The re-fit is reduced to 500 draws across two chains.

#md # ```@raw html
#md # <details><summary>Re-fit the joint under the faster clock rate (reduced)</summary>
#md # ```

## The faster early-epidemic clock dates the common ancestor about 17 days
## more recently, so the bound on the outbreak age sits that many days
## closer to the cut-off, with a tighter spread from its narrower interval.
clock_alt_offset = value(Date("2026-04-11") - Date("2026-03-25"))
tmrca_days_alt = obs.tmrca_days - clock_alt_offset

chn_joint_fast_clock = refit_joint_variant(
    tmrca_days = tmrca_days_alt, tmrca_days_sd = 9.0)

posterior_C_fast_clock = vec(Array(chn_joint_fast_clock[:C_T]))
T_baseline_clock = vec(Array(chn_joint[:T]))
T_fast_clock = vec(Array(chn_joint_fast_clock[:T]))

#md # ```@raw html
#md # </details>
#md # ```

# The infection count to date under the two clock rates, side by side.

#md # ```@raw html
#md # <details><summary>Clock-rate infection-count table</summary>
#md # ```

clock_sensitivity_C_table = streams_table(
    "baseline clock" => posterior_C_joint,
    "faster clock" => posterior_C_fast_clock);

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_C_table #hide

#md # ```@raw html
#md # <details><summary>Clock-rate infection-count density plot</summary>
#md # ```

clock_sensitivity_C_fig = plot_cumulative_cases(
    "baseline clock" => posterior_C_joint,
    "faster clock" => posterior_C_fast_clock;
    scenarios = []);

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

clock_sensitivity_T_table = streams_table(
    "baseline clock" => T_baseline_clock,
    "faster clock" => T_fast_clock;
    digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_T_table #hide

#md # ```@raw html
#md # <details><summary>Clock-rate outbreak-age density plot</summary>
#md # ```

clock_sensitivity_T_fig = plot_density_overlay(
    "baseline clock" => T_baseline_clock,
    "faster clock" => T_fast_clock;
    xlabel = "Outbreak age (days before cut-off)",
    title = "Posterior outbreak age by clock rate");

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
# posterior draws, and a copy of the input `observations.toml` so the
# exact data that produced each result is recorded alongside it.

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
posterior_draws = DataFrame(
    r = vec(Array(chn_joint[:r])),
    r0 = vec(Array(chn_joint[:r0])),
    doubling_time = vec(Array(chn_joint[:doubling_time])),
    T = vec(Array(chn_joint[:T])),
    R_T = vec(Array(chn_joint[:R_T])),
    CFR = vec(Array(chn_joint[:CFR])),
    p_drc = vec(Array(chn_joint[:p_drc])),
    p_uganda = vec(Array(chn_joint[:p_uganda])),
    C_T = vec(Array(chn_joint[:C_T]))
)[1:10:end, :]
CSV.write(joinpath(output_dir, "posterior_draws.csv"), posterior_draws);

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
