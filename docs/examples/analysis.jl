#md # ```@eval
#md # using BVDOutbreakSize, Markdown
#md # readme = read(joinpath(pkgdir(BVDOutbreakSize), "README.md"), String)
#md # body = strip(match(r"^(.*?)<!-- SHARED:END -->"s, readme).captures[1])
#md # body = replace(body,
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
# **See:**
# [current outbreak size](@ref "Summary") ·
# [comparison with McCabe et al.](@ref "Comparison with McCabe et al.") ·
# [how the data streams compare](@ref "How the data streams compare") ·
# [limitations](@ref "Limitations").
#
# ## What we do differently from McCabe et al.
#
# We reimplement the McCabe et al. [mccabe2026](@cite) report as a single
# Bayesian model fitted jointly to every data stream. The points below
# summarise how it differs from the report; the Methods section carries
# the full treatment and the limitations are listed separately.
#
#md # ```@raw html
#md # <details><summary>Expand: differences from the report</summary>
#md # ```
#md #
# **Latent process and parameters**
#
# - *Infections as the latent quantity.* We model the cumulative infection
#   count and map it through an incubation period to onsets and reported
#   cases, recovering the case count McCabe et al. estimate for comparison.
# - *Joint posterior, not 15 scenario estimates.* The growth rate, case
#   fatality ratio (CFR), delays, traveller volume and surveillance
#   dispersion are given priors and sampled jointly, where McCabe et al.
#   fix each and sweep scenarios.
#
# **Delays and convolutions**
#
# - *Explicit infection→detection delay for exports.* Exports follow a
#   delay convolution of the infection trajectory rather than a fixed
#   detection window, capturing pre-symptomatic travel. The window form is
#   kept for the comparison.
# - *Exact deaths convolution.* We evaluate the onset-to-death convolution
#   in exact closed form, where McCabe et al. use a large-time
#   approximation.
#
# **Likelihoods and data streams**
#
# - *Over-dispersed likelihoods.* Deaths and reported cases use a
#   negative-binomial likelihood with a shared surveillance dispersion,
#   where McCabe et al. use Poisson for the deaths and do not model the
#   reported cases.
# - *Ascertainment extension.* A reporting-fraction hyperprior gives a
#   joint posterior over the reported suspected-case count alongside the
#   deaths and exports.
# - *Suspected, confirmed and samples-received streams.* Three streams
#   absent from McCabe et al., tying the reported counts and the
#   laboratory pipeline to the same latent outbreak.
# - *Per-vintage fit of the DRC streams.* The suspected and confirmed
#   series are fitted across successive situation reports rather than as
#   single totals, so the growth rate is informed by the reported
#   trajectory.
#
# **Extensions**
#
# - *No-onward-transmission counterfactual.* A lower bound on the eventual
#   death toll if all onward transmission stopped at the cut-off.
# - *Posterior-predictive forecasts.* A one-week-ahead projection of each
#   stream, plus a retrospective refit of the original report's data to
#   the current cut-off.
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
# - *Fitted only to aggregate reported counts.* The data are a handful of
#   summary figures from press and situation reports: total suspected
#   cases and deaths in the DRC, and cases, with one death, detected in
#   Uganda. There is no line list and no temporal information, so no onset
#   dates, epidemic curve or per-case data. The model also has no
#   knowledge of conditions on the ground such as case definitions,
#   testing capacity, affected areas, interventions or reporting
#   completeness.
#   Every estimate is a model-based extrapolation from sparse summary
#   statistics under strong assumptions, not a measurement.
# - *Prior-driven inference where data is scarce.* A dozen suspected
#   exports, ~$10^2$ deaths, and a single reported-case total give little
#   information about $\tau$, $m$, the surveillance dispersion or the
#   reporting fraction individually, so the posteriors track their priors
#   closely.
# - *Per-sitrep increments are not clean new incidence.* Fitting the new
#   cases and deaths in each sitrep treats the between-vintage change as
#   fresh counts, but later sitreps likely backfill earlier cases and add
#   newly-reporting health zones, and ascertainment likely rose as the
#   response scaled up. The increments therefore mix true incidence with
#   backfill and changing detection, which the ascertainment fraction does
#   not absorb. The most recent increment is also not corrected for
#   right-truncation, so it is exposed as is, with the same caveat as the
#   latest cumulative total.
# - *The reporting format has changed.* From SitRep 013 on 27 May, INSP
#   began revising the suspected-case line list downward as suspects were
#   confirmed or ruled out, and from SitRep 014 on 28 May it stopped
#   publishing a suspected-death headline. The cut-off is 28 May, SitRep
#   014; the confirmed and laboratory streams run to it, but the
#   suspected case and death streams are frozen at their last clean value
#   on 26 May, the SitRep 012 revised re-issue. The same reclassification may
#   already be acting within the fitted window, so the reported counts may
#   carry more uncertainty than the cut-off implies.
# - *Confirmed-cases stream rests on weak priors.* The non-BVD background
#   rate is identified from the suspected/confirmed contrast rather than
#   external surveillance data, the severe-first share curve has no
#   per-sample data for this outbreak, and PCR sensitivity is taken from
#   earlier validation studies.
# - *Data conflict not explored in detail.* We combine four data streams
#   jointly but have not systematically checked whether they conflict, for
#   instance whether the exports and the deaths streams imply different
#   outbreak sizes. Characterising data-source properties and conflict is
#   part of the modelling workflow we otherwise follow
#   [abbott_workflow](@cite); a fuller treatment is left for future work.
#
# **Model assumptions and design**
#
# - *Inherits McCabe et al.'s epidemiological assumptions and core
#   structures.* Exponential growth from a single zoonotic seed; a case
#   trajectory treated as a deterministic function of the latent state,
#   where only the observation counts carry sampling noise via Poisson or
#   NegBinomial, rather than a stochastic incidence process; the
#   cumulative-case and deaths convolution structure for Method 2; the
#   geographic-spread and detection-window structure for Method 1; and no
#   spatial structure beyond the Ituri and Nord Kivu split.
# - *Constant growth rate assumed to hold throughout.* A single
#   exponential rate $r$ governs the whole trajectory, including the
#   per-vintage fit and the one-week-ahead forecast. The first WHO and
#   Africa CDC reports almost certainly triggered a response such as case
#   finding, isolation and safe burials, which would bend the real curve
#   away from constant growth, but the model has no compartment for depletion
#   of susceptibles or intervention effects, so an early-window growth
#   rate is projected forward unchanged. Estimates should be read as "if
#   early growth continues", not as a prediction net of the response.
# - *Onset-to-death delay based on Isiro 2012.* A single-outbreak fit;
#   the delay reporting follows [charniga2024](@cite) but cross-outbreak
#   heterogeneity is unmodelled. The baseline uses the full Rosello
#   onset-to-death distribution, as in McCabe et al. The
#   [delay sensitivity](#Delay-sensitivity) section refits the joint model
#   with the community-only delay, covering the $n = 5$ cases who died
#   without admission as weak evidence of a shorter delay, to show how much
#   the outbreak-size estimate leans on this assumption.
# - *Genetic seeding bound depends on a fixed clock rate.* The time to the
#   most recent common ancestor (TMRCA) is dated under an external Ebola
#   clock rate, and the sampled tree is almost entirely from Bunia. The
#   [clock-rate sensitivity](#Clock-rate-sensitivity) section refits under
#   the faster early-epidemic rate to show how much the timing, growth and
#   outbreak-size estimates move.
# - *Export delay proxied by the DRC reporting delay.* The default model
#   replaces McCabe et al.'s lumped detection window with an
#   infection→detection delay convolution, with the at-risk clock running
#   from infection to capture pre-symptomatic travel. Because the imports
#   are dated by their Uganda admittance dates, with import #1 admitted
#   11 May, the relevant step abroad is onset-to-admittance, a
#   care-seeking delay rather than an administrative reporting lag. We have
#   no Uganda onset-to-admittance data and the few exports cannot identify
#   their own delay, so we **assume it equals the DRC onset-to-report
#   delay** and reinterpret that draw accordingly. The imports were
#   travellers actively seeking care, and import #1 went from admission to
#   death in three days, so their true onset-to-admittance delay is
#   plausibly *shorter*; a longer borrowed delay pushes implied infection
#   times earlier and biases the export-implied outbreak size upward.
# - *Exports fit at their Uganda detection dates, stopped at the last
#   import.* The three imports enter as a dated daily series at their
#   detection dates, with import #1 admitted 11 May, #2 confirmed 16 May
#   and #3 announced 23 May, fit through the infection→detection delay as a
#   per-day Poisson process, with the pre-detection survival term carrying
#   the earliest-detection timing bound that McCabe et al. add separately.
#   Both travel-gated streams, exports and deaths-among-exports, run only
#   to the most recent import on 23 May, not the 28 May cut-off, because
#   they assume a constant per-capita travel rate while cross-border
#   movement likely shifts over the outbreak and the most recent days are
#   right-truncated, so the trailing days carry no informative zero.
#   Import #3's 23 May date is a public announcement, so
#   it carries reporting-pipeline lag.
# - *Exports treated as DRC importations only.* The exports likelihood
#   conditions on the three WHO-confirmed travel-related cases in Uganda
#   and excludes the two domestic contacts, a driver and a healthcare
#   worker linked to the first import, so it no longer conflates onward
#   transmission with importation. It still relies on the source
#   classification of each case; with only three counts, reclassifying any
#   one shifts the implied outbreak size and ascertainment.
# - *Exports assume one-way travel.* An infected traveller is counted as a
#   Ugandan importation once they cross the border, with no return leg, so
#   anyone who travels to Uganda and back before detection is counted in
#   error. The travel input is a one-directional points-of-entry flow
#   from McCabe et al. Table 3, while IOM DTM flow monitoring on the same
#   border shows roughly balanced movement dominated by routine economic
#   round trips, so the at-risk crossings are over-counted. Only the
#   product of Uganda ascertainment and travel rate is identified, so a
#   roughly constant round-trip fraction is absorbed into the fitted
#   ascertainment; what it cannot absorb is a systematic shift in movement,
#   which is why both travel-gated streams stop at the last import.
# - *Ascertainment partially pooled, not separately identified.* Uganda's
#   exported-case ascertainment $p_{\text{Uganda}}$ and DRC's
#   reported-case ascertainment $p_{\text{DRC}}$ share a logit-scale
#   hyperprior. With a handful of suspected exports and one export death,
#   the Uganda fraction is weakly identified and leans on the pooled mean
#   and the DRC side.
# - *Observation streams share one case pool.* The four streams are
#   modelled as conditionally independent given the latent cumulative
#   incidence, but they observe overlapping individuals: exported cases
#   are a subset of all cases and may also be DRC-reported, and expected
#   DRC deaths are computed over all incidence including those who
#   travelled. Ignoring this overlap can double-count evidence and
#   understate uncertainty. The effect is small here because the Uganda
#   counts are small.
# - *Deaths-among-exports is an approximate construction.* The expected
#   count weights the exported at-risk person-time by the death CDF
#   $F_d(T-s)$ rather than convolving the death delay against the
#   exported-case incidence, treating the cohort present at time $s$ as if
#   infected at $s$. A more direct construction would convolve the death
#   delay against the exported-case incidence. The death timing is keyed
#   to infection, so deaths are not timed one incubation period too early,
#   at the cost of a slight double-count of the incubation period shared
#   between the detection and death delays.
# - *Selection bias in deaths-among-exports.* The likelihood assumes
#   Uganda's surveillance retains detected exports through to any
#   subsequent death. If the system loses cases that progress to death,
#   the observed count is biased downward and the constraint it places on
#   $T$ is overstated.
# - *Constant forwarded fraction.* A single $\tau_{\text{forward}}$
#   under-fits the late jump in samples received; a time-varying forwarded
#   fraction is a natural extension.
#
# **Implementation**
#
# - *LLM-driven reimplementation.* The model code, priors, convolution
#   implementation and analysis were drafted by a language model from the
#   published McCabe et al. [mccabe2026](@cite) report and the companion
#   delay reanalysis, then reviewed and revised. Not independently
#   replicated against the authors' code.
# - *First-vintage positivity peak.* The model slightly under-shoots the
#   first vintage's positivity peak.
#md #
#md # ```@raw html
#md # </details>
#md # ```
#
#md # ```@raw html
#md # <details><summary>Load packages and seed the RNG</summary>
#md # ```

using Turing
using Turing: to_submodel, @varname
using Distributions
using StatsFuns: logit, logistic
using DataFrames: DataFrame, rename
import CSV
using Random
using Markdown
using Dates: Date, Day, value
using BVDOutbreakSize
using BVDOutbreakSize: integrate_cumulative,
                       integrate_exports_deaths, delay_convolution
import CairoMakie

## Render figures at higher resolution so they stay crisp in the docs.
CairoMakie.activate!(type = "png", px_per_unit = 3)

Random.seed!(20260518)

## Opt-in TensorBoard tracing of the model fits: set the `BVD_TRACE_DIR`
## environment variable to a directory and every fit below streams its
## draws (warmup included) to a per-fit run under it. Off by default, so
## normal builds are unchanged. Inspect with `tensorboard --logdir <dir>`.
## Loading TensorBoardLogger activates the package's tracing extension.
haskey(ENV, "BVD_TRACE_DIR") && @eval using TensorBoardLogger
function trace_kw(label)
    haskey(ENV, "BVD_TRACE_DIR") ?
    (; callback = tensorboard_callback(
            joinpath(ENV["BVD_TRACE_DIR"], label)),
        warmup = true) : NamedTuple()
end

#md # ```@raw html
#md # </details>
#md # ```

# ## Methods
#
# ### Data
#
# The analysis uses a handful of aggregate counts. The DRC suspected
# cases, suspected deaths and laboratory-confirmed cases are the
# national cumulative totals from the INSP situation reports
# [insp_sitrep_2026](@cite), read from the report PDFs (archived by
# INRB-UMIE [inrb_umie_2026](@cite)). We draw straight from the sitreps
# rather than the published per-zone CSVs because the regional (health
# zone) breakdown is inconsistent with the national totals: the zone
# sums omit cases not yet attributed to a zone, understating the count.
# The figures were read from the PDFs with an LLM agent and
# independently re-scanned by a second agent; the values, with their
# per-vintage sources and caveats, are recorded in
# `data/insp_sitrep_scanned.csv` and `data/observations.toml`. The
# 18-22 May points predate that reporting format and use WHO AFRO Weekly
# External Situation Report 01 [who_afro_sitrep01_2026](@cite). The
# Uganda export-case counts and deaths, the first-export detection date
# and the dated export death come from WHO Disease Outbreak News DON602
# [who_don_2026_602](@cite); the cross-border traveller volume and
# source population from the McCabe et al. report [mccabe2026](@cite).
# The first table lists each figure as of the cut-off (the suspected
# counts are unconfirmed); the source population is fixed and the
# traveller volume is given a normal prior around the McCabe et al.
# figure. The three DRC streams are additionally resolved by sitrep
# vintage and fitted as between-vintage increments, shown in the second
# table.

#md # ```@raw html
#md # <details><summary>Loading observations and building the data table</summary>
#md # ```

obs = load_observations()
observations_table = DataFrame(
    field = [
        "exported_cases",
        "exports_deaths",
        "total_deaths",
        "reported_cases",
        "daily_outbound_travellers",
        "daily_outbound_travellers_sd",
        "source_population"
    ],
    value = [
        obs.exported_cases,
        obs.exports_deaths,
        obs.total_deaths,
        obs.reported_cases,
        obs.daily_outbound_travellers,
        obs.daily_outbound_travellers_sd,
        obs.source_population
    ],
    source = [
        obs.sources.exported_cases,
        obs.sources.exports_deaths,
        obs.sources.total_deaths,
        obs.sources.reported_cases,
        obs.sources.daily_outbound_travellers,
        obs.sources.daily_outbound_travellers_sd,
        obs.sources.source_population
    ]
);

const ITURI_POPULATION = obs.source_population
const ITURI_DAILY_TRAVEL = obs.daily_outbound_travellers
const EXPORTED_CASES = obs.exported_cases
const EXPORTS_DEATHS = obs.exports_deaths

#md # ```@raw html
#md # </details>
#md # ```

observations_table #hide

# The per-vintage cumulative history of the DRC sitrep streams, the
# national totals at each INSP situation-report date. The joint model
# fits the between-vintage increments of these series (a single vintage
# reduces to the cut-off total). The suspected streams are frozen at
# 26 May and the confirmed streams run to the 28 May cut-off, so the
# table is aligned on the confirmed-case dates and the frozen suspected
# columns are blank (`missing`) for 27-28 May. The 23-26 May points are
# the report totals; 18-22 May use the WHO AFRO / early-report baseline.
# See `data/observations.toml` and `data/insp_sitrep_scanned.csv` for
# the per-stream sources.

#md # ```@raw html
#md # <details><summary>Building the per-vintage time-series table</summary>
#md # ```

## Align every stream on the confirmed-case dates (the longest history,
## running to the cut-off) and pad with `missing` where a frozen
## suspected stream has no vintage for that date.
_align(h, dates) = map(
    d -> begin
        i = findfirst(==(d), h.dates)
        i === nothing ? missing : h.values[i]
    end, dates)
vintage_dates = obs.confirmed_case_history.dates
vintage_table = DataFrame(
    sitrep_date = vintage_dates,
    suspected_cases = _align(obs.reported_case_history, vintage_dates),
    confirmed_cases = obs.confirmed_case_history.values,
    suspected_deaths = _align(obs.death_history, vintage_dates),
    confirmed_deaths = _align(obs.confirmed_death_history, vintage_dates)
);

#md # ```@raw html
#md # </details>
#md # ```

vintage_table #hide

# ### Model
#
# #### Model overview
#
# We model a single outbreak seeded by one zoonotic introduction that
# then grows exponentially, so the cumulative number of people ever
# infected by outbreak age $s$ is $C(s) = \exp(r s)$, set by a growth
# rate $r$ (equivalently a doubling time). We never observe infections
# directly. Instead each available data stream observes a thinned,
# delayed or transformed view of that same latent incidence curve.
# Reported cases in the DRC are an ascertained fraction of the
# cumulative cases. Suspected deaths in the DRC are the case fatality
# ratio applied to past incidence, convolved with the onset-to-death
# delay. Cases exported to Uganda are the fraction of recent cases that
# crossed the border, set by the travel rate and a detection window.
# Deaths among the exported cases are those cases weighted by their
# probability of having died by now.
#
# Fitting all four streams together gives the posterior for the latent
# cumulative case count $C(T)$ at the report date $T$ — the quantity we
# care about — while sharing the growth, delay, fatality and
# ascertainment parameters across the streams that depend on them.
#
# In implementation terms, the model is assembled from small reusable
# Turing [ge2018turing](@cite) submodels of three kinds, defined in the
# order they appear below: building-block submodels, observation submodels
# and composers.
#
# The table below shows which building-block parameters feed each
# observation submodel. The *confirmed & received* column covers the
# laboratory pipeline, which observes the confirmed (PCR-positive) cases
# and the samples received from the suspect pool:
#
# | Parameter | Exports | Deaths | Cases | Confirmed & received | Export deaths (time-resolved) | First export-detection timing | Genetic seeding |
# |---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
# | Growth $C(s) = e^{rs}$ | ● | ● | ● | ● | ● | ● | ● |
# | Incubation period | ● | ● | ● | ● | ● | ● |  |
# | Onset-to-death delay |  | ● |  |  | ● |  |  |
# | Case-fatality ratio |  | ● |  |  | ● |  |  |
# | Onset-to-report delay | ● |  | ● | ● | ● | ● |  |
# | Report-to-lab delay |  |  |  | ● |  |  |  |
# | PCR sensitivity $s$ |  |  |  | ● |  |  |  |
# | PCR specificity |  |  |  | ● |  |  |  |
# | Severe-first share $q_0, q_\infty, \text{decay}$ |  |  |  | ● |  |  |  |
# | Forwarded fraction $\tau_{\text{forward}}$ |  |  |  | ● |  |  |  |
# | Background rate $\lambda_{\text{bg}}$ |  |  | ● | ● |  |  |  |
# | Surveillance dispersion |  | ● | ● | ● |  |  |  |
# | Ascertainment | ● |  | ● | ● | ● | ● |  |
# | Traveller volume | ● |  |  |  | ● | ● |  |
#
# The model components, in the order they appear below:
#
# 1. **Building-block submodels** — one per parameter family, listed in
#    the table above. Each samples its own priors and returns a small
#    NamedTuple of values, introducing only the maths for its own
#    parameters.
# 2. **Observation submodels** — exports, deaths, cases, the
#    laboratory pipeline (confirmed cases and samples received),
#    deaths-among-exports (time-resolved), and the first
#    export-detection timing. Each takes the growth state as input,
#    introduces the
#    forward integral it needs and the likelihood, and ties one data
#    stream to the latent $C(T)$.
# 3. **Composers** — one per analysis: the five single-stream fits, a
#    two-stream reimplementation of the report's methods (exports and
#    deaths), and the full joint fit. Each samples the building blocks
#    and the relevant observation
#    submodels. A composer conditionally includes only the likelihoods
#    for the streams it uses, so a single-stream fit never instantiates
#    the other observation submodels.
#
# Beyond the model itself:
#
# 4. **Inference** — prior predictive, the four No-U-Turn Sampler
#    (NUTS) fits, posterior
#    summaries, posterior-predictive plots.
# 5. **Counterfactual, forecast and sense check** — a
#    no-onward-transmission lower bound on cumulative deaths, a
#    one-week-ahead forecast, and a `Turing.fix`-pinned reproduction
#    of Method 2 main scenario via the exports-and-deaths
#    composer.

#md # ```@setup main
#md # using BVDOutbreakSize, CodeTracking, Revise
#md # ```

# #### Building-block submodels
#
# Each building-block submodel introduces only the mathematical objects
# and priors for one parameter family; the likelihoods and forward
# integrals enter later, in the observation submodels that use them.
#
# The implementation approach taken here is based on the hantavirus
# modelling project [hantavirus_2026](@cite): Mooncake
# [mooncake_jl](@cite) automatic differentiation, the
# Integrals [integrals_jl](@cite) quadrature helpers, a NaN-safe
# NegBinomial, FlexiChains, and PairPlots [pairplots_jl](@cite) with
# AlgebraOfGraphics [danisch2021makie](@cite) for the figures.

# ##### Growth — exponential
#
# The outbreak is seeded $T$ days ago by a single zoonotic case and
# grows exponentially with doubling time $\tau$, giving the cumulative-
# incidence trajectory
#
# ```math
# C(s) = \exp(r\,s), \qquad r = \frac{\log 2}{\tau}, \tag{1}
# ```
#
# so that the cumulative case count at the cut-off is $C(T) = 2^m$
# with $m = T/\tau$ the number of doublings since seeding. McCabe et al.'s
# primary assumption is the doubling time, which they vary over a
# sensitivity sweep of 7 / 14 / 21 days; each choice of doubling time
# implies a growth rate $r = \log 2/\tau$. We place the prior on that
# implied growth rate $r$ directly, centred at the main-scenario doubling
# time (14 d, so $r = \log 2/14$) with log-SD 0.4:
#
# ```math
# r \sim \mathrm{LogNormal}(\log(\log 2 / 14),\ 0.4),
# \qquad
# \tau = \frac{\log 2}{r}. \tag{2}
# ```
#
# Because $r = \log 2/\tau$ is a reciprocal, putting this LogNormal on $r$
# is exactly equivalent to a $\mathrm{LogNormal}(\log 14, 0.4)$ prior on
# the doubling time: a reciprocal negates and shifts the log-scale mean
# but preserves the log-scale SD 0.4, so the implied doubling-time prior
# is unchanged (95% interval roughly $(6, 31)$ d, spanning the full
# sweep), as is every derived quantity. Only the sampled coordinate
# differs.
#
# Rather than sampling $T$ directly (ridge-correlated with $r$ through
# $C(T) = \exp(r T)$), the model samples the *doubling count*
# $m = T/\tau$. Then $C(T) = 2^m$ is near-orthogonal to $r$. We give $m$
# a truncated-Normal prior with SD 3, centred on a base assumption that
# advances with the cut-off date:
#
# ```math
# m \sim \mathrm{Normal}(m_0,\ 3)\ \text{on}\ (0, \infty),
# \qquad
# m_0 = 9 + \frac{\text{cut-off} - \text{18 May 2026}}{14}. \tag{3}
# ```
#
# The base assumption is McCabe et al.'s first report (18 May 2026): its
# Method 2 central scenario of 501 cases is a doubling count
# $m = \log_2 501 \approx 9$. Each day that the cut-off runs past that
# report adds a fraction of a doubling at the central 14-day doubling
# time, so the size prior stays centred on the plausible outbreak as the
# data are refreshed rather than being fixed at the report-date value.
# The doubling time $\tau$, the elapsed time $T = m\cdot\tau$ and $C(T)$
# are exposed as deterministics so they appear in posterior tables and
# pair plots.

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

# ##### Incubation period
#
# The growth trajectory above counts *infections*. People are only seen
# once they develop symptoms, so the onsets that drive every observed
# stream are the infections shifted forward by the incubation period.
# Treating the trajectory as onsets directly would date the outbreak
# wrongly and understate how many people are already infected. We
# therefore sample the incubation period as a delay, like every other
# delay in the model, and convolve it with the infection trajectory to
# get onsets; the per-stream reporting delays then act on those onsets.
#
# The incubation period cannot be fitted from the BDBV line list, which
# has no exposure dates [bdbv_linelist_analysis_2026](@cite), so we take
# the Bundibugyo estimate from the 2007 Uganda outbreak: a mean of 6.3
# days, 95% CI 5.2-7.3 [macneil2010](@cite). We put the prior on the
# mean and the coefficient of variation so MacNeil et al.'s reported
# uncertainty on the mean is carried through directly, and recover the
# gamma shape and scale from them. MacNeil et al. give no interval on the
# spread, so the coefficient of variation has a weakly-informative prior
# chosen to keep individual incubation times within their observed 2-20
# day range.

#md # ```@raw html
#md # <details><summary>Submodel: incubation_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.incubation_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Genetic seeding bound
#
# A BEAST time tree of the first ten sequenced genomes
# [virological2026](@cite) places the TMRCA, the age of the oldest
# internal node of the tree, at a mean of 25 March 2026, with a 95% HPD
# interval of about $\pm 30$ days.
# The temporal sampling range is too short to estimate the clock, so it
# is fixed. The source analysis considers two literature rates for the
# 2013–2016 West African Ebola epidemic [holmes2016](@cite): a
# $1.2\times10^{-3}$ substitutions/site/year rate across all public
# data, and a faster $1.9\times10^{-3}$ early-epidemic rate that dates
# the TMRCA more recently. We use the $1.2\times10^{-3}$ rate in the
# main analysis and the $1.9\times10^{-3}$ rate in the
# [clock-rate sensitivity](#Clock-rate-sensitivity).
# This is a lower bound on the seeding time $T$. Adding sequences, or
# more geographically representative ones, can only push the TMRCA
# earlier, never later, as the sampled tree is almost entirely from Bunia.
# Combining the genetic TMRCA with the other data streams as a seeding
# bound follows a suggestion of N. Ferguson [ferguson2026](@cite).
# We parameterise the bound as an uncertain threshold
# $B \sim \mathrm{Normal}(g, \sigma)$, where
# $g = t_{\mathrm{cut}} - t_{\mathrm{TMRCA}}$ is the data cut-off date
# minus the reported TMRCA date, so it tracks the cut-off rather than a
# fixed offset. We require $T \ge B$, leaving $T$ free above it.
# Marginalising over $B$ gives a soft one-sided likelihood,
#
# ```math
# p_\text{gen}(T) = \Pr[B \le T] = \Phi\!\left(\frac{T - g}{\sigma}\right),
# \qquad \sigma = 15\ \text{d}. \tag{3a}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: genetic_seeding_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.genetic_seeding_model(100.0, 50.0; tmrca_days_sd = 15.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# Observing $g$ (`tmrca_days`) at the upper censoring point of
# $\mathrm{Normal}(T, \sigma)$ contributes the log probability of the
# censored upper tail, which is the soft bound of Eq. (3a):
#
# ```math
# \Pr[\mathrm{Normal}(T, \sigma) \ge g]
#   = 1 - \Phi\!\left(\frac{g - T}{\sigma}\right)
#   = \Phi\!\left(\frac{T - g}{\sigma}\right).
# ```

# ##### Onset-to-death delay
#
# Following McCabe et al., we assume the symptom-onset-to-death delay is
# gamma distributed with shape $\alpha$ and scale $\theta$, with density
# $f$ and CDF $F_d$:
#
# ```math
# \text{delay} \sim \mathrm{Gamma}(\alpha,\ \theta). \tag{4}
# ```
#
# The McCabe et al. report uses the point estimate of
# [rosello2015](@cite). We instead use the companion Bayesian reanalysis
# of the same Isiro line list [bdbv_linelist_analysis_2026](@cite),
# which re-estimates the delay with uncertainty. The two prior means are
# that reanalysis' posterior estimates of the gamma shape and scale, and
# the two standard deviations are set so each prior reproduces the
# reanalysis' 95% credible interval on that parameter. The published
# uncertainty therefore enters the fit directly rather than collapsing
# onto a single point estimate:
#
# ```math
# \alpha \sim \mathrm{Normal}^{+}(4.3,\ 1.22), \qquad
# \theta \sim \mathrm{Normal}^{+}(2.6,\ 0.82). \tag{5}
# ```
#
# This implies a prior mean onset-to-death delay of
# $\alpha\,\theta \approx 11$ days. The delay estimation in that
# reanalysis follows the recommendations of [charniga2024](@cite).

#md # ```@raw html
#md # <details><summary>Submodel: delay_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.delay_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Onset-to-report delay
#
# Suspected and laboratory-confirmed counts see the growth curve
# through one delay: an onset-to-report delay $f_{\text{rep}}$ from
# symptom onset to the case appearing on the suspected line list. The
# laboratory pipeline reads from the same reported backlog rather than
# adding a second delay; the lab selects which received samples
# to run, as in the confirmed-case likelihood below.
#
# The delay is Gamma-distributed, with shape and scale given
# truncated-Normal priors in the same way as the onset-to-death delay
# above. The prior is taken from the companion BDBV linelist reanalysis
# of the 2012 Isiro outbreak [bdbv_linelist_analysis_2026](@cite), whose
# Gamma-family posterior on the onset-to-notification delay has median
# around $11$ days, loosened slightly to allow for 2026-specific
# deviations:
#
# ```math
# \alpha_{\text{rep}} \sim \mathrm{Normal}^{+}(2.5,\ 1.0), \qquad
# \theta_{\text{rep}} \sim \mathrm{Normal}^{+}(4.5,\ 1.5), \tag{5a}
# ```
#
# a prior mean onset-to-report delay of about $11$ days. The shape and
# scale are truncated at $0.1$ to keep the Gamma well-defined under the
# sampler.

#md # ```@raw html
#md # <details><summary>Submodel: report_delay_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.report_delay_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: test_selection_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.test_selection_model()), "\n```"))
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

# ##### Case-fatality ratio
#
# The US Centers for Disease Control and Prevention (CDC) summary for
# the two previous BVD outbreaks is $55$ deaths in $169$ cases,
# $\approx 33\%$
# ([CDC outbreak history](https://www.cdc.gov/ebola/outbreaks/index.html)),
# with confidence bands spanning
# roughly $26$-$40\%$. The companion Bundibugyo virus (BDBV) reanalysis
# reports a baseline of $0.47$, $95\%$ CrI $0.31$-$0.65$, for
# non-healthcare-worker (non-HCW) confirmed cases. The prior on the
# case-fatality ratio is
#
# ```math
# \mathrm{CFR} \sim \mathrm{Beta}(6.6,\ 13.4), \tag{6}
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

# ##### Detection window
#
# $w$ is the mean time during which a case is still infectious and
# detectable abroad, the incubation plus onset-to-detection time. This
# rectangular window is the McCabe et al. mechanism, kept for the
# comparison through `imperial_only_model`; the default export model
# instead uses the onset-to-detection delay convolution described next.
# The prior is based on the detection windows McCabe et al. sweep in
# their Method 1 scenarios of 10, 15 and 20 days. It is centred on their
# central 15-day value with an SD wide enough to cover the 10–20 day
# range.
#
# ```math
# w \sim \mathrm{Normal}^{+}(15,\ 5)\ \text{days}. \tag{7}
# ```

#md # ```@raw html
#md # <details><summary>Submodel: detection_window_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.detection_window_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Infection→detection delay
#
# The default export model replaces McCabe's fixed detection window with
# an explicit infection→detection delay convolution.
# Exports are travel-gated, so the at-risk clock starts at infection. A
# person travels during incubation, before symptoms, and is detected
# abroad only after the full infection→detection delay has elapsed.
# The at-risk prevalence is the difference between cumulative
# infections and the infections that have already completed
# infection→detection, and integrating this prevalence over the per-day
# per-capita travel rate gives the expected detected exports.
# The window form is recovered exactly as the infection→detection delay
# collapses to a point mass, so this generalises the McCabe assumption.
# The infection→detection delay is the incubation period convolved with
# the DRC onset-to-report delay $f_{\text{rep}}$ (`report_delay_model`),
# moment-matched to a single Gamma via `combined_delay`. Entering
# surveillance as a suspected case is taken to be the same process abroad
# as in the DRC, and the export count is a single datum that cannot
# identify its own delay, so the incubation and reported-cases streams
# pin it.
# Because the incubation period sits inside this delay, the export
# likelihoods carry no separate incubation rescale, and the survival is 1
# at infection, unlike a flat onset rescale.
# Its mean is the incubation mean ($\approx 6.3$ days) plus the
# onset-to-report mean ($\approx 11.25$ days), $\approx 17.5$ days.
# The moment-match is an approximation: the sum of two Gammas is not
# Gamma, but matching the first two moments is accurate for the convex
# survival integrals here.
# The `imperial_only_model` composer keeps the rectangular detection
# window for comparison.

# ##### Daily traveller volume
#
# The number of people crossing from the source area to Uganda each day
# sets the travel rate in the exports likelihood. We treat it as an
# estimated quantity rather than a fixed input. McCabe et al. Table 3
# records mean weekly passenger counts across seven points of entry; the
# Ituri-side daily total of $1871$ is a sample mean across roughly
# $15$-$21$ point-of-entry-weeks. We use a Normal prior centred on
# $1871$ with SD $200$, $\approx 10\%$ CV, truncated at zero, covering
# point-of-entry variation and the sitrep sampling uncertainty. The
# source population is kept fixed at the census value.

#md # ```@raw html
#md # <details><summary>Submodel: traveller volume</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.traveller_volume_model()), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Surveillance dispersion
#
# We assume the surveillance counts are reported with negative binomial
# observation error around their expected value, with a single shared
# dispersion $k$ for the suspected-death, reported-case and
# confirmed-case streams. Under the mean-$\mu$ / dispersion-$k$
# parameterisation a count $Y$ has
#
# ```math
# Y \sim \mathrm{NegBinomial}(\mu,\ k), \qquad
# \mathrm{Var}(Y) = \mu + \frac{\mu^2}{k}. \tag{8}
# ```
#
# The dispersion captures passive-surveillance noise, namely
# under-reporting that varies by district, weekend reporting effects, and
# batched updates, not transmission heterogeneity.
# We judge this noise to be substantial, so a priori we expect
# meaningful overdispersion rather than near-Poisson counts.
# Following the Stan prior-choice recommendations
# [stan_prior_choice](@cite), the dispersion is sampled on the
# $1/\sqrt{k}$ scale, which behaves like a standard deviation, with a
# weakly-informative prior centred on that expected overdispersion,
#
# ```math
# 1/\sqrt{k} \sim \mathrm{Normal}^{+}(0.6,\ 0.2), \tag{9}
# ```
#
# giving $k$ a prior median near $3$ with a 90% range of about $1$-$14$.
# Because the prior allows near-Poisson counts, $k$ itself ranges over
# many orders of magnitude, so the pair plots and summary table show
# dispersion on both the sampled $1/\sqrt{k}$ scale, which is easier to
# read, and the more familiar $k$ scale.
# This extends the McCabe et al. report, which uses a Poisson likelihood
# for the Method 2 deaths and does not model the reported case counts at
# all. The negative binomial adds overdispersion to absorb
# passive-surveillance noise.

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

# ##### Ascertainment — partial pooling between DRC and Uganda
#
# Two surveillance systems detect cases: DRC passive surveillance (the
# reported suspected-case count) and Uganda's point-of-entry / hospital
# surveillance (the exported-case count). Each captures only a fraction
# of the true cases passing through it, and each fraction is informed
# by essentially a single aggregate data point, the one reported-case
# total and the one export count, so neither is well identified on its
# own. We therefore centre the prior on an assumed reporting fraction
# of 25% and partially pool the two fractions so they share strength.
# Treating them as identical would conflate two different systems,
# while treating them as independent would leave the Uganda fraction
# almost wholly prior-driven.
#
# Both ascertainment fractions $p_{\text{DRC}}$ and $p_{\text{Uganda}}$
# share a logit-scale hyperprior with mean $\mu$ and SD $\tau$:
#
# ```math
# \mu \sim \mathrm{Normal}(\mathrm{logit}(0.25),\ 1),
# \qquad
# \tau \sim \mathrm{Normal}^{+}(0,\ 0.5), \tag{10}
# ```
#
# ```math
# \mathrm{logit}(p_{\text{DRC}}) \sim
#     \mathrm{Normal}(\mu,\ \tau),
# \qquad
# \mathrm{logit}(p_{\text{Uganda}}) \sim
#     \mathrm{Normal}(\mu,\ \tau), \tag{11}
# ```
#
# with $p = \mathrm{logistic}(\mathrm{logit}\,p)$. Here $\tau$ is the
# pooling strength: small $\tau$ pulls the two fractions together towards
# the shared-fraction limit, large $\tau$ lets them move independently
# towards the separate-fraction limit. The cases likelihood uses
# $p_{\text{DRC}}$; the two Uganda-side likelihoods use
# $p_{\text{Uganda}}$.
#
# We sample this prior in its non-centred form: draw offsets
# $z_{\text{DRC}}, z_{\text{Uganda}} \sim \mathrm{Normal}(0, 1)$ and set
# $\mathrm{logit}(p) = \mu + \tau z$. This is the same prior but avoids
# the funnel geometry of the centred form, which gave hundreds of
# divergent transitions.

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

# We model both the deaths and the reported cases with a negative
# binomial, so we define one small constructor for it and share it.
# It is parameterised by mean $\mu$ and dispersion $k$, so the variance
# is given by equation (8), with NaN / Inf-safe clamping on the
# success probability so extreme NUTS proposals during warmup do not
# trip the distribution domain check. It is used by the deaths and
# cases observation submodels below.

#md # ```@raw html
#md # <details><summary>Function: safe_nbinomial</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.safe_nbinomial(1.0, 1.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# #### Observation submodels
#
# Each observation submodel takes the growth state as input,
# introduces the forward integral it needs, and ties one data stream
# to the latent $C(T)$. The forward integrals are solved numerically:
# the at-risk person-time integral for exports, the gamma convolution
# for deaths, and the deaths-among-exports convolution. Each submodel
# refers back to the parameters defined in equations (1)-(11).

# ##### Exports — Method 1 (geographic spread)
#
# Each case in the source population travels to Uganda on any given
# day with probability
# $q = \text{daily travellers}/\text{source population}$,
# treating cases as exchangeable with the general population. A case
# is *detection-eligible* for $w$ days from infection (equation (7)).
# For a case infected at outbreak age $s \leq T$, the accumulated
# probability of being detected in Uganda by $T$ is
#
# ```math
# P(\text{detected by } T \mid \text{infected at } s)
#     = q \cdot \min(T - s,\ w). \tag{12}
# ```
#
# Splitting at $s = T - w$ (full window elapsed before $T-w$, partial
# window after) and summing across incidence $i(s)$ gives the full
# export integral
#
# ```math
# \mathbb{E}[\text{exports by }T]
#     = q \cdot \Bigl[ w \cdot C(T-w)
#          + \int_{T-w}^{T} i(s) \, (T - s) \, ds \Bigr], \tag{13}
# ```
#
# which integration by parts collapses to the at-risk person-time form
# using the cumulative-incidence trajectory $C(s)$ of
# equation (1):
#
# ```math
# \mathbb{E}[\text{exports by }T] = q \cdot \int_{T-w}^{T} C(s)\, ds. \tag{14}
# ```
#
# For exponential growth this evaluates to
# $q\cdot(C(T) - C(T-w))/r$. We evaluate equation (14) numerically so
# the same form works for any growth parameterisation, scale by the
# Uganda ascertainment fraction $p_{\text{Uganda}}$ (equation (11)),
# and apply a Poisson likelihood:
#
# ```math
# \mu_e = p_{\text{Uganda}} \cdot q \cdot \int_{T-w}^{T} C(s)\, ds,
# \qquad
# Y_{\text{exports}} \sim \mathrm{Poisson}(\mu_e). \tag{15}
# ```
#
# !!! note "Comparison with McCabe et al. / Imai 2020"
#     McCabe et al. (and [imai2020](@cite) before them) use the
#     small-$rw$ simplification $\mu_e \approx q\cdot w\cdot C(T)$, the
#     limit of equation (14) as $r \to 0$.
#     For BVD's prior range $rw \in 0.33 - 2.0$ the simplification
#     under-estimates $C(T)$ by roughly $15$-$57\%$. We use a Poisson
#     likelihood for the detected exports; at the small detection
#     probability here ($p \approx q\cdot w \approx 6\cdot 10^{-3}$) it
#     is indistinguishable from a binomial detection model.

#md # ```@raw html
#md # <details><summary>Submodel: exports_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_model(1, nothing, 0.25)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: exports_delay_model (default)</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_delay_model(1, nothing, 0.25, BVDOutbreakSize.Gamma(2.5, 4.5))), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Deaths — Method 2 (back-calculation from deaths)
#
# Expected cumulative deaths by $T$ from a single seeding case is the
# CFR-weighted convolution of the cumulative-incidence trajectory
# $C(s)$ (equation (1)) with the onset-to-death density $f$
# (equation (4)):
#
# ```math
# \mathbb{E}[D(T)] = \mathrm{CFR} \cdot
#     \int_0^T e^{r s}\, f(T - s)\, ds. \tag{16}
# ```
#
# For a gamma delay this integral has an exact closed form carrying a
# $\gamma(\alpha, (\beta + r)T)/\Gamma(\alpha)$ factor from the finite
# upper limit. McCabe et al. use the large-$T$ simplification
# $D(T) \approx \mathrm{CFR}\cdot C(T)\cdot(1 + r/\beta)^{-\alpha}$,
# valid for $T \gtrsim 12/(\beta+r)$, which drops that factor. We
# evaluate equation (16) numerically instead. Its Gamma CDF is
# differentiated under reverse-mode AD by a hand-written rule, with
# the shape-parameter derivative of the regularized incomplete gamma
# function from a Kummer series following the Stan Math Library
# [carpenter2015stanmath](@cite). The observed deaths follow the
# NegBinomial likelihood of equation (8) with the dispersion $k$ of
# equation (9), shared with the cases likelihood:
#
# ```math
# Y_{\text{deaths}} \sim \mathrm{NegBinomial}(\mathbb{E}[D(T)],\ k). \tag{17}
# ```
#
# The DRC death count is *suspected* deaths, so it may include both
# missed BVD deaths under-reported and non-BVD deaths over-reported
# under the suspected case definition. A multiplicative drift factor
# $p_{\text{deaths}}$ absorbs that drift, with a tight prior centred
# at unity:
#
# ```math
# p_{\text{deaths}} \sim \mathrm{Normal}(1.0,\ 0.05),\ \text{truncated at 0}. \tag{17a}
# ```
#
# We do not see a single death total but a run of INSP sitreps, each
# reporting the deaths recorded since the last. Mapping each sitrep
# date to elapsed time $s_v = T - (d_{\text{as of}} - d_v)$, the new
# deaths in sitrep $v$ have mean
#
# ```math
# \mu_v^{\text{deaths}} = p_{\text{deaths}}\,
#     \bigl(\mathbb{E}[D(s_v)] - \mathbb{E}[D(s_{v-1})]\bigr),
# \quad
# \Delta Y_v^{\text{deaths}} \sim \mathrm{NegBinomial}(\mu_v^{\text{deaths}}, k), \tag{17b}
# ```
#
# one term per sitrep sharing $k$, with the first running from
# $s_0 = 0$ so a single death total recovers equation (17).
# The suspected-case and confirmed streams below are modelled the same way.

#md # ```@raw html
#md # <details><summary>Submodel: deaths_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.deaths_model(
#md #         Int[], nothing, 1.0, Float64[])), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Reported cases
#
# $C(T)$ is the latent cumulative *infection* count (equation (1)), not a
# count of reported, tested or confirmed cases. The BVD-driven suspected
# term acts on the *onset* trajectory, infections mapped through the
# incubation period, written as $e^{rs}$ below in the same shorthand as
# the deaths convolution (equation (16)). Suspected reports also include
# test-negative referrals such as malaria or other febrile illness, whose
# rate we assume is set by background prevalence and surveillance
# intensity, not by BVD growth. We therefore model the suspected stream
# as the sum of a BVD-driven component and a non-BVD background that
# accrues with elapsed surveillance time:
#
# ```math
# \mu_{\text{BVD}} = p_{\text{DRC}}
#     \int_0^T e^{r s}\, f_{\text{rep}}(T - s)\, ds, \qquad
# \mu_{\text{bg}} = \lambda_{\text{bg}}\, T, \tag{18}
# ```
#
# ```math
# \mu_{\text{cases}} = \mu_{\text{BVD}} + \mu_{\text{bg}}, \qquad
# Y_{\text{cases}} \sim \mathrm{NegBinomial}(\mu_{\text{cases}},\ k). \tag{19}
# ```
#
# $\lambda_{\text{bg}}$ is the expected non-BVD suspected reports per
# day. The implied per-suspected positivity at the cut-off is
# $\pi = \mu_{\text{BVD}} / \mu_{\text{cases}}$, distinct from the
# per-test positivity used by the lab pipeline below. Half-Normal
# prior:
#
# ```math
# \lambda_{\text{bg}} \sim \mathrm{Normal}^{+}(0,\ 1)\ \text{per day}. \tag{20}
# ```
#
# The prior is informative because $\lambda_{\text{bg}}$ is degenerate
# with outbreak size: the per-bin reported mean is
# $p_{\text{DRC}}\,\Delta\mu_{\text{BVD}} + \lambda_{\text{bg}}\,\Delta t$,
# so a diffuse prior lets the background absorb arbitrarily many
# suspected cases and resolve at the high end where the deaths and
# exports streams pin $C(T)$. The background must not explain more
# suspected cases than were ever reported, so the SD of $1$ is chosen
# from the observed suspected total $1077$ at the $26$ May cut-off over
# the elapsed window $T \approx 132$ days. It puts the median background
# at $\approx 0.67$ per day, $\approx 89$ cases or $\approx 8\%$ of
# observed, and the $95\%$ prior bound at $\approx 2.0$ per day,
# $\approx 259$ cases or $\approx 24\%$ of observed. A wider SD
# ($\approx 1.5$) leaves a second posterior mode in which
# $\lambda_{\text{bg}}$ runs to $\approx 8$ per day and the background
# explains the majority of suspected cases; SD $1$ keeps the fit in the
# regime where the BVD trajectory drives the suspected total.
# $\lambda_{\text{bg}}$ is identified from the suspected/confirmed
# contrast and from the samples-received series, which conditions on the
# combined BVD and background backlog. The dispersion $k$ (equation (9))
# is shared with the deaths and confirmed likelihoods.
#
# $\mu_{\text{bg}} = \lambda_{\text{bg}}\, T$ assumes the non-BVD
# background rate is constant in time and independent of the outbreak.
# With a single cumulative suspected count we can only identify a
# cumulative rate, so we keep the constant-rate parameterisation.
#
# As for the deaths, we model each sitrep's new suspected cases rather
# than only the latest cumulative total. The DRC ascertainment fraction
# $p_{\text{DRC}}$ (equation (11)) is applied to every sitrep's
# increment. Writing the unit-ascertainment BVD-suspected cumulative as
# $\mu_{\text{BVD},0}(s) = \int_0^s e^{r u} f_{\text{rep}}(s - u)\,du$,
# the new suspected cases in sitrep $v$ have mean
#
# ```math
# \mu_v^{\text{rep}}
#   = p_{\text{DRC}}\,\bigl(\mu_{\text{BVD},0}(s_v)
#     - \mu_{\text{BVD},0}(s_{v-1})\bigr)
#     + \lambda_{\text{bg}}\,(s_v - s_{v-1}),
# \quad
# \Delta Y_v^{\text{rep}} \sim
#     \mathrm{NegBinomial}(\mu_v^{\text{rep}}, k). \tag{20}
# ```
#
# ##### Per-sitrep conditional predictive
#
# The per-sitrep increment likelihoods (equations (17b) and (20)) give a
# matching per-sitrep diagnostic.
# Let a stream have observed cumulative counts $y_1, \dots, y_n$ at the
# sitrep edges.
# The conditional one-step-ahead predictive of the cumulative at sitrep
# $v$ conditions on the *observed* previous cumulative and predicts only
# the new increment,
#
# ```math
# \hat{y}_v = y_{v-1} + \Delta_v, \qquad y_0 = 0, \tag{20a}
# ```
#
# where $\Delta_v$ is the posterior-predictive between-vintage increment
# for bin $v$, a draw from $\mathrm{NegBinomial}(\mu_v, k)$, the same
# increment the model samples in predictive mode.
# This is the "filtered" one-step-ahead predictive. It differs from
# reconstructing the cumulative as the running sum of the modelled
# increments $\sum_{u \le v} \Delta_u$, whose errors compound across
# sitreps because every step builds on earlier modelled increments
# rather than on what was observed. Conditioning each step on the
# observed previous cumulative isolates the model's one-step increment
# prediction, so a mismatch at one sitrep is attributed to that sitrep
# alone and does not propagate down the series.

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
#md # <details><summary>Submodel: reported_cases_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.reported_cases_model(
#md #         Int[], nothing, 1.0, Float64[], Float64[])), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Confirmed cases and samples received
#
# The laboratory pipeline observes two per-vintage series: the confirmed
# (PCR-positive) cases and the samples received from the suspected pool.
# Testing is selection from an accumulated backlog. Every suspected
# sample received stays eligible, and the lab draws which to run, so
# there is no report-to-lab delay. The confirmed series reads from the
# same reported backlog as the suspected stream.
#
# Each vintage's confirmed count is a Binomial draw on its observed
# number of samples analysed, with the analysed count $A_v$ a known
# denominator from the sitrep rather than a modelled quantity:
#
# ```math
# C_v \sim \mathrm{Binomial}(A_v,\ p_{\text{pos},v}). \tag{21}
# ```
#
# Per-test positivity mixes true and false positives. With PCR
# sensitivity $s$, specificity $\text{spec}$ and $q_v$ the BVD share of
# the analysed batch at vintage $v$,
#
# ```math
# p_{\text{pos},v} = s\, q_v + (1 - \text{spec})\,(1 - q_v). \tag{22}
# ```
#
# The lab tests the most-likely-BVD cases first, the obvious severe
# cluster, so the BVD share of the tested pool starts high and relaxes
# to a baseline as testing widens to the broad suspect pool:
#
# ```math
# q(c) = q_\infty + (q_0 - q_\infty)\, e^{-c / \text{decay}},
# \qquad c = \max(t - t_{\text{report}},\ 0), \tag{23}
# ```
#
# with $c$ the time since surveillance onset. $q_0$ (near 1) is the early
# severe-cluster BVD fraction, $q_\infty$ the broad-pool baseline, and
# $\text{decay}$ the timescale over which the share relaxes. The early
# vintages, where the BVD share is high, inform the sensitivity, and the
# plateau positivity $s\, q_\infty + (1 - \text{spec})(1 - q_\infty)$ holds
# the later vintages.
#
# The samples-received series conditions the fraction of suspects
# forwarded to the lab. The cumulative suspect backlog at each vintage is
# the BVD-suspected term of equation (18) plus the non-BVD background,
# $N_{\text{susp},v} = \mu_{\text{BVD}}(t_v) + \lambda_{\text{bg}}\,t_v$.
# The received count is a fraction $\tau_{\text{forward}}$ of that
# backlog,
#
# ```math
# R_v \sim \mathrm{NegBinomial}(\tau_{\text{forward}}\,
#     N_{\text{susp},v},\ k), \tag{24}
# ```
#
# which pins $\tau_{\text{forward}}$ directly from received-versus-
# suspected.
#
# The priors are a Beta on the PCR sensitivity, a high Beta on the
# specificity, a near-1 Beta on $q_0$, weakly-informative priors on
# $q_\infty$ and the decay timescale, and a Beta on the forwarded
# fraction:
#
# ```math
# s \sim \mathrm{Beta}(6,\ 2), \qquad
# \text{spec} \sim \mathrm{Beta}(50,\ 1.5), \qquad
# \tau_{\text{forward}} \sim \mathrm{Beta}(5,\ 2), \tag{25}
# ```
#
# ```math
# q_0 \sim \mathrm{Beta}(20,\ 1.5), \qquad
# q_\infty \sim \mathrm{Beta}(6,\ 6), \qquad
# \text{decay} \sim \mathrm{Normal}^{+}(0,\ 10)\ \text{days}. \tag{25a}
# ```
#
# Confirmation runs on the altona RealStar Filovirus Screen RT-PCR at
# INRB. The rapid Cepheid GeneXpert Ebola assay is
# Zaire-ebolavirus-specific and does not reliably detect Bundibugyo virus
# [cepheid_xpert_ebola_ifu](@cite). The RealStar kit detects Bundibugyo
# virus at $11$-$67$ RNA copies per reaction in laboratory and field
# evaluation [rieger2016](@cite) but has no published field sensitivity
# for this variant, and early low-viral-load specimens and field handling
# lower real-world detection. The $\mathrm{Beta}(6, 2)$ prior, mean
# $0.75$ and $95\%$ interval $\sim 0.39$-$0.97$, keeps good analytical
# sensitivity plausible while carrying substantial downside mass for those
# field losses. The specificity prior is tight and high, mean
# $\approx 0.97$, reflecting the assay's strong analytical specificity and
# keeping the false-positive term from absorbing the positivity signal.
#
# The forwarded fraction $\tau_{\text{forward}}$ has a
# $\mathrm{Beta}(5, 2)$ prior, mean $0.71$, expressing that a majority
# but not all of the suspected backlog is forwarded; the samples-received
# series pulls it to the value the data imply. The $q_0$ prior sits near
# 1, mean $\approx 0.93$; $q_\infty$ is centred at mean $0.5$ near the
# cut-off positivity-implied share; the decay prior has median
# $\approx 6.7$ days, spanning the lab window. The forwarded fraction and
# the share curve carry no outbreak-specific data beyond the received and
# confirmed series, so they are weakly informative.

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
#md # <details><summary>Submodel: confirmed_cases_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.confirmed_cases_model(
#md #         Int[], Int[], Union{Missing, Int}[], missing, nothing, 1.0,
#md #         Float64[], 0.0, 0.7, nothing, Float64[], 1.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Deaths among exports
#
# Uganda reports a single death among its detected exports. That count
# is informative. If the exports happened long ago, we would expect
# more deaths among them by now under the onset-to-death gamma, so the
# observed deaths-among-exports bound how recently the exports occurred
# and help constrain $T$ and $C(T)$. The expected count reuses the
# at-risk export integrand of equation (14) but weights each case by its
# probability of having died by $T$, the onset-to-death CDF
# $F_d(T - s)$ (equation (4)), then scales by the CFR, the travel rate
# $q$ and the Uganda ascertainment fraction $p_{\text{Uganda}}$:
#
# ```math
# \mathbb{E}[D_{\text{Uganda}}]
#     = \mathrm{CFR} \cdot p_{\text{Uganda}} \cdot q
#       \cdot \int_{T-w}^{T} C(s)\, F_d(T - s)\, ds. \tag{19}
# ```
#
# Equation (19) is evaluated numerically, writing $F_d$ as the inner
# integral of the density $f$ so the whole expression differentiates
# through $f$ alone, since the reverse-mode AD backend does not support
# the gamma CDF shape-parameter derivative. The detection window $w$ and
# daily traveller volume are shared with the exports likelihood so the
# two Uganda-side observations use the same person-time.
#
# Equation (19) is the McCabe-window form, kept for the comparison. In
# the default delay mechanism the top-hat window is replaced by the
# infection→detection survival $\overline{F}_{\text{det}}(T-s)$, and the
# death CDF $F_d$ is the infection→death delay, incubation $\oplus$
# onset-to-death, so both clocks run from infection in step with $C(s)$.
# Detection and death share the same onset, so incubation enters both
# delays, a slight accepted double-count of the shared incubation period.
#
# Uganda's export deaths are point-of-entry or hospital-detected, so
# their dates are recorded directly and carry information beyond the
# total count: a death seen early bounds how old the outbreak can be. We
# model the detected export deaths as an inhomogeneous Poisson process
# with cumulative intensity $\mathbb{E}[D_{\text{Uganda}}(t)]$
# (equation (19), at any elapsed time $t$) and use its time-resolved
# likelihood. Split $[0, T]$ at the earliest dated death, offset
# $\Delta_1$ before the cut-off at elapsed time $T-\Delta_1$. Before it
# no export death was seen, contributing one continuous survival term
# over $[0, T-\Delta_1]$. From that day to the cut-off each day $d$
# carries a Poisson count of the export deaths that day, with bin mean
# $\mu_d$,
#
# ```math
# \log L = \sum_d y_d \log \mu_d - \mathbb{E}[D_{\text{Uganda}}(T)],
# \quad
# \mu_d = \mathbb{E}[D_{\text{Uganda}}(t_d)]
#         - \mathbb{E}[D_{\text{Uganda}}(t_{d-1})]. \tag{20}
# ```
#
# This is the same per-bin construction as the DRC streams, differencing
# the cumulative expected count between successive dates as in equation
# (20b), with one addition: the continuous survival weight over the long
# pre-death zero stretch before the first dated death. That term
# collapses the run of pre-death zeros into a single weight, so there is
# no per-day vector of zeros, while the recent window is resolved per
# day. It is used here because the exact export-death dates are known;
# the DRC streams have no such fixed zero stretch and are differenced
# only across the reported sitrep dates. The number of daily bins is
# fixed by the earliest death's date, so the likelihood is well defined
# even though $T$ is latent. With one death the series is a single count,
# and it takes more dated deaths as they are reported.
#
# !!! note "Selection-bias caveat"
#     This assumes Uganda's surveillance retains detected exports
#     through to any subsequent death. If the system instead loses
#     cases that progress to death, the observed deaths-among-exports
#     count is selection-biased downward and the constraint it places
#     on $T$ is overstated.

#md # ```@raw html
#md # <details><summary>Submodel: exports_deaths_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_deaths_model( Int[], nothing, 0.33, nothing, 0.25; window = 15.0, daily_travellers = 1871.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: exports_deaths_delay_model (default)</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_deaths_delay_model( Int[], nothing, 0.33, nothing, 0.25; f_det = BVDOutbreakSize.Gamma(2.0, 4.5), daily_travellers = 1871.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### First export detection — timing survival term
#
# The same logic applies to the *first detected export case* (Uganda's
# first hospital admission), using the at-risk export person-time
# intensity $\mathbb{E}[\text{exports}(t)]$ (equation (15)) in place of
# the export-death intensity. With $\Delta$ the offset from the first
# admission date to the cut-off and $t_1 = T - \Delta$,
#
# ```math
# \Pr(\text{no export detected before } t_1)
#     = \exp\!\bigl(-\mathbb{E}[\text{exports}(t_1)]\bigr). \tag{22}
# ```
#
# As with the export-death term, this is one-sided and only marginally
# constrains the posterior because the Uganda detections sit only days
# before the cut-off.

#md # ```@raw html
#md # <details><summary>Submodel: exports_detection_timing_model</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_detection_timing_model( nothing, 0.25; delta = missing, window = 15.0, daily_travellers = 1871.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Submodel: exports_detection_timing_delay_model (default)</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_detection_timing_delay_model( nothing, 0.25; delta = missing, f_det = BVDOutbreakSize.Gamma(2.0, 4.5), daily_travellers = 1871.0)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# #### Composers
#
# These composers combine the building blocks into the full model for
# each analysis. McCabe et al. invert a
# deterministic summary formula at fixed nuisance parameters; here we
# sample all of them — growth, delay, CFR, detection
# window, traveller volume, dispersion, ascertainment — and condition
# on the observed counts. Each composer includes only the likelihoods
# for the streams it carries, so a single-stream composer never
# instantiates the other observation submodels.
# Of the observation streams, the geographic-spread exports reproduce
# McCabe et al.'s Method 1 and the back-calculation from deaths their
# Method 2; the reported-cases ascertainment, the deaths-among-exports,
# the export-detection-timing, and the genetic seeding terms are
# extensions with no counterpart in McCabe et al.
#
# The joint composer samples a single dispersion scale $k$ and passes it
# to both the deaths and cases likelihoods, so they share one passive-
# surveillance noise scale. It samples one pooled set of
# ascertainment fractions, threading $p_{\text{DRC}}$ to the cases
# likelihood and $p_{\text{Uganda}}$ to the two Uganda-side likelihoods.
# The window $w$ and daily traveller volume sampled by the exports
# likelihood are reused by the deaths-among-exports likelihood so the
# two share person-time.
#
# We write single-stream composers for the four count-based streams
# only.
# The export-detection-timing and genetic seeding terms constrain the
# outbreak start $T$ rather than the size directly and are weakly
# identified in isolation, so we do not fit them on their own. The
# joint composer still conditions on both.

# ##### Exports-only fit — Method 1 analogue

#md # ```@raw html
#md # <details><summary>Composer: exports-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_only_model(1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Deaths-only fit — Method 2 analogue

#md # ```@raw html
#md # <details><summary>Composer: deaths-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.deaths_only_model(1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Cases-only fit — ascertainment extension
#md # ```@raw html
#md # <details><summary>Composer: cases-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.cases_only_model(1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Deaths-among-exports-only fit
#md # ```@raw html
#md # <details><summary>Composer: exports-deaths-only fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.exports_deaths_only_model(Int[])), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### Joint fit — full posterior over $C(T)$ from all data streams
#
# The joint composer is the same generative process conditioned on all four
# observed streams simultaneously. Each stream argument may be passed
# as `missing` to drop it; the matching likelihood is then not
# instantiated. Passing all streams as `missing` turns the composer into
# a generator for the prior- and posterior-predictive checks.

#md # ```@raw html
#md # <details><summary>Composer: joint fit</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.bvd_joint(1, [1], [1];
#md #         reported_offsets = [0])), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ##### McCabe et al. reimplementation — exports and deaths only
#
# McCabe et al.'s joint configuration uses exactly two data sources: the
# geographic-spread exports (Method 1) and the back-calculation from
# deaths (Method 2). It has no reported-cases ascertainment model and
# no deaths-among-exports likelihood. This composer wraps just
# those two observation submodels, so the sense-check can fix the model
# to the McCabe et al. joint configuration. Either stream may
# be `missing`; passing `missing` for exports recovers a pure
# deaths-only Method 2 fit.
#
# McCabe et al. estimate a case count and do not separate infection from
# symptom onset, so this reimplementation fits no incubation period and
# its trajectory is the case count directly. Our joint model instead
# treats the trajectory as infections and recovers cases downstream
# through the incubation period. To compare on the same footing we
# report the case count from both, our joint $C(T)$ mapped to onsets and
# the reimplementation's $C(T)$ taken as cases. The infection count is
# reported only for the joint model.

#md # ```@raw html
#md # <details><summary>Composer: report reimplementation</summary>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, CodeTracking, Markdown
#md # Markdown.parse(string("```julia\n",
#md #     (@code_string BVDOutbreakSize.imperial_only_model(1, 1)), "\n```"))
#md # ```

#md # ```@raw html
#md # </details>
#md # ```

# ### Model fitting and evaluation
#
# #### Prior predictive check
#
# Before any observation is taken into account, what does the joint
# prior imply about replicated exports, deaths, reported cases and
# deaths among exports? Draws from the prior over the unobserved data
# should bracket the observed counts.

#md # ```@raw html
#md # <details><summary>Sample the joint prior</summary>
#md # ```

## Per-vintage observation arguments for the joint model. Each DRC
## stream is fitted as between-vintage increments of its cumulative
## sitrep history (the first increment is the cumulative at the first
## vintage), so a stream with a single vintage reduces to the
## cumulative single-total likelihood. `joint_obs` packages the
## increment vectors and offsets; `observe = false` swaps the counts
## for `missing` so the same model structure generates predictive draws.
function _increments(v)
    d = similar(v, Int)
    prev = 0
    for i in eachindex(v)
        d[i] = v[i] - prev
        prev = v[i]
    end
    return d
end

function joint_obs(o; observe = true)
    _stream(h,
        s) = h === missing ?
             (Union{Missing, Int}[observe ? s : missing], [0]) :
             (observe ? _increments(h.values) :
              fill(missing, length(h.values)), h.offsets)
    rep, rep_off = _stream(o.reported_case_history, o.reported_cases)
    dth, dth_off = _stream(o.death_history, o.total_deaths)
    ## Confirmed cases enter PER VINTAGE over the lab vintages that carry an
    ## analysed denominator (`tests_analysed_history`, the 23-28 May
    ## sitreps with denominators 211/295/295/403/648/755). Each cumulative
    ## confirmed count is observed as a Binomial on its analysed count, so
    ## the joint conditions on the per-vintage positivity trajectory
    ## (exhaustion + modelled specificity in `confirmed_cases_model`) rather
    ## than a single cumulative total. Confirmed counts are aligned to the
    ## analysed offsets. When no per-vintage analysed history is available,
    ## fall back to the single cumulative total at the cut-off.
    have_conf = o.confirmed_case_history !== missing ||
                o.confirmed_cases !== missing
    have_pervintage = have_conf &&
                      o.confirmed_case_history !== missing &&
                      o.tests_analysed_history !== missing
    if have_pervintage
        sa = o.tests_analysed_history
        ch = o.confirmed_case_history
        idx = [findfirst(==(off), ch.offsets) for off in sa.offsets]
        any(isnothing, idx) &&
            error("confirmed history missing an analysed-vintage offset")
        ## Merge any vintage that analysed no new samples (ΔA = 0, e.g. the
        ## 25 May stall) into the next, so every window has a positive
        ## denominator. Keep the first vintage and any with more cumulative
        ## analysed than the previous kept one, then difference to the
        ## between-vintage increments the confirmed model conditions on.
        keep = [i == 1 || sa.values[i] > sa.values[i - 1]
                for i in eachindex(sa.values)]
        aoff = collect(sa.offsets)[keep]
        analysed = Union{Missing, Int}[_increments(sa.values[keep])...]
        ccum = [ch.values[i] for i in idx][keep]
        conf = observe ? Union{Missing, Int}[_increments(ccum)...] :
               fill(missing, length(ccum))
        conf_off = aoff
        ## Tests-received increments (`Cumul échantillons reçus`), aligned
        ## to the kept offsets, condition τ_forward via the received-count
        ## NegBinomial.
        if o.tests_received_history !== missing
            sr = o.tests_received_history
            ridx = [findfirst(==(off), sr.offsets) for off in aoff]
            received = any(isnothing, ridx) ? Union{Missing, Int}[] :
                       (observe ?
                        Union{Missing, Int}[_increments([sr.values[i]
                                                                             for i in ridx])...] :
                        fill(missing, length(ridx)))
        else
            received = Union{Missing, Int}[]
        end
    else
        conf_total = o.confirmed_cases !== missing ? o.confirmed_cases :
                     o.confirmed_case_history === missing ? missing :
                     o.confirmed_case_history.values[end]
        conf,
        conf_off = have_conf ?
                   (Union{Missing, Int}[observe ? conf_total : missing],
            [0]) : (Union{Missing, Int}[], Int[])
        analysed = Union{Missing, Int}[]
        received = Union{Missing, Int}[]
    end
    ## Travel-gated export streams stop at the most recent reported import
    ## to Uganda. Cross-border movement patterns likely shift over the
    ## outbreak and the days after the last import are right-truncated by
    ## reporting lag, so the trailing zeros are dropped from both the
    ## export-case and deaths-among-exports series and each stream is run
    ## only up to that date (see `exports_daily_delay_model`).
    ec_full = o.exported_cases_daily
    last_import = isempty(ec_full) ? nothing : findlast(!=(0), ec_full)
    export_last_offset = last_import === nothing ? 0 :
                         length(ec_full) - last_import
    _truncate(v) = v[1:max(length(v) - export_last_offset, 0)]
    ecases = isempty(ec_full) ? ec_full :
             (observe ? _truncate(ec_full) :
              fill(missing, length(_truncate(ec_full))))
    ed_trunc = _truncate(o.export_deaths_daily)
    edaily = observe ? ed_trunc : fill(missing, length(ed_trunc))
    return (deaths = dth, reported = rep, export_deaths = edaily,
        kw = (; reported_offsets = rep_off, death_offsets = dth_off,
            confirmed_cases = conf, confirmed_offsets = conf_off,
            samples_analysed = analysed,
            samples_received = received,
            exported_cases_daily = ecases,
            export_last_offset = export_last_offset,
            tests_analysed = observe ? o.cumulative_tests_analysed :
                             missing, tests_offset = 0))
end

## Dummy non-missing confirmed/tested counts instantiate the laboratory
## submodel so its priors (s, spec, the severe-first selection q0/qinf/
## decay_scale, plus the derived per-test positivity) appear in the prior
## chain for the lab-pipeline pair plot. Under `Prior()` the likelihood is
## discarded, so the placeholder values do not influence the sampled priors.
prior_args = joint_obs(obs; observe = false)
prior_chn = sample(
    bvd_joint(missing, prior_args.deaths, prior_args.reported,
        prior_args.export_deaths; prior_args.kw...,
        confirmed_q_random_effect = confirmed_q_re_model),
    Prior(), 2_000; progress = false);

prior_C_table = summary_table(prior_chn, [:cumulative_cases]; digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

prior_C_table #hide

# Pair plot of the prior over the latent quantities, for spotting prior
# correlations before any data has been seen.

#md # ```@raw html
#md # <details><summary>Prior pair plot</summary>
#md # ```

prior_pair_fig = plot_pair(prior_chn,
    [:r, :τ, :m, :cumulative_cases, :CFR, :α_rep, :inv_sqrt_k, :k,
        :p_drc, :p_uganda, :τ_logit,
        :λ_bg, :τ_forward, :s_test, :spec_test,
        :q0, :qinf, :decay_scale,
        :positivity, :p_positive]);

#md # ```@raw html
#md # </details>
#md # ```

prior_pair_fig #hide

# #### Fitting the models
#
# NUTS [hoffman2014nuts](@cite) with Mooncake [mooncake_jl](@cite)
# reverse-mode automatic differentiation, four chains, 1000 post-warmup
# draws each, with a target acceptance probability of 0.95. Chains
# initialise from the prior to keep the sampler away from the boundary
# of $r$ and $m$. We fit the joint model and the four single-stream
# models so the per-stream posteriors over $C(T)$ can be compared with
# the joint.

#md # ```@raw html
#md # <details><summary>Run the joint and per-stream NUTS fits</summary>
#md # ```

genetic_seeding = T -> genetic_seeding_model(T, obs.genetic_tmrca_days;
    tmrca_days_sd = obs.genetic_tmrca_days_sd)

## Growth submodel whose doubling-count prior centre advances with the
## cut-off date (base value at McCabe's first report; see
## `m_prior_centre`), so every fit uses the size prior appropriate to its
## own `as_of_date`.
function growth_for(o)
    exponential_growth_model(
        m_prior = truncated(Normal(m_prior_centre(o.as_of_date), 3.0);
        lower = 0))
end
growth_now = growth_for(obs)

## Per-vintage tested-BVD-share random effect for the confirmed stream:
## the 28 May per-window positivity (0.48, 0.05, 0.15, 0.02, 0.79) is
## non-monotone, so a partially-pooled logit-scale offset lets each
## confirmed vintage fit its own positivity while the assay sensitivity
## stays fixed (see `confirmed_q_re_model`).
fit_args = joint_obs(obs)
chn_joint = nuts_sample(
    bvd_joint(obs.exported_cases, fit_args.deaths, fit_args.reported,
        fit_args.export_deaths; fit_args.kw...,
        growth = growth_now,
        first_export_detection_delta = obs.first_export_detection_delta,
        report_onset_offset = report_onset_offset(obs.as_of_date),
        confirmed_q_random_effect = confirmed_q_re_model,
        genetic = genetic_seeding); trace_kw("joint")...);

chn_exports = nuts_sample(
    exports_only_model(obs.exported_cases; growth = growth_now);
    trace_kw("exports")...);
chn_deaths = nuts_sample(
    deaths_only_model(obs.total_deaths; growth = growth_now);
    trace_kw("deaths")...);
chn_cases = nuts_sample(
    cases_only_model(obs.reported_cases; growth = growth_now);
    trace_kw("cases")...);
chn_confirmed = nuts_sample(
    confirmed_only_model(obs.confirmed_cases, obs.cumulative_tests_analysed,
        obs.tests_received_history.values[end]; growth = growth_now);
    trace_kw("confirmed")...);
chn_exports_deaths = nuts_sample(
    exports_deaths_only_model(obs.export_deaths_daily; growth = growth_now);
    trace_kw("exports_deaths")...);

posterior_C_joint = vec(Array(chn_joint[:cumulative_cases]));
posterior_C_infections_joint = vec(Array(chn_joint[:cumulative_infections]));
posterior_C_exports = vec(Array(chn_exports[:cumulative_cases]));
posterior_C_deaths = vec(Array(chn_deaths[:cumulative_cases]));
posterior_C_cases = vec(Array(chn_cases[:cumulative_cases]));
posterior_C_confirmed = vec(Array(chn_confirmed[:cumulative_cases]));
posterior_C_exports_deaths = vec(Array(chn_exports_deaths[:cumulative_cases]));

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
    "confirmed (DRC)" => chn_confirmed) #hide

#md # ```@raw html
#md # </details>
#md # ```

# ### Additional analyses
#
# On top of the main analysis we run four supporting analyses.
#
# #### Counterfactual: no-onward-transmission lower bound
#
# Suppose every onward transmission stopped at the report date $T$. The
# cohort already infected by $T$ still carries future expected deaths in
# the onset-to-death tail: a case infected at outbreak age $s$ has died
# by the report date with probability $F_d(T - s)$ (equation (4)), so a
# fraction
# $1 - F_d(T - s)$ of its CFR-weighted contribution has not yet been
# observed. Integrating against the incidence
# $i(s) = r\cdot\exp(r\cdot s)$ from
# equation (1) gives the additional future expected deaths
#
# ```math
# \Delta D = \mathrm{CFR} \cdot \int_0^T r\,\exp(r\,s)
#            \,\bigl(1 - F_d(T - s)\bigr)\,ds, \tag{23}
# ```
#
# and a lower bound on the cumulative-death endpoint of
# $D(T) + \Delta D$, evaluated per posterior draw.
#
# #### One-week-ahead forecast
#
# If the fitted model is taken at face value, what should the next
# week's situation reports show? We continue the fitted exponential
# growth seven days past the report date $T$ and apply the same
# observation models to forecast the cumulative reported cases (DRC),
# deaths (DRC) and exports (Uganda) by $T + 7$. We also forecast the
# *new* counts expected over the coming week, the cumulative at $T + 7$
# minus the count already observed. These are posterior-predictive: each
# draw yields a replicated integer count, so the intervals carry both
# parameter and observation uncertainty.
#
# This assumes growth continues unchanged for the week, with no
# interventions and no saturation, so it is a "no-change" projection.
#
# #### Delay sensitivity
#
# The deaths back-calculation (equation (16)) depends on the onset-to-
# death delay. The baseline fit builds the gamma shape $\alpha$ and
# scale
# $\theta$ from the all-deaths Isiro mixture (equation (5)). The companion
# reanalysis [bdbv_linelist_analysis_2026](@cite) also reports a
# *community-only* pathway, the $n = 5$ cases who died without hospital
# admission, with a shorter but far more uncertain delay: a shape of
# about $5.6$ ($95\%$ CrI $1.0$-$25.9$) and a scale of about $1.4$
# ($95\%$ CrI $0.3$-$9.5$). A shorter delay means deaths appear sooner
# after infection, so a given death count back-calculates to a smaller
# $C(T)$.
#
# We refit the joint model once with the onset-to-death delay priors
# rebuilt from the community-only pathway, building truncated-Normal
# priors from those credible intervals exactly as the baseline delay
# priors (equation (5)) are constructed:
#
# ```math
# \alpha \sim \mathrm{Normal}^{+}(5.6,\ 6.35), \qquad
# \theta \sim \mathrm{Normal}^{+}(1.4,\ 2.35). \tag{24}
# ```
#
# The comparison shows how sensitive the
# outbreak-size estimate is to the delay assumption.
#
# !!! warning "Sensitivity only, not a preferred estimate"
#     The community-only delay is fitted from $n = 5$ deaths, so the
#     evidence is weak and the priors in equation (24) are very wide.
#     This section is included to probe sensitivity, not as a preferred
#     alternative to the baseline.
#
# #### Report reproduction and validation
#
# How does our joint posterior sit against what McCabe et al. reported,
# and how much of any difference is the method rather than the newer
# data? The report itself was revised: the 18 May version
# [mccabe2026](@cite) used $88$ deaths and a central CFR of $30\%$,
# while the 20 May version [mccabe2026update](@cite) used $131$ deaths
# and corrected the central CFR to $33\%$. We work through both in
# sequence so the effect of the report's own deaths-plus-CFR correction
# is visible separately from our newer data and joint method. For each
# version we fit our full joint model to that version's data snapshot
# (`data/report-snapshot.toml` for 18 May, `report-snapshot-20may.toml`
# for 20 May), so the only difference between a version's fit
# and our headline fit is the data.
#
# McCabe et al. Method 2 reports Poisson intervals (no overdispersion,
# $k \to \infty$). We reproduce it by fixing the exports-and-deaths
# composer to their Method 2 central assumptions and dropping exports.
# $1/\sqrt{k}$ is fixed to a small positive value ($k \approx 10^6$,
# Poisson-like) because exactly $0$ gives $k \approx 4.5\times10^{15}$,
# where the NegBinomial saturates and reverse-mode AD returns NaN
# gradients.
#
# As a sense check we ask whether our machinery recovers McCabe et
# al.'s Method 2 headline when given their inputs. The 18 May reported
# Method 2 central estimate is $501$ cases and the 20 May one is $678$
# cases. Each reproduction drops exports so only the deaths likelihood
# is instantiated, conditions on that version's deaths ($88$ for 18 May,
# $131$ for 20 May), and `Turing.fix`-pins the Method 2 main-scenario
# values (the growth rate to the $\tau = 14$ d scenario via
# $r = \log 2/14$, $\mathrm{CFR} = 30\%$ for 18 May and $33\%$ for
# 20 May, $\alpha = 4.42$, $\beta = 0.388$/d), with the deaths
# NegBinomial made Poisson-like. The only sampled latent is $m$, the
# number of doublings since seeding ($C(T) = 2^m$). A close match
# confirms our back-calculation matches the report. Any remaining gap
# to our headline estimate is method (joint fit, exact convolution,
# sampled nuisance parameters) and newer data.
#
# This sense check covers the deaths (Method 2) side. The exports
# (Method 1) side differs by construction: we use the exact cumulative
# integral $q\int_{T-w}^{T} C(s)\,ds$ rather than the small-$rw$
# simplification $q\,w\,C(T)$, so our exports-implied size is expected
# to sit above a Method 1 reproduction (by the $15$-$57\%$ noted
# earlier) rather than match it.

# ## Results
#
# ### Summary
#
# For the response the question that matters is how many people have
# already been infected. The reported counts capture only part of the
# outbreak, and planning for beds, contacts and vaccine needs depends
# on the true total. The numbers below are our current best estimate of
# that total, computed from the joint posterior and refreshed on every
# build. For each headline number we give the equal-tailed 30%, 60% and
# 90% credible intervals. The same intervals appear in the tables below.

#md # ```@raw html
#md # <details><summary>Compute the headline ranges</summary>
#md # ```

summary_ranges = let
    med(x) = quantile(x, 0.5)
    iqr(x) = quantile(x, 0.75) - quantile(x, 0.25)
    ## Posterior-minus-prior shift in units of the parameter's prior
    ## IQR, reusing the prior draws so nothing is respecified here.
    shift(post, prior) = round((med(post) - med(prior)) / iqr(prior); digits = 2)
    ints_i(s) = string(
        "30% ", round(Int, s.lo30), "–", round(Int, s.hi30),
        ", 60% ", round(Int, s.lo60), "–", round(Int, s.hi60),
        ", 90% ", round(Int, s.lo90), "–", round(Int, s.hi90))
    ints_f(s,
        d) = string(
        "30% ", round(s.lo30; digits = d), "–", round(s.hi30; digits = d),
        ", 60% ", round(s.lo60; digits = d), "–", round(s.hi60; digits = d),
        ", 90% ", round(s.lo90; digits = d), "–", round(s.hi90; digits = d))
    ## Seeding-start dates derived from the T posterior. Higher T means
    ## earlier seeding, so the start-date range flips the lo/hi labels.
    start_from(t) = Date(obs.as_of_date) - Day(round(Int, t))
    ints_d(s) = string(
        "30% ", start_from(s.hi30), "–", start_from(s.lo30),
        ", 60% ", start_from(s.hi60), "–", start_from(s.lo60),
        ", 90% ", start_from(s.hi90), "–", start_from(s.lo90))

    C = posterior_C_joint
    Td = vec(Array(chn_joint[:T]))
    τd = vec(Array(chn_joint[:τ]))
    rd = vec(Array(chn_joint[:r]))
    sC = posterior_summary(C)
    sT = posterior_summary(Td)
    sτ = posterior_summary(τd)
    sr = posterior_summary(rd)
    f_lo = round(sC.lo90 / obs.reported_cases; digits = 1)
    f_hi = round(sC.hi90 / obs.reported_cases; digits = 1)

    moves = ["cumulative case load" => shift(C, vec(Array(
            prior_chn[:cumulative_cases]))),
        "time since seeding" => shift(Td, vec(Array(prior_chn[:T]))),
        "doubling time" => shift(τd, vec(Array(prior_chn[:τ])))]
    biggest = argmax(p -> abs(p.second), moves)

    Markdown.parse("""
    - **Current cumulative case load:** we estimate $(ints_i(sC)) cases,
      combining all four data streams (reported and as-yet-unreported).
    - That is roughly $(f_lo)–$(f_hi)× the $(obs.reported_cases) cases
      reported to date, so most infections are not yet reported. This
      multiplier is one over the DRC reporting fraction; see
      [what the reporting fraction means](#Joint-model-estimates).
    - **Time since seeding:** we estimate $(ints_i(sT)) days, placing
      the start of sustained transmission at $(ints_d(sT)).
    - **Doubling time and growth rate:** we estimate a doubling time of
      $(ints_f(sτ, 1)) days, and an implied growth rate of
      $(ints_f(sr, 3)) per day.
    - **Shift from priors:** how far the data has moved each estimate
      from its prior, measured in prior interquartile ranges (IQRs) — a
      value of 1 means the posterior median sits one prior IQR from the
      prior median, 0 means unchanged, and the sign gives the direction.
      The fit moves the cumulative case load by $(moves[1].second), the
      time since seeding by $(moves[2].second) and the doubling time by
      $(moves[3].second); the largest move is in the $(biggest.first).
    """)
end;

#md # ```@raw html
#md # </details>
#md # ```

summary_ranges #hide

# **Why our estimate is higher than McCabe et al.**
# Our central estimate sits above the McCabe et al. [mccabe2026](@cite)
# report for three reasons. We fit more data streams and average over the
# uncertainty in the parameters the report fixes one scenario at a time.
# We use more recent counts. And we fit each stream to its full run of
# situation reports rather than a single total, so the growth rate is
# informed by the shape of the reported trajectory.
# See [what we do differently](#What-we-do-differently-from-McCabe-et-al.),
# the [comparison with McCabe et al.](#Comparison-with-McCabe-et-al.) and
# the [limitations](#Limitations) for the detail behind this.

# ### Joint model estimates
#
# Our main result is an estimate of how far the outbreak has spread by
# the data cut-off, obtained by fitting all four data streams together:
# the cases exported to Uganda, the suspected deaths in the DRC, the
# reported cases in the DRC (with an ascertainment component) and the
# deaths among exported cases in Uganda.
#
# We report two cumulative counts, both as of the data cut-off. The
# headline quantity is the cumulative *infection* count, everyone
# infected by the cut-off whether or not they have yet developed
# symptoms. The cumulative *case* count is the subset whose symptoms have
# appeared by the cut-off. It is smaller than the infection count because
# some of the most recently infected are still incubating. It is the
# quantity McCabe et al. estimate, so we carry it for the comparison
# below. We report each as a credible-interval table and a posterior
# density, infections first.

#md # ```@raw html
#md # <details><summary>Cumulative infection count summary table</summary>
#md # ```

cumulative_infections_summary = summary_table(
    chn_joint, [:cumulative_infections]; digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_infections_summary #hide

#md # ```@raw html
#md # <details><summary>Cumulative case count summary table</summary>
#md # ```

cumulative_cases_summary = summary_table(
    chn_joint, [:cumulative_cases]; digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_cases_summary #hide

# Against the 210 laboratory-confirmed cases recorded by the cut-off,
# the estimated infections give the under-ascertainment multiplier,
# infections per confirmed case. We report it alongside the case
# fatality ratio as headline quantities.

#md # ```@raw html
#md # <details><summary>Headline multiplier and CFR summary table</summary>
#md # ```

headline_summary = let
    mult = posterior_C_infections_joint ./ obs.confirmed_cases
    cfr = vec(Array(chn_joint[:CFR]))
    df = DataFrame(quantity = String[],
        lower_90 = Float64[], lower_60 = Float64[], lower_30 = Float64[],
        upper_30 = Float64[], upper_60 = Float64[], upper_90 = Float64[])
    for (
        name, xs) in ["Infections per confirmed case" => mult,
        "Case fatality ratio" => cfr]
        s = posterior_summary(xs)
        push!(df,
            (name, round(s.lo90; digits = 2), round(s.lo60; digits = 2),
                round(s.lo30; digits = 2), round(s.hi30; digits = 2),
                round(s.hi60; digits = 2), round(s.hi90; digits = 2)))
    end
    rename(df,
        ["Quantity", "Lower 90%", "Lower 60%", "Lower 30%",
            "Upper 30%", "Upper 60%", "Upper 90%"])
end;

#md # ```@raw html
#md # </details>
#md # ```

headline_summary #hide

#md # ```@raw html
#md # <details><summary>Cumulative infection and case count densities</summary>
#md # ```

joint_infections_density_fig = plot_cumulative_cases(
    "infections" => posterior_C_infections_joint; scenarios = []);

joint_density_fig = plot_cumulative_cases(
    "cases" => posterior_C_joint; scenarios = []);

#md # ```@raw html
#md # </details>
#md # ```

joint_infections_density_fig #hide

joint_density_fig #hide

# Reading these together:
#
# - **Infections to date** is the headline count, $C(T) = \exp(r T)$:
#   everyone infected by the data cut-off.
# - **Cases (symptom onset to date)** is the smaller subset who have
#   already developed symptoms by the cut-off, $C(T)$ scaled by the share
#   past the incubation period. The gap between the two is the recently
#   infected still incubating, and it is this case count that lines up
#   with McCabe et al. and the reported surveillance counts.
#
# The cumulative infection count $C(T) = \exp(r T)$ is set jointly by the
# doubling time $\tau$, equivalently the growth rate $r = \log 2/\tau$,
# and the time since seeding $T$. Read as a calendar date, $T$ places the
# start of sustained transmission at the report date minus $T$ days.
# The left panel below shows the posterior for that start date. The
# right panel shows the joint $(\tau, T)$ posterior, which is positively
# correlated: slower growth (larger $\tau$) needs a longer elapsed $T$
# to reach the same observed counts.

#md # ```@raw html
#md # <details><summary>Outbreak start date and (τ, T) posterior</summary>
#md # ```

start_date_fig = plot_start_date_pair(chn_joint;
    as_of_date = obs.as_of_date);

#md # ```@raw html
#md # </details>
#md # ```

start_date_fig #hide

# The full posterior summary table reports equal-tailed 30%, 60% and
# 90% credible intervals on the key joint-fit parameters: growth rate
# $r$, doubling count $m$, days since seeding $T$, CFR, the
# DRC and Uganda ascertainment fractions $p_{\text{DRC}}$ and $p_{\text{Uganda}}$, the
# pooling SD $\tau_{\text{logit}}$, the surveillance dispersion on both
# the sampled $1/\sqrt{k}$ scale and the more familiar $k$ scale, the
# laboratory-pipeline parameters — onset-to-report delay shape and scale
# $\alpha_{\text{rep}}$, $\theta_{\text{rep}}$, PCR sensitivity $s$,
# specificity $\text{spec}$, the forwarded fraction
# $\tau_{\text{forward}}$, the severe-first share parameters $q_0$,
# $q_\infty$ and $\text{decay}$, the non-BVD background rate
# $\lambda_{\text{bg}}$, the per-suspected positivity $\pi$ and the
# per-test positivity $p_{\text{pos}}$ — and cumulative cases $C(T)$.

#md # ```@raw html
#md # <details><summary>Joint posterior summary table</summary>
#md # ```

joint_summary = summary_table(chn_joint,
    [:r, :τ, :m, :T, :CFR, :p_drc, :p_uganda, :τ_logit,
        :inv_sqrt_k, :k, :α_rep, :θ_rep,
        :s_test, :spec_test, :τ_forward, :λ_bg, :q0, :qinf, :decay_scale,
        :positivity, :p_positive, :q_cutoff, :q_baseline,
        :cumulative_infections, :cumulative_cases]; digits = 2);

#md # ```@raw html
#md # </details>
#md # ```

joint_summary #hide

# The posterior pair plot shows the joint distribution of the key
# parameters, with the prior overlaid so the data's contribution to
# each marginal is visible.

#md # ```@raw html
#md # <details><summary>Posterior pair plot (prior overlaid)</summary>
#md # ```

posterior_pair_fig = plot_pair(chn_joint,
    [:r, :τ, :m, :cumulative_cases, :CFR, :α_rep, :inv_sqrt_k, :k,
        :p_drc, :p_uganda, :τ_logit];
    prior = prior_chn);

#md # ```@raw html
#md # </details>
#md # ```

posterior_pair_fig #hide

# ### Reporting process
#
# The streams above all observe the same latent infections: infections
# become onsets after the incubation period, onsets enter the suspected
# line list after the onset-to-report delay, and the lab selects which of
# the received samples to confirm. The observed suspected count mixes an
# *accepted* (BVD-attributable) part $\mu_{\text{BVD}}$ and a non-BVD
# background $\mu_{\text{bg}}$ (equation (18)). Only the total is
# observed. The accepted share is the per-suspected positivity
# $\pi = \mu_{\text{BVD}} / \mu_{\text{cases}}$ in the summary table, so
# $\pi$ times the suspected total recovers the unobserved accepted-BVD
# count. This differs again from the lab-confirmed count, since only the
# tested samples enter the confirmed series.
#
# A pair plot covers the laboratory-pipeline parameters that the
# confirmed and samples-received streams add: the onset-to-report delay
# shape and scale $\alpha_{\text{rep}}$, $\theta_{\text{rep}}$, PCR
# sensitivity $s$, specificity $\text{spec}$, the forwarded fraction
# $\tau_{\text{forward}}$, the severe-first share parameters $q_0$,
# $q_\infty$ and $\text{decay}$, and the non-BVD background rate
# $\lambda_{\text{bg}}$, against cumulative cases $C(T)$. The prior
# is overlaid so the contribution of the lab observations to each
# marginal is visible. The delay and sensitivity priors are only weakly
# updated, while $\tau_{\text{forward}}$ and $\lambda_{\text{bg}}$ are
# informed by the received and confirmed series.

#md # ```@raw html
#md # <details><summary>Laboratory-pipeline pair plot (prior overlaid)</summary>
#md # ```

lab_pair_fig = plot_pair(chn_joint,
    [:α_rep, :θ_rep, :s_test, :spec_test, :τ_forward,
        :λ_bg, :q0, :qinf, :decay_scale, :cumulative_cases];
    prior = prior_chn);

#md # ```@raw html
#md # </details>
#md # ```

lab_pair_fig #hide

# The DRC reporting fraction $p_{\text{DRC}}$ is the share of true cases
# that reach the reported suspected-case count. The reported total
# therefore scales up to the cumulative case load by about
# $1/p_{\text{DRC}}$, the multiplier quoted in the
# [summary](@ref "Summary"): a reporting fraction near $0.25$ implies a
# roughly fourfold gap between reported and true cases. The pair plot
# above shows its posterior against the prior. How far below one the
# fraction sits is what sets that scaling.
#
# The confirmed-case stream enters the joint fit per vintage, each
# vintage a Binomial on its observed analysed denominator (equation
# (21)). The joint therefore conditions on the positivity trajectory
# across the lab vintages, with the severe-first share carrying the fall
# in positivity from the early severe cluster to the broad-pool baseline
# (equation (23)). The samples-received series sets the forwarded
# fraction $\tau_{\text{forward}}$, and the per-vintage positivity sets
# the sensitivity and the share curve.
#
# A posterior predictive check draws replicated observations from the
# fitted joint model and compares them to the observed counts. If the
# fit is reasonable the observed value (red line) sits inside the bulk
# of its replicate distribution. The four panels are the four data
# streams: exported cases and deaths among exports in Uganda, and deaths
# and reported cases in the DRC.

#md # ```@raw html
#md # <details><summary>Joint posterior predictive plot</summary>
#md # ```

pp_args = joint_obs(obs; observe = false)
pp_joint = predict(
    bvd_joint(missing, pp_args.deaths, pp_args.reported,
        pp_args.export_deaths; pp_args.kw...,
        pre_start_deaths = missing,
        pre_detection_exports = missing,
        first_export_detection_delta = obs.first_export_detection_delta,
        report_onset_offset = report_onset_offset(obs.as_of_date),
        confirmed_q_random_effect = confirmed_q_re_model,
        genetic = genetic_seeding),
    chn_joint);
## Exports are now a dated per-day series (to the last reported import),
## so sum each draw's days for the cumulative-total posterior predictive,
## matching the observed total.
pp_exports = vec(sum.(pp_joint[@varname(exported_cases_daily)]));
## All DRC streams are now per-vintage increment vectors (deaths, reported,
## confirmed positives, samples received); sum each draw's bins for the
## cumulative-total posterior predictive. Export deaths are a per-day series
## summed the same way.
pp_deaths = vec(sum.(pp_joint[@varname(total_deaths)]));
pp_cases = vec(sum.(pp_joint[@varname(reported_cases)]));
pp_confirmed = vec(sum.(pp_joint[@varname(confirmed_cases)]));
pp_tests = vec(sum.(pp_joint[@varname(samples_received)]));
pp_exports_deaths = vec(sum.(pp_joint[@varname(export_deaths_daily)]));

joint_ppc_fig = plot_posterior_predictive(
    pp_exports, pp_deaths,
    obs.exported_cases, obs.total_deaths;
    pp_cases = pp_cases,
    obs_cases = obs.reported_cases,
    pp_exports_deaths = pp_exports_deaths,
    obs_exports_deaths = obs.exports_deaths,
    pp_tests = pp_tests,
    obs_tests = obs.tests_received_history.values[end],
    pp_confirmed = pp_confirmed,
    obs_confirmed = obs.confirmed_cases);

#md # ```@raw html
#md # </details>
#md # ```

joint_ppc_fig #hide

# ### Conditional one-step-ahead predictive across the sitrep series
#
# The panels above collapse each DRC stream to a single total. Because
# the suspected-case and death streams are fitted per vintage, we can
# also check the fit one sitrep at a time. For each vintage we condition
# on what was actually observed at the previous sitrep and predict only
# the new between-vintage increment, $\hat{y}_v = y_{v-1} + \Delta_v$
# with $y_0 = 0$, defined in the [observation model](@ref "Per-sitrep
# conditional predictive"). Each step carries the full posterior
# uncertainty of the new increment while conditioning on the observed
# previous cumulative, so the diagnostic isolates the model's per-step
# increment prediction and its errors do not compound across the series.
# The observed cumulative counts are overlaid as points. If the fit is
# reasonable each point sits inside the predictive ribbon for its
# vintage.

#md # ```@raw html
#md # <details><summary>Per-vintage conditional one-step-ahead predictive</summary>
#md # ```

## Confirmed cases span the merged lab vintages (the 25 May stall folded
## into 26 May); align the observed cumulative to the kept offsets.
conf_keep = [i == 1 ||
             obs.tests_analysed_history.values[i] >
             obs.tests_analysed_history.values[i - 1]
             for i in eachindex(obs.tests_analysed_history.values)]
conf_cidx = [findfirst(==(off), obs.confirmed_case_history.offsets)
             for off in obs.tests_analysed_history.offsets]
vintage_ppc_fig = plot_vintage_conditional_ppc([
    (; title = "Suspected cases",
        dates = obs.reported_case_history.dates,
        replicates = collect(pp_joint[@varname(reported_cases)]),
        observed = obs.reported_case_history.values, colour = :steelblue),
    (; title = "Confirmed cases",
        dates = obs.tests_analysed_history.dates[conf_keep],
        replicates = collect(pp_joint[@varname(confirmed_cases)]),
        observed = [obs.confirmed_case_history.values[i]
                    for i in conf_cidx][conf_keep], colour = :seagreen),
    (; title = "Suspected deaths",
        dates = obs.death_history.dates,
        replicates = collect(pp_joint[@varname(total_deaths)]),
        observed = obs.death_history.values, colour = :firebrick)]);

#md # ```@raw html
#md # </details>
#md # ```

vintage_ppc_fig #hide

# ### Counterfactual: lower bound under no further transmission
#
# The lower bound on cumulative deaths if transmission stopped at the
# report date: still-expected and projected-total deaths per draw.

#md # ```@raw html
#md # <details><summary>Project no-onward deaths and summarise</summary>
#md # ```

no_onward = predict_no_onward_deaths(chn_joint; obs_deaths = obs.total_deaths);

no_onward_table = streams_table(
    "no-onward total" => no_onward.total_projected;
    digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

no_onward_table #hide

# Two panels. On the left, the *still expected* deaths $\Delta D$, future
# deaths in cases already infected by $T$ net of those already observed.
# On the right, the *projected total* $D(T) + \Delta D$ with the observed
# death count marked as a dashed black rule.

#md # ```@raw html
#md # <details><summary>No-onward projected-deaths plot</summary>
#md # ```

no_onward_fig = plot_no_onward_deaths(no_onward; obs_deaths = obs.total_deaths);

#md # ```@raw html
#md # </details>
#md # ```

no_onward_fig #hide

# ### One-week-ahead forecast
#
# The seven-day no-change projection: cumulative and new expected counts
# per stream by $T + 7$. Exports are dropped from this projection.
# Cross-border travel is unlikely to be continuing at its baseline rate,
# so the forward travel rate the export forecast relies on no longer
# holds.

#md # ```@raw html
#md # <details><summary>Generate the one-week-ahead forecast</summary>
#md # ```

forecast = forecast_reported(chn_joint;
    horizon = 7,
    daily_travellers = ITURI_DAILY_TRAVEL,
    source_population = ITURI_POPULATION,
    obs_cases = obs.reported_cases,
    obs_deaths = obs.total_deaths,
    obs_confirmed = obs.confirmed_cases,
    obs_tests = obs.tests_received_history.values[end],
    obs_analysed = obs.cumulative_tests_analysed,
    forecast_exports = false,
    report_onset_offset = report_onset_offset(obs.as_of_date));
forecast_summary = forecast_table(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_summary #hide

# New counts expected over the coming week, by stream:

#md # ```@raw html
#md # <details><summary>One-week-ahead forecast plot</summary>
#md # ```

forecast_fig = plot_forecast(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_fig #hide

# ### Forecast validation against later data
#
# The one-week-ahead forecast above projects from the current fit, so it
# cannot yet be checked. We instead validate the same machinery
# retrospectively: take our joint fit to the original McCabe et al.
# report's data, project each posterior draw forward to the current data
# cut-off, and compare the predicted cumulative counts against the counts
# observed since.

#md # ```@raw html
#md # <details><summary>Fit the joint model to the original report's data and forecast it forward</summary>
#md # ```

obs_report = load_observations(
    joinpath(pkgdir(BVDOutbreakSize), "data", "report-snapshot.toml"));

report_args = joint_obs(obs_report)
chn_joint_report = nuts_sample(
    bvd_joint(obs_report.exported_cases, report_args.deaths,
    report_args.reported, report_args.export_deaths; report_args.kw...,
    growth = growth_for(obs_report),
    first_export_detection_delta =
    obs_report.first_export_detection_delta));
posterior_C_joint_report = vec(Array(chn_joint_report[:cumulative_cases]));

validation_horizon = value(Date(obs.as_of_date) - Date(obs_report.as_of_date))

## The report-snapshot fit predates the laboratory streams (no
## `confirmed_cases` or `cumulative_tests_analysed` in its observation
## toml), so the chain does not carry the lab-pipeline parameters and
## the validation forecast covers the three streams that were
## available at the snapshot date only.
forecast_validation = forecast_reported(chn_joint_report;
    horizon = validation_horizon,
    daily_travellers = ITURI_DAILY_TRAVEL,
    source_population = ITURI_POPULATION,
    obs_cases = obs_report.reported_cases,
    obs_deaths = obs_report.total_deaths,
    obs_exports = obs_report.exported_cases);

forecast_validation_table = forecast_vs_truth(forecast_validation;
    cases = obs.reported_cases,
    deaths = obs.total_deaths,
    exports = obs.exported_cases);

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@eval
#md # using BVDOutbreakSize, Markdown
#md # using Dates: Date, value
#md # rep = load_observations(joinpath(pkgdir(BVDOutbreakSize), "data",
#md #     "report-snapshot.toml")).as_of_date
#md # cur = load_observations().as_of_date
#md # h = value(Date(cur) - Date(rep))
#md # Markdown.parse("Projecting the original-report fit $(h) days " *
#md #     "forward to the current ($(cur)) data, against the counts " *
#md #     "observed by then:")
#md # ```

forecast_validation_table #hide

# The top row shows the cumulative forecast per stream and the bottom row
# the new counts over the horizon, mirroring the one-week-ahead forecast.
# Each panel shades the 90% predictive interval and draws the
# later-observed count as a dashed rule, so coverage can be read off
# directly. The confirmed and samples-received streams are absent here.
# The validation refits the original report snapshot, which predates any
# laboratory data, so that fit carries no lab parameters to project
# forward.

#md # ```@raw html
#md # <details><summary>Forecast-validation plot</summary>
#md # ```

forecast_validation_fig = plot_forecast_vs_truth(forecast_validation;
    cases = obs.reported_cases,
    deaths = obs.total_deaths,
    exports = obs.exported_cases,
    baseline_cases = obs_report.reported_cases,
    baseline_deaths = obs_report.total_deaths,
    baseline_exports = obs_report.exported_cases);

#md # ```@raw html
#md # </details>
#md # ```

forecast_validation_fig #hide

# The check above scores only the endpoint at the current cut-off. Since
# The INSP reports give a cumulative count at every sitrep date, so we can
# score the whole horizon: project the same original-report fit forward
# to each vintage date between the report and now, and compare against
# the count observed at that date. Each row is one stream at one date,
# with the 90% predictive interval and whether the observed cumulative
# fell inside it, so forecast coverage can be read across the horizon
# rather than at a single point.

#md # ```@raw html
#md # <details><summary>Score the report fit against the observed daily trajectory</summary>
#md # ```

forecast_trajectory_table = forecast_vs_truth_trajectory(chn_joint_report;
    dates = obs.reported_case_history.dates,
    cases = obs.reported_case_history.values,
    deaths = obs.death_history.values,
    snapshot_date = obs_report.as_of_date,
    daily_travellers = ITURI_DAILY_TRAVEL,
    source_population = ITURI_POPULATION,
    baseline_cases = obs_report.reported_cases,
    baseline_deaths = obs_report.total_deaths);

#md # ```@raw html
#md # </details>
#md # ```

forecast_trajectory_table #hide

# ### Delay sensitivity
#
# Refit under the community-only onset-to-death delay: the baseline and
# refitted $C(T)$ posteriors side by side.

#md # ```@raw html
#md # <details><summary>Refit the joint model with the community-only delay</summary>
#md # ```

community_delay = delay_model(;
    alpha_prior = truncated(Normal(5.6, (25.9 - 1.0) / 3.92); lower = 0),
    theta_prior = truncated(Normal(1.4, (9.5 - 0.3) / 3.92); lower = 0))

chn_joint_community = nuts_sample(
    bvd_joint(obs.exported_cases, fit_args.deaths, fit_args.reported,
    fit_args.export_deaths; fit_args.kw...,
    growth = growth_now,
    first_export_detection_delta = obs.first_export_detection_delta,
    genetic = genetic_seeding,
    deaths = (total_deaths,
        growth_state,
        k,
        t_edges;
        kwargs...) -> deaths_model(total_deaths, growth_state, k, t_edges;
        delay = community_delay, kwargs...)));

posterior_C_community = vec(Array(chn_joint_community[:cumulative_cases]));

#md # ```@raw html
#md # </details>
#md # ```

# Fit diagnostics for the community-only delay refit.

#md # ```@raw html
#md # <details><summary>Fit diagnostics</summary>
#md # ```

diagnostics_table( #hide
    "joint (community-only delay)" => chn_joint_community) #hide

#md # ```@raw html
#md # </details>
#md # ```

# Baseline versus community-only delay, side by side:

#md # ```@raw html
#md # <details><summary>Delay-sensitivity C_T table</summary>
#md # ```

delay_sensitivity_table = streams_table(
    "joint (baseline delay)" => posterior_C_joint,
    "joint (community-only delay)" => posterior_C_community);

#md # ```@raw html
#md # </details>
#md # ```

delay_sensitivity_table #hide

# Overlaid posterior densities of $C(T)$ under the two delay
# assumptions:

#md # ```@raw html
#md # <details><summary>Delay-sensitivity C_T density plot</summary>
#md # ```

delay_sensitivity_fig = plot_cumulative_cases(
    "baseline delay" => posterior_C_joint,
    "community-only delay" => posterior_C_community;
    scenarios = []);

#md # ```@raw html
#md # </details>
#md # ```

delay_sensitivity_fig #hide

# ### Clock-rate sensitivity
#
# The main analysis fixes the molecular clock to the
# $1.2\times10^{-3}$ substitutions/site/year rate (see the genetic
# seeding bound above). The source analysis also reports a faster
# $1.9\times10^{-3}$ early-epidemic rate, without favouring either
# [virological2026](@cite). Under it the TMRCA is dated more recently.
# We refit the joint model under that alternative bound and compare the
# impact on the outbreak size $C(T)$, the seeding time $T$ and the
# growth rate $r$.

#md # ```@raw html
#md # <details><summary>Refit the joint model under the 1.9e-3 clock</summary>
#md # ```

genetic_seeding_clock19 = T -> genetic_seeding_model(T,
    obs.genetic_tmrca_alt_days;
    tmrca_days_sd = obs.genetic_tmrca_alt_days_sd)

chn_joint_clock19 = nuts_sample(
    bvd_joint(obs.exported_cases, fit_args.deaths, fit_args.reported,
    fit_args.export_deaths; fit_args.kw...,
    growth = growth_now,
    first_export_detection_delta = obs.first_export_detection_delta,
    genetic = genetic_seeding_clock19));

posterior_C_clock19 = vec(Array(chn_joint_clock19[:cumulative_cases]));

#md # ```@raw html
#md # </details>
#md # ```

# Fit diagnostics for the 1.9e-3 clock-rate refit.

#md # ```@raw html
#md # <details><summary>Fit diagnostics</summary>
#md # ```

diagnostics_table( #hide
    "joint (1.9e-3 clock)" => chn_joint_clock19) #hide

#md # ```@raw html
#md # </details>
#md # ```

# Each quantity is shown as a side-by-side table followed by overlaid
# posterior densities under the two clock rates.
#
# Outbreak size $C(T)$:

#md # ```@raw html
#md # <details><summary>Clock-rate C_T table</summary>
#md # ```

clock_sensitivity_C_table = streams_table(
    "joint (1.2e-3 clock)" => posterior_C_joint,
    "joint (1.9e-3 clock)" => posterior_C_clock19);

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_C_table #hide

#md # ```@raw html
#md # <details><summary>Clock-rate C_T density plot</summary>
#md # ```

clock_sensitivity_C_fig = plot_cumulative_cases(
    "1.2e-3 clock" => posterior_C_joint,
    "1.9e-3 clock" => posterior_C_clock19;
    scenarios = []);

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_C_fig #hide

# Seeding time $T$ (days before the cut-off); a more recent TMRCA
# permits later seeding:

#md # ```@raw html
#md # <details><summary>Clock-rate seeding-time table</summary>
#md # ```

T_clock12 = vec(Array(chn_joint[:T]));
T_clock19 = vec(Array(chn_joint_clock19[:T]));

clock_sensitivity_T_table = streams_table(
    "joint (1.2e-3 clock)" => T_clock12,
    "joint (1.9e-3 clock)" => T_clock19;
    digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_T_table #hide

#md # ```@raw html
#md # <details><summary>Clock-rate seeding-time density plot</summary>
#md # ```

clock_sensitivity_T_fig = plot_density_overlay(
    "1.2e-3 clock" => T_clock12,
    "1.9e-3 clock" => T_clock19;
    xlabel = "Seeding time T (days before cut-off)",
    title = "Posterior seeding time by clock rate");

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_T_fig #hide

# Growth rate $r$ (per day); later seeding implies faster growth to
# reach the same observed counts:

#md # ```@raw html
#md # <details><summary>Clock-rate growth-rate table</summary>
#md # ```

r_clock12 = vec(Array(chn_joint[:r]));
r_clock19 = vec(Array(chn_joint_clock19[:r]));

clock_sensitivity_r_table = streams_table(
    "joint (1.2e-3 clock)" => r_clock12,
    "joint (1.9e-3 clock)" => r_clock19;
    digits = 3);

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_r_table #hide

#md # ```@raw html
#md # <details><summary>Clock-rate growth-rate density plot</summary>
#md # ```

clock_sensitivity_r_fig = plot_density_overlay(
    "1.2e-3 clock" => r_clock12,
    "1.9e-3 clock" => r_clock19;
    xlabel = "Growth rate r (per day)",
    title = "Posterior growth rate by clock rate");

#md # ```@raw html
#md # </details>
#md # ```

clock_sensitivity_r_fig #hide

# ### How the data streams compare
#
# Each data stream constrains the latent outbreak size differently.
# The table below puts the posteriors over $C(T)$ side by side, the
# five single-stream fits and the joint, to show what each stream buys
# on its own and what the joint combination adds.
# The single-stream fits cover the four count-based streams and the
# laboratory pipeline (confirmed and samples-received fit together). The
# joint additionally conditions on the export-detection-timing and
# genetic seeding terms, which constrain $T$ rather than the size and
# so are not fit in isolation.

#md # ```@raw html
#md # <details><summary>Per-stream C_T table</summary>
#md # ```

streams_C_table = streams_table(
    "exports (cases)" => posterior_C_exports,
    "exports (deaths)" => posterior_C_exports_deaths,
    "deaths (DRC)" => posterior_C_deaths,
    "cases (DRC)" => posterior_C_cases,
    "confirmed (DRC)" => posterior_C_confirmed,
    "joint" => posterior_C_joint);

#md # ```@raw html
#md # </details>
#md # ```

streams_C_table #hide

# The grid below has replicates from the per-stream fits on the
# top row and the joint fit on the bottom row, comparable column-wise
# so it is easy to see what each per-stream fit constrains and how the
# joint combination shifts the predictives. The confirmed and
# samples-received columns come from the laboratory-pipeline fit, which
# conditions on both lab observations together.

#md # ```@raw html
#md # <details><summary>Per-stream vs joint posterior predictive grid</summary>
#md # ```

pp_exports_only = vec(Array(predict(
    exports_only_model(missing), chn_exports)[:exported_cases]));
pp_deaths_only = vec(Array(predict(
    deaths_only_model(missing), chn_deaths)[:total_deaths]));
pp_cases_only = vec(Array(predict(
    cases_only_model(missing), chn_cases)[:reported_cases]));
pp_exports_deaths_only = vec(sum.(predict(
    exports_deaths_only_model(fill(missing, length(obs.export_deaths_daily));
        pre_start_deaths = missing),
    chn_exports_deaths)[@varname(export_deaths_daily)]));
## Laboratory-pipeline fit: predict the confirmed Binomial and the
## received-count NegBinomial from the confirmed-only posterior for the
## individual row of the grid. Both are single-vintage vectors here.
pp_confirmed_only_chn = predict(
    confirmed_only_model(missing, obs.cumulative_tests_analysed),
    chn_confirmed);
pp_confirmed_only = vec(sum.(
    pp_confirmed_only_chn[@varname(confirmed_cases)]));
pp_tests_only = vec(sum.(
    pp_confirmed_only_chn[@varname(samples_received)]));

ppc_grid_fig = plot_posterior_predictive_grid(;
    individual = (; exports = pp_exports_only,
        exports_deaths = pp_exports_deaths_only,
        deaths = pp_deaths_only,
        cases = pp_cases_only,
        tests = pp_tests_only,
        confirmed = pp_confirmed_only),
    joint = (; exports = pp_exports,
        exports_deaths = pp_exports_deaths,
        deaths = pp_deaths,
        cases = pp_cases,
        tests = pp_tests,
        confirmed = pp_confirmed),
    observed = (; exports = obs.exported_cases,
        exports_deaths = obs.exports_deaths,
        deaths = obs.total_deaths,
        cases = obs.reported_cases,
        tests = obs.tests_received_history.values[end],
        confirmed = obs.confirmed_cases)
);

#md # ```@raw html
#md # </details>
#md # ```

ppc_grid_fig #hide

# Overlaid posterior densities of $C(T)$ from the five fits:

#md # ```@raw html
#md # <details><summary>Overlaid C_T density plot</summary>
#md # ```

## Clip the x-axis to keep the bulk of every density legible: the
## exports-deaths fit has a very heavy upper tail (95% reaching several
## thousand), which otherwise stretches the axis and compresses the
## other curves. The cap is the widest 95% upper across the remaining
## fits, so only the exports-deaths tail is truncated.
density_xmax = 1.1 * maximum(quantile(v, 0.95)
for v in (
    posterior_C_exports, posterior_C_deaths, posterior_C_cases,
    posterior_C_confirmed, posterior_C_joint))
cumulative_density_fig = plot_cumulative_cases(
    "exports (cases)" => posterior_C_exports,
    "exports (deaths)" => posterior_C_exports_deaths,
    "deaths (DRC)" => posterior_C_deaths,
    "cases (DRC)" => posterior_C_cases,
    "confirmed (DRC)" => posterior_C_confirmed,
    "joint" => posterior_C_joint;
    scenarios = [], xmax = density_xmax);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_density_fig #hide

# ### Comparison with McCabe et al.
#
# Our joint fit against the McCabe et al. estimates and our Method 2
# reproduction: point estimates with 95% intervals. We step through
# both report versions in turn — the 18 May version [mccabe2026](@cite)
# and the 20 May version [mccabe2026update](@cite) — and end with our
# joint fit to the current data.

#md # ```@raw html
#md # <details><summary>Fit our model to each report version's data, and run the Method 2 reproductions</summary>
#md # ```

## `obs_report` and `chn_joint_report` (the 18 May report's data) are
## loaded and fitted earlier, in the forecast-validation section.
obs_report_20may = load_observations(
    joinpath(pkgdir(BVDOutbreakSize), "data",
    "report-snapshot-20may.toml"));

report_20may_args = joint_obs(obs_report_20may)
chn_joint_report_20may = nuts_sample(
    bvd_joint(obs_report_20may.exported_cases,
    report_20may_args.deaths, report_20may_args.reported,
    report_20may_args.export_deaths; report_20may_args.kw...,
    growth = growth_for(obs_report_20may),
    first_export_detection_delta =
    obs_report_20may.first_export_detection_delta));
posterior_C_joint_report_20may = vec(Array(chn_joint_report_20may[:cumulative_cases]));

imperial_fixed = Turing.fix(
    imperial_only_model(missing, 88),       # exports missing → pure Method 2
    ## Pin the τ = 14 d scenario through the growth rate r = log(2)/14.
    (r = log(2) / 14, CFR = 0.30, α = 4.42, θ = 1/0.388,
        inv_sqrt_k = 1e-3)
)
chn_imperial = nuts_sample(imperial_fixed);
posterior_C_imperial = vec(Array(chn_imperial[:cumulative_cases]));

imperial_fixed_20may = Turing.fix(
    imperial_only_model(missing, 131),      # exports missing → pure Method 2
    ## Pin the τ = 14 d scenario through the growth rate r = log(2)/14.
    (r = log(2) / 14, CFR = 0.33, α = 4.42, θ = 1/0.388,
        inv_sqrt_k = 1e-3)
)
chn_imperial_20may = nuts_sample(imperial_fixed_20may);
posterior_C_imperial_20may = vec(Array(chn_imperial_20may[:cumulative_cases]));

#md # ```@raw html
#md # </details>
#md # ```

# The plot places each estimate of $C(T)$ on one axis: the central
# estimate as a point, the 95% interval as a bar. Rows are grouped by
# report version. For each version the first two rows are McCabe et
# al.'s published Method 1 and Method 2 headline scenarios with their
# reported 95% confidence intervals, followed by our Method 2
# reproduction and our joint fit to that version's data, both as 95%
# credible intervals. The final row is our joint fit to the current
# data.

#md # ```@raw html
#md # <details><summary>Build the comparison</summary>
#md # ```

## Our model rows use 95% equal-tailed credible intervals, matching the
## 95% confidence intervals McCabe et al. report for both methods.
function ci95(xs)
    (round(Int, quantile(xs, 0.5)),
        round(Int, quantile(xs, 0.025)),
        round(Int, quantile(xs, 0.975)))
end

joint_ci = ci95(posterior_C_joint)
joint_report_ci = ci95(posterior_C_joint_report)
joint_report_20may_ci = ci95(posterior_C_joint_report_20may)
imperial_ci = ci95(posterior_C_imperial)
imperial_20may_ci = ci95(posterior_C_imperial_20may)

comparison_rows = [
    ("18 May: McCabe Method 1 (Ituri, w=15 d)", 313, 39, 870),
    ("18 May: McCabe Method 2 (τ=14 d, CFR 30%)", 501, 402, 612),
    ("18 May: Our Method 2 reproduction", imperial_ci...),
    ("18 May: Our joint (report data)", joint_report_ci...),
    ("20 May: McCabe Method 1 (Ituri, w=15 d)", 313, 39, 870),
    ("20 May: McCabe Method 2 (τ=14 d, CFR 33%)", 678, 568, 800),
    ("20 May: Our Method 2 reproduction", imperial_20may_ci...),
    ("20 May: Our joint (report data)", joint_report_20may_ci...),
    ("Our joint (current data)", joint_ci...)
]

comparison_fig = plot_estimate_comparison(comparison_rows);

#md # ```@raw html
#md # </details>
#md # ```

comparison_fig #hide

# Fit diagnostics for the two report-data joint fits and the two Method
# 2 reproductions.

#md # ```@raw html
#md # <details><summary>Fit diagnostics</summary>
#md # ```

diagnostics_table( #hide
    "joint (18 May report)" => chn_joint_report, #hide
    "joint (20 May report)" => chn_joint_report_20may, #hide
    "Method 2 reproduction (18 May)" => chn_imperial, #hide
    "Method 2 reproduction (20 May)" => chn_imperial_20may) #hide

#md # ```@raw html
#md # </details>
#md # ```

# The same comparison as a table, with a column for the report version:

#md # ```@raw html
#md # <details><summary>Comparison table</summary>
#md # ```

comparison_version = [
    "18 May", "18 May", "18 May", "18 May",
    "20 May", "20 May", "20 May", "20 May",
    "current"
]

main_comparison = DataFrame(
    "Report version" => comparison_version,
    "Source" => [r[1] for r in comparison_rows],
    "Central estimate" => [r[2] for r in comparison_rows],
    "Lower 95%" => [r[3] for r in comparison_rows],
    "Upper 95%" => [r[4] for r in comparison_rows]
);

#md # ```@raw html
#md # </details>
#md # ```

main_comparison #hide

# Both Method 2 reproductions land on McCabe et al.'s reported Method 2
# central estimates: the 18 May reproduction against their reported
# $501$ ($88$ deaths, CFR $30\%$) and the 20 May reproduction against
# their reported $678$ ($131$ deaths, CFR $33\%$).

#md # ```@raw html
#md # <details><summary>Reproductions vs McCabe et al. Method 2</summary>
#md # ```

imperial_sense_check = let
    rep, lo, hi = imperial_ci
    rep2, lo2, hi2 = imperial_20may_ci
    delta = round(100 * (rep - 501) / 501; digits = 1)
    delta2 = round(100 * (rep2 - 678) / 678; digits = 1)
    Markdown.parse("""
    18 May reproduction: **$(rep) cases** (95% CrI $(lo)–$(hi)) against
    McCabe et al.'s reported **501** — a difference of $(delta)%.

    20 May reproduction: **$(rep2) cases** (95% CrI $(lo2)–$(hi2))
    against McCabe et al.'s reported **678** — a difference of
    $(delta2)%.
    """)
end;

#md # ```@raw html
#md # </details>
#md # ```

imperial_sense_check #hide

# Joint posterior coverage of all 15 published McCabe et al. scenarios
# — for each scenario, the narrowest joint credible interval that
# contains it:

#md # ```@raw html
#md # <details><summary>Joint coverage table</summary>
#md # ```

coverage_table = comparison_table(posterior_C_joint);

#md # ```@raw html
#md # </details>
#md # ```

coverage_table #hide

# The joint $C(T)$ density with the 15 published scenario point
# estimates overlaid as faint dashed rules, for our current-data fit
# and our fits to each report version's data:

#md # ```@raw html
#md # <details><summary>Joint C_T density with published scenarios</summary>
#md # ```

imperial_density_fig = plot_cumulative_cases(
    "joint (current data)" => posterior_C_joint,
    "joint (18 May report)" => posterior_C_joint_report,
    "joint (20 May report)" => posterior_C_joint_report_20may);

#md # ```@raw html
#md # </details>
#md # ```

imperial_density_fig #hide

# ## Saving results
#
# The tables above are written to an `output/` directory at the repo
# root so they can be archived and shared. On every push to `main` a
# GitHub Actions workflow regenerates these files and publishes them
# as a GitHub Release, downloadable from the repository's releases
# page (<https://github.com/epiforecasts/BVDOutbreakSize/releases>).
# The release bundles the four summary tables, a thinned set of
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

CSV.write(joinpath(output_dir, "posterior_summary.csv"), joint_summary)
CSV.write(joinpath(output_dir, "cumulative_cases_by_stream.csv"),
    streams_C_table)
CSV.write(joinpath(output_dir, "imperial_comparison.csv"), main_comparison)
CSV.write(joinpath(output_dir, "scenario_coverage.csv"), coverage_table)
CSV.write(joinpath(output_dir, "forecast_validation.csv"),
    forecast_validation_table)

## Copy the input data so the release records what produced these
## results.
cp(joinpath(pkgdir(BVDOutbreakSize), "data", "observations.toml"),
    joinpath(output_dir, "observations.toml"); force = true)

## Thinned posterior draws of the key joint parameters (every 10th
## draw) so downstream users can recompute their own summaries.
posterior_draws = DataFrame(
    τ = vec(Array(chn_joint[:τ])),
    r = vec(Array(chn_joint[:r])),
    m = vec(Array(chn_joint[:m])),
    T = vec(Array(chn_joint[:T])),
    CFR = vec(Array(chn_joint[:CFR])),
    p_drc = vec(Array(chn_joint[:p_drc])),
    p_uganda = vec(Array(chn_joint[:p_uganda])),
    cumulative_cases = vec(Array(chn_joint[:cumulative_cases]))
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
