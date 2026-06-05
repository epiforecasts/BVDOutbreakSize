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
# - *Discrete-time renewal model.* The whole model runs on a daily grid.
#   Infections follow the discrete renewal equation $I_t = R_t \sum_{s
#   \ge 1} I_{t-s} g_s$, where $g$ is the discretised generation-interval
#   PMF, and every delay is applied as a discrete convolution. McCabe et
#   al. [mccabe2026](@cite) use continuous-time closed forms.
# - *Time-varying reproduction number.* $R_t$ follows a weekly Gaussian
#   random walk on the log scale, interpolated within weeks, with a
#   logistic intervention ramp at the first WHO situation report. McCabe
#   et al. use one constant exponential growth rate.
# - *Joint posterior, not 15 scenario estimates.* The reproduction
#   number, case-fatality ratio, all delays, traveller volume and
#   surveillance dispersion have priors and are sampled together. McCabe
#   et al. fix each and sweep.
# - *Delays sampled from priors and discretised.* Generation interval,
#   incubation period, onset-to-death, onset-to-report,
#   onset-to-confirmation and onset-to-detection-abroad each get a prior
#   centred on published Ebola estimates, discretised with double
#   interval censoring [charniga2024](@cite). No delay is fixed.
# - *Euler–Lotka seeding.* The seeding window grows exponentially at the
#   rate implied by the initial reproduction number and generation
#   interval, so infections start smoothly rather than from a single
#   seed.
# - *Per-vintage time-series fitting.* The DRC streams (suspected cases,
#   confirmed cases, deaths) are fitted as cumulative series of
#   between-vintage increments across successive sitreps, which sharpens
#   $R_t$. McCabe et al. condition on a single cumulative total.
# - *Ascertainment extension.* The DRC and Uganda reporting fractions
#   share a logit-scale hyperprior, giving a joint posterior over
#   ascertainment alongside outbreak size. Not in McCabe et al.
# - *Comparison against published scenarios.* The model's expected
#   cumulative confirmed cases are compared against all 15 McCabe et al.
#   scenario estimates with a coverage table, while $C_T$ (the latent
#   infection count, summed from the renewal trajectory) is reported
#   separately.
# - *No-onward-transmission counterfactual and one-week-ahead
#   forecasts.* Future expected deaths from infections already seeded,
#   and a posterior-predictive projection of each stream. Neither is in
#   McCabe et al.
#
# ## Limitations
#
# - *Fitted only to aggregate reported counts.* The data are a handful
#   of summary figures from press and situation reports: suspected cases
#   and deaths in the DRC, and cases (with one death) in Uganda. There
#   is no line list, and no information on case definitions, testing
#   capacity or reporting completeness. Every estimate is a model-based
#   extrapolation under strong assumptions, not a measurement.
# - *Prior-driven where data is scarce.* The sitrep trajectory pins down
#   $R_t$, but a few totals say little about the delays, surveillance
#   dispersion or reporting fraction on their own, so those posteriors
#   track their priors.
# - *Per-sitrep increments are not clean new incidence.* Later sitreps
#   likely backfill earlier cases and add newly-reporting health zones,
#   and ascertainment probably rose over the window.
# - *Inherits McCabe et al.'s epidemiological assumptions.* A single
#   zoonotic seed, an assumed generation interval, no spatial structure
#   beyond the Ituri / Nord Kivu split, and no depletion of
#   susceptibles. The onset-to-death delay is anchored on Isiro 2012 and
#   the genetic seeding bound on an external clock rate, neither
#   propagating cross-outbreak or clock uncertainty.
# - *Intervention ramp is weakly identified.* With only a few sitreps
#   straddling it, the ramp effect and the pre-ramp reproduction number
#   are not well separated.
# - *Ascertainment partially pooled, not separately identified.* The
#   DRC and Uganda reporting fractions share a hyperprior; with a
#   handful of exports the Uganda fraction leans on the DRC side.
# - *Streams share one case pool.* They are fitted as conditionally
#   independent given latent incidence but observe overlapping people,
#   which can understate uncertainty. We have not checked whether the
#   streams imply conflicting outbreak sizes.
# - *LLM-driven reimplementation.* The model code, priors and analysis
#   were drafted by a language model from the McCabe et al.
#   [mccabe2026](@cite) report and the companion delay reanalysis, then
#   reviewed and revised. Not independently replicated against the
#   authors' code.
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
# The analysis uses a handful of aggregate counts. The DRC suspected
# cases, suspected deaths and laboratory-confirmed cases are the
# national cumulative totals from the INSP situation reports
# [insp_sitrep_2026](@cite), read from the report PDFs (archived by
# INRB-UMIE [inrb_umie_2026](@cite)). We draw straight from the sitreps
# rather than the published per-zone CSVs because the regional (health
# zone) breakdown is inconsistent with the national totals: the zone
# sums omit cases not yet attributed to a zone, understating the count.
# The Uganda export-case counts and deaths come from WHO Disease Outbreak
# News DON602 [who_don_2026_602](@cite); the cross-border traveller
# volume and source population from McCabe et al. [mccabe2026](@cite).
# The first table lists each figure as of the cut-off; the source
# population is fixed and the traveller volume is given a Normal prior
# around the McCabe et al. figure. The three DRC streams are
# additionally resolved by sitrep vintage and fitted as between-vintage
# increments, shown in the second table.

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
        "confirmed_cases",
        "daily_outbound_travellers (prior mean)",
        "daily_outbound_travellers_sd (prior SD)",
        "source_population"
    ],
    value = [
        obs.exported_cases,
        obs.exports_deaths,
        obs.total_deaths,
        obs.reported_cases,
        obs.confirmed_cases,
        ITURI_DAILY_TRAVEL,
        ITURI_DAILY_TRAVEL_SD,
        ITURI_POPULATION
    ]
);

#md # ```@raw html
#md # </details>
#md # ```

observations_table #hide

# The per-vintage cumulative history of the three DRC sitrep streams,
# the national totals at each INSP situation-report date. The joint
# model fits the between-vintage increments of these series (a single
# vintage reduces to the cut-off total). See `data/observations.toml`
# for the per-stream sources.

#md # ```@raw html
#md # <details><summary>Building the per-vintage time-series table</summary>
#md # ```

vintage_table = let
    dh = obs.deaths_history
    rh = obs.reported_history
    ch = obs.confirmed_history
    nm = maximum([length(dh.days), length(rh.days), length(ch.days)])
    _pad(v) = vcat(v, fill(missing, nm - length(v)))
    DataFrame(
        deaths_day = _pad(dh.days),
        deaths_count = _pad(dh.counts),
        reported_day = _pad(rh.days),
        reported_count = _pad(rh.counts),
        confirmed_day = _pad(ch.days),
        confirmed_count = _pad(ch.counts)
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
# infections directly: each available data stream observes a thinned,
# delayed or transformed view of the same latent incidence.
#
# The model is assembled from small reusable Turing [ge2018turing](@cite)
# submodels. Building-block submodels own the maths and priors for one
# parameter family. Observation submodels assemble those blocks,
# apply the convolution chain and the likelihood, and tie one stream to
# the latent incidence. Composers combine the observation submodels for
# the per-stream and joint fits.
#
# The table below shows which parameters feed each observation submodel.
# The *received* column is the received-specimen volume, the laboratory
# stream fitted as a count; the *confirmed* positives are scored as a
# Binomial of the observed analysed denominator with a free per-window
# positivity, so they are deliberately decoupled from the infection
# process (the outbreak size is pinned by the deaths and exports). The
# *conf. deaths* column thins the suspected deaths by the case composition:
#
# | Parameter | Exports | Deaths | Cases | Received | Confirmed | Conf. deaths | Export deaths |
# |---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
# | Reproduction number $R_t$ | ● | ● | ● | ● |  | ● | ● |
# | Generation interval | ● | ● | ● | ● |  | ● | ● |
# | Incubation period | ● | ● | ● | ● |  | ● | ● |
# | Seed $I_0$ | ● | ● | ● | ● |  | ● | ● |
# | Onset-to-death delay |  | ● |  |  |  |  | ● |
# | Case-fatality ratio |  | ● |  |  |  |  | ● |
# | Onset-to-report delay |  |  | ● | ● |  | ● |  |
# | Receipt delay |  |  |  | ● |  |  |  |
# | Onset-to-detection delay | ● |  |  |  |  |  |  |
# | Laboratory positivity (per-window RE) |  |  |  |  | ● |  |  |
# | Confirmed-death enrichment $m_{\text{death}}$ |  |  |  |  |  | ● |  |
# | Testing fraction $\tau_{\text{test}}$ |  |  |  | ● |  |  |  |
# | Background rate $\lambda_{\text{bg}}$ |  |  | ● | ● |  | ● |  |
# | Surveillance dispersion |  | ● | ● | ● |  |  |  |
# | Ascertainment | ● |  | ● | ● |  | ● | ● |
# | Traveller volume | ● |  |  |  |  |  | ● |

#md # ```@setup main
#md # using BVDOutbreakSize, CodeTracking, Revise
#md # ```

# #### Building-block submodels
#
# The implementation uses Mooncake [mooncake_jl](@cite) reverse-mode
# automatic differentiation, CensoredDistributions for delay
# discretisation, FlexiChains for chain handling, and PairPlots
# [pairplots_jl](@cite) with AlgebraOfGraphics
# [danisch2021makie](@cite) for the figures.
# Each building-block submodel owns the maths and priors for one parameter
# family; its source is shown in the collapsible block beneath its prose.
#
# ##### Reproduction number — weekly random walk with intervention ramp
#
# Knots are placed at weekly intervals (day 1, 8, 15, …, $n$).
# The knot values follow a Gaussian random walk on the log scale in
# non-centred cumulative-sum form:
#
# ```math
# \log R_0 \sim \mathrm{Normal}(\log 1.3,\ 0.4), \qquad
# \sigma_{\text{rw}} \sim \mathrm{Normal}^{+}(0,\ 0.2), \tag{2}
# ```
#
# ```math
# \log R_k = \log R_0 + \sigma_{\text{rw}}
#            \sum_{j=1}^{k} z_j, \quad
# z_j \sim \mathrm{Normal}(0, 1). \tag{3}
# ```
#
# Daily $\log R_t$ is the piecewise-linear interpolation between knots.
# An intervention at the first WHO situation report adds a sampled
# effect shaped by a logistic ramp over seven days:
#
# ```math
# \log R_t \mathrel{+}= \delta \cdot
#     \mathrm{logistic}\!\left(\frac{t - t_{\text{bp}}}{7}\right),
# \qquad
# \delta \sim \mathrm{Normal}(0,\ 0.5). \tag{4}
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

# ##### Generation interval and incubation period
#
# The generation-interval PMF $g$ is sampled from a prior centred on
# the Ebola virus disease serial interval as a generation-time proxy
# (mean 15.3 d, SD 9.3 d; WHO Ebola Response Team 2014, NEJM), then
# discretised with double interval censoring [charniga2024](@cite). The
# prior is centred on those published moments, with an assumed
# weakly-informative spread (the SDs on the mean and SD priors are our
# choice, not from the source). The lag-0 bin is dropped and the
# remainder renormalised so an infectee is always strictly later than
# its infector:
#
# ```math
# \mu_g \sim \mathrm{Normal}^{+}(15.3,\ 3.0), \qquad
# \sigma_g \sim \mathrm{Normal}^{+}(9.3,\ 2.0). \tag{5}
# ```
#
# The incubation period is similarly discretised with a prior centred on
# the Bundibugyo virus incubation estimate from the 2007 Uganda outbreak
# (mean 6.3 d, 95% CI 5.2-7.3, $n = 24$; [macneil2010](@cite)). The
# line-list reanalysis cannot fit incubation, as the line list has no
# exposure dates, so it recommends this estimate instead. The mean prior
# reproduces MacNeil et al.'s 95% CI; the spread prior is a
# weakly-informative modelling choice, as they report no interval on the
# SD:
#
# ```math
# \mu_{\text{inc}} \sim \mathrm{Normal}^{+}(6.3,\ 0.54), \qquad
# \sigma_{\text{inc}} \sim \mathrm{Normal}^{+}(3.5,\ 0.8). \tag{6}
# ```
#
# All LogNormal parameters are recovered by moment-matching from the
# sampled mean and SD. Every delay in the model shares the same
# double-interval-censored discretisation [charniga2024](@cite); the
# generation interval wraps it to drop the lag-0 bin.

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

# ##### Seeding — Euler–Lotka implied growth
#
# The seeding window (length $L$ = generation-interval support) is
# filled with exponential growth at the rate $r_0$ implied by the
# initial reproduction number $R_0 = R_t[1]$ and generation interval
# $g$ via the Euler–Lotka relation
#
# ```math
# R_0 \sum_{s=1}^{L} g_s \, e^{-r_0 s} = 1, \tag{7}
# ```
#
# solved by Newton iteration. The seed count on the last seeding day
# has a prior centred on a single introduction:
#
# ```math
# I_0 \sim \mathrm{Normal}^{+}(1.0,\ 1.0). \tag{8}
# ```

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

# ##### Generating infection process and onset staging
#
# Infections $I_t$ are produced by equation (1) for $t > L$, with
# $I_{1:L}$ from the seeding window. Cumulative infections are the
# running sum $C_t = \sum_{s=1}^{t} I_s$; the cut-off cumulative is
# $C_T = C_n$. The renewal model's current growth rate $r$ and doubling
# time are derived from the day-over-day log-ratio at the cut-off:
#
# ```math
# r = \log I_n - \log I_{n-1}, \qquad
# \tau_{1/2} = \log(2) / r. \tag{9}
# ```
#
# Infections are convolved with the incubation PMF to produce daily
# symptom-onset incidence, which every downstream stream then consumes.
# The incubation delay is an injected delay submodel, defaulting to the
# Bundibugyo incubation prior of equation (6).

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

# ##### Onset-to-death delay
#
# The McCabe et al. report uses the point estimate of
# [rosello2015](@cite). We instead anchor the onset-to-death delay on the
# companion Bayesian reanalysis of the same Isiro 2012 BDBV line list
# [bdbv_linelist_analysis_2026](@cite), which re-estimates the delay with
# uncertainty. The renewal samples the delay by its mean and SD rather
# than a Gamma shape and scale. The prior means are centred on the
# reanalysis' posterior mean (11.2 d) and SD (5.4 d), with an assumed
# weakly-informative spread on each so the fit reproduces the reanalysis
# uncertainty rather than collapsing onto a single point estimate:
#
# ```math
# \mu_d \sim \mathrm{Normal}^{+}(11.2,\ 2.0), \qquad
# \sigma_d \sim \mathrm{Normal}^{+}(5.4,\ 1.5). \tag{10}
# ```
#
# The sampled mean and SD are moment-matched to a LogNormal and
# discretised with double interval censoring [charniga2024](@cite), as
# for every delay above. The delay estimation in that reanalysis follows
# the same recommendations. The submodel source is shown with the deaths
# observation submodel below, where the delay is injected.
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
# non-healthcare-worker (non-HCW) confirmed cases. The prior on the
# case-fatality ratio is
#
# ```math
# \mathrm{CFR} \sim \mathrm{Beta}(6.6,\ 13.4), \tag{11}
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

# ##### Daily traveller volume
#
# The number of people crossing from the source area to Uganda each day
# sets the travel rate in the exports likelihood. We treat it as an
# estimated quantity rather than a fixed input. McCabe et al. Table 3
# records mean weekly passenger counts across seven points of entry; the
# Ituri-side daily total of $1871$ is a sample mean across roughly
# $15$-$21$ point-of-entry-weeks. We use a Normal prior centred on
# $1871$ with SD $200$ ($\approx 10\%$ CV), truncated at zero, covering
# point-of-entry variation and the sitrep sampling uncertainty; the
# source population is kept fixed (census).
#
# ```math
# N_{\text{travel}} \sim \mathrm{Normal}^{+}(1871,\ 200). \tag{12}
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

# ##### Surveillance dispersion
#
# Passive-surveillance counts are modelled with negative-binomial
# observation error with a single shared dispersion $k$ for the DRC
# suspected deaths, reported cases and confirmed cases. Following Stan
# prior-choice recommendations [stan_prior_choice](@cite), the
# dispersion is sampled on the $1/\sqrt{k}$ scale:
#
# ```math
# 1/\sqrt{k} \sim \mathrm{Normal}^{+}(0.6,\ 0.2). \tag{13}
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

# ##### Ascertainment — partial pooling between DRC and Uganda
#
# Two surveillance systems detect cases: DRC passive community
# surveillance (the reported suspected-case count) and Uganda's
# point-of-entry / hospital surveillance (the exported-case count). Each
# captures only a fraction of the true cases passing through it, and each
# fraction is informed by essentially a single aggregate data point. The
# two ascertainment fractions $p_{\text{DRC}}$ and $p_{\text{Uganda}}$
# share a logit-scale hyperprior with mean $\mu$ and pooling strength
# $\tau$, centred on an assumed reporting fraction of $25\%$:
#
# ```math
# \mu \sim \mathrm{Normal}(\mathrm{logit}(0.25),\ 1),
# \qquad
# \tau \sim \mathrm{Normal}^{+}(0,\ 0.5), \tag{14}
# ```
#
# ```math
# \mathrm{logit}(p_{\text{DRC}}) \sim \mathrm{Normal}(\mu,\ \tau),
# \qquad
# \mathrm{logit}(p_{\text{Uganda}}) \sim \mathrm{Normal}(\mu,\ \tau). \tag{15}
# ```
#
# The prior is sampled in non-centred form to avoid the funnel geometry.
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

# ##### Genetic seeding bound
#
# A BEAST time tree of the first ten sequenced genomes
# [virological2026](@cite) places the TMRCA, the age of the oldest
# internal node of the tree, at a mean of 25 March 2026. The temporal
# sampling range is too short to estimate the molecular clock, so it is
# fixed to the $1.2\times10^{-3}$ substitutions/site/year rate of the
# 2013-2016 West African Ebola epidemic [holmes2016](@cite).
# The TMRCA is a lower bound on the seeding time $T$. Adding sequences,
# or more geographically representative ones, can only push it earlier,
# never later, as the sampled tree is almost entirely from Bunia.
# Using the genetic TMRCA as a one-sided seeding bound rather than a
# point estimate follows a suggestion of N. Ferguson [ferguson2026](@cite).
#
# We treat the TMRCA as a right-censored, noisy reading of the seeding
# time. Writing $g = t_{\text{cut}} - t_{\text{TMRCA}}$ for the seeding
# age implied by the reported TMRCA date, so that $g$ tracks the cut-off
# rather than a fixed offset, the bound contributes the probability that
# a $\mathrm{Normal}(T, \sigma)$ reading falls at or beyond $g$,
#
# ```math
# p_\text{gen}(T) = \Pr[\mathrm{Normal}(T, \sigma) \ge g]
#   = \Phi\!\left(\frac{T - g}{\sigma}\right),
# \qquad \sigma = 15\ \text{d}. \tag{16}
# ```
#
# The bound is one-sided. It penalises outbreak ages younger than the
# TMRCA but leaves $T$ free above it, with the clock fixed and no
# propagation of cross-outbreak or clock uncertainty.

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

# ##### Laboratory priors — testing fraction, background, positivity and enrichment
#
# The laboratory pipeline adds building-block priors. The testing fraction
# $\tau_{\text{test}}$ is the share of suspected cases routed to the lab,
# with a $\mathrm{Beta}(5, 2)$ prior (mean $\approx 0.71$); the non-BVD
# background rate $\lambda_{\text{bg}}$ is the expected non-BVD suspected
# reports per day, with an informative half-normal prior identified from
# the suspected/confirmed contrast (a diffuse prior lets the background
# absorb the whole suspected stream). The per-window laboratory positivity
# is a partially-pooled logit-normal random effect with baseline
# $q_\mu$ centred on the cut-off cumulative positivity ($210 / 755
# \approx 0.28$) and pooling SD $\sigma_q$. The confirmed-death enrichment
# $m_{\text{death}}$ scales the case composition on the odds scale. The
# report-to-laboratory receipt delay is centred on a short turnaround with
# a heavy right tail for specimen shipment, with no per-sample anchor.
#
# ```math
# \tau_{\text{test}} \sim \mathrm{Beta}(5,\ 2), \qquad
# \lambda_{\text{bg}} \sim \mathrm{Normal}^{+}(0,\ 1)\ \text{per day},
# ```
#
# ```math
# q_\mu \sim \mathrm{Normal}(\operatorname{logit}(0.28),\ 0.7), \qquad
# \sigma_q \sim \mathrm{Normal}^{+}(0,\ 1), \qquad
# m_{\text{death}} \sim \mathrm{LogNormal}(0,\ 1). \tag{17}
# ```

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

# #### Observation submodels
#
# Each observation submodel takes the shared daily onset incidence,
# convolves it with a sampled onset-to-event delay PMF, scales by the
# relevant ascertainment, CFR or positivity factor, and reads the
# modelled cumulative count off the daily series at each vintage day.
# Likelihoods score the between-vintage increments.
#
# ##### Exports — geographic spread
#
# Each onset-day case travels to Uganda with daily rate
# $q = N_{\text{travel}} / N_{\text{pop}}$. The expected detected
# exports sum the onset-to-detection daily series:
#
# ```math
# \mu_e = p_{\text{Uganda}} \cdot q
#         \cdot \sum_{t=1}^{n} \mathrm{onsets}_t \cdot f_{\text{det}}(n - t),
# \qquad
# Y_{\text{exports}} \sim \mathrm{Poisson}(\mu_e). \tag{18}
# ```
#
# The onset-to-detection delay is centred on the Ebola
# onset-to-hospitalisation delay (mean 5.0 d, SD 4.7 d; WHO Ebola
# Response Team 2014, NEJM), again with an assumed weakly-informative
# spread rather than uncertainty taken from the source:
#
# ```math
# \mu_{\text{det}} \sim \mathrm{Normal}^{+}(5.0,\ 2.0), \qquad
# \sigma_{\text{det}} \sim \mathrm{Normal}^{+}(4.7,\ 1.5). \tag{19}
# ```
#
# The exports stream uses this onset-to-detection delay in place of a
# fixed detection window: rather than a single mean window over which a
# case is detectable abroad, exports are the export-onset series convolved
# with a sampled onset-to-detection delay PMF, so the timing of detection
# is carried explicitly.

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

# ##### Deaths — back-calculation
#
# Expected cumulative deaths at the cut-off are the CFR-weighted
# discrete convolution of onsets with the onset-to-death PMF (equation
# (10)). The model conditions on the between-vintage increment at each
# sitrep date with a NegBinomial likelihood sharing $k$. The death
# history ends at the cut-off, so the cut-off total is the final
# increment's cumulative and is not scored separately. Writing
# $\mu_i = \mathrm{CFR} \cdot \mathrm{conv}(\mathrm{onsets},\, f_d)$
# cumulated to sitrep day $d_i$, the increment at vintage $i$ is
#
# ```math
# Y_{\text{deaths},i} - Y_{\text{deaths},i-1} \sim
#     \mathrm{NegBinomial}(\mu_i - \mu_{i-1},\ k),
#     \qquad \mu_0 = 0. \tag{20}
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

# ##### Reported cases
#
# Reported suspected cases are the sum of two parts.
# A BVD-driven component, onsets convolved with the onset-to-report
# delay (mean 4.5 d, SD 3.6 d) and scaled by the DRC ascertainment
# $p_{\text{DRC}}$, plus an additive non-BVD background accruing at
# $\lambda_{\text{bg}}$ per day, so a suspected case need not be a true
# BVD infection. The increments are
# scored per vintage with a NegBinomial sharing $k$:
#
# ```math
# Y_{\text{cases}} \sim \mathrm{NegBinomial}\!\bigl(p_{\text{DRC}}
#     \cdot \mathrm{conv}(\mathrm{onsets},\, f_{\text{rep}})[n]
#     + \lambda_{\text{bg}}\, n,\ k\bigr). \tag{21}
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

# ##### Laboratory pipeline: received specimens and confirmed cases
#
# The laboratory pipeline fits two streams. The received-specimen volume
# is the suspected daily pipeline ($p_{\text{DRC}}$-scaled BVD
# onset-to-report signal plus the non-BVD background $\lambda_{\text{bg}}$)
# carried through a receipt delay $f_{\text{rec}}$ and thinned by the
# testing fraction $\tau_{\text{test}}$, scored per vintage with the
# shared $k$:
#
# ```math
# Y_{\text{received}} \sim \mathrm{NegBinomial}\!\bigl(\tau_{\text{test}}
#     \cdot \mathrm{conv}(p_{\text{DRC}}\, \mathrm{BVD}_{\text{rep}}
#     + \lambda_{\text{bg}},\, f_{\text{rec}}),\ k\bigr). \tag{22}
# ```
#
# The confirmed positives are scored as a Binomial of the *observed*
# specimens-analysed denominator $A_v$ in each laboratory window, with a
# partially-pooled per-window positivity $p_{\text{pos},v} =
# \operatorname{logistic}(q_\mu + \sigma_q z_v)$. Conditioning on the
# observed denominator rather than passing the confirmed count through the
# multiplicative ascertainment ridge $p_{\text{DRC}}\, s_{\text{test}}\,
# \tau_{\text{test}}$ decouples the confirmed counts from the outbreak
# size, which the deaths and exports streams pin instead:
#
# ```math
# C_v \sim \mathrm{Binomial}(A_v,\ p_{\text{pos},v}). \tag{23}
# ```
#
# The early confirmed vintages (18-23 May) have no per-vintage analysed
# denominator, so they are scored as NegativeBinomial counts against the
# modelled laboratory volume $V_v$ with the same partially-pooled
# positivity, extending the use of the confirmed data to where no
# laboratory denominator is observed:
#
# ```math
# C_v^{\text{early}} \sim \mathrm{NegBinomial}(p_{\text{pos},v}\, V_v,\ k).
# \tag{23b}
# ```

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
# Confirmed deaths are a thinning of the observed suspected deaths. The
# confirmation probability is the suspected-case BVD composition
# $q_{\text{susp}} = p_{\text{DRC}}\,\mathrm{BVD} /
# (p_{\text{DRC}}\,\mathrm{BVD} + \lambda_{\text{bg}})$ enriched on the
# odds scale by $m_{\text{death}}$, so it stays in $(0, 1)$ without a hard
# clamp and a confirmed-death observation informs the background and
# ascertainment:
#
# ```math
# Y_{\text{conf-deaths}} \sim \mathrm{Binomial}\bigl(D_{\text{susp}},\
#     \operatorname{logistic}(\operatorname{logit}(q_{\text{susp}})
#     + \log m_{\text{death}})\bigr). \tag{24}
# ```

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

# ##### Deaths among exports
#
# The expected deaths among detected exports reuse the export-onset
# series from the exports submodel, convolving it with the onset-to-death
# PMF and scaling by the CFR. A Poisson likelihood is used because the
# Uganda death count is small:
#
# ```math
# Y_{\text{exp-deaths}} \sim \mathrm{Poisson}(\mathrm{CFR}
#     \cdot \mathrm{conv}(\mathrm{export\_onsets},\, f_d)[n]). \tag{25}
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

# #### Composers
#
# Composers combine building blocks into the full model for each
# analysis. Each observation stream argument may be `missing` to drop
# it, so the same composer structure generates prior- and
# posterior-predictive draws.
#
# The joint composer runs the generating infection process once and
# routes the shared onsets into all five observation submodels. It
# samples a single dispersion $k$ and the pooled ascertainment fractions,
# threading $p_{\text{DRC}}$ to the cases and laboratory likelihoods and
# $p_{\text{Uganda}}$ to the two Uganda-side likelihoods.
#
# We write single-stream composers for each of the five count-based
# streams: exports-only, deaths-only, cases-only, confirmed-only and
# exports-deaths-only.

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
    [:C_T, :R_T, :r, :doubling_time, :T, :CFR, :k,
        :p_drc, :p_uganda]);

#md # ```@raw html
#md # </details>
#md # ```

prior_pair_fig #hide

# #### Fitting the models
#
# NUTS [hoffman2014nuts](@cite) with Mooncake [mooncake_jl](@cite)
# reverse-mode automatic differentiation, four chains, 1000
# post-warmup draws each, with a target acceptance probability of 0.95.
# Chains initialise from the prior to keep the sampler away from the
# boundary of the renewal recursion. We fit the joint model and the
# five single-stream models so the per-stream posteriors over $C_T$ can
# be compared with the joint.

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
    breakpoint = _BREAKPOINT,
    background_re = true,
    confirmed_positivity_link = :composition,
    genetic = genetic_seeding_model,
    tmrca_days = obs.tmrca_days));

chn_exports = nuts_sample(
    exports_only_model(obs.n, obs.exported_cases;
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

# ## Results
#
# ### Summary
#
# The question that matters for the response is how many people have
# already been infected.
# The confirmed counts capture only part of the outbreak, and planning
# for beds, contacts and vaccine needs depends on the true total.
# The numbers below are our current best estimate of that total, from
# the joint posterior and refreshed on every build.
# Each headline number is given as equal-tailed 30%, 60% and 90%
# credible intervals, the same intervals used in the tables below.

#md # ```@raw html
#md # <details><summary>Compute the headline ranges</summary>
#md # ```

summary_ranges = let
    C = posterior_C_joint
    Td = vec(Array(chn_joint[:T]))
    rd = vec(Array(chn_joint[:r]))
    dt = vec(Array(chn_joint[:doubling_time]))
    sC = posterior_summary(C)
    sT = posterior_summary(Td)
    sr = posterior_summary(rd)
    sdt = posterior_summary(dt)

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
    - **Cumulative infections \$C_T\$:** the posterior is $(ints_i(sC))
      infections, the latent pool behind every stream.
    - Against the $(obs.confirmed_cases) laboratory-confirmed cases by the
      cut-off that is roughly $(f_lo)–$(f_hi)× as many infections, so
      confirmed cases capture only a small share of the outbreak.
    - **Confirmed-case fit:** the model expects $(ints_i(ec)) confirmed
      cases by the cut-off, against $(obs.confirmed_cases) observed.
    - **Time since seeding:** the posterior is $(ints_i(sT)) days, placing
      the start of sustained transmission at $(ints_d(sT)).
    - **Growth rate and doubling time:** the current growth rate is
      $(ints_f(sr, 3)) per day.
      The implied doubling time is $(ints_f(sdt, 1)) days.
    """)
end;

#md # ```@raw html
#md # </details>
#md # ```

summary_ranges #hide

# **Why our estimate may differ from McCabe et al.**
# Our estimate fits all streams jointly, samples the nuisance parameters
# that McCabe et al. vary in scenario sweeps, and uses a time-varying
# reproduction number constrained by the sitrep trajectory.
# See [what we do differently](#What-we-do-differently-from-McCabe-et-al.),
# the [comparison with McCabe et al.](#Comparison-with-McCabe-et-al.) and
# the [limitations](#Limitations) for the detail behind this.

# ### Joint model estimates
#
# Our main result is the current cumulative case load, reported and
# unreported, at the report date.
# It is the joint posterior over the cumulative case count $C_T$, fitting
# all five data streams together: cases exported to Uganda, suspected
# deaths in the DRC, reported cases in the DRC (with an ascertainment
# component), laboratory-confirmed cases in the DRC, and deaths among
# exported cases in Uganda.
# We show it first as a credible-interval table and then as a posterior
# density.

#md # ```@raw html
#md # <details><summary>Cumulative case count summary table</summary>
#md # ```

cumulative_cases_summary = summary_table(
    chn_joint, [:C_T]; digits = 0);

#md # ```@raw html
#md # </details>
#md # ```

cumulative_cases_summary #hide

#md # ```@raw html
#md # <details><summary>Cumulative case count density</summary>
#md # ```

## The joint C_T posterior has a long upper tail (the slow-growth,
## late-seeding draws), so clip the x-axis to the 97.5th percentile to
## keep the bulk of the density legible rather than rendering the full
## heavy tail.
joint_density_fig = plot_cumulative_cases(
    "joint" => posterior_C_joint; scenarios = [],
    xmax = quantile(posterior_C_joint, 0.975));

#md # ```@raw html
#md # </details>
#md # ```

joint_density_fig #hide

# The cumulative case count $C_T$ is set by the reproduction number
# trajectory and the seeding time $T$. Read as a calendar date, $T$
# places the start of sustained transmission at the cut-off date minus
# $T$ days. The left panel below shows the posterior for that start
# date; the right panel shows the joint $(T, \tau_{1/2})$ posterior.

#md # ```@raw html
#md # <details><summary>Outbreak start date and seeding-time posterior</summary>
#md # ```

start_date_fig = plot_start_date_pair(chn_joint;
    as_of_date = string(obs.cutoff));

#md # ```@raw html
#md # </details>
#md # ```

start_date_fig #hide

# The full posterior summary table reports equal-tailed 30%, 60% and
# 90% credible intervals on the key joint-fit parameters.

#md # ```@raw html
#md # <details><summary>Joint posterior summary table</summary>
#md # ```

joint_summary = summary_table(chn_joint,
    [:r, :r0, :doubling_time, :T, :R_T, :CFR,
        :p_drc, :p_uganda, :k, :C_T]; digits = 2);

#md # ```@raw html
#md # </details>
#md # ```

joint_summary #hide

# ### Reproduction number over time
#
# The model fits a weekly random-walk reproduction number, so we can show
# the full daily $R_t$ trajectory rather than only its cut-off value
# $R_T$. The saved chain stores the sampled random-walk parameters
# (`rt_state.log_R0`, `rt_state.sigma_rw`, the innovation vector
# `rt_state.z` and `rt_state.intervention_effect`), so the daily $R_t$ is
# reconstructed per draw by mirroring the model's `rt_walk_model`. The
# median and 50%/90% ribbons are shown only over the established-outbreak
# window, the days on or after each draw's cumulative infections first
# reach one. The earlier seeding window (shaded) is prior-driven and is
# not an $R_t$ of an established epidemic, so it is left unplotted. The
# WHO-response breakpoint (red dashed) and the data cut-off (grey dashed)
# are marked.

#md # ```@raw html
#md # <details><summary>Reproduction-number trajectory</summary>
#md # ```

rt_fig = plot_rt(chn_joint;
    n = obs.n, breakpoint = _BREAKPOINT,
    rt_start = clamp(obs.n - round(Int, obs.tmrca_days), 1, obs.n),
    as_of_date = string(obs.cutoff), seeding = obs.seeding);

#md # ```@raw html
#md # </details>
#md # ```

rt_fig #hide

# The WHO response is fitted as a sampled multiplicative shift on $R_t$
# applied over a logistic ramp at the first situation report. The
# intervention effect is constrained to be non-positive, so it can only
# reduce transmission. The table reports its posterior on the
# multiplicative $\exp(\text{effect})$ scale: a value below one is the
# factor by which the response lowers $R_t$ once the ramp completes.

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

# The laboratory pipeline fits the received-specimen volume through the
# testing fraction `tau_test` and a receipt delay, and scores the confirmed
# positives as a Binomial of the observed specimens-analysed denominator
# with a partially-pooled per-window positivity. The implied per-suspected
# and per-test positivity and the non-BVD background rate `lambda_bg` are
# surfaced for comparison with the sitrep figures. The confirmed deaths are
# a thinning of the suspected deaths: the BVD composition among suspected
# (`death_composition`) enriched on the odds scale by `m_death` gives the
# death-confirmation probability (`death_confirmation`).

#md # ```@raw html
#md # <details><summary>Laboratory-pipeline posterior summary table</summary>
#md # ```

lab_summary = summary_table(chn_joint,
    [:tau_test, :lambda_bg, :suspected_positivity, :test_positivity,
        :expected_confirmed_T, :expected_received_T,
        :m_death, :death_composition, :death_confirmation,
        :expected_confirmed_deaths_T];
    digits = 3);

#md # ```@raw html
#md # </details>
#md # ```

lab_summary #hide

# The posterior pair plot shows the joint distribution of the key
# parameters, with the prior overlaid so the data's contribution to
# each marginal is visible.

#md # ```@raw html
#md # <details><summary>Posterior pair plot (prior overlaid)</summary>
#md # ```

posterior_pair_fig = plot_pair(chn_joint,
    [:C_T, :R_T, :r, :doubling_time, :T, :CFR, :k,
        :p_drc, :p_uganda];
    prior = prior_chn);

#md # ```@raw html
#md # </details>
#md # ```

posterior_pair_fig #hide

# A pair plot of the laboratory-pipeline parameters and the derived
# per-test positivity shows their joint posterior and any trade-offs, with
# the prior overlaid so the laboratory observations' contribution to each
# marginal is visible.

#md # ```@raw html
#md # <details><summary>Laboratory-pipeline pair plot (prior overlaid)</summary>
#md # ```

lab_pair_fig = plot_pair(chn_joint,
    [:tau_test, :lambda_bg, :test_positivity, :death_confirmation];
    prior = prior_chn);

#md # ```@raw html
#md # </details>
#md # ```

lab_pair_fig #hide

# A posterior predictive check draws replicated observations from the
# fitted joint model and compares them to the observed counts. The five
# dated DRC streams are real per-vintage observations, so each replicated
# cumulative trajectory is shown across the situation-report dates with the
# observed series overlaid: suspected cases, confirmed cases, suspected
# deaths, confirmed deaths and specimens received. The single-count Uganda
# export and export-death streams keep their scalar predictive check.

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
                  for (e, p, l) in
                      zip(vec(_conf_early), vec(_conf_obs), vec(_conf_late))],
    observed = [_confirmed_at(d) for d in _conf_window_days],
    colour = :goldenrod);

## Confirmed deaths are now a per-vintage stream (the series grows 17→64
## over 26 May-3 June), scored as increments of the modelled confirmed-death
## trajectory, so they get the same cumulative conditional check.
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

# The Uganda export and export-death streams are single cut-off counts,
# checked as scalar posterior predictives.

#md # ```@raw html
#md # <details><summary>Scalar posterior predictive plot</summary>
#md # ```

## The replicated counts are nested under their submodel prefix; match the
## full prefixed varname so the deterministic `expected_*_T` quantities are
## not picked up by a loose substring.
function _scalar_replicates(pp, name)
    key = first(k for k in keys(pp) if occursin(name, string(k)))
    return vec(Array(pp[key]))
end;

pp_exports = _scalar_replicates(pp_joint, "exports_state.exported_cases");
pp_exports_deaths = _scalar_replicates(
    pp_joint, "exports_deaths_state.exports_deaths");

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
# The lower bound on cumulative deaths if transmission stopped at the
# report date: every infection present by the cut-off still dies with
# probability CFR, so the committed future deaths are
# $\Delta D = \mathrm{CFR} \cdot C_T - \mathbb{E}[D_T]$.

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

# ### One-week-ahead forecast
#
# The seven-day no-change projection: cumulative and new expected counts
# per stream by $T + 7$, for the four DRC streams (suspected reported
# cases, suspected deaths, laboratory-confirmed cases and confirmed
# deaths). This continues the current growth rate (no interventions, no
# saturation) and carries both parameter and observation uncertainty.
# Exports are not forecast: cross-border travel is unlikely to be
# continuing at its baseline rate, so the forward travel rate the export
# model relies on no longer holds.

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

# New counts expected over the coming week, by stream:

#md # ```@raw html
#md # <details><summary>One-week-ahead forecast plot</summary>
#md # ```

forecast_fig = plot_forecast(forecast);

#md # ```@raw html
#md # </details>
#md # ```

forecast_fig #hide

# ### How the data streams compare
#
# Each data stream constrains the latent outbreak size differently.
# The table below puts the posteriors over $C_T$ side by side, the five
# single-stream fits and the joint, to show what each stream implies on
# its own and what the joint adds.

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

# Overlaid posterior densities of $C_T$ from the five fits:

#md # ```@raw html
#md # <details><summary>Overlaid C_T density plot</summary>
#md # ```

## Clip the x-axis so the exports-deaths heavy tail does not
## compress the other curves.
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
# The headline comparison is the model's estimate of the cumulative
# reported-case count, the ascertained quantity McCabe et al. report,
# against their 15 published scenario estimates. The McCabe et al.
# scenarios are cumulative reported *cases*, so they are not directly
# comparable to $C_T$, which counts latent *infections*; we compare like
# for like against the model's expected cumulative reported cases.
# For each scenario the table gives the narrowest joint credible
# interval that contains it, so coverage reads off directly.
#
# This renewal model has diverged from McCabe et al.'s exact
# construction: it fits a time-varying reproduction number on a renewal
# process rather than the report's closed-form exponential growth, and
# scores five jointly-fitted streams rather than their two-source Method 1
# and Method 2 sweep. So this is not an exact replication of either McCabe
# et al. method, and the per-scenario coverage below should be read as a
# like-for-like check on the reported-case count rather than a
# reproduction of their estimator. The single-stream
# [per-stream comparison](@ref "How the data streams compare") keeps each
# stream as close to McCabe et al.'s assumptions as the renewal
# parameterisation allows, for the most direct stream-by-stream contrast.

#md # ```@raw html
#md # <details><summary>Joint coverage table</summary>
#md # ```

posterior_reports_joint = vec(Array(chn_joint[:expected_reports_T]));
coverage_table = comparison_table(posterior_reports_joint);

#md # ```@raw html
#md # </details>
#md # ```

coverage_table #hide

# The joint expected-reported-cases density with the 15 published
# scenario estimates overlaid as faint dashed rules (both are
# reported-case counts, so the overlay is like for like):

#md # ```@raw html
#md # <details><summary>Joint reported-cases density with published scenarios</summary>
#md # ```

imperial_density_fig = plot_cumulative_cases(
    "joint (current data)" => posterior_reports_joint);

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
CSV.write(joinpath(output_dir, "scenario_coverage.csv"), coverage_table)

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
