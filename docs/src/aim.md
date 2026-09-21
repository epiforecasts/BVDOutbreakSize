# Aim and origins

The aim is a transparent estimate of how far the outbreak has already grown,
and of how much each published data stream contributes to that estimate.
Most infections are not yet reported, so the current size has to be inferred
from the surveillance data that are available.
What the estimate can and cannot support is set out in the
[limitations](limitations.md).

## Origins of this work

This work began as a replication of the [mccabe2026](@citet) report.
It has since evolved into a real-time joint Bayesian estimate of the current outbreak size.
The model is a discrete-time renewal process with a time-varying reproduction number, fitted to more of the available data streams than the original.
The [methods](analysis.md#Methods) carry the full treatment, and the [comparison with McCabe et al.](@ref "Comparison with McCabe et al.") sets the current estimates against theirs.

**Latent process and parameters**

- *Discrete-time meta-population renewal model.* The whole model runs on a daily grid.
  Each province follows the discrete renewal equation $I_{p,t} = R_{p,t} \sum_{s \ge 1} I_{p,t-s} g_s$, where $g$ is the discretised generation-interval PMF, and the provinces are coupled by importation.
  National incidence is their sum, and every delay is applied as a discrete convolution.
  [mccabe2026](@citet) use continuous-time closed forms.
- *Time-varying reproduction number.* The trend $R^{\text{trend}}_t$ the provinces pool toward is held flat at the established $R_0$ until the first WHO situation report (18 May 2026).
  It then follows a weekly Gaussian random walk on the log scale, interpolated within weeks.
  A logistic outbreak-response ramp of about three weeks starts from that report.
  Each province's reproduction number is that trend plus a mean-reverting deviation, and the deviations sum to zero.
  McCabe et al. use one constant exponential growth rate.
- *Joint posterior rather than scenario estimates.* The reproduction number, case-fatality ratio, all delays, traveller volume and surveillance dispersion have priors and are sampled together.
  [mccabe2026](@citet) fix each and report a set of scenarios.
- *Two-phase seeding with a wide, genetically-floored outbreak age.* A single import grows through an unobserved cryptic exponential phase before the renewal process takes over.
  Growth follows the rate the genetic estimate informs, reaching a magnitude set by a prior on the number of cryptic generations.
  The established reproduction number is derived forward from that growth rate.
  The genetic time to the most recent common ancestor floors the cryptic duration from below.
  McCabe et al. fix the start from a single seed.

**Delays and convolutions**

- *Delays re-estimated with uncertainty.* [mccabe2026](@citet) take the onset-to-death delay from the Isiro 2012 point estimate of [rosello2015](@citet).
  We instead use a Bayesian reanalysis of the same line list [bdbv_linelist_analysis_2026](@cite) that re-estimates the delay with uncertainty.
  We sample every other delay (generation interval, incubation period, onset-to-report, onset-to-confirmation and onset-to-hospitalisation abroad) from a prior centred on published Ebola estimates.
  Each is discretised with double interval censoring [charniga2024](@cite), so the delay uncertainty propagates.

**Likelihoods and data streams**

- *More streams fitted.* [mccabe2026](@citet) fit the Uganda export cases and deaths.
  We add the DRC suspected cases, the laboratory-confirmed cases, the confirmed deaths and the deaths among the Uganda exports.
- *Per-vintage time-series fitting.* The DRC streams are fitted on the incidence scale, as the between-vintage increments across successive sitreps (the first vintage being the cumulative count to that date).
  This sharpens $R_t$.
  McCabe et al. condition on a single cumulative total.
- *Ascertainment estimated.* We jointly estimate the outbreak size and the fraction of cases each surveillance system reports.
  McCabe et al. have no ascertainment component.
- *Comparison against published scenarios.* The model is set beside the [mccabe2026](@citet) scenario estimates as an external sense-check, matched in time at the cut-off each scenario was computed.
  The cumulative infection count, the running sum of the daily infections, is the headline quantity reported separately.
  A forward projection from a frozen fit is also set against the [chamla2026](@citet) confirmed-case projection and the data observed since.

**Extensions**

- *No-onward-transmission counterfactual and one-week-ahead forecasts.* Future expected deaths from infections already seeded, and a posterior-predictive projection of each stream.
