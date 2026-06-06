# News

Release notes for BVDOutbreakSize.
Major versions of the report are kept as
[GitHub Releases](https://github.com/epiforecasts/BVDOutbreakSize/releases);
each push to `main` also republishes the rendered analysis and the
`output/` artifacts.

## v1.3.0

### Data

- Advanced the cut-off to 3 June 2026 and switched the DRC streams to the
  INSP national cumulative totals from the situation-report PDFs
  (SitReps 009-020), rather than the per-zone CSVs whose zone sums drop
  cases not yet attributed to a zone.
- Added per-sitrep-vintage confirmed cases, confirmed deaths and
  laboratory throughput (samples received and analysed) to
  `data/observations.toml`, alongside the suspected cases and deaths.

### Modelling

- The joint model now fits four DRC streams per sitrep vintage (suspected
  cases, suspected deaths, laboratory-confirmed cases and
  laboratory-confirmed deaths) by conditioning on the between-vintage
  increments, alongside the Uganda exports and export deaths. A
  single-vintage stream reduces to the cumulative likelihood.
- Confirmed cases are fitted through a laboratory-throughput queue.
  Suspects enter a received backlog after a report-to-receipt delay prior,
  a capacity-limited drain (a log random walk on daily capacity, centred
  on the external ~150/day figure) sets the samples analysed, and the new
  positives in each window are a Binomial on the samples newly analysed
  (`ΔC ~ Binomial(ΔA, p_pos)`). Windows with no published analysed count
  use the queue's expected throughput and the Poisson-thinned marginal
  `ΔC ~ Poisson(μ_A · p_pos)`, so no free per-window denominator is
  introduced.
- Test positivity is severity-first: early specimens skew toward severe
  presentations and relax toward the latent case composition as cumulative
  analysed volume accrues (`severity_enrichment_model`).
- Suspected cases are a BVD onset-to-report convolution plus an additive
  non-BVD background rate `λ_bg`.
- Confirmed deaths are a laboratory process matching the cases: suspected
  deaths carry a non-BVD background `λ_bg_death`, and confirmed deaths are
  positives among the forwarded death specimens, sharing the case-lab PCR
  sensitivity and specificity (issue #193).
- Ascertainment is split into independent DRC and Uganda reporting
  fractions, with the DRC centre at 0.75: under active case finding most
  cases reach the suspected count, so the under-counting that matters is at
  laboratory confirmation.
- The growth prior is recentred on the molecular clock,
  `r ~ LogNormal(log(log(2)/20), 0.15)` (20-day doubling time from the BDBV
  phylodynamic reanalysis, [cuomodannenburg2026](@cite)), and the
  doubling-count base advances with the cut-off.
- Exports and export deaths are timed from infection via an
  infection→detection delay convolution (incubation ⊕ onset-to-report)
  rather than a rectangular detection window; both reduce to the McCabe et
  al. window as the delay collapses to a point mass, which stays available
  for comparison.
- The headline estimand is cumulative infections (`2^m`), with the
  under-ascertainment multiplier anchored on the laboratory-confirmed
  cases.
- Noted that the fatality ratio is applied per infection: with no
  asymptomatic fraction and no case-ascertainment on the death denominator
  it multiplies the latent infection trajectory directly, so it coincides
  with the infection-fatality ratio (IFR); the conventional CFR label is
  kept.

### Outputs

- Posterior summary table and a laboratory-pipeline pair plot over the
  report and lab delays, PCR sensitivity, positivity and background rate.
- Posterior-predictive panels for the confirmed-case and confirmed-death
  streams in the per-stream-versus-joint grid.
- Recast the forecast around the four trusted quantities (infections, true
  BVD deaths, confirmed cases, confirmed deaths) over two horizons: a
  one-week-ahead `forecast_reported` and a counterfactual-year
  `predict_committed` (committed totals under no onward transmission). The
  untrusted suspected cases, suspected deaths and tests-analysed streams are
  dropped from the forecast.
- Restored the forecast validation as a last-week-vs-now out-of-sample
  check on the two observable targets: fit the joint to the data through 28
  May (via the new `load_observations` `as_of_override` truncation),
  forecast six days, and score the predicted confirmed cases and confirmed
  deaths against the 3 June counts (`forecast_vs_truth`,
  `forecast_vs_truth_trajectory`, `plot_forecast_vs_truth`). Latent
  infections and all-BVD deaths are not directly validatable.
- The per-stream C(T) overlay fits the exports streams only up to their last
  data (the at-risk integral runs to the last import) while still reporting
  the implied C(T) at the cut-off.
- `plot_vintage_conditional_ppc`: a conditional one-step-ahead predictive
  across the sitrep series, each vintage conditioning on the observed
  previous cumulative and predicting only the new increment, replacing the
  unconditional version whose running sum compounded errors.

### Documentation

- Surfaced the delay priors as equations, clarified that the latent pool
  is the true-case count rather than the tested or confirmed count, and
  added limitations on the constant-growth assumption and on per-sitrep
  increments mixing incidence with backfill.

### Infrastructure

- Added streaming progress to `nuts_sample` via an optional `callback`:
  `progress_callback` writes a dependency-free file stream and
  `tensorboard_callback` streams step statistics to TensorBoard.
- Added optional Enzyme reverse-mode AD alongside the default Mooncake
  backend, with matching gradients across every model.

## v1.2.0

### Modelling

- Improved the comparison to the McCabe et al. report by making sure that 95% credible intervals are being compared and reordering it.
- Added a custom chain rule for `SpecialFunctions.gamma_inc`. This allows us to differentiate through the analytical solution to the gamma convolution integral.

### Data

- Moved the cut-off to 23 May 2026 and switched the DRC source from
  the WHO AFRO joint sitrep to the situation reports of the Institut
  National de Santé Publique (INSP), transcribed by
  [INRB-UMIE/Ebola_DRC_2026](https://github.com/INRB-UMIE/Ebola_DRC_2026).
  The INSP series gives a per-zone, per-sitrep daily vintage trajectory
  (suspected and confirmed; this analysis uses suspected). Cumulative
  counts at 23 May: 905 suspected DRC cases, 220 suspected DRC deaths,
  across the 12 reporting health zones. The 18 May INSP vintage (516
  cases, 131 deaths) matches the WHO joint sitrep 01 total exactly.
- Updated Uganda to three travel-related imports with one death,
  reflecting the third import announced on 23 May 2026 (woman from DRC
  who travelled Arua to Entebbe to Kampala; tested positive on
  follow-up). Two further Uganda-confirmed cases announced the same
  day (a driver and a healthcare worker) are domestic contacts of the
  first import and are excluded from `exported_cases` because the
  model treats Uganda as imports only.
- Added a `reported_case_history` block in `data/observations.toml`
  with eight INSP sitrep vintages (14 May to 23 May 2026), ready for
  the cumulative-trajectory likelihood once it merges.

### Infrastructure

- Moved the submodels out of the analysis file and into the supporting package. Instead we now print these in the analysis.
- Added additional package infrastructure including `Aqua.jl` and `Jet.jl`.
- Streamlined the package unit tests.

## v1.1.0

### Modelling

- Bound the seeding time `T` from below with a soft prior on the
  genetic time to the most recent common ancestor (TMRCA), following a
  suggestion from Neil Ferguson to combine the genetic signal with the
  other data streams as a seeding bound.
- Switched the export deaths to a daily (time-resolved binned) Poisson
  process: a continuous survival weight for the no-death stretch before
  the first dated death, then a per-day Poisson from that day to the
  cut-off.
- Bound `T` with export-death timing through that survival weight, and
  with case-export detection timing through a first-export-detection
  survival term on the Uganda admission date. Dates supplied in
  `data/observations.toml`.
- Death-convolution quadrature adapted to the sampled delay scale.
- Added a clock-rate sensitivity: refit the joint model under the
  faster 1.9e-3 early-epidemic TMRCA estimate and compare the impact on
  outbreak size, seeding time and growth rate against the 1.2e-3
  baseline.
- Sped up the deaths-among-exports likelihood: precompute the
  onset-to-death CDF once and reuse it across bin edges
  (`ExportDeathDelay`), replacing the per-node nested quadrature.
- Removed hardcoded death and case constants that diverged from the
  observations in `data/observations.toml`.
- Added a forecast validation: fit the joint model to the original
  report's data, project it forward to the current cut-off, and compare
  the predicted cumulative and new counts per stream against the counts
  observed since, as a table and a 2×3 coverage plot.

### Data

- Updated to the McCabe et al. 20 May 2026 report, comparing both
  report versions.
- Sourced the genetic TMRCA seeding bound from the BEAST temporal-tree
  estimate in the 2026-05-21
  [virological.org](https://virological.org/t/initial-genomes-from-may-2026-bundibugyo-virus-disease-outbreak-in-the-democratic-republic-of-the-congo-and-uganda/1032)
  update (mean 2026-03-25, 95% HPD 2026-02-20 to 2026-04-20, at the
  1.2e-3 EBOV clock rate this analysis assumes).

### Infrastructure

- Dropped MCMCChains for FlexiChains and prepared for registry
  release.
- CI docs preview PR comments and version-bump automation.

### Docs

- Added a scope note to the README and analysis report framing the
  work as an external view built on our understanding of real-time
  infectious disease dynamics, and inviting feedback, reuse and
  adaptation.
- Surfaced results from the README and analysis landing page, added
  stable and dev docs badges.
- Plotting and labelling fixes: surveillance dispersion on the 1/√k
  scale, predictive histograms labelled as frequency, and coarser
  (four-weekly) start-date axis ticks so the labels stay readable.
- Reworked the headline summary to report the credible intervals as
  sentences rather than leading with a median, defined the prior-IQR
  shift, and explained the reported-case scaling in terms of the DRC
  reporting fraction with a link to the pair plot.
- Replaced the model-structure diagram with a parameter-to-observation
  table.
- Culled promotional register in the analysis report.

## v1.0.0

First release.
A joint Bayesian re-analysis of the McCabe et al. report that fits all
data streams together in a single Turing model over the latent
cumulative case count.

- Conditions on the exported cases and DRC deaths the report uses,
  plus reported DRC cases (with an ascertainment component) and deaths
  among exported cases.
- Adds a no-onward-transmission projected-deaths counterfactual, a
  one-week-ahead forecast of newly reported cases, deaths and exports,
  and an onset-to-death delay sensitivity analysis.
- Replaces the deaths-convolution and small-growth-rate exports
  closed-form approximations with their exact forms.
- Maths-first analysis page with code folded behind dropdowns and a
  diagram of the model build-up.
- Compares against a joint reimplementation of the report's approach
  and its original published estimates.
