# News

Release notes for BVDOutbreakSize.
Major versions of the report are kept as
[GitHub Releases](https://github.com/epiforecasts/BVDOutbreakSize/releases);
each push to `main` also republishes the rendered analysis and the
`output/` artifacts.

## v1.3.0

### Data

- Added `confirmed_case_history` and `death_history` blocks alongside
  `reported_case_history` in `data/observations.toml`, per INSP sitrep
  vintage, consumed by the per-vintage likelihoods.
- Added the laboratory observations from the sitrep section IV.3
  LABORATOIRE (`cumulative_tests_analysed`, `confirmed_cases`).
- Advanced the cut-off to 26 May 2026 and switched the DRC streams to the
  national cumulative totals read from the INSP situation-report PDFs
  (SitReps 009-012), rather than the per-zone CSVs whose zone sums drop
  cases not yet attributed to a zone. Figures were read by a
  language-model agent and independently re-scanned; recorded in
  `data/insp_sitrep_scanned.csv`. Cut-off (26 May): 1077 suspected cases,
  238 suspected deaths, 121 confirmed, 403 samples analysed; the 23 May
  suspected-death total uses the SitRep 009 zone-row sum (220), the
  headline (119) being a data-entry error.

### Modelling

- Reparameterised growth to sample the growth rate `r` directly
  (`LogNormal(log(log(2)/14), 0.4)`, McCabe et al.'s primary assumption),
  with the doubling time `τ = log(2)/r` and `m`, `T`, `C(T)` still
  exposed. This is the exact reciprocal pushforward of the old `τ` prior,
  so the implied priors are unchanged.
- Recentred the doubling-count prior on a base that advances with the
  cut-off: `m ~ Normal(m_prior_centre(as_of), 3)`, centre
  `9 + (cut-off − 18 May)/14` doublings, based on McCabe et al.'s
  first-report Method 2 central (501 cases, `log2 ≈ 9`) and the 14-day
  doubling time, so it tracks data refreshes.
- Added a partially-pooled DRC and Uganda ascertainment extension: a
  logit-scale reporting fraction applied to the latent incidence, fitting
  the reported suspected-case count alongside the deaths and exports.
- Added a laboratory pipeline coupling the cumulative tests-analysed and
  confirmed-case streams to the latent incidence, introducing a testing
  fraction, PCR sensitivity and a report-to-confirmation (lab-turnaround)
  delay, with right-truncation of the tested observation handled by the
  lab-delay CDF.
- Rewrote the suspected-cases stream as a BVD-driven onset-to-report
  convolution plus an additive non-BVD background rate, exposing the
  implied per-suspected positivity as a derived quantity.
- Fit the DRC suspected-case, laboratory-confirmed and suspected-death
  streams per sitrep vintage: `bvd_joint` conditions on the
  between-vintage increments rather than a single cut-off total, and a
  single-vintage stream reduces exactly to the cumulative likelihood,
  recovering the McCabe et al. configuration. Confirmed cases enter as
  per-vintage NegBinomial increments with per-test positivity a derived
  quantity, and each stream carries its own vintage offsets so a lagging
  stream is not assumed to run to the cut-off.
- Added `confirmed_only_model`, a single-stream composer that fits the
  laboratory pipeline in isolation for the per-stream comparison.
- Added `forecast_vs_truth_trajectory`: scores the retrospective forecast
  against the observed cumulative at every sitrep date across the horizon,
  not just the endpoint.
- Cut quadratures from the time-varying convolution: the laboratory
  background tested-volume integral now uses a closed-form gamma-CDF
  integral (`_gamma_cdf_integral`) instead of a per-draw quadrature, with
  an analytic reverse-mode rule, speeding up the lab-pipeline likelihood
  without changing the model.

### Outputs

- Posterior summary table and a laboratory-pipeline pair plot covering the
  report and lab delays, PCR sensitivity, testing fraction, background
  rate and the per-suspected and per-test positivity.
- Posterior-predictive panels for the confirmed and tests-analysed
  streams, included in the per-stream-versus-joint grid and the
  one-week-ahead forecast; the laboratory streams and the per-vintage
  time-series table also appear in the data table.
- `plot_vintage_ppc`: a posterior-predictive-across-the-sitrep-series
  figure that reconstructs the cumulative replicate at each vintage and
  overlays the observed trajectory, checking the fit against the whole
  series rather than only the latest total.

### Documentation

- Surfaced the onset-to-report and report-to-lab delay priors as
  equations; the onset-to-death prior means are the BDBV reanalysis
  estimates (about an 11-day mean) with standard deviations reproducing
  its 95% credible intervals.
- Clarified that the latent cumulative count is the true-case pool, not
  the tested or confirmed count, and framed the testing-fraction prior as
  weakly informative with no outbreak-specific data.
- Distributed the per-vintage increment maths into each submodel section.
- Cited the INSP situation reports and the INRB-UMIE archive.
- Added limitations on the constant exponential growth-rate assumption
  holding beyond the report period, and on per-sitrep increments mixing
  true incidence with backfill and rising ascertainment.

### Infrastructure

- Added optional Enzyme reverse-mode AD, selected with `enzyme_adtype()`,
  alongside the default Mooncake backend. Gradients match Mooncake across
  every model including the full joint and fitting runs at the same speed,
  so Mooncake stays the default.
- Fixed a posterior-predictive grid regression under AlgebraOfGraphics
  0.12 and widened the AoG compat bound to include 0.12; bumped the
  `softprops/action-gh-release` Action to v3.

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
