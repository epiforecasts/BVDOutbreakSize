# News

Release notes for BVDOutbreakSize.
Major versions of the report are kept as
[GitHub Releases](https://github.com/epiforecasts/BVDOutbreakSize/releases);
each push to `main` also republishes the rendered analysis and the
`output/` artifacts.

## v1.6.0

Changes since v1.5.0.

### Model

- The isolation BVD treatment stay (admission → outcome length of stay) is now
  grounded on our BDBV line-list reanalysis (`bdbv-linelist-analysis`
  submodule; doubly-censored Gamma fits with their posterior uncertainty, not
  the source's raw means): the fitted admission→death (7.6 d, 95% CrI 5.6–10.4)
  and admission→discharge (7.7 d, 95% CrI 4.8–13.8) stays, mixed at the
  admitted-case fatality (22/37 ≈ 0.59), give a ≈7.6 d mean, 6.0 d SD stay,
  replacing the previous weakly-informative ~12 d prior. Because the isolation
  occupancy reads beds back into a suspect inflow as `occupancy ≈ admission
  fraction · admissions · (stay + 1)`, the shorter, evidence-based stay relaxes
  the downward pull the treatment-bed stream placed on the inferred outbreak
  size.
- The reproduction-number random walk now starts two weeks BEFORE the first
  situation report (the new `RT_WALK_LEAD = 14` lead, exposed as the
  `bvd_joint` keyword `rt_walk_lead`) rather than exactly at it, so `R_t` is
  free to move over the fortnight of transmission leading up to the first
  report instead of being held flat at `R0` right to it. The walk start is
  floored at the renewal start so it never precedes the seeded trajectory, and
  the `plot_rt` reconstruction uses the same earlier knot grid.
- Added a supply-limited isolation/treatment-bed stream ("Patients en
  isolement"), the renewal analogue of the convolution secondary-observation
  model of EpiNow2.
  Bed occupancy has been supply-driven (demand has outstripped supply), so
  the model fits a latent bed DEMAND — the suspect inflow carried through a
  length-of-stay survival (BVD cases with a sampled treatment stay, non-BVD
  suspects leaving after a sampled rule-out stay) — and the occupancy is
  that demand soft-capped at the bed capacity,
  `occupancy = C(t)·(1−exp(−D/C(t)))`, so the admitted fraction is
  availability-driven. The capacity `C(t)` is a random walk
  (`bed_capacity_walk_model`), so the ceiling tracks the beds being added and
  can be projected forward; it is pinned by the implied bed count (reported
  occupancy / "Taux d'occupation" rate) on the days a rate is published.
  The stream exposes the bed demand, occupancy, capacity, shortfall and
  utilisation, and carries its own observation dispersion. Added
  `convolve_survival`, the `treatment_admission_model` observation submodel,
  the `isolation_admission_model`, `bed_capacity_model` and
  `bed_capacity_walk_model` priors and the `treatment_only_model`
  single-stream composer (resolves #265).
  A key limitation is that this is a single national model: it cannot
  represent local bed saturation (Ituri at 93.9% occupancy on 13 June against
  Sud-Kivu 21.9%), which is the level the supply constraint operates at, so
  the national shortfall understates the local unmet need. The renewal model
  does not carry per-province inflow, so it cannot be split to the
  per-province bed model the constraint needs.
- Gave the non-BVD isolation rule-out stay its own sampled length-of-stay
  (`ruleout_los`) in `treatment_admission_model`, separate from the
  report-to-receipt laboratory delay it previously shared. Sharing one delay
  let the daily-observed occupancy stream, where the delay mean scales the
  non-BVD bed demand, collapse the lab-turnaround delay that the testing,
  composition and confirmed-death streams use only as a time shift. The
  occupancy now identifies the rule-out stay on its own clock, and the
  lab-turnaround delay is set by those streams plus its prior. Exposes
  `isolation_ruleout_los_mean`.
- Added a recovered-among-confirmed stream ("cumul guéris"), the
  secondary-observation incidence analogue: survivors among the modelled
  daily confirmed cases, scaled by a recovery proportion and convolved with a
  confirmation-to-recovery delay, with its own observation dispersion. The
  recovery proportion is grounded on the case-fatality ratio (a recovered
  case is one that did not die) with a log-odds adjustment for the confirmed
  population, rather than estimated from scratch.
  The confirmed model now exposes one daily confirmed-case series
  (`confirmed_daily`) that both the recovered stream and the cumulative-
  confirmed trajectory reuse. Added the `recovered_model` submodel and the
  `recovery_probability_model` prior.
- Added the `exports_joint_only_model` composer, which fits the Uganda export
  cases and deaths together over the one travel-gated at-risk prevalence. The
  single-stream comparison in the walkthrough now shows one joint "exports"
  fit instead of separate export-case and export-death fits.
- The one-week-ahead forecast now also projects the isolation/treatment beds
  — both the bed demand a week ahead (need under unconstrained supply) and the
  supply-limited occupancy it produces, whose gap is the projected bed
  shortfall — and the cumulative recovered total, each replicated with its
  own dispersion. Added `plot_forecast_beds`, which shows the projected bed
  need against the supply-limited occupancy and the shortfall in the
  walkthrough's forecast section.
- The forecast-versus-frozen validation now also scores the isolation beds:
  the frozen one-week-back fit conditions on the isolation occupancy, and the
  projected bed occupancy is compared against the beds actually held a week
  later (`forecast_vs_truth` gains an `isolation` argument, and
  `plot_forecast_beds_vs_truth` plots the projected occupancy against the
  observed beds). The bed check is weak at a one-week-back freeze because the
  reported occupancy rate starts only on 9 June, so the capacity rides its
  random walk back to the freeze date.
- Replaced the per-vintage step background random effect with a smooth daily
  lognormal random walk (`background_walk_model`), fixing the joint convergence
  failure on the current data (the per-vintage effect opened a second
  "background explains the majority of suspected cases" posterior mode that
  split the chains; the published dev joint had max R-hat 2.5).
  The background is now per day with no reporting-vintage steps, gated to begin
  a report-to-receipt lead before the first suspected-case report (it does not
  exist before surveillance began), shared across the suspected-case and
  suspected-death streams, with a half-normal baseline and a tight random-walk
  innovation SD.
  This also removes the per-vintage step in the modelled cumulative-death
  trajectory.
- Widened the non-BVD background level prior so the laboratory positivity
  identifies it.
  The previous prior, calibrated for the un-gated full-grid background, forced
  the background to a tiny fraction of the suspect pool, inconsistent with the
  observed per-test positivity (210/755 ≈ 0.28).
  The suspect pool is now inferred to be a minority BVD, raising the non-BVD
  background substantially and lowering the cumulative-infection estimate
  (C_T) accordingly, with a wider credible interval.
- Gated the laboratory analysed-specimen capacity to the testing onset, so no
  specimens are modelled as analysed before testing existed.
  This removes a large over-prediction of the early confirmed-case counts (the
  modelled cumulative confirmed previously ran several-fold above the observed
  early values).
- Redesigned the death pathway.
  Suspected deaths carry a death ascertainment `p_death` (the death analogue
  of the case ascertainment, with an informative prior centred high) and a
  non-BVD death background that applies a background CFR (`cfr_bg`) to the
  suspected-case background and lags it by the onset-to-death delay, so a
  background death follows its background case. The death background tracks the
  identified case background rather than a second free, outbreak-size-
  degenerate rate, and inherits the case background's smooth gated daily shape,
  so the modelled cumulative-death trajectory is smooth.
  Added the `death_ascertainment_model` and `background_cfr_model` priors.
- Rebuilt the confirmed-death stream as a laboratory pipeline mirroring the
  confirmed cases.
  The death analysed volume scales the modelled case analysed volume at the
  per-day suspected death-to-case ratio, times a testing-intensity scaling
  (`LogNormal(0, 0.25)`, centred on one), so death testing follows the
  laboratory's realised throughput; the death-to-case ratio carries the
  suspect-pool severity and the suspected-death level. The volume is scored
  through a death-pool composition positivity
  `p = s·q_death + (1−spec)(1−q_death)`, with `q_death` the BVD share of the
  suspected deaths from the death series' own components. The case volume
  carries the laboratory capacity onset, so the death volume inherits it and no
  deaths are confirmed before testing began.
  The joint exposes the `death_ascertainment`, `background_cfr`,
  `death_testing_scaling`, `tau_death` and `death_composition` deterministics
  and drops `m_death`.

### Data

- Added the daily "Patients en isolement" occupancy for 1-11 June (SitReps
  018-028) as a structured `patients_isolated` column and the
  `[isolation_history]` manifest block.
  The fitted series begins 1 June where the all-patients column definition
  is stable; the narrower suspects-only count in SitReps 016-017 is a
  different quantity and is excluded.
  Corrected the SitRep 020 note (the PDF headline occupancy is 233, not the
  173 the note claimed).
- Added the implied bed-capacity series (occupancy / reported "Taux
  d'occupation" rate ≈ 400-452 beds, 9-13 June) as the `[bed_capacity_history]`
  manifest block, which pins the bed capacity in the supply-limited isolation
  model.
- Added the cumulative "cumul guéris" recovered-among-confirmed total for
  6-11 June (SitReps 023-028) as a structured `cumul_recovered` column and
  the `[recovered_history]` manifest block.
- Added the peer-reviewed McCabe et al. Lancet Infectious Diseases publication
  (online first 9 June 2026, DOI 10.1016/S1473-3099(26)00299-9) as a third
  scenario vintage in `REPORT_SCENARIOS_CI`, with inputs as of 27 May 2026
  (1031 DRC cases, 240 deaths, three Uganda imports).
  Both methods now vary the epidemic doubling time (7/10/14 d); the
  back-calculation assumes 30% of deaths are attributable to Ebola.
  The published paper swaps the method numbers relative to the Imperial
  reports, which the noted convention reconciles.
  Recorded the matching frozen-data snapshot in
  `data/report-snapshot-27may.toml`.

### Analysis

- The walkthrough adds posterior-predictive panels for the isolation
  occupancy and recovered streams, a single-stream "in isolation" fit for the
  isolation occupancy, and surfaces the isolation length-of-stay and
  confirmation-to-recovery delays in the observation-delay table and pair
  plot, and the admission proportion, recovery probability and the two new
  per-stream dispersions in the surveillance-parameter table.
- Cite EpiNow2 [epinow2](@cite) for the convolution-and-scaling
  secondary-observation analogy, and fix the `epinow2` bibliography entry so
  the documentation build no longer warns about a missing field.
- The McCabe et al. scenario comparison now carries a third vintage, the
  27 May 2026 Lancet publication, plotted beside our renewal estimate on
  27 May, with a frozen re-fit at the 27 May cut-off added to the
  frozen-fit outbreak-size table.

### Outputs

- Added the latent symptom onsets (the "symptomatic cases" outcome) to the
  shared posterior outputs.
  `posterior_draws.csv` gains a `cumulative_onsets_T` column, the cumulative
  symptom onsets by the cut-off per draw (the onset analogue of `C_T`), and a
  new `onsets_over_time.csv` records the daily new and cumulative onset
  trajectory over time with 30/60/90% credible intervals, the curve plotted
  to show the outbreak grow.
  Exposed through `onsets_over_time`.

## v1.5.0

Changes since v1.4.0.

### Model

- Confirmed deaths now carry the report-to-receipt laboratory delay, so the
  laboratory-confirmed-death series lags the death event rather than tracking
  it instantaneously and the confirmed case and death streams pay a consistent
  laboratory delay.
- Added an optional daily new-suspect inflow stream ("nouveaux cas suspects
  du jour") to the suspected-case likelihood.
  The post-26 May per-day counts are scored against the modelled daily
  suspected series at each report day, continuing the suspected signal where
  the frozen cumulative series stops, on days disjoint from it (#222).
- Added the deaths analogue, an optional daily new suspected-death inflow
  stream ("cas suspects du jour N (M deces)") to the suspected-death
  likelihood.
  The post-26 May per-day counts are scored against the modelled daily
  suspected-death series at each report day, continuing the suspected-death
  signal where the frozen cumulative series stops, on days disjoint from it,
  and a matching "New suspected deaths/day" posterior-predictive panel is
  added alongside the new-suspects-per-day panel.
- Collapsed the laboratory pipeline onto a single suspected-to-analysed
  volume, fit to the specimens-analysed series through one report-to-analysed
  delay.
  The received stream is still recorded but no longer fitted, and the
  post-cut-off 24-hour analysed volume is now fit directly.
- Late reporting windows, where the cumulative national analysed denominator
  stops, are scored in one submodel.
  A day with a published 24-hour analysed count anchors its positivity as a
  binomial on that count, and the remaining days are scored against the
  modelled laboratory volume.
  These are a reporting-format change rather than data blackouts, so the
  earlier "dark window" framing is dropped.
- Added the delay-corrected confirmed case-fatality ratio, the Nishiura et al.
  (2009) real-time correction computed per posterior draw on the modelled
  confirmed trajectory and sampled confirmation-to-death delay.
  The denominator shrinks from all confirmed cases to those expected to have
  had a fatal outcome resolve by the cut-off, debiasing the naive confirmed
  ratio.

### Forecast

- The one-week-ahead forecast and its validation now target the
  laboratory-confirmed case and confirmed death streams.
  The suspected reported cases and deaths are no longer published, so they no
  longer serve as forecast targets or as the last-week-versus-now comparison.

### Report

- Added a confirmed case-fatality ratio section, setting the delay-corrected
  confirmed CFR against the structural infection-based CFR and the naive
  confirmed ratio, with a comparison table and posterior-density plot.
- The estimate-evolution figure now draws each release as a discrete per-fit
  estimate with nested 30/60/90% intervals, read from
  `data/released_estimates.csv` rather than a hand-maintained literal.
  Frozen renewal re-fits are restricted to the integral-era release cut-offs,
  since renewal-era releases already are renewal fits.
  A new `scripts/refresh_releases.jl` pulls the per-release estimates from the
  tagged results releases.

### Data

- Captured the daily new-suspect counts (SitReps 021-024, 4-7 June) as a
  structured `new_daily_suspects` column in the scanned situation-report CSV
  and a `suspected_daily_history` block in the observation manifest.
- Captured the daily new suspected-death counts ("cas suspects du jour N (M
  deces)", SitReps 024-032, 7-15 June) as a structured
  `new_daily_suspected_deaths` column in the scanned situation-report CSV and a
  `suspected_daily_deaths_history` block in the observation manifest, the
  deaths analogue of the daily new-suspect inflow.
- Extended the confirmed case and death series to SitRep 025 (8 June).
- Added the trusted-day 24-hour analysed laboratory counts (1, 4-7 June) as a
  `tests_analysed_daily_history` block to anchor late-window positivity.

## v1.4.0

The methods switch flagged in v1.3.0: the continuous-time, fixed-growth-rate
model is replaced by a discrete-time renewal model that is simpler and avoids
the single-stream-versus-joint size tension of issue #212.
This is a substantial revision; the changes below are relative to v1.3.0.

### Model

- Replaced the integral exponential-growth model with a discrete-time renewal
  process on a daily grid.
  Infections follow the renewal equation under a time-varying reproduction
  number (a weekly log-scale random walk with an intervention ramp), and every
  observed stream sits downstream of latent onsets through its own sampled,
  discretised delay.
- The prior is placed on the growth rate (the molecular-clock doubling time)
  and the first reproduction number is derived forward through Euler–Lotka.
  The generation interval is a Gamma with shape and scale taken from the cited
  source and its reported uncertainty.
- The onset-to-event delays are taken from a Bayesian reanalysis of the 2012
  Isiro line list on their natural Gamma parameters, with one onset-to-
  admission delay serving both suspected-case reporting and export detection
  and onset-to-death the convolution of two atomic components.
- Two-phase seeding: a single import grows through an unobserved cryptic
  exponential phase to the renewal start, with the outbreak age bounded by the
  genetic time to the most recent common ancestor.
- Confirmed positivity is tied to the suspect-pool composition through an assay
  sensitivity and specificity, and exports are travel-gated from infection and
  scored on their dated detection days.
- The DRC streams are fitted on the incidence scale, as the between-vintage
  increments across successive situation reports (the first vintage being the
  cumulative count to that date).

### Report

The report was rebuilt around the renewal model; the analyses carried over
from v1.3.0 (the one-week-ahead forecast and its validation, the
no-onward-transmission counterfactual, the delay and clock-rate sensitivity
analyses, and the McCabe et al. comparison) were re-implemented for the new
model rather than added here.

- Restructured the methods in generative order (infections, epidemiological
  processes, observation models, the joint model) with the model maths given
  explicitly.
- Reworked the figures (reproduction number with credible ribbons and sampled
  trajectories, cumulative infections, onsets and deaths, outbreak size by data
  stream, and estimate evolution across releases).
- The McCabe et al. scenario comparison now carries their reported 95%
  confidence intervals.
- The one-week-ahead forecast and its last-week-versus-now validation now
  target the laboratory-confirmed cases and confirmed deaths. The suspected
  reported cases and deaths are no longer reported, so they are dropped as
  forecast targets.

### Data

- Advanced the cut-off to 7 June 2026 (SitRep 024).
  The laboratory-confirmed streams run to the cut-off while the suspected
  streams stay frozen at their 26 May values.

## v1.3.0

!!! warning "Final release of this model formulation"

    This is the last planned release of the continuous-time,
    fixed-growth-rate model. The laboratory and testing observation
    model has outgrown the available data, and the joint fit now implies
    a larger outbreak than any single data stream does on its own
    (issue #212). We are replacing this model with a discrete-time
    renewal model, which is simpler and avoids these problems, in a
    follow-up release that will note the methods switch. Treat the
    estimates here as provisional.

Changes since v1.2.0.

### Data

- Advanced the model cut-off to 28 May 2026 and switched the DRC streams
  to the INSP national cumulative totals read from the situation-report
  PDFs, rather than the per-zone CSVs whose zone sums drop cases not yet
  attributed to a zone. The suspected streams are frozen at their 26 May
  values, after which INSP stopped publishing a national suspected total;
  the confirmed and laboratory streams run to 28 May.
- Added per-sitrep-vintage confirmed cases, confirmed deaths and
  laboratory throughput (samples received and analysed) to
  `data/observations.toml`, alongside the suspected cases and deaths.
- Extended the observed series through 5 June 2026 to validate the
  forecast out of sample.

### Modelling

- The joint model now fits four DRC streams per sitrep vintage (suspected
  cases, suspected deaths, laboratory-confirmed cases and
  laboratory-confirmed deaths) by conditioning on the between-vintage
  increments, alongside the Uganda exports and export deaths. A
  single-vintage stream reduces to the cumulative likelihood.
- Confirmed cases are fitted through a laboratory-throughput queue:
  suspects enter a received backlog after a report-to-receipt delay, a
  capacity-limited drain sets the samples analysed, and the new positives
  in each window are a Binomial on the samples newly analysed. Windows
  with no published analysed count fall back to the queue's expected
  throughput, so no free per-window denominator is introduced.
- Test positivity is severity-first: early specimens skew toward severe
  presentations and relax toward the latent case composition as analysed
  volume accrues.
- Suspected cases and deaths are BVD onset-to-report convolutions plus
  additive non-BVD background rates; confirmed deaths share the case-lab
  PCR sensitivity and specificity.
- Split ascertainment into independent DRC and Uganda reporting fractions
  (DRC centre 0.75), and recentred the growth prior on the molecular-clock
  20-day doubling time ([cuomodannenburg2026](@cite)).
- Exports and export deaths are timed from infection via an
  infection→detection delay convolution rather than a rectangular
  detection window, reducing to the McCabe et al. window as the delay
  collapses to a point mass.
- The headline estimand is cumulative infections (`2^m`), with the
  under-ascertainment multiplier anchored on the laboratory-confirmed
  cases.

### Outputs

- Posterior summary table, a laboratory-pipeline pair plot, and
  posterior-predictive panels for the confirmed-case and confirmed-death
  streams in the per-stream-versus-joint grid.
- Recast the forecast around the four trusted quantities (infections,
  true BVD deaths, confirmed cases, confirmed deaths) over a
  one-week-ahead and a counterfactual-year horizon, dropping the
  untrusted suspected and tests-analysed streams.
- Restored the forecast validation as a last-week-vs-now out-of-sample
  check: fit the joint through 28 May, forecast forward, and score the
  predicted confirmed cases and deaths against the observed counts.
- Added a conditional one-step-ahead predictive across the sitrep series,
  each vintage predicting only its new increment.

### Documentation

- Surfaced the delay priors as equations, clarified that the latent pool
  is the true-case count rather than the tested or confirmed count, and
  added limitations on the constant-growth assumption and on per-sitrep
  increments mixing incidence with backfill.

### Infrastructure

- Added streaming progress to `nuts_sample` via an optional callback
  (a dependency-free file stream or TensorBoard), and optional Enzyme
  reverse-mode AD alongside the default Mooncake backend.

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
