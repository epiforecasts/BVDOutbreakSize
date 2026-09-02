# News

Release notes for BVDOutbreakSize.
Major versions of the report are kept as
[GitHub Releases](https://github.com/epiforecasts/BVDOutbreakSize/releases);
each push to `main` also republishes the rendered analysis and the
`output/` artifacts.

## v1.16.0

Changes since v1.15.0

### Report

- The cross-release forecast scores, the forecast-versus-observed validation tables and the situation-report data tables stay tables on the page.
A `DataFrame` is html-showable, so each one went into the page as a raw HTML block, and Documenter compiles a regex from every raw block's own text, which fails once a block passes PCRE's compiled-pattern limit.
The scoring tables grow with every release and crossed it after v1.15.0, taking the docs build down with them.
The stop-gap printed each table as a fixed-width block of text instead, which built but is not a table.
They now go out as markdown tables, which carry no such limit.
Float columns are rounded rather than printed to their full binary expansion, and numeric columns are right-aligned.

### Data

- Advanced the model cut-off from SitRep 102 (24 August) to SitRep 108 (30 August), six vintages.
Confirmed cases reach 6100 and confirmed deaths 2950.
Two health zones enter the case table at SitRep 104, Biena and Manguredjipa in Nord-Kivu, taking the affected count from 58 to 60 of 151.
Bas-Uélé prints its first bed-capacity figure at SitRep 106.
The daily suspected case and death series and the seven treatment-centre streams stay frozen, as they have since the reports dropped the tables that carried them.
- SitRep 107's isolation occupancy is excluded rather than recorded.
Its page-1 tile reconciles exactly to the province prose, but only because Nord-Kivu is absent from that section altogether, having carried 236 the day before and 282 the day after.
Bed capacity for the same date is kept, on the split already established at SitRep 085.
- The onset-curve digitiser reads SitRep 108.
That vintage embeds its chart at a size no neighbour uses, which anti-aliased the axis border and the weekly ticks too light for the mask the digitiser reads dark pixels with.
A looser mask now runs only where the strict one finds nothing, so every earlier vintage still digitises to the committed values.

### Dependencies

- Turing 0.47 is allowed alongside 0.45 and 0.46, in the package and in the Enzyme test environment.

## v1.15.0

Changes since v1.14.0

### Fixes

Five defects in the one-week-ahead forecast are corrected here.
All five affected the headline forecast in the published report, not only the evaluation, so skill numbers from earlier releases were measured against a forecaster carrying them.
That is why 90% coverage sat at 1.00 at nearly every stream and horizon against a nominal 0.90.

- Each observed count stream now projects from its own cumulative trajectory.
`cumulative_reports` and `cumulative_deaths_total` were named by the forecast and by the stream registry but defined in no model, so reported cases, suspected deaths and confirmed deaths all fell back to inverting the cumulative total under exponential growth.
That inversion collapses towards zero once the fitted growth rate is at or below zero, and on one cached fit gave a median 0.55 confirmed deaths per day against an observed 17.
- The confirmed-death forecast is no longer capped at the cumulative suspected deaths.
The model does impose a thinning, but per day and on the latent pool, before positivity scales it down.
The forecast capped a different quantity against a reported headline that froze at 246 on 26 May while confirmed deaths passed 246 on 19 June, clamping the forecast below its own origin and flooring the new count at zero in most draws.
- Each stream is replicated through its own dispersion rather than the population mean.
The per-stream forecaster already used the right one, so the two disagreed for the same model on the same stream.
- The observation replicate accumulates day by day rather than in one draw on the horizon total.
The dispersion is fitted against single-day counts and short vintage increments, so carrying it on the summed mean inflated the overdispersion term by roughly the length of the horizon.
- The reproduction number continues its fitted walk past the cut-off.
The forecast repeated the last sampled innovation as a fixed daily slope for the whole horizon, so the spread in log-`R_t` grew with the horizon rather than with its square root, putting the four-week reproduction number between 0.245 and 14.9.
It now draws fresh weekly innovations at the fitted step scale and interpolates between them, as the walk is built in the first place.

The forecast validation is corrected with them.

- A retrospective harmonisation now comes out of the cumulative truth as well as the new-count truth.
A projection cannot contain an administrative reattachment, so leaving it in the truth scored the forecast against something it could not produce.
This is a no-op at the current cut-off and bites whenever a break day falls inside a validation week.

### Report

- The joint posterior predictive checks are split by whether a stream is still reported, cumulative then per-vintage within each group, and the vintage axes tick weekly rather than crowding as the outbreak runs on.
- The forecast is validated only against streams the situation reports still update.
Reported cases and suspected deaths stopped on 26 May, so their new-count truth was a guaranteed zero and the figure was scoring a forecast against an unmoving series, which the release scoring already withheld.
They are still drawn, in their own figure, as a projection rather than a validation.
- The sensitivity page said the confirmed new-count rows keep any retrospective harmonisation step.
The code has subtracted it since the break-day correction was consolidated.

### Data

- Advanced the model cut-off from SitRep 100 (22 August) to SitRep 102 (24 August).
Both are clean days, with the cumulative case and death increments equal to the printed 24h gross.
Bas-Uélé prints its first traveller count at SitRep 102.

### Infrastructure

- The docs-preview cleanup collects previews that land after the pull request closes.
The close-triggered job assumed the preview existed by then, but the docs build is fanned across runners and publishes hours later, so it found nothing and exited green every time, leaving two hundred preview directories against three open pull requests.
A nightly sweep now drops every preview directory with no open pull request behind it, and refuses to run when the pull-request listing fails.

### Known issues

- The release comparison mixes two forecast constructions.
Every row in it is reconstructed by running each release tag's own code, so the historical rows keep the defects corrected above and cannot be regenerated without rewriting what those tags would have produced.
Read a change in forecast width across this release as a change in the code rather than in the outbreak.
- The automatic version increment is wedged by a stale branch and skips silently on every push, so `Project.toml` does not advance on its own (#607).

## v1.14.0

Changes since v1.13.2

### Data

- Advanced the model cut-off from SitRep 089 (11 August) to SitRep 100
  (22 August). The reports returned to their full-length format at SitRep
  090 after five vintages of the shorter "MVEBDB" brief, restoring the
  laboratory, point-of-entry and continuity-of-care sections. The
  confirmed-case and confirmed-death series, recoveries, isolation
  occupancy, bed capacity and the 24h analysed volume all advance with
  them, and the symptom-onset curve is digitised for each new vintage that
  carries a fresh figure. The daily new-suspected-case count stays frozen
  at 5 August and the treatment-centre patient-movement series at 2 August.
- A sixth province, Bas-Uélé, records its first confirmed case at SitRep
  090 (ZS Buta). No isolation, occupancy or laboratory figures are printed
  for it in any vintage, so it contributes zero to every per-province
  stream, following the convention already in place for a non-reporting
  province.
- The 17 August point was recorded from the INRB-UMIE mirror alone while
  INSP had published nothing beyond SitRep 093, verified against its posts
  and media endpoints rather than inferred from a failed fetch. Publication
  has since resumed and every stream is scanned from the reports again.
  `scripts/check_new_sitreps.jl` reports how far the mirror leads insp.cd,
  so a mirror-only point stays visible rather than silently absorbed.
- `candidate_signals.csv` gains a `province` column and holds one row per
  signal per vintage per province, replacing the national figure that was
  a hand-made sum over whichever provinces happened to print a count. A
  province the report is silent on now gets no row, which distinguishes it
  from a printed zero. `province = Ensemble` marks a national figure the
  report prints itself. Coverage changes are now readable from the data
  rather than only from the prose in `source_note`: SitRep 069, for
  example, gives Ituri and Nord-Kivu for `eds_death_alerts` but Ituri
  alone for `eds_investigations_performed`.
- Corrected SitRep 089's `eds_investigations_performed`, whose national
  value of 107 included the 15 CTE EDS its own note recorded as excluded.
  The per-province rows carry community-alert EDS only and sum to 92.

### Fixes

- The still-expected-deaths density no longer extends below zero. The draws
  were never negative; the kernel density spread mass past the smallest
  draw, putting a tail on the impossible side of the bound. The axis now
  stops at the bound, as the count and CFR panels already did. The same
  applies to the projected-total panel, which had been drawn crossing left
  of the observed-deaths rule.
- The forecast reproduction-number density stops at zero, the same bound
  already applied to the still-expected-deaths and projected-total
  densities. A posterior sitting close to zero picks up a tail on the
  impossible side from the kernel.
- The posterior outbreak age on the tree-prior sensitivity figure stops at
  zero for the same reason. An age in days cannot be negative, and that
  figure renders in every published build.
- The doubling-time interval is the image of the growth rate's interval
  rather than the quantiles of its own draws, which bounded nothing once
  the posterior for the growth rate spanned zero and reported a range of
  -306 to 301 days. The row now runs from the fastest decline, through the
  zero-growth pole, to the fastest growth, and is ordered by growth rate
  rather than by value. The narrative summary carried the same fault and is
  fixed with it.
- The joint model's score tables show the joint model alone. Every fit was
  rendered as a row, so each stream's individual fit appeared twice, once
  in the headline table and again in the section built for it. A frozen
  fit, which is the joint model refit at an earlier cut-off, was also read
  as an individual fit rather than a joint one.
- The persistence baseline is withheld where its own window is not covered
  by the vintage's reporting, rather than scored from a degenerate centre.
  An uncovered window let the centre saturate at the whole cumulative to
  the made date, or collapse to a point mass at zero, and either flattered
  the fits it was compared against. A test now pins that the baseline sees
  nothing after the date it is made on, and fails under three deliberate
  leaks.
- The in-report forecast validation subtracts a retrospective harmonisation
  from its new-count truth, as the release scoring already did. A frozen
  cut-off sitting before a listed break day and a current cut-off after it
  put records into the confirmed truth that were never notified that week.
  The cumulative rows are unchanged, the harmonisation being genuinely part
  of the reported total.

### Infrastructure

- The macOS test cell runs the platform-sensitive items only. It was the
  one cell skipping nothing, so it re-ran the sampling fits and the quality
  checks on the slowest runner, never finished, and was cancelled at
  GitHub's six-hour ceiling on every run. It now completes in about twelve
  minutes. The job also carries a timeout below that ceiling, so a stall
  reports as a failure rather than as a cancellation hours later.
- The automatic version increment deletes its branch and warns when GitHub
  refuses the pull request it opens, rather than leaving the branch behind.
  A leftover branch tripped the increment's own already-exists guard on
  every later push, so the first failure was loud and the rest were silent.
- Turing's compat bound admits 0.46.

### Documentation

- Rewrote the reduced-data-streams warning in the README and the summary
  page. The situation reports returned to their full format on 12 August,
  so the note now covers only what stayed frozen: the daily
  new-suspected-case count and the treatment-centre patient-movement
  series.
- Rewrote the abstract around the situation reports and the data they
  publish rather than a list of streams, in the wording the analysis
  report itself uses. It now covers the digitised symptom-onset curve and
  the recoveries, which the earlier text left out, and says that each
  release projects each DRC stream a week ahead and scores those
  forecasts against later data and a persistence baseline.
- A newly spotted surveillance signal now goes to an issue rather than a
  row in `data/README.md`'s scan table. The table is a short index, and
  rows added to it on a data update had to be reverted; the values still
  accumulate in `candidate_signals.csv` on the same pull request either
  way.
- The docs deploy job serialises on the `gh-pages` ref it writes to. The
  workflow's own concurrency group is keyed on the source branch, so two
  pull requests publishing previews at once raced and the loser's push was
  rejected.
- Gave one display block each to the tested BVD share and the positivity,
  the cumulative infections and the completed-detection-delay term, and the
  daily export intensity and its running sum, finishing the
  one-quantity-per-block pass. The export-death prevalence gets its own
  block rather than an inline definition.
- Split the death-pool BVD share and the assay positivity built on it into
  a display block each, matching the confirmed-case pipeline.
- Cut the symptom-onset reporting-delay methods section and its results
  text to the length and shape the other observation modules carry. The
  model definition, the priors and the limitations all remain. The
  identifiability discussion is condensed, and the prose runs one sentence
  per line like the rest of the document.

## v1.13.2

Changes since v1.13.1

### Data

- Advanced the model cut-off from SitRep 082 (4 August) to SitRep 089
  (11 August): added SitReps 083-089, advancing the confirmed-case and
  confirmed-death series and the recovered, isolation-occupancy and 24h
  laboratory-analysed streams, and digitising the symptom-onset curve
  vintages restored from SitRep 087 (the 084-086 brief-format reports carry
  a different, notification-week chart and so contribute no onset
  snapshot).
- From SitRep 084 (6 August) the INSP switched to a shorter "MVEBDB" brief
  format that no longer publishes the daily new-suspected-case stream
  (`suspected_daily_history` freezes at SitRep 083) or the treatment-centre
  patient-movement streams (frozen at SitRep 080), so the fit now depends
  on just the confirmed-case, confirmed-death, recovered and isolation
  streams rather than the full surveillance suite. The README and report
  carry a warning to this effect; see the inclusion rules in
  `data/README.md`.

## v1.13.1

Changes since v1.13.0

### Data

- Advanced the model cut-off from SitRep 081 (3 August) to SitRep 082
  (4 August): added SitRep 082 and the missing 3 August (SitRep 081)
  confirmed-case / confirmed-death cumulative points that the previous
  data update left at 2 August, and fixed the SitRep 081 onset-curve
  digitisation (the x-axis tick detection had undercounted that figure
  to −54%; the fix brings it inside the noise band) so both the 3 and
  4 August onset snapshots are now included.

## v1.13.0

Changes since v1.12.0

### Data

- Advanced the model cut-off from SitRep 079 (1 August) to SitRep 081
  (3 August): added SitRep 080 (2 August) and SitRep 081 (3 August) (#557).

### Model

- Widened the growth-rate `r` prior (eq 9) to
  `LogNormal(log(log 2 / 11.7), 0.40)`, a doubling time of 5.3--25.6 d at
  95%. The centre still matches the BEAST X reanalysis of 139 BDBV
  genomes [mbalaplacide2026](@cite), which reports an Exponential-growth
  doubling time of 11.7 d (95% HPD 6.8--17.5). That HPD is conditional on
  a single-rate coalescent. This is the assumption the Mongbwalu field
  epidemiology contradicts [kupferschmidt2026](@cite), the same evidence
  behind the `m` change in v1.12.0. An independent reanalysis of the
  earlier genomes puts the doubling time at 15.2--24.5 d
  [cuomodannenburg2026](@cite), which the old spread largely excluded.
- Corrected the growth-rate citation in the `rt_walk_model` docstring,
  which still named a generic 20 d molecular-clock estimate the model no
  longer uses.
- Fitted the digitised symptom-onset reporting triangle as a stream
  (`onset_reporting_model`), scored on its between-vintage increments
  through a discrete reporting-delay hazard that is nonparametric over the
  delay and drifts over calendar time. The loader collapses reprinted
  vintage figures by exact value equality, filters to the manifest cut-off,
  and builds between-vintage increment cells over a trailing 28-day window.
  Fitted by `bvd_joint` alongside every other stream and by a new
  `onsets_only_model` single-stream composer.
- Gave the onset stream an explicit ascertainment. The reporting hazard is
  normalised to reach one, so it carries the delay shape alone, and
  ascertainment is a separate level anchored on the confirmed pipeline's
  `p_drc · τ_test · positivity` averaged over the delay PMF, with a
  logit-scale offset and a slow walk over onset date. `onsets_only_model`
  falls back to a constant 0.15 anchor.
- Added a symptom-onset nowcast and forecast (`forecast_onsets`), splitting
  the coming week's new onset reports into the part arising from onsets that
  have already happened but are not yet reported (`onset_reports_backfill`)
  and the part from onsets that have not yet happened
  (`onset_reports_future`). `cumulative_onsets` moved from `bvd_joint` to
  the shared `_latent` submodel so every composer carries it. The reported
  increment is scored across releases as the `onset reports` stream.

### Fixes

- Widened the onset calendar-walk prior so it can follow reporting drift
  (#530).
- Bound the observation-scale slack below at one, preventing a fitted scale
  below the measurement floor from outvoting other streams in the joint fit.
- Fall back to four degrees of freedom on a degenerate Student-t nu argument
  instead of producing a Cauchy likelihood with no mean or variance.

### Documentation and infrastructure

- The analysis prose stated the confirmation-process sensitivity prior as
  `Beta(10, 1.76)` where `test_sensitivity_model` samples `Beta(38, 2)`, and
  gave the wrong rationale for it. Both the number and the rationale are
  corrected (#548).
- The analysis prose stated a 15-day TMRCA censoring SD where the model
  uses 16. Corrected to match the code (#548).
- `plot_rt`, `plot_cumulative_trajectories` and the single-stream overlay
  (`plot_stream_trajectories`) each drew a 50%-width credible band while an
  inline comment and docstring described it as part of a 30/60/90% ribbon
  trio. All three now compute a genuine 30/60/90% ribbon and drop the
  unused median-line description, so the figure and its documentation
  agree.
- `docs/src/contributing.md` documented an abstract-marker mechanism
  (`<!-- ABSTRACT:START/END -->`) that does not exist in `README.md` (the
  real marker is `<!-- SHARED:END -->`), eleven integrator, submodel and
  composer function names that do not exist in `src/`, two integration
  constants that do not exist either, and a surveillance dispersion prior of
  `truncated(Normal(0, 1); lower = 0)` where the code uses
  `truncated(Normal(0.6, 0.2); lower = 0)`. Rewrote the repository layout,
  running-and-testing and model-architecture sections against the current
  code and `Taskfile.yml`.
- Reflowed every write-up page, comment and docstring to the repository's
  prose rules: one sentence per line, shorter sentences, no run-ons, no
  restated points, no development history, and no code identifiers inside
  narrative prose.
- Removed emphasis-caps and issue and pull request references from `src/`,
  `test/`, `scripts/` and the report pages. That history now lives here
  instead.
- Added a carve-out to the no-issue-numbers rule in `contributing.md`. An
  issue number may stay where it is the provenance record for a still-open
  decision, rather than an account of what changed.
- Removed `predict_committed` and its three helpers from
  `src/counterfactual.jl`. They called `delay_convolution` and
  `DEATH_INTEGRAL_ALG`, neither of which is defined anywhere in `src/`,
  left behind by the rename to `convolve_delay` in `renewal.jl`. Nothing
  referenced them, so removing them changes nothing (#545).
- Removed eight unused one-off scripts. `scripts/Project.toml` gained a
  missing `Base64` dependency, and `scripts/README.md` now documents which
  Julia project each script needs.
- Removed the `notes/` folder, including the prose style guide the rules
  were derived from, which now lives in `contributing.md`. Its remaining
  open items are filed as issues #544 to #549.
- Added `AGENTS.md`, pointing at `README.md`, `contributing.md` and
  `scripts/README.md`, so a session starts with the rules in context.

## v1.12.0

Changes since v1.11.0

### Model

- Widened and shifted the prior on the cryptic-phase doubling count `m`,
  from `truncated(Normal(3, 3); lower = 0)` to
  `truncated(Normal(5, 4); lower = 0)`. Field epidemiology in Mongbwalu
  traced a sustained transmission chain to a death on 25 January 2026,
  with 500+ suspected cases between mid-January and mid-May
  [kupferschmidt2026](@cite), and the genetic TMRCA
  [mbalaplacide2026](@cite) is a lower bound on the outbreak age
  consistent with an origin that early. The new centre places the implied
  prior on the outbreak origin at the end of January (at the central
  11.7-day doubling, `m = 5` gives a cryptic duration of ~58.5 d, i.e. an
  origin around 30 January), with the wider SD letting the data and the
  genetic seeding bound pull it earlier or later (#533). The main fit's `m`
  prior is now set directly in `exponential_growth_model`; the backfill's
  advancing centre (`m_prior_centre`) is unchanged and tracked separately
  by #534.
- Corrected the growth-rate `r` prior prose (eq 9) to match the code and
  its actual source: the BEAST X reanalysis of 139 BDBV genomes
  [mbalaplacide2026](@cite) reports an Exponential-growth doubling time
  of 11.7 d (95% HPD 6.8--17.5), so the prior is
  `LogNormal(log(log 2 / 11.7), 0.28)`. The text had previously described
  an older, generic 20 d estimate that the model no longer uses.

## v1.11.0

Changes since v1.10.0

- Updated data processing and added digitised onset data
- Updated forecast evalaution but this remains highly experimental
- Added break days for confirmed cases and deaths due to a data harmoisation adding backdated data.

## v1.10.0

Changes since v1.9.0.

### Data

- Advanced the model cut-off from SitRep 055 (8 July) to **SitRep 064
  (17 July 2026)** across four batches: SitReps 056–057 (9–10 July),
  SitRep 058 (11 July), SitReps 059–061 (12–14 July, from the INSP source),
  and SitReps 062–064 (15–17 July, stepping over the unpublished 063).
  `as_of_date` → `2026-07-17`. Every fitted stream (confirmed cases/deaths,
  daily new-suspects, occupancy, bed capacity, recovered, 24h analysed,
  Tableau 6 treatment flows) is extended through each batch.
- Captured symptom-onset epidemic curves from the DHIS2 line list published
  in the analytique SitRep figures (059–062, 064). Digitised via
  `scripts/digitize_onset_curve.jl` (Julia) and a byte-identical Python port,
  and recorded as `data/onset_curve_scanned.csv`. Fitted as a stream, see
  Model below.
- Re-ran the confirmed/suspect in-care occupancy split on the resumed
  Tableau 7 data (SitReps 052–055), extending
  `treatment_confirmed_incare_history` and
  `treatment_suspect_incare_history` through 8 July (#413, revisits #373).
- Added `scripts/check_new_sitreps.jl` to detect stale manifests by
  comparing published INSP SitReps against the latest recorded vintage.

### Model

- Made the digitised symptom-onset reporting triangle a fitted stream
  (`onset_reporting_model`), scored on its between-vintage increments
  through a discrete reporting-delay hazard that is nonparametric over the
  delay and drifts over calendar time. Fitted by `bvd_joint` alongside
  every other stream and by a new `onsets_only_model` single-stream
  composer. The snapshot figure is a posterior predictive, so its band
  carries the measurement error the likelihood gives a digitised bar.

- Gave that stream an explicit ascertainment. The reporting hazard is
  normalised to reach one, so it carries the delay shape alone, and
  ascertainment is a separate level. It is anchored on the confirmed
  pipeline's own `p_drc · τ_test · positivity` averaged over the delay
  distribution, with a logit-scale offset and a slow walk over onset date.
  Not yet validated. The ascertainment walk and the reproduction-number
  walk share the onset axis and are both least constrained late, so $R_t$
  over the final fortnight and $C_T$ need reporting either side of this
  change.
- Added a symptom-onset nowcast and forecast (`forecast_onsets`),
  splitting the coming week's new onset reports into the part arising from
  onsets that have already happened but are not yet reported and the part
  from onsets that have not yet happened. The latent onset trajectory
  `cumulative_onsets` moved from `bvd_joint` to the shared `_latent`
  submodel so every composer carries it. The reported increment is scored
  across releases as the `onset reports` stream; the triangle's cumulative
  level is deliberately not scored, since every vintage rereads the whole
  figure and its ≈4% per-scan level error revises the total both ways.

- Added a results figure reading the triangle along the onset date: the
  latest digitised bar for each onset date against the posterior predictive
  for that bar and against the modelled onsets. The nowcast and forecast
  figure separates the latent components from the replicate the next
  vintage will print, since only the replicate is scored.

- De-boxed the anonymous `map(do)` closures on `bvd_joint`'s default
  log-density path (`confirmed_positivity_link = :composition`): extracted
  the composition-positivity loop to a plain `composition_positivity()`
  function and replaced three occupancy/split `map(1:n) do t …` blocks with
  elementwise broadcasts. Pure type-stability tidy-up; Mooncake gradient
  is bit-identical (same `hash(g)`). Clears the **layer-1 prerequisite**
  for the opt-in Enzyme reverse-mode backend — Enzyme cannot construct a
  shadow for anonymous closures that capture conditionally-scoped variables
  (they get boxed in `Base.RefValue`), and after this fix reaches LLVM's
  `nodecayed_phis!` pass (next blocker tracked in #445; #446).

### Fixes

- Restored the renamed clock-sensitivity chain reference (`chn_joint_fast_clock`
  → `chn_joint_exp_growth_clock`) in analysis diagnostics after a merge
  inadvertently reverted it to the old, undefined variable. This had broken the
  `Documenter → Render analysis` step on `main` when `BVD_RUN_SENSITIVITY=true`
  (#418).

### Documentation and infrastructure

- Added `data/README.md` documenting the two data sources, the direct-INSP
  fetch workflow via the WordPress REST API, the per-SitRep field checklist,
  and the drift-check protocol.
- Added Julia (`scripts/digitize_onset_curve.jl`) and Python
  (`scripts/digitize_onset_curve.py`) digitisation scripts for the symptom-onset
  epidemic curves, producing byte-identical CSVs. Python deps managed via `uv`
  (PEP 723).

### Dependencies

- Widened the Turing compat from `0.45` to `0.45, 0.46`, allowing the
  upstream package evolution while keeping 0.45 pinned as the working default
  (#415, #416).
- Updated compat bounds for ADTypes (1.22.2), CairoMakie (0.15.13),
  Distributions (0.25.129), LogDensityProblems (2.2.0), Integrals (4, 5.4),
  StatsFuns (2.2.0), SHA (0.7.0), Aqua (0.8.16), TestItemRunner (1.1.5),
  CensoredDistributions (0.2.22), and TestItems (1.0.0).

## v1.9.0

Changes since v1.8.0.

### Data

- Updated data including movement flows for those in isolation.

### Model

- Updated the molecular-clock time estimate to use the outbreak-specific
  BEAST X analysis (mbalaplacide2026, 139 BDBV genomes from 16 health
  zones, ~1.1E-3 subs/site/year). The genetic TMRCA baseline moves from
  2026-03-25 (SD 15, fixed 1.2E-3 EBOV rate) to **2026-03-15** (SD 16,
  95 HPD 09 Feb--12 Apr) under the Skygrid coalescent prior. Added a
  tree-prior sensitivity comparing the Exponential growth estimate
  (2026-03-08, SD 16, 95 HPD 01 Feb--05 Apr).
- Updated the growth-rate prior to the outbreak-specific doubling time
  from the same analysis: centre moves from 20 d (Cuomo-Dannenburg &
  Ghafari) to **11.7 d** (95 HPD 6.8--17.5, mbalaplacide2026), with the
  log-SD widened from 0.15 to 0.3 to match the wider credible interval.

### Report and forecasts

- Saved the one- to four-week-ahead forecasts as a release asset
  (`forecast.csv`), plus the one-week-back validation forecast
  (`forecast_validation.csv`), at each results release. Past forecasts were
  shown in the report but never stored, so they had to be reconstructed by
  re-running each past release's own code on that release's data; these are
  published separately as a backfill release. Each release now records the
  forecast it made, so it can later be scored against what is observed. This
  underpins the new cross-release forecast scoring (CRPS, log-scale CRPS,
  coverage, bias, and relative skill against a persistence baseline) shown
  on the sensitivity page.

### Documentation

- Updated all prose references to the genetic bound and growth-rate prior
  to cite the new virological.org report (v1045, mbalaplacide2026) and
  the outbreak-specific rate and doubling time.
- Added a tree-prior sensitivity section comparing the Skygrid and
  Exponential growth TMRCA estimates.

### Performance

- Halved the default post-warmup draws in `nuts_sample` from 2000
  (2 chains x 1000) to 1000 total (2 chains x 500), and moved the docs
  fit registry (`build_fit_specs`, `fit_key`, `fit_content_hash`) to the
  same 500 x 2 setting. This roughly halves the sampling wall-clock at the
  cost of some effective sample size.
- Trimmed the default NUTS warmup cap in `nuts_sample` from
  `min(250, samples ÷ 2)` to `min(200, samples ÷ 2)`, so at the new
  `samples = 500` default each fit runs 200 adaptation steps rather than
  250. Shortens warmup by a further ~20% per fit.

## v1.8.0

Changes since v1.7.0.

### Model

- Credited the repeat-control confirmation process in the confirmation
  sensitivity prior. Rule-out is investigative rather than a single negative
  PCR, so the effective sensitivity is higher than one assay draw (two controls
  give about 0.98). The headline `test_sensitivity_model` prior moves from the
  single-assay `Beta(10, 1.76)` (mean 0.85) to `Beta(38, 2)` (mean 0.95) on the
  confirmed and confirmed-deaths streams. The outbreak-size estimate is robust
  because the sensitivity enters the multiplicative ascertainment ridge
  (`p_drc · s_test · τ_test`); the ascertainment posterior re-centres (resolves
  the retesting/rule-out part of #374).
- Grounded the confirmed onset-to-sample delay on the NEJM DRC 2026 cohort
  (Akilimali et al.) as a standard part of the joint model. The onset-to-report
  and report-to-receipt legs already convolve to onset-to-sample for confirmed
  cases. The convolution's mean (the sum of the report and receipt leg means)
  and median (Wilson-Hilferty of the summed leg variances) are fitted to the
  reported 7.4 d and 4.8 d as soft Normal observations, with the reported 95%
  credible intervals as their SDs. This grounds the otherwise-unidentified
  laboratory-turnaround delay directly from the cohort's own uncertainty and
  adds no latent parameter. On by default in the joint fit; the single-stream
  and isolation fits lack the laboratory receipt leg and so do not carry it
  (resolves #359).

- Scored the confirmed-case laboratory positives as an overdispersed
  `BetaBinomial` of the observed analysed denominator instead of a plain
  `Binomial`. A plain `Binomial` on denominators of several hundred
  specimens gave posterior-predictive intervals far too tight, so the
  confirmed stream was systematically under-covered: the smooth pooled /
  composition-linked per-window positivity does not capture the day-to-day
  laboratory batching and within-window positivity heterogeneity the
  confirmed counts carry. A single intra-window overdispersion `ρ`
  (`confirmed_overdispersion_model`, weakly-informative `Beta(1, 24)`)
  inflates each window's variance to `n·p·(1 − p)·(1 + (n − 1)·ρ)`,
  identified across the laboratory windows, and recovers the `Binomial` as
  `ρ → 0`. The mean structure and the composition link that identifies the
  background `λ_bg` are unchanged.

### Documentation and infrastructure

- Excluded the published `released_estimates.csv` overlay from the fit
  content hash, so the render jobs reuse the matrix fits instead of
  refitting every model. The render step rewrites that overlay before
  rendering, which changed the data-tree digest and busted every fit
  cache key; the overlay feeds only the estimate-evolution figure and is
  not a fit input. `tree_sha256` and `content_hash` gained an exclude for
  non-input data files, and genuine data changes still refit.

## v1.7.0

Changes since v1.6.0.

### Data

- Added the situation-report `Tableau 6` treatment-centre patient-movement
  flows (CTE/CT/CI) as optional daily streams: admissions, in-care deaths,
  rule-outs and absconded patients (13–23 June). Each stream is resilient: an
  empty history is a no-op, so the model degrades to the occupancy backbone
  where a flow is not reported. Advanced the data through situation report
  046 (29 June).

### Model

- Reworked the isolation submodel into a treatment-centre flow model that
  fits the Tableau 6 flows alongside occupancy. The bed length-of-stay is an
  outcome mixture, with the death and recovery branches weighted by an in-care
  case-fatality `CFR_iso = logistic(logit(CFR) + β_iso)` — a reported modifier
  on the infection case-fatality, identified by the in-care death flow rather
  than estimated independently. The daily discharge flows are scored as
  optional negative-binomial streams (resolves #338).
- Added a manual, opt-in occupancy reclassification-break offset. Break days
  are listed explicitly in `[occupancy_break_dates]`, replacing the removed
  threshold detector that flagged too many days. A level step is fitted into
  the modelled occupancy mean at each listed day, so the fit tracks a known
  between-report measurement-basis discontinuity without bending Rt to chase
  it. Each step is centred on zero, so the fit partitions it into reporting
  artifact vs real demand. The 19 June DHIS2 reclassification (occupancy 416 →
  361, confirmed by the au-lit start-of-day stock) is listed; the joint fit
  with the isolation stream had been bending Rt up and down to chase this and
  the later missing-SitRep steps, which no single-stream fit shows.

### Report and forecasts

- Added one-week-ahead forecasts for the treatment-centre admissions, in-care
  deaths and rule-outs.
- Added a per-stream calibration plot (with the calibration table kept in a
  collapsible block), and posterior-predictive panels for the four flow
  streams. The treatment-centre flow methods section was rewritten, and the
  flow streams added to the data-overview table.
- Added a comparison of the confirmed-case projection against Chamla et al.
  as a second external comparator, forward-projected from a dedicated frozen
  fit at their 8 June confirmed-case calibration anchor. This carries the
  confirmed-case testing history, replacing the poorly-identified 27 May proxy
  that had effectively no testing data. The 8 June fit also shows as a vintage
  in the estimate-evolution overlay (resolves #340, #349).
- Refreshed the released-estimate evolution overlay to v1.6.0 and refresh it
  automatically in continuous integration before each documentation deploy.
  Dropped the per-release current-model re-fits from the estimate-evolution
  plot, keeping the matched-McCabe cut-offs and the one-week-back validation
  fit (resolves #341).
- Tightened the reproduction-number plot y-axis to 1.2 times the 90% upper
  bound so the credible band is legible (resolves #342).
- Added a posterior correlation heatmap across the key estimates (outbreak
  size, reproduction number, outbreak age, CFR, ascertainment, background,
  fraction tested and each stream's expected total) and a pairs plot of the
  per-stream modelled totals against each other and the observed value, so
  the size-versus-ascertainment trade-off and per-stream over/undershoot are
  visible in one place (resolves #346).

### Performance

- Halved the NUTS warmup: `nuts_sample` now defaults `n_adapts` to
  `min(250, samples ÷ 2)` (250 adaptation steps at the standard 1000
  draws) instead of Turing's `min(1000, samples ÷ 2)` (500), cutting the
  discarded warmup iterations on every report fit. Pass `n_adapts`
  explicitly to override.
- Centred the per-stream pooled negative-binomial dispersion
  (`pooled_dispersion_model`) instead of drawing it non-centred, and made
  centred the default. The surveillance streams are data-rich, so the
  non-centred pooling funnelled as `τ → 0`; the centred form is an exact
  reparameterisation that removes the funnel (resolves #352).
- Put the bed-capacity baseline `C0` on the log scale (`LogNormal(log 450,
  0.42)` in `bed_capacity_model` and `bed_capacity_walk_model`) instead of a
  truncated normal, so the whole capacity `C(t) = C0·exp(walk)` is
  log-consistent with no hard boundary, improving the worst-mixing capacity
  block (resolves #358).

### Fixes

- Relaxed the `tau_death` test assertion to non-negativity: the joint exposes
  the realised cut-off death-testing intensity (analysed over suspected, a
  diagnostic computed independently of the death volume), which is not a
  probability and can exceed one in a backlog regime.
- Widened the per-stream dispersion pooling-SD prior `τ` in
  `pooled_dispersion_model` from `HalfNormal(0.3)` to `HalfNormal(0.6)`,
  resolving a prior-data conflict where the posterior `τ` sat entirely above
  the old prior's tail because the stream dispersions genuinely span ~9×
  (resolves #336).
- Fixed the one-week-ahead forecast to project new counts over the horizon
  and add them to the cut-off cumulative, instead of scaling the cumulative
  stock by `exp(r·horizon)`, which made a below-one reproduction number
  imply an impossible shrinking cumulative in the Chamla comparison
  (`cases_cum`, `deaths_cum`, `confirmed_cum`, `confirmed_deaths_cum` and
  `recovered_cum`; resolves #351).

### Documentation and infrastructure

- Split the report into two literate pages rendered from a shared setup: an
  analysis page (methods, results, one-week-ahead forecast) and a sensitivity
  page (forecast validation, per-stream outbreak size, estimate evolution,
  McCabe and Chamla comparisons, delay and clock sensitivity), so the deploy
  no longer renders one 4.3k-line file in a single job (resolves #364).
- Fanned the documentation build into a `list → fit → render → combine` CI
  grid: one content-addressed, cached NUTS fit per matrix job, the two pages
  rendered in parallel from the cached chains, and a combine job that deploys
  and publishes, dropping the critical path from all fits serialised behind
  one runner to roughly the slowest single fit.
- Pinned the shared workspace-root `Manifest.toml` as the artifact the docs
  grid uploads and restores, since `docs/Project.toml` is a workspace member
  and resolves the root manifest rather than `docs/Manifest.toml`, so the
  `fit`, `render` and `combine` jobs resolve the same package set (resolves
  #368).
- Set `include-matrix: false` on the `render` job's env-cache step so the two
  render jobs restore the precompiled depot shared by `list` and `fit` instead
  of keying on the matrix `page`, which missed the shared cache and made each
  render re-precompile the whole stack (~530 deps) before rendering.

## v1.6.0

Changes since v1.5.0.

### Model

- The isolation BVD treatment length-of-stay uses the BDBV line-list
  admission-to-death delay as its prior.
- The reproduction-number random walk starts a month before the first
  situation report (`RT_WALK_LEAD = 28`, exposed as the `bvd_joint` keyword
  `rt_walk_lead`), so `R_t` can move over the weeks of transmission
  leading up to that report. The walk start is floored at the renewal start,
  and the `plot_rt` reconstruction uses the same knot grid.
- Added a supply-limited isolation/treatment-bed stream ("Patients en
  isolement"), the renewal analogue of the convolution secondary-observation
  model of EpiNow2.
  Bed occupancy may be supply-driven (demand can outstrip supply), so the model
  fits a latent bed demand, the suspect inflow carried through a length-of-stay
  survival (BVD cases with a sampled treatment stay, non-BVD suspects leaving
  after a sampled rule-out stay), right-censored at an effective bed capacity
  `ρ·C(t)` (a censored negative binomial). The capacity `C(t)` is a random walk
  (`bed_capacity_walk_model`) that tracks the beds being added and can be
  projected forward, pinned by the implied bed count (reported occupancy /
  "Taux d'occupation" rate) on the days a rate is published.
  The stream exposes the bed demand, occupancy, capacity, shortfall and
  utilisation, and carries its own observation dispersion. Added
  `convolve_survival`, the `treatment_admission_model` observation submodel,
  the `isolation_admission_model`, `bed_capacity_model` and
  `bed_capacity_walk_model` priors and the `treatment_only_model`
  single-stream composer (resolves #265).
  This is a single national model, so it cannot represent local bed saturation
  (Ituri at 93.9% occupancy on 13 June against Sud-Kivu 21.9%); the national
  shortfall understates the local unmet need.
- Gave the non-BVD isolation rule-out stay its own sampled length-of-stay
  (`ruleout_los`) in `treatment_admission_model`, separate from the
  report-to-receipt laboratory delay, so the occupancy identifies the rule-out
  stay on its own clock and the lab-turnaround delay is set by the testing,
  composition and confirmed-death streams. Exposes
  `isolation_ruleout_los_mean`.
- Added a recovered-among-confirmed stream ("cumul guéris"), the
  secondary-observation incidence analogue: survivors among the modelled
  daily confirmed cases, scaled by a recovery proportion and convolved with a
  confirmation-to-recovery delay, with its own observation dispersion. The
  recovery proportion is grounded on the case-fatality ratio (a recovered
  case is one that did not die) with a log-odds adjustment for the confirmed
  population, rather than estimated independently.
  The confirmed model exposes one daily confirmed-case series
  (`confirmed_daily`) that both the recovered stream and the cumulative-
  confirmed trajectory reuse. Added the `recovered_model` submodel and the
  `recovery_probability_model` prior.
- Added the `exports_joint_only_model` composer, which fits the Uganda export
  cases and deaths together over the one travel-gated at-risk prevalence. The
  single-stream comparison in the walkthrough now shows one joint "exports"
  fit instead of separate export-case and export-death fits.
- The one-week-ahead forecast also projects the isolation/treatment beds: the
  bed demand a week ahead (need under unconstrained supply) and the
  supply-limited occupancy it produces, whose gap is the projected bed
  shortfall, and the cumulative recovered total, each replicated with its
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
  lognormal random walk (`background_walk_model`): a per-day background with no
  reporting-vintage steps, gated to begin a report-to-receipt lead before the
  first suspected-case report, shared across the suspected-case and
  suspected-death streams, with a half-normal baseline and a tight random-walk
  innovation SD.
  This also removes the per-vintage step in the modelled cumulative-death
  trajectory.
- Widened the non-BVD background level prior so the laboratory positivity
  (210/755 ≈ 0.28 positive) identifies it.
  The suspect pool is inferred to be a minority BVD, which lowers the
  cumulative-infection estimate (`C_T`) with a wider credible interval.
- Gated the laboratory analysed-specimen capacity to the testing onset, so no
  specimens are modelled as analysed before testing existed.
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
- Quantified the per-vintage posterior-predictive checks with a per-stream
  calibration table (`stream_calibration`): the mean forecast bias and the
  empirical 50%/90% interval coverage of each stream's one-step-ahead
  conditional predictive, so the streams the joint fit reproduces less well
  can be read off rather than eyeballed. Added the `bias_sample` scoring
  helper (resolves #269).

### Documentation

- Added a one-page [Summary dashboard](@ref) for readers with limited time:
  the headline estimates as prose and tables alongside the reproduction
  number, infections-over-time and modelled-versus-observed reported-case
  figures. It reuses the artifacts written by the analysis build rather than
  re-fitting, so it refreshes whenever the data updates.

### Outputs

- Added the latent symptom onsets (the "symptomatic cases" outcome) to the
  shared posterior outputs.
  `posterior_draws.csv` gains a `cumulative_onsets_T` column, the cumulative
  symptom onsets by the cut-off per draw (the onset analogue of `C_T`), and a
  new `onsets_over_time.csv` records the daily new and cumulative onset
  trajectory over time with 30/60/90% credible intervals.
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
