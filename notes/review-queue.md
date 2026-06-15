# Report review queue (renewal branch)

Status of every item raised in the live report review.

## Done (committed)
- Removed Cuomo-Dannenburg doubling-time aside from the Rt subsection.
- Moved derived-R0 (1.31-1.55) note to the growth subsection.
- Reworded step-size prior: Rt unlikely to change >~10%/week (2 SDs).
- Generation interval: "We assume g is a Gamma"; cut serial-interval caveat.
- "From that seed we assume…"; m prior centred on three doublings N+(3,3).
- Discretisation maths -> double-interval (primary-event) censoring integral.
- "confirming cases via laboratory testing"; de-garbled sensitivity sentence.
- Cut renewal-start repetition in the infection-process subsection.
- Dropped r0 from the infection summary table (kept current-growth r).
- Added Rt-walk sigma + intervention effect to the infection pair plot.
- Removed the export pair plot, then the whole export-parameters section.
- plot_rt: y-axis capped at 1.2x the upper 90% band.
- Lab-testing-model report written (notes/lab-testing-model-report.md).
- Maths-presentation audit roadmap (notes/maths-presentation-fixes.md).

## In flight (background agents)
- Delays -> natural-parameter Gamma priors from the bdbv-linelist submodule
  (onset-to-report, onset-to-hospitalisation, onset-to-death convolution);
  incubation + receipt stay. Worktree.
- Plots cleanup: remove median/central lines EVERYWHERE; fix estimate-evolution
  to plot cumulative infections over time (not a constant ribbon), align the
  release dates, add dotted per-release vertical lines. Worktree.

## Needs a decision (blocks lab work)
- Lab model bundle: A = composition-link headline (lab identifies λ_bg, size
  ~5000, strong assumption) vs B = free-positivity headline (positivity follows
  the test data, λ_bg weak, size ~2280 deaths/exports-pinned) vs hybrid. See
  notes/lab-testing-model-report.md. Mutually exclusive for the headline.

## Queued model changes (then one joint refit)
- Lab do-anyway (either bundle): lag the suspect-pool composition by the
  receipt delay (Q2 bug); switch on the suspected-death background λ_bg_death
  and align the death lab model to the chosen positivity philosophy; make the
  tested fraction τ_test time-varying (Q1).
- Export deaths: apply the CFR to the exported cases BEFORE ascertainment
  (deaths would be reported), not the detected/ascertained exports. Correct my
  earlier "100% ascertainment" report wording too.
- McCabe comparison: show our modelled cumulative cases AT the McCabe report
  dates (18/20 May) read off our time series, not the current cut-off C_T.
- Delay sensitivity: verify the section uses the corrected delays after the
  natural-params change.

## Queued presentation (after model changes)
- Apply the 9 maths-presentation fixes (sentence-maths style, split the NegBin
  likelihoods over lines, display the composition equations).
- Cleaner chain names: turn off submodel prefixing / FlexiChains rename
  (e.g. rt_state.sigma_rw -> sigma_rw) for tables and pair plots.

## Final
- One joint refit, then rebuild + host the docs.
