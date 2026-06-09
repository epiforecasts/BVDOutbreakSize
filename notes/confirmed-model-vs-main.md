# Confirmed cases and confirmed deaths: renewal vs main

How the renewal branch models the laboratory-confirmed case and death
streams, and how that differs from `main`. Written for review alongside the
experimental "main-like confirmed deaths" PR (issue #225).

## Shared starting point

Both branches link laboratory positivity to the composition of the suspect
pool rather than fitting a free positivity curve. The tested-positive
probability is

```
p = s·q + (1 − spec)·(1 − q)
```

where `q` is the BVD share of the relevant suspect pool, upsampled by a
decaying severity enrichment (more severe, more Ebola-like suspects are
tested first), `s` is the assay sensitivity and `spec` the specificity. The
false-positive term `(1 − spec)(1 − q)` is what lets the positivity data
identify the non-BVD background `λ_bg` instead of leaving it free. Both
branches sample `s` and `spec`.

## Confirmed cases

The two branches share the composition-link positivity above. They differ
in how the analysed-specimen volume — the denominator the confirmed counts
are drawn from — is handled, especially for the post-28-May "dark" vintages
where INSP published confirmed counts but no national analysed denominator.

- **main** wraps the positivity in an explicit laboratory-throughput model
  (#190 + #200): a capacity-limited processing queue (centre ≈150 tests per
  day), an optional epi-exclusion fraction, a forward delay, and a
  per-window positivity random effect. The dark vintages are fitted *through*
  this modelled capacity — the lab drains a received backlog at a modelled
  rate, so the unobserved denominators are inferred from the throughput
  dynamics.

- **renewal** is deliberately lighter. Received specimens are the testing
  fraction times the receipt-delayed suspect pipeline; the observed windows
  are a Binomial on the *real published analysed denominator*, and the dark
  windows are a NegBinomial on the *modelled* received volume times the
  positivity. There is no capacity queue, no forward delay, no per-window
  positivity random effect.

The renewal choice is for stability, not simplicity for its own sake. The
constant-capacity / FIFO and per-vintage queue variants basin-split in the
renewal joint (recorded R-hat ≈2.1, ESS ≈2), whereas binomial-on-analysed
plus the composition link converges cleanly (R-hat ≈1.006). The dark-window
denominator is handled by a zero-degrees-of-freedom deterministic volume
extrapolation rather than a stochastic capacity walk, which also converges.

## Confirmed deaths

This is the larger divergence and the subject of the experimental PR.

- **main** (#191, folded into #200) models confirmed deaths as a genuine
  death-side laboratory process. The positivity
  `p_pos_death = s·q_death + (1 − spec)(1 − q_death)` is built on the
  *death-pool* composition `q_death = bvd_deaths / (bvd_deaths + bg_death)`,
  with the **same assay `s`/`spec` shared from the case stream**, and the
  expected confirmed-death increment is `τ_death · p_pos_death · ΔN_death`, a
  forwarded fraction of the suspect-death *backlog* increment. Scored
  per-edge NegBinomial.

- **renewal** is lighter. The confirmation probability is
  `p_death_conf = logistic(logit(q_susp) + log(m_death))`: it takes the
  **case-pool** composition `q_susp` (BVD share of suspected *cases*) and
  applies a free odds-enrichment `m_death`, then thins the modelled
  suspected-*death* daily trajectory directly. There is no death-pool
  composition, no assay `s`/`spec` on the death stream, and no backlog
  forwarding. Scored per-vintage NegBinomial.

The renewal rationale (now corrected in the docstring): the death-pool
composition `q_death` collapses to 1 whenever the suspect-death background is
off, which is renewal's default, so a death-side composition link would be
degenerate; the case background `λ_bg` is always on, so `q_susp` stays
informative, and the enrichment anchors there. The original docstring also
argued the renewal case lab carried no `s`/`spec` to share — that reason has
lapsed, because the renewal confirmed-case lab now *does* sample `s`/`spec`
through the composition link. So the main obstacle to porting main's death
lab process is now only the `q_death` degeneracy, which needs the
suspect-death background turned on. That is exactly what the experimental
PR tests.

## Notes found while auditing

- **Stale docstrings fixed.** Six model docstrings referenced superseded
  facts (the `s`/`spec` rationale above; several values frozen at the 28 May
  / 3 June cut-off rather than 5 June). All corrected to be cut-off-agnostic.
- **Dead integral-lineage code.** `src/counterfactual.jl`'s
  `predict_committed` and its `_committed_*` / `τ_forward` / `τ_death`
  helpers are leftover from the integral model and are not called anywhere
  in the renewal source. This is why main's "committed-totals" counterfactual
  was not carried over: the underlying machinery is non-functional on this
  branch. Worth removing or rebuilding as a follow-up.

## Implications

Renewal's confirmed-case lab is a deliberate, well-converging
simplification of main's throughput model and is unlikely to want changing.
Renewal's confirmed-death model is the weaker of the two: it borrows the
case composition and a free enrichment rather than modelling the death
specimens. Now that `s`/`spec` are available to share, porting main's
death-side lab process is feasible if the suspect-death background is
enabled — the experimental PR measures whether that improves the fit or
reintroduces the degeneracy.
