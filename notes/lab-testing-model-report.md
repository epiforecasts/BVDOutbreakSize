# Laboratory testing model: decision report (renewal branch)

## Current model: exact generative structure

All file:line references are in the `renewal` worktree.

**Shared inputs.** `bvd_joint` (`src/models/joint.jl:248`) runs the renewal
once, stages onsets, then samples the shared lab/positivity parameters inside
`reported_cases_model` via `test_positivity_model` (`src/models/priors.jl:450`):
the scalar non-BVD background `λ_bg ~ Normal⁺(0,1)` and the tested fraction
`τ_test ~ Beta(5,2)`. The suspected case stream (`src/models/observations.jl:293`)
computes the unit-ascertainment BVD onset-to-report series `bvd_reports_daily`
(line 310) and the per-day background `bg_daily` (lines 322-331, constant
`λ_bg` by default). These two daily series, plus `p_drc` and `τ_test`, are
handed to `confirmed_cases_model` (joint.jl:329-334).

**Received-specimen volume (lab denominator driver).** In
`confirmed_cases_model` (observations.jl:556):
- `suspected_daily = p_drc·bvd_reports_daily + bg_daily` (line 579).
- `received_daily = τ_test · convolve_delay(suspected_daily, receipt_pmf)`
  (line 580), with the receipt delay sampled by `lab_delay_model`
  (priors.jl:564, mean ≈ 4.5 d).
- The between-vintage increments of `received_daily` are fit against
  `tests_received_history` as NegBinomial counts (lines 591-593). This
  identifies `τ_test` and the receipt delay.

So the modelled testing volume is **entirely a deterministic function of the
ascertained suspected-case pipeline** (`p_drc·bvd + λ_bg`) thinned by a scalar
`τ_test`. There is no latent testing-volume process. Note `data.jl:99-101`:
the **analysed** series is the Binomial denominator; the **received** series is
what `received_daily` is fitted to.

**Confirmed positives, three window groups** (`confirmed_positivity_windows`,
observations.jl:419; consumed at 595-709):
- *Observed* windows: confirmed increments scored as
  `Binomial(obs_analysed, p_pos)` against the **observed analysed denominator**
  (line 688).
- *Early* windows (before first lab date) and *late* windows (after last lab
  date): no observed denominator, so scored as
  `NegBinomial(p_pos · modelled_volume, k)` where the volume is binned from
  `received_daily` (lines 677-682, 699-709). This is already "positivity ×
  modelled received volume".

**Positivity link** (lines 619-672). Default `:composition`:
- Suspect-pool composition over each window
  `φ = (p_drc·bvd)_window / ((p_drc·bvd)_window + bg_window)`
  (lines 636, 650-655), built from the **same suspected-stream daily series**.
- Severity enrichment `δ_i = δ0 · exp(−c_window[i]/decay)` on a
  cumulative-analysed-volume clock (lines 640-646, 656;
  `severity_enrichment_model` priors.jl:635).
- Tested BVD share `q = logistic(logit(φ) + δ_i)` (line 660).
- Positivity `p_pos = s·q + (1−spec)(1−q)` (line 663; `s` priors.jl:532,
  `spec` priors.jl:549).

The false-positive term `(1−spec)(1−q)` is what makes confirmed counts respond
to the non-BVD share and thereby identifies `λ_bg`. Alternative `:free` link
(line 670) uses `confirmed_positivity_model` (priors.jl:590): a
partially-pooled per-window random-effect positivity decoupled from `λ_bg`.

**Confirmed deaths** (`confirmed_deaths_model`, observations.jl:957).
Composition-linked odds enrichment:
`q_susp = p_drc·Σbvd / (p_drc·Σbvd + Σbg_daily)` (line 969, a **single scalar**
over the whole grid, not windowed),
`p_death_conf = logistic(logit(q_susp) + log m_death)` (line 971, `m_death`
priors.jl:668), then `confirmed_death_daily = p_death_conf · deaths_daily`
scored as NegBinomial increments (lines 980-985). There is no
received→tested→positive pipeline for deaths and no shared `s`/`spec` here.

**Background deaths.** A suspected-death background `λ_bg_death` exists
(`death_background_model`, priors.jl:410; wired in `deaths_model`,
observations.jl:218-256), **but it is off by default in the joint**:
`bvd_joint` passes neither `death_background` nor a `background_re` unless
`background_re=true` (joint.jl:312-325, default `false`). So the default deaths
stream is pure BVD with `λ_bg_death = 0` (observations.jl:252-256), which
collapses any death-pool composition to 1, which is precisely why the death
stream is grounded on the case `q_susp` instead (docstring
observations.jl:937-955).

---

## Q1. Testing scale-up / latent testing volume V_t

**Now.** Testing volume is
`received_daily = τ_test · convolve(p_drc·bvd + λ_bg, receipt)`
(observations.jl:580). It is rigidly tied to ascertained incidence: volume can
only ramp as modelled suspected cases ramp, scaled by a single constant
`τ_test`. A real lab scale-up (capacity coming online, a testing campaign) that
does not track incidence cannot be represented. The objection is correct: the
denominator is functionally pinned to ascertained cases.

**Options.**
- **(1a) Latent volume RW.** Replace `received_daily` with its own latent
  log-random-walk `V_t`, fit to `tests_received_history` directly, decoupled
  from incidence. Positives in dark windows become `p_pos · V_t`.
- **(1b) Time-varying fraction τ_test,t.** Keep the incidence-driven structure
  but make `τ_test` a smooth RW/RE instead of a scalar, so the *fraction*
  tested scales up over the window. Cheaper, fewer new latents, still anchored
  to incidence shape.
- **(1c) Keep as-is** (constant `τ_test`).

**Trade-offs.** (1a) captures genuine ramp-up but adds a full latent series
identified only by the received-specimen vintages (few, noisy). It also
**severs the composition link's volume clock**: the severity decay
`exp(−c/decay)` (line 656) reads `c` off `received_daily`; if volume is a free
latent, the decay clock floats, and more importantly the dark-window positives
`p_pos·V_t` no longer constrain incidence at all, weakening the lab stream's
contribution to `C_T`. (1b) is the smallest change that answers the objection
while keeping volume incidence-shaped: it lets the tested fraction climb
without inventing volume from nothing, and `λ_bg`/composition identification
survives because the suspect pool still drives the shape. (1c) is the status
quo.

**Recommendation.** Prefer **(1b)** — a time-varying `τ_test,t` (tight RW or
low-dimensional ramp) fit to the received series — as the principled, low-risk
answer to "testing scaling up". Reserve (1a) for a sensitivity analysis only,
because a free `V_t` decouples the lab stream from incidence and removes most
of its identifying power. Either way, fit the **received** series for volume
and keep the **analysed** series as the Binomial denominator where observed,
since those are the real quantities.

---

## Q2. Suspect-pool composition timing (the suspected bug)

**Finding: this is right.** In `:composition` mode the pool composition `φ` is
built from `bvd_window = bin_increments(p_drc·bvd_reports_daily, window_days)`
and `bg_window = bin_increments(bg_daily, window_days)`
(observations.jl:636-637, 650-655). Both `bvd_reports_daily` and `bg_daily` are
indexed at **report/onset time** — they are the suspected-case daily series,
*not* carried through the receipt delay. The only series that passes through
the receipt delay is `received_daily` (line 580), and that is used only for the
volume and the decay clock, never for `φ`.

So a specimen received and tested on day t has its BVD share computed from the
suspect-pool composition *as it was at report time t*, with no lag for the
~4.5-day receipt delay. During a growing outbreak the BVD share rises over
time, so evaluating `φ` un-lagged overstates the BVD fraction of what is
actually on the bench, biasing positivity high and `λ_bg` low. This is a
genuine timing inconsistency: volume is lagged, composition is not.

**Options.**
- **(2a) Lag both numerator and denominator by the receipt delay.** Form
  `received_bvd_daily = τ_test·convolve(p_drc·bvd, receipt)` and
  `received_bg_daily = τ_test·convolve(bg_daily, receipt)`, then
  `φ_window = bin(received_bvd)/bin(received_bvd+received_bg)`. The ratio
  `τ_test` cancels, so `φ` becomes the delayed (received-time) composition.
  Cheap and exactly consistent with the volume.
- **(2b)** Leave un-lagged (status quo) and document the approximation.

**Recommendation.** Implement **(2a)**. It is the small correct fix, reuses the
receipt PMF already sampled, keeps the `λ_bg` identification mechanism intact,
and removes a directional bias on `λ_bg`. This is the cleanest standalone
change and should be done regardless of the larger Q3 decision.

---

## Q3. Positivity as a random effect vs the composition link (the crux)

**Now.** Default `:composition` derives `p_pos` from `φ` (suspect-pool BVD
share) + severity enrichment + assay `s`/`spec` (observations.jl:619-668). The
`(1−spec)(1−q)` false-positive term is the only thing tying confirmed counts to
the non-BVD share, hence to `λ_bg`. The `:free` link (line 670, priors.jl:590)
already implements exactly the described alternative: a partially-pooled
per-window positivity `p_pos` fit to the test-positivity data, with dark-window
positives = `p_pos · volume`. So option (b) is **already in the codebase as a
switch**.

**The identifiability consequence (per memory + docstrings).** The joint
`C(T) ≈ 4900-5200` exceeds every single-stream fit (≤ 2280) *because* the
composition link pins `λ_bg` low, attributing more of the suspected pool to BVD.
Switching to `:free` positivity:
- **Gains:** a positivity curve that tracks the noisy per-vintage data
  directly; cleaner mixing (the `:free` RE historically converges,
  R-hat ≈ 1.006); removes the structural assumption that the *only* explanation
  for less-than-100% positivity is the suspect-pool composition.
- **Loses:** `λ_bg` identification. With `:free`, positivity absorbs the
  non-BVD share into its own curve, `λ_bg` floats up toward its prior, and
  `C(T)` is then pinned only by deaths + exports (the single-stream regime).
  The framing "positivity × received, false-positives = (1−spec)×received" is
  essentially `:free` with a fixed specificity floor — it does not constrain
  `λ_bg` because `p_t` and `(1−spec)` are not linked to the suspect-pool
  composition.

**Options.**
- **(3a) Keep `:composition` as headline** (with the Q2 lag fix). `λ_bg`/`C(T)`
  identified by the lab data, at the cost of a strong structural assumption.
- **(3b) Switch to `:free` positivity** as headline. Positivity follows the
  data; `λ_bg` reverts to weakly-identified; `C(T)` set by deaths + exports;
  report a smaller, single-stream-consistent size.
- **(3c) Hybrid:** `:free`-style positivity RE for the *fit to test-positivity
  data*, but retain a weak composition-anchored prior on the baseline `q_mu` so
  positivity is data-driven yet still leans on `φ`. Softer `λ_bg`
  identification than (3a), less prior-driven than pure `:free`.

**Recommendation.** This is a modelling-philosophy choice, not a bug. If the
goal is *the lab data should drive positivity and we should not over-claim
`λ_bg` identification*, choose **(3b)** and accept the smaller,
deaths/exports-pinned `C(T)`. If the goal is *use the lab signal to bound the
non-BVD background and the joint size*, keep **(3a)**. Given the long history of
the composition link being the thing that inflates `C(T)` beyond every single
stream, leading with **(3b)** as the headline and keeping `:composition` as the
labelled sensitivity analysis makes the `λ_bg`-pinning assumption explicit
rather than load-bearing. **(3a) and (3b) are mutually exclusive** for the
headline.

---

## Q4. Background deaths and death/case lab alignment

**Now.** Confirmed deaths are grounded on the *case* composition `q_susp` (a
single scalar, observations.jl:969) enriched by `m_death`, scoring
`p_death_conf · deaths_daily` (line 980). There is no death-side
received→tested→positive pipeline, no death-side `s`/`spec`, and the death
pool's *own* composition is deliberately not used because the default
suspected-death background `λ_bg_death = 0` makes it collapse to 1 (docstring
937-955). The background-death term `λ_bg_death` exists (priors.jl:410) but is
**off by default** in the joint (joint.jl:312-325).

**Assessment of alignment.** Two requested things are partly present, partly
missing:
- *Symmetric background deaths:* the machinery exists
  (`death_background_model`, `deaths_model` background branches) but is not
  wired on by default. Turning it on is a one-line composer change, and it is
  the prerequisite for any death-side composition.
- *Aligned lab models:* currently asymmetric. Cases have a windowed `φ` +
  `s·q+(1−spec)(1−q)` + received-volume pipeline; deaths have a scalar `q_susp`
  + odds enrichment, no volume, no assay terms.

**Options.**
- **(4a) Switch on `λ_bg_death`, then give deaths a death-pool composition
  link** `q_death = bvd_deaths/(bvd_deaths + bg_death)`, with the shared
  `s`/`spec`, and (if death-specimen volumes exist) a received→analysed
  denominator like cases. Fully symmetric.
- **(4b) Switch on `λ_bg_death` only**, keep the death-confirmation grounded on
  `q_susp` but enrich by the *true* death composition once the background is
  non-degenerate. Partial alignment.
- **(4c) Status quo** — `q_susp`-grounded thinning, no background deaths.

**Trade-offs.** (4a) is the principled symmetric model but needs death-specimen
analysed/received data to have a denominator; if those data do not exist
(deaths are few, post-mortem swabbing rare), the death lab model is
under-identified and (4a) buys little over (4b). Turning on `λ_bg_death`
(4a/4b) does **not** by itself fix the death level (the death model is
secondary and over-prediction is upstream), but it does make the death stream
consistent with the case stream's treatment of suspected counts as
non-pure-BVD.

**Recommendation.** Switch on the **suspected-death background `λ_bg_death`** so
the deaths stream mirrors the cases stream's "suspected ≠ true BVD" structure
(it is already built, just off). Then move to **(4b)**: keep the
death-confirmation tied to composition but use the now-non-degenerate
death-pool share, sharing the case stream's `s`/`spec`. Only go to full
**(4a)** (received→analysed death pipeline) if death-specimen volume data
exist. Crucially, the death lab link must use the **same positivity philosophy
as the cases** — so it should follow whatever is chosen in Q3, not diverge.

---

## Recommended package (coherent sets)

These changes are not all independent. Two consistent bundles:

**Bundle A — keep the joint-size identification (composition headline):**
- Q2 fix (lag `φ` by the receipt delay) — do this regardless.
- Q1: time-varying `τ_test,t` (1b), not a free `V_t`.
- Q3: keep `:composition` (3a) as headline, `:free` as sensitivity.
- Q4: switch on `λ_bg_death`, align deaths to the *composition* link (4b/4a)
  with shared `s`/`spec`.
- Net: `λ_bg` and `C(T)` stay lab-identified (≈ 5000); symmetric, internally
  consistent; strong structural assumption stated explicitly.

**Bundle B — let the data drive positivity (free headline):**
- Q2 fix still applies.
- Q1: time-varying `τ_test,t` (1b); a free `V_t` (1a) is acceptable here since
  `λ_bg` is no longer lab-pinned anyway.
- Q3: switch to `:free` positivity (3b) as headline; positivity follows the
  test-positivity data; `C(T)` pinned by deaths + exports (smaller,
  single-stream-consistent).
- Q4: switch on `λ_bg_death`, align deaths to the *free* positivity philosophy
  (death positivity as its own RE on death-pool data), not the composition
  link.

**Mutually exclusive:** Q3 (3a) vs (3b) — composition-link `λ_bg`
identification cannot coexist as the headline with a free positivity RE; pick
one. The Q4 death alignment must follow whichever Q3 choice is made. Q1 (1a
free volume) is only coherent inside Bundle B.

**Do-anyway:** the Q2 receipt-delay lag on `φ` and switching on `λ_bg_death`
are improvements under either bundle.
