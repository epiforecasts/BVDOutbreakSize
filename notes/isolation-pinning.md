# Why the isolation stream pins the joint outbreak size low

## Observation

In recent dev fits the joint posterior outbreak size `C_T` (cumulative
infections at the cut-off) sits well below what most single-stream fits imply
on their own. From the released analysis (per-stream `C_T`, 90% CrI):

| stream          | 90% CrI of C_T |
|-----------------|----------------|
| exports         | 915 – 25 733   |
| deaths (DRC)    | 4 866 – 19 289 |
| cases (DRC)     | 5 867 – 17 429 |
| confirmed (DRC) | 7 412 – 28 773 |
| **isolation**   | **1 548 – 12 869** |
| **joint**       | **1 723 – 4 112** |

The deaths / cases / confirmed streams individually want a *large* outbreak
(several thousand to tens of thousands). The **isolation** stream is the one
single-stream fit whose lower bound anchors near the joint, and the joint
ends up close to the bottom of that range. So when isolation is added to the
mix it pulls `C_T` down hard — the user's "pinned low" observation.

## Mechanism

The isolation/treatment-bed stream (`treatment_admission_model` in
`src/models/observations.jl`) fits the daily occupied-bed count
("Patients en isolement") as a length-of-stay survival of the *admitted*
suspect inflow:

- BVD demand  = `convolve_survival(p_iso_bvd · p_drc · bvd_reports_daily, bvd_los.pmf)`
- background  = `convolve_survival(p_iso · bg_daily, ruleout_los.pmf)`
- occupancy   = `min(demand, usable_frac · C)` (censored saturation), observed
  beds ~ NegBinomial(occupancy).

By Little's law (stated in the `isolation_admission_model` docstring):

```
mean occupancy  ≈  p_iso · admissions · (E[LOS] + 1)
```

The occupied-bed count is *observed* (~376 beds at the cut-off). For a fixed
observed occupancy, a **high admission fraction × long stay** means each
occupied bed corresponds to *few* suspect admissions, hence a *small* implied
suspect inflow and a *small* outbreak. So the isolation stream maps the bed
count into `C_T` through the product `p_iso · E[LOS]`, and that product is the
lever.

## Hypothesis (what to test)

The low pin is driven by the **admission fraction `p_iso` (`Beta(2,2)`, mean
0.5) × BVD length-of-stay (~12 d) product** being effectively *high* and only
*weakly identified*. `p_iso` and `E[LOS]` are explicitly confounded for the
occupancy *level* (only their product is pinned by the level; the shape/lag of
the daily curve adds weak extra information). So the implied outbreak size
rides on the prior centres of these two quantities, not on the data, and those
centres happen to sit high — squeezing `C_T` down.

The released isolation posteriors are consistent with this:

- `isolation_admission` (`p_iso`): 90% CrI **0.19 – 0.54**, median ~0.29 — the
  data pulls it *down* from the prior mean of 0.5 but it stays substantial and
  the interval is wide (weakly identified, as predicted).
- `isolation BVD length-of-stay mean`: 90% CrI **5.4 – 19.3 d**, median ~12 d —
  essentially the prior; the occupancy does not sharpen it much.
- `isolation non-BVD rule-out stay mean`: 90% CrI **1.4 – 6.7 d**.
- bed capacity ~ 410 – 462, expected occupancy ~ 371 – 434, bed shortfall
  ~ 0 (censoring inactive at current utilisation).

A high admission-rate × ~12 d stay is exactly the configuration that turns
~376 beds into a small suspect count.

## The experiment (`scripts/experiment_isolation_sensitivity.jl`)

Two complementary tests, both reading posterior `C_T`:

**(A) Joint leave-one-out.** Fit the full joint with isolation *included* vs
*excluded* (excluded by passing an empty `isolation_history`, so the treatment
submodel samples from its prior and contributes no likelihood). If isolation
pins `C_T` low, removing it should **raise** the joint `C_T` toward what the
other streams want. This is the direct demonstration of the pull.

**(B) Isolation single-stream assumption variants.** Refit the isolation
single-stream (`treatment_only_model`) under:

- *baseline* — shipped priors (`p_iso ~ Beta(2,2)`, BVD LOS ~12 d);
- *low admission* — inject `admission = isolation_admission_model(p_prior = Beta(2,6))`
  (mean 0.25 instead of 0.5);
- *short BVD LOS* — re-centre the BVD stay on ~6 d instead of ~12 d.

Both alternatives *reduce* the `p_iso · E[LOS]` product, so each should
**raise** the single-stream `C_T`. Whichever variant moves `C_T` most is the
assumption that most pins it low. The hypothesis predicts both move it up, and
quantifies their relative weight (admission fraction vs length-of-stay).

The script prints a comparison table of `C_T` medians / 90% / 60% intervals
across all five settings and saves a density overlay to
`logs/isolation_sensitivity_C_T.png`.

### How to read the output

- **(A):** "joint, isolation EXCLUDED" median **higher** than "...INCLUDED"
  confirms isolation pins the joint low; the printed delta quantifies the pull
  in infections.
- **(B):** if *low admission* and/or *short BVD LOS* raise the single-stream
  `C_T` substantially, that assumption is a driver. Compare the two deltas to
  rank admission fraction vs length-of-stay. If the single-stream baseline
  `C_T` is already low and the variants lift it toward the deaths/cases level,
  the prior centres — not the bed data — were holding it down.

## Judgement on reasonableness

The mechanism is structurally sound: occupancy genuinely *is* admissions ×
stay, and that is the correct way to read a prevalence (stock) stream into an
incidence (flow) quantity. The concern is **identifiability, not the model
form**. Because only the *product* `p_iso · E[LOS]` is pinned by the occupancy
level, the implied outbreak size is effectively *set by the priors* on those
two quantities — and `Beta(2,2)` (mean 0.5 admission) combined with a ~12 d
BVD stay is a fairly aggressive, high-retention assumption. If the true
admission fraction is lower (many suspects are never bedded) or the effective
BVD stay shorter (faster turnover/transfer/death), the bed count implies a
*larger* outbreak, and the joint should not be dragged so far below the
deaths/cases/confirmed consensus.

So the pin is **plausible but fragile**: it is reasonable *if* you believe
roughly half of ascertained suspects occupy a bed for ~12 days, and
*unreasonable to lean on* otherwise, because the data barely constrains that
product. The recommendation is to (i) confirm the size of the pull via test
(A), (ii) check whether a more defensible, lower admission prior (e.g.
`Beta(2,6)`) materially relaxes it via test (B), and (iii) if so, consider a
tighter, better-justified prior on `p_iso` or `E[LOS]` (or down-weighting the
isolation stream's influence on `C_T`) rather than letting a weakly-identified
product anchor the headline outbreak size.

This note and the script are a *diagnostic*; they do not modify the model.
