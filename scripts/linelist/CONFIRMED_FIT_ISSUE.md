# `confirmed` does not sample on line-list data

Written 2026-08-20, against the data cut-off `LINELIST_AS_OF=2026-08-10`.

This is a write-up for someone who did not run the fits and needs enough to judge the diagnosis and propose an alternative. It records what was measured, what was ruled out and how, and what is still open.

`confirmed` has been removed from the comparison grid as a result. `onsets` carries the comparison and samples acceptably on both sources. `fit_single.jl` will still run `confirmed` on request.

## Symptom

`confirmed_only_model` fitted to the line-list confirmed-case history produces chains that do not move. Released sampler settings throughout: 500 draws, 2 chains, from `docs/fits/registry.jl`.

| | max R-hat | min bulk ESS | divergences |
|---|---|---|---|
| `confirmed`, sitrep | 1.015 | 164.4 | 35 |
| `confirmed`, linelist_known | 2.630 | 2.4 | 0 |

A bulk ESS of 2.4 out of 1000 draws with zero divergent transitions is the important pairing. Divergences indicate a sampler fighting difficult curvature; zero divergences alongside an ESS of 2.4 indicates a sampler that is not exploring at all.

For scale, this repository records its own healthy figures for `confirmed_only_model` at the same 500 x 2 settings in `src/data.jl:207`: roughly 477 to 522 bulk ESS with 20 to 22 divergences, and a documented pathological case at 94 divergences and a bulk ESS of 15. The line-list fit is far below the pathological case on ESS.

## What the chains are doing

Splitting the cached chain by chain, on `dispersion_state.inv_sqrt_k`:

| Run | chain | median | min | max | draws < 1e-8 |
|---|---|---|---|---|---|
| linelist_known | 1 | 5.319e-14 | 5.304e-14 | 5.334e-14 | 100% |
| linelist_known | 2 | 1.327 | 1.323 | 1.346 | 0% |
| sitrep | 1 | 0.4423 | 0.3484 | 0.5567 | 0% |
| sitrep | 2 | 0.4382 | 0.3364 | 0.6469 | 0% |

Both line-list chains are frozen, at different points. Chain 1 varies by 3e-17 across 500 draws. Chain 2 varies by 0.02. The situation-report chains range over roughly 0.34 to 0.65 and mix with each other.

Median log-posterior and a few scalars per chain:

| Run | chain | lp | `growth_state.C_T` | R0 | `k` |
|---|---|---|---|---|---|
| linelist_known | 1 | -2191.4 | 44 | 3.228 | 4.503599627313107e15 |
| linelist_known | 2 | -11442.1 | 2095 | 2.149 | 0.568 |
| sitrep | 1 | -860.5 | 17 | 2.031 | 5.112 |
| sitrep | 2 | -860.2 | 18 | 2.072 | 5.207 |

Worst R-hat by parameter on the line-list fit, top of the list:

```
Extra(:step_size)                      11.489
Extra(:nom_step_size)                  11.489
Parameter(dispersion_state.k)           2.630
Parameter(dispersion_state.inv_sqrt_k)  2.218
Parameter(gi_state.θ)                   2.125
Parameter(growth_state.m)               2.125
Parameter(growth_state.C_T)             2.125
Parameter(rt_state.log_R0)              2.125
... every remaining latent parameter    2.125
```

The identical 2.125 across every latent parameter is consistent with one chain being near-constant rather than with genuine disagreement about a particular quantity. The two chains adapted to different step sizes, reported in the run log as initial epsilon 0.00625 and 0.025.

## Proposed mechanism

`surveillance_dispersion_model` (`src/models/joint.jl`) is:

```julia
@model function surveillance_dispersion_model(;
        inv_sqrt_k_prior = truncated(Normal(0.6, 0.2); lower = 0))
    inv_sqrt_k ~ inv_sqrt_k_prior
    k := 1.0 / (inv_sqrt_k^2 + eps(typeof(inv_sqrt_k)))
    return (; k, inv_sqrt_k)
end
```

`k` is the negative-binomial dispersion; `k` to infinity is the Poisson limit, reached as `inv_sqrt_k` goes to zero, which is the truncation boundary of the prior. The `+ eps` guard caps `k` at `1/eps(Float64)` = 4.503599627370496e15. The observed median `k` in line-list chain 1 is 4.503599627313107e15, i.e. the chain is sitting on that boundary.

The reading is that the line-list confirmed series carries less dispersion than the negative-binomial observation model can express, so the likelihood pushes `inv_sqrt_k` to its boundary and it pins there.

The candidate reason is the construction of the series. `confirmed_case_history` in the line-list manifest has 35 vintages against the situation reports' 82, and the known-by indexing counts each case at the snapshot that first held it, so increments are multi-day sums. Aggregating records into snapshots averages out the day-to-day reporting noise the observation model expects to see.

Increments in the line-list series run 14 to 210 with no zeros and no non-monotonic steps.

## Ruled out

**Harmonisation break days.** `confirmed_break_dates` is deleted from the line-list manifest by `manifest.jl:236`, with the reasoning recorded there — the line list carries its own revisions inside the case records. Verified absent from the substituted manifest and present in the baseline. This mattered because `src/data.jl:207` documents a break day declared without a backlog to absorb as the cause of a 94-divergence, ESS-15 fit, which is the nearest known pathology in this codebase.

**Chain count and initialisation.** A four-chain run was started at the same settings. At the point it was stopped, two chains sat at log-posterior around -2190 and two around -11000, the same split as the two-chain run. More chains would supply more stuck chains, not a fix. Initialisation determines where a chain freezes rather than whether it explores, so pathfinder-style initialisation is not expected to help on its own — though it has not been tried, and the `pathfinder-init` branch exists.

**Two genuine posterior modes.** Rejected because neither chain is exploring. A multimodality story requires at least one chain to be sampling within a mode; here the within-chain ranges are 3e-17 and 0.02.

## Open, and worth a second opinion

- Whether the underdispersion reading is right, or whether the boundary collapse has another cause. Nobody has yet compared the empirical dispersion of the two confirmed series directly against what the model expects, which would settle it.
- Whether the fix belongs in the data construction (a denser known-by series, daily rather than 35 snapshots, if the line list supports it) or in the model (a prior on `inv_sqrt_k` bounded away from zero, or an observation model that admits underdispersion). The first keeps the released model intact and is preferred if the line list can supply it.
- Why the `sitrep` arm also sits well below this repository's own benchmark for the same fit — 164 bulk ESS against a documented 477 to 522. That arm converges, but not comfortably, and it is not obvious why. It may be unrelated; it may be the same effect in milder form.

## A second, separate problem with `confirmed`

Independent of the sampling failure, `confirmed_only_model` conditions on the laboratory volumes as well as the confirmed-case history. `tests_analysed_history` and `tests_analysed_daily_history` are outside `LINELIST_BLOCKS`, so they are byte-identical in both manifests: the fit pairs line-list confirmed counts against situation-report test denominators in the positivity likelihood.

At the cut-off the numerator differs by about 13%, 3,862 against 4,449, with the denominator held fixed. All 63 days with a daily test count were checked and none implies a positivity above 100% under either source, so this is not on its own a broken likelihood. But a difference from `confirmed` would be a difference between confirmed-case series measured against one shared laboratory denominator, not between two independent constructions of the outbreak.

`onsets` has no equivalent problem: it conditions on the reporting triangle alone, which the line list supplies in full.

## Reproducing

```bash
cd ~/Documents/Github/BVDOutbreakSize
export LINELIST_INPUT_DIR=$PWD/ignore/linelist/inputs
export LINELIST_RELEASED_MANIFEST=$PWD/data/observations.toml
export LINELIST_AS_OF=2026-08-10

julia -t 2 --project=docs scripts/linelist/fit_single.jl confirmed \
  --data=linelist_known --delays=repo
```

About 35 minutes. Writes `linelist_confirmed_linelist_known_repo_diagnostics.csv` alongside the trajectory, carrying the R-hat, bulk ESS and divergence count.

A cached chain is reused, so re-running the command above after a fit only re-exports. `--refit` forces a new fit. `--chains=N` changes the chain count and is recorded in the output tag so a diagnostic run cannot be mistaken for a comparison arm.

Per-chain inspection of a cached chain, which is how the tables above were produced, needs `repair_chain_keys` from `docs/fits/cache.jl` and `BVDOutbreakSize.FlexiChains`; the chain is VarName-keyed, so index it with `Symbol("dispersion_state.inv_sqrt_k")` rather than a `FlexiChains.Parameter`.
