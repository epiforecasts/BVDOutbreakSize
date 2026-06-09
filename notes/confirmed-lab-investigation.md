# Confirmed-cases lab model: why the posterior is high + prior audit

Findings from the read-only investigation. File:line refs in the renewal tree.

## Why the confirmed-ONLY C_T is high and very wide
- `confirmed_only_model` (joint.jl:107-136) runs the suspected stream as a
  pure predictive generator (empty history + missing cut-off, line 126-128),
  so λ_bg, τ_test and p_drc are prior draws, constrained only through the
  confirmed/received likelihood.
- It passes no `tmrca_days`, so there is no genetic seeding bound and the Rt
  walk runs free from day 1; C_T = 2^m grown forward is nearly prior-driven
  (m ~ Normal⁺(3,3), very wide).
- The only confirmed term that touches the infection level is the early/late
  window `NegBinomial(p_pos · modelled_volume)` (observations.jl:686-721); the
  observed-window Binomial conditions on the OBSERVED analysed denominator and
  is decoupled from the size by construction.
- National analysed denominators stop on 28 May, so the 305 confirmed cases
  over 29 May-6 June (210->515) are scored against the MODELLED received
  volume `τ_test·(p_drc·bvd + λ_bg)`. The cheapest way to reproduce them is to
  raise the BVD infection level and the terminal Rt, so C_T inflates and runs
  wide. This is structural, not a coding bug.

## Confirmed > suspected is a date/freeze artefact (not a pathology)
- Data load unfrozen to 6 June; suspected streams are frozen at 26 May and
  simply absent after (observations.toml:36-39); confirmed runs to 6 June.
- The vintage PPC plots confirmed-to-6-June against suspected-frozen-at-26-May
  on different date grids (analysis.jl ~2173-2226). Within one fit confirmed is
  a thinning of the same pipeline (observations.jl:579) and stays <= suspected
  on identical dates. Fix: plot both on the same grid (or freeze confirmed to
  26 May) to confirm.

## Lab prior audit (recommendations)
- **decay_scale ~ Normal⁺(0,200)** (priors.jl:637): FAR too diffuse. A large
  draw keeps the severity enrichment δ0 from decaying, so the tested BVD share
  q stays near 1, positivity stops tracking the composition φ, and the lab data
  stop identifying λ_bg. TIGHTEN substantially — highest-value change.
- **τ_test ~ Beta(5,2)** (priors.jl:452): centred too high (~0.71) and directly
  scales the late-window modelled volume. RECENTRE LOWER / tighten (e.g.
  Beta(2,3), mean 0.4).
- **δ0 ~ Normal⁺(1.5,0.75)** (priors.jl:636): reinforces the above; tighten
  alongside decay_scale.
- s ~ Beta(10,1.76), spec ~ Beta(60,2), m_death ~ LogNormal(0,1), λ_bg ~
  Normal⁺(0,1): reasonable, leave. λ_bg is simply unidentified in the
  confirmed-only fit by construction.

## Wild lab pair plot (analysis.jl ~2093)
Weak identification, not a broken sampler: marginals (λ_bg, τ_test,
death_confirmation) sit near the overlaid prior, and there is a
τ_test–λ_bg–p_drc ridge because the received volume constrains only their
product. Tightening decay_scale (restoring the composition link) and τ_test is
the highest-value fix.
