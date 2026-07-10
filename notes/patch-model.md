# Patch (meta-population) model for the 2026 DRC BVD outbreak

## Motivation

The current model treats the DRC outbreak as a single well-mixed population.
In reality, transmission is spatially structured: the outbreak began in Ituri
Province and has since spread to Nord-Kivu and (to a much lesser extent)
Sud-Kivu. The INSP sitreps publish per-province spatial tables (Tableau 1 &
Tableau 6) for confirmed cases, confirmed deaths, isolation occupancy and
laboratory throughput. A patch model lets us:

1. **Fit per-province data streams directly** — confirmed cases by province,
   lab analysed by province, isolation occupancy by province — instead of
   only the national aggregates.
2. **Estimate province-specific reproduction numbers** that share strength
   through a hierarchical national `Rt`.
3. **Model between-province spread** via an importation term, capturing how
   the outbreak moves from the Ituri epicentre into Nord-Kivu and Sud-Kivu.
4. **Keep fitting the national aggregates** by summing the patch trajectories,
   so the model still talks to the headline numbers.

## Patch structure

Three patches matching the INSP spatial-table provinces:

| Patch | Province | Population | Role |
|-------|----------|------------|------|
| 1     | Ituri    | ~4.4M      | Epicentre (first cases, largest burden) |
| 2     | Nord-Kivu | ~6.6M      | Secondary spread |
| 3     | Sud-Kivu  | ~5.8M      | Low-level incursions |

Note: The national model's source population `ITURI_POPULATION` (4.39M) and
travel volume `ITURI_DAILY_TRAVEL` (1871) are Ituri-specific. The patch model
will need population estimates for each province, and a travel/mobility matrix
for the importation kernel.

## Latent process

### Per-patch renewal with importation

For each patch `p` on day `t`:

```
I_{p,t} = R_{p,t} · Σ_s I_{p,t−s} · g_s  +  import_{p,t}
```

where:

- `g_s` is the **shared** generation-interval PMF (same gi for all patches,
  sampled once from the existing [`generation_interval_model`](@ref)).
- `R_{p,t}` is the **patch-specific** reproduction number.
- `import_{p,t}` is imported infections arriving in patch `p` from other
  patches on day `t`.

### Seeding

Each patch has its own cryptically-growing seed:

- Ituri gets the existing `seed_model` (I0 ~ N⁺(0.1, 0.1), exponential growth
  at the sampled `r`). The doubling-count prior `m ~ N⁺(3, 3)` applies to the
  Ituri patch only — Ituri is the origin patch.
- Nord-Kivu and Sud-Kivu are seeded later, at or after the first confirmed case
  appears in that province. The seeding time for each secondary patch is a
  sampled offset from the Ituri renewal start, with a prior that reflects the
  date of the first confirmed case in that province.

A simpler alternative for the first pass: all three patches share the same
renewal-start day but the secondary patches initialise from a much smaller
seed (e.g. `I0 ~ N⁺(0.01, 0.01)`) so their early trajectory is near-zero and
grows only through importation from Ituri.

### Reproduction numbers

```
log R_{p,t} = log R_national(t) + δ_p(t)
```

where:

- `log R_national(t)` is the **national** log-Rt trajectory from the existing
  `rt_walk_model` — a weekly-knot random walk with a sigmoid intervention
  ramp at the first sitrep date.
- `δ_p(t)` is a **patch modifier** on log-Rt, allowing province-specific
  deviations from the national trend.

**Option A — constant modifier (simpler, better-identified):**

```
δ_p ~ Normal(0, σ_region)
```

The modifier is constant over time, so the patch Rt inherits the national
shape (same growth/decay pattern) but scaled up or down. This captures
persistent differences in transmission (e.g. higher density in Bunia,
different response effectiveness).

**Option B — time-varying modifier (richer, harder to identify):**

```
δ_p(t) interpolated from weekly knots with a random walk
δ_p ~ Normal(δ_p(t−1), σ_patch_rw)
```

This lets the province-specific trajectory diverge from the national shape,
at the cost of P×K additional knot parameters. Likely over-parameterised for
the data available.

**Recommendation**: Start with Option A, the constant modifier. If the
per-province data are informative enough, extend to Option B.

**Option C — independent multivariate-normal random walk (implemented):**

The current implementation uses `patch_rt_mvwalk_model`, which replaces
the national Rt walk + patch modifiers with a single multivariate-normal
random walk on joint log-Rt for all patches at weekly knots:

```
log R_1(t_k)       \nlog R_2(t_k)        ∼ MVN( log R_1(t_{k-1}), Σ )
...                 \nlog R_p(t_k)       /
```

The covariance matrix Σ is decomposed as `diag(σ_patch) · Ω · diag(σ_patch)`
where `σ_patch[p]` is the patch-specific step SD (sampled) and Ω is a
correlation matrix with an LKJ(2.0) prior. This lets the walk learn
cross-patch correlation from the data — patches whose log-Rt co-move
(e.g. both reacting to the same national intervention) share strength
through Ω. An `intervention_effect` with a half-normal prior
`N⁺(0, 0.3)` is applied via a sigmoid ramp at the breakpoint, allowing
a post-intervention decline in Rt across all patches.

Per-patch initial log-Rt `log_R0_patch[p]` is centred on `log_R0_base`
with patch-specific variation `σ_R0 · z0[p]`.

### Importation kernel

```
import_{p,t} = ε · Σ_{q ≠ p} ( I_{q,t−1} · m_{q→p} )
```

where:

- `m_{q→p}` is the daily **per-capita travel rate** from province `q` to
  province `p`, derived from known mobility data (road networks, PoE flows).
- `ε` is the **importation intensity** — the fraction of travellers who are
  infectious and successfully establish a secondary infection. Sampled with a
  strong prior (e.g. `ε ~ Beta(1, 100)` — rare events).

For the first pass, a simplified isotropic kernel:

```
import_{p,t} = ε · Σ_{q ≠ p} I_{q,t−1} · (w_q / Σ w_q)
```

where `w_q` weights each source patch by its population (more people = more
travellers). This avoids needing a full O-D matrix.

Even simpler: one-directional importation from Ituri only, since Ituri is the
epicentre and the Nord-Kivu / Sud-Kivu cases are overwhelmingly Ituri-linked:

```
import_{NK,t} = ε · I_IT,t−1 · (travel_volume / pop_IT)
import_{SK,t} = ε · I_IT,t−1 · (travel_volume / pop_IT) · r_SK
```

where `r_SK < 1` reflects much lower travel volume to South Kivu.

## Data streams

### National-level (sum of patches)

These are the existing national data streams, fitted to the sum of patch
trajectories. No change to the observation model structure.

| Stream | Patch sum |
|--------|-----------|
| Reported (suspected) cases | Σ_p onsets_p (through ascertainment) |
| Suspected deaths | Σ_p deaths_p (through CFR) |
| Confirmed cases (national) | Σ_p confirmed_p (through lab) |
| Confirmed deaths (national) | Σ_p confirmed_deaths_p |
| Exports to Uganda | Ituri-only (exports are Ituri → Uganda) |

### Province-level (individual patches)

New observation models that fit per-province data from the spatial tables.

| Stream | Patch | Observation model |
|--------|-------|-------------------|
| Confirmed cases by province | Each p | Per-province cumulative confirmed (spatial table) fitted through the same lab-positivity model but with province-specific positivity |
| Lab analysed by province | Each p | Per-province 24h analysed counts (section 4.3) fitted as Binomial denominator against the province's new confirmed |
| Isolation occupancy by province | Each p | Per-province isolation (Tableau 6 Fin J) fitted through the length-of-stay model |
| Treatment flows by province | Each p | Per-province admissions / deaths / ruleouts (Tableau 6) |

## Implementation plan

### 1. Add `patch_infections()` to `renewal.jl`

A new multi-patch renewal function that takes a matrix of patch-Rts
(`n_patches × n_days`) and returns a matrix of patch infections
(`n_patches × n_days`), with an importation kernel:

```julia
function patch_infections(Rt_matrix, g, seeds_matrix, importation_kernel)
```

- Each row is one patch's daily infection trajectory.
- Between-patch importation: `import_{p,t} = ε · Σ_q K_{p,q} · I_{q,t−1}`
  where `K` is the importation kernel.
- The same generation interval `g` is shared.

AD-transparent under Mooncake.

### 2. Add `patch_rt_model()` to `models/priors.jl`

A hierarchical Rt model that samples the national Rt walk and per-patch
modifiers, combining them into a `n_patches × n_days` Rt matrix:

```julia
@model function patch_rt_model(n, n_patches; ...)
    # National Rt walk (existing)
    # Per-patch constant modifiers δ_p ~ Normal(0, σ_region)
    # Rt_matrix[p, t] = Rt_national[t] * exp(δ_p)
end
```

### 3. Add `patch_infection_model()` to `models/priors.jl`

The latent infection process for the patch model. Builds per-patch
seeds, runs the multi-patch renewal, and computes per-patch onsets:

```julia
@model function patch_infection_model(n, n_patches, breakpoint; ...)
    # Generation interval (shared, sampled once)
    # Patch Rt matrix from patch_rt_model
    # Per-patch seeds
    # Multi-patch renewal via patch_infections()
    # Per-patch onsets via onset_incidence_model per patch
end
```

### 4. Add `_patch_latent()` and `bvd_patch_joint()` to `models/joint.jl`

The joint composer that:
1. Runs `patch_infection_model` for all patches
2. Sums patch onsets → national-level streams (existing observation models)
3. Routes individual patch onsets → per-province observation models

### 5. Add per-province observation models to `models/observations.jl`

New observation submodels that consume a per-patch trajectory:

- `confirmed_cases_patch_model` — per-province confirmed cases from spatial
  tables (Tableau 1 in sitreps)
- `isolation_patch_model` — per-province isolation occupancy (Tableau 6
  "Patients en isolement Fin J")
- `lab_analysed_patch_model` — per-province daily lab analysed counts
  (section 4.3)

These mirror the national observation models but consume a single patch's
onset/confirmed trajectory.

### 6. Data loading

The per-province data lives in `data/insp_sitrep_scanned.csv` (the LLM-scanned
sitrep values) and in the INRB-UMIE per-province daily CSVs. A new loading
function `load_patch_observations()` in `data.jl` would parse these alongside
the national data.

For the first exploration pass, the per-province spatial table figures can
be hard-coded as constants or loaded directly from the scanned CSV, matching
the dates of the national sitreps.

## Identifiability considerations

1. **The importation rate `ε` and the seed sizes for secondary patches**
   are both informed by the lag between Ituri's growth and Nord-Kivu's.
   If Nord-Kivu cases only appear weeks after Ituri cases, `ε` is small
   and/or the NK seed is delayed.

2. **Constant patch modifiers `δ_p` vs. time-varying national `Rt`**:
   The national `Rt` already captures the overall growth/decline shape.
   The constant modifier captures a persistent offset (e.g. NK's Rt
   is consistently 10% lower than the national average). These are
   identified by the relative case counts — higher case counts in Ituri
   than NK at the same time imply Rt_IT > Rt_NK.

3. **Importation vs. local seeding in secondary patches**: If Nord-Kivu
   grows faster than importation from Ituri can explain, the model will
   infer a larger NK seed or a higher NK patch modifier. The prior on
   `ε` (rare events) and the strong link to travel volume help separate
   these.

## Next steps

1. [x] Write this design document
2. [x] Implement `patch_infections()` in `renewal.jl`
3. [x] Implement `patch_rt_model()` in `models/priors.jl`
4. [x] Implement `patch_infection_model()` in `models/priors.jl`
5. [x] Add `_patch_latent()` and `bvd_patch_joint()` to `models/joint.jl`
6. [x] Fix wiring: treatment model call, duplicate `:=` deterministics
7. [x] Wire up exports in `BVDOutbreakSize.jl`
8. [x] Verify prior predictive (patch_infection_model + bvd_patch_joint)
9. [x] Add per-patch outbreak summaries (`patch_summary_table`)

## Per-location forecasts (next priority)

Once the patch model runs, add per-patch forecast capability:

### 9. Add `forecast_patch()` to `forecast.jl`

A per-patch analogue of `forecast_reported` that projects each patch's
latent infections, onsets, and deaths over the forecast horizon. The
per-patch reproduction number continues evolving (same terminal drift
as the national walk, plus the per-patch modifier). Key differences:
- Reads `patch_state.infections_matrix[p, :]`, `patch_state.onsets_matrix[p, :]`
  from the chain for per-patch daily incidence
- Uses the evolving national Rt (same `_evolving_rates`) plus the constant
  `δ_patch[p]` modifier to get per-patch future Rt
- Returns a `DataFrame` with per-patch columns

### 10. Add `forecast_vs_truth_patch()` to `forecast.jl`

Validates per-patch forecasts against observed per-province data. Uses the
same pattern as `forecast_vs_truth` but checks against spatial-table data
(confirmed cases by province, isolation by province).

### 11. Data: per-province observations for forecast validation

The per-province spatial-table data (confirmed cases by province from
Tableau 1, isolation by province from Tableau 6 Fin J) needs to be loaded
for the forecast-vs-truth comparison. This can come from the existing
INRB-UMIE per-province CSVs or the scanned data.
