# BVDOutbreakSize benchmarks

Per-component AD gradient benchmarks.
Each benchmark times one unconstrained log-density evaluation, and one gradient of it per AD backend, for a single model component.

## Why per component

The sampler's cost is the gradient of `bvd_joint`, but on its own that is not a useful benchmark.
One gradient is about 14 ms over 76 parameters, behind a cold compile of roughly 18 minutes, most of it type inference rather than Mooncake.
A full-joint number alone says nothing about where the time went.

The components are the units the joint is built from.
The observation submodels in `src/models/observations.jl` are one unit each, evaluated on a fixed prior draw of the latent trajectory.
The single-stream composers in `src/models/joint.jl` are the same submodels with the shared infection and onset process attached, so a composer minus the `latent` baseline is the marginal cost that stream's likelihood adds.
That difference is what guides optimisation: earlier profiling found `onset_reporting` and `treatment_flow` together were about 44% of the three-patch gradient and the spatial structure about 18%.

The joint is a component too, off by default so a local run of the components stays quick.
`BVD_BENCH_JOINT=true` adds it; expect the compile cost above.
The benchmark workflow sets it, because this comparison is the only place the joint's gradient is exercised: the test suite asserts the components differentiate and leaves the joint to this and to the fits the docs build runs.

## Running

```bash
task benchmark                        # results.json in the working directory
task benchmark -- out.json            # somewhere else
BVD_BENCH_JOINT=true task benchmark   # plus the full joint
BVD_BENCH_ENZYME=true task benchmark  # plus the Enzyme backend

# Compare two saved runs the way CI does
task benchmark-compare -- pr.json main.json comment.md
```

## Structure

```
Log density/
  Latent/     latent, patch_infection_model (uncoupled, coupled)
  Submodel/   reported_cases, confirmed_cases, deaths, exports,
              treatment_flow, onset_reporting, province_composition
  Composer/   exports_only, deaths_only, cases_only, confirmed_only,
              treatment_only, onsets_only
AD gradients/
  <same groups>/<component>/<backend>
```

Backends are Mooncake, the package default, and Enzyme, the opt-in backend from the package's Enzyme extension.
Enzyme is off unless `BVD_BENCH_ENZYME=true` is set.
Both backends over every component do not fit a CI run: that sweep was cancelled at the 90 minute cap, where Mooncake alone finished in 39 minutes on the same cold runner.

## Known-broken pairs

Every (component, backend) pair is smoke-tested before it is registered: the gradient is taken once and checked finite and non-trivial.
A pair that throws or comes back degenerate is skipped with a line on stderr, so the suite runs unattended rather than aborting on a backend that cannot compile a component.

Enzyme cannot differentiate `bvd_joint` at all, for two stacked reasons: boxed `map(do)` closures in the model, and then an upstream `nodecayed_phis!` LLVM bug once those are removed.
That is [issue #445](https://github.com/epiforecasts/BVDOutbreakSize/issues/445), and it is the standing reason for the smoke test.

Because a pair is registered only when its smoke test passes, a pair that disappears between two revisions is a component that stopped differentiating.
`compare.jl` calls that out separately from the timing tables.

## Shared fixtures

The component list lives in `test/ad_fixtures.jl`, not here, and this suite includes that file.
`test/test_ad_gradients.jl` asserts every component differentiates under Mooncake from the same list, so the benchmarked surface and the tested surface cannot drift apart.
Add a model there and it is both timed and asserted.

## Relationship to `scripts/bench_*.jl`

`scripts/bench_convolve.jl` and `scripts/bench_discretise.jl` stay where they are.
They are one-off diagnostics that time a superseded implementation against the current one to record why the current one was kept, so their comparison arms are deliberately dead code with a fixed answer.
Folding them in would mean carrying those reimplementations in the benchmark environment to report a ratio that never moves.
They also measure pure helpers below the component level this suite reports.

## CI

`.github/workflows/benchmark.yaml` runs the suite on pull requests that touch `src/`, `ext/`, `benchmark/` or `test/ad_fixtures.jl`, once per revision in its own job, and posts a single comparison comment.
It does not run on pushes to `main` and records no history.
