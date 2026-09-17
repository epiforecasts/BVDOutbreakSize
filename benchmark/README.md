# BVDOutbreakSize benchmarks

Per-component AD gradient benchmarks.
Each benchmark times one unconstrained log-density evaluation, and one gradient of it per AD backend, for a single model component.

## Why per component

The sampler's cost is the gradient of `bvd_joint`, but that is not a useful benchmark.
One gradient is about 14 ms over 76 parameters, behind a cold compile of roughly 18 minutes, most of it type inference rather than Mooncake.
A full-joint benchmark would cost more per CI run than the rest of the test suite and would report a single number that says nothing about where the time went.

The components are the units the joint is built from.
The observation submodels in `src/models/observations.jl` are one unit each, evaluated on a fixed prior draw of the latent trajectory.
The single-stream composers in `src/models/joint.jl` are the same submodels with the shared infection and onset process attached, so a composer minus the `latent` baseline is the marginal cost that stream's likelihood adds.
That difference is what guides optimisation: earlier profiling found `onset_reporting` and `treatment_flow` together were about 44% of the three-patch gradient and the spatial structure about 18%.

The joint is available as a component but off by default.
`BVD_BENCH_JOINT=true` adds it; expect the compile cost above.

## Why both revisions on one machine

A comparison is only as good as the machine underneath it.
Each revision used to be timed in its own CI job, which made every reported ratio a division of one hosted runner's speed by another's.
GitHub's pool is heterogeneous by about a factor of two and the difference landed whole on the ratio.

That is not a theoretical worry.
Across nine comparisons, `province_composition_model` ranged from 971 ns to 1.90 μs on code that barely changed.
Four of those pull requests touched no differentiated source at all: a formatter swap, a scoring change, a dependency pin and a forecast fix.
They reported every one of the sixteen log-density benchmarks moving together, by 1.37×, 0.61×, 0.78× and 0.98× respectively.
The direction was not even consistent, which is why it read as signal rather than as noise.

[AirspeedVelocity](https://github.com/MilesCranmer/AirspeedVelocity.jl) removes that at the root.
Its `benchpkg` obtains both revisions itself and benchmarks them in one job on one machine, so there is no second runner to divide by.

## What the comparison does and does not tell you

AirspeedVelocity reports, per benchmark, each revision's median with an interquartile range, and a ratio with the error propagated from those two ranges.
The ratio is `main / PR`, so **above 1 means the pull request is faster**, which is the opposite of the convention the old comment used.

It has no neutral band and no noise threshold.
It runs each revision once, in its own process, one after the other, with no interleaving and no repeated rounds.
The two arms are therefore separated by the twenty minutes it takes to compile the second one's gradients, and slow drift over that gap lands on the ratio with nothing to flag it.
The `±` on a ratio is the only uncertainty signal, and it is within-trial dispersion rather than run-to-run drift.

So the large, systematic, cross-machine error is gone and a smaller within-machine one is not measured.
Read a ratio whose `±` overlaps 1 as unresolved.

## The suite is frozen at `main`

`benchpkg` resolves `benchmark/benchmarks.jl`, and everything it includes, once from the revision named by `--bench-on`, then runs that single suite definition against each revision's `src/`.
CI leaves `--bench-on` at its default, the repository's default branch.

Two consequences, both deliberate.

A change to this directory or to `test/ad_fixtures.jl` is not exercised by its own pull request.
It takes effect on the next one.

A pull request that adds a model component does not get that component benchmarked, and one that renames a model cannot break the baseline arm by making the fixtures unresolvable there.
Graceful degradation was preferred to measuring a new component in the pull request that introduces it.

## Running

```bash
task benchmark                        # results.json in the working directory
task benchmark -- out.json            # somewhere else
BVD_BENCH_JOINT=true task benchmark   # plus the full joint
BVD_BENCH_ENZYME=true task benchmark  # plus the Enzyme backend
```

`task benchmark` times one revision, which is what you want while profiling a change in place.
There is no task that reproduces the CI comparison, because that comparison is AirspeedVelocity's.
To run it by hand, install `benchpkg` and point it at two revisions:

```bash
julia -e 'import Pkg; Pkg.add("AirspeedVelocity"); Pkg.build("AirspeedVelocity")'
~/.julia/bin/benchpkg BVDOutbreakSize --path=. --bench-on=main \
  --rev=main,HEAD --output-dir=results
~/.julia/bin/benchpkgtable BVDOutbreakSize --rev=main,HEAD \
  --input-dir=results --ratio --mode=time,memory
```

That invocation is not exercised in CI, which drives the action rather than the CLI, so treat it as a starting point.
Do not pass `--bench-on=dirty`: with a local path it also discards `benchmark/Project.toml`, leaving the benchmark environment without its dependencies.

## Benchmark parameters

`benchmarks.jl` sets `BenchmarkTools.DEFAULT_PARAMETERS` before it builds the suite, because AirspeedVelocity calls `run(SUITE)` with no arguments and a benchmark's own parameters are the only place left to set a budget.
The budget is one second per benchmark rather than the default five.
`gctrial` and `gcsample` are off: both exist to keep a mean or a median honest, this suite reports a minimum, and three full collections per benchmark were most of an 11 minute pass.

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

AirspeedVelocity flattens this to one row per leaf, joining the group keys with `/`, so the comment carries 32 rows folded into a `<details>` block per mode.

Backends are Mooncake, the package default, and Enzyme, the opt-in backend from the package's Enzyme extension.
Enzyme is off unless `BVD_BENCH_ENZYME=true` is set.
Both backends over every component do not fit a CI run: that sweep was cancelled at the 90 minute cap, where Mooncake alone finished in 39 minutes on the same cold runner.

## Known-broken pairs

Every (component, backend) pair is smoke-tested before it is registered: the gradient is taken once and checked finite and non-trivial.
A pair that throws or comes back degenerate is skipped with a line on stderr, so the suite runs unattended rather than aborting on a backend that cannot compile a component.

Enzyme cannot differentiate `bvd_joint` at all, for two stacked reasons: boxed `map(do)` closures in the model, and then an upstream `nodecayed_phis!` LLVM bug once those are removed.
That is [issue #445](https://github.com/epiforecasts/BVDOutbreakSize/issues/445), and it is the standing reason for the smoke test.

The suite is built once per revision, so the smoke test runs against each revision's own `src/`.
A pair that stops differentiating therefore still shows up, as a row present for one revision and blank for the other.
AirspeedVelocity leaves the ratio cell empty rather than calling it out, so that case needs reading for rather than jumping out.

## Shared fixtures

The component list lives in `test/ad_fixtures.jl`, not here, and this suite includes that file.
`test/test_ad_gradients.jl` asserts every component differentiates under Mooncake from the same list, so the benchmarked surface and the tested surface cannot drift apart.
Add a model there and it is both timed and asserted.

`benchmarks.jl` reaches it through `@__DIR__`, which under `benchpkg` resolves inside the package checkout of the `--bench-on` revision, so the include finds a real file rather than a stripped temp directory.
The `[sources]` entry in `Project.toml` naming `BVDOutbreakSize` is stripped by `benchpkg`, which supplies the package at an explicit revision instead.
Do not add a `[sources]` entry for anything else: only `benchmark/Project.toml` is copied into the temp environment, so a relative path would be resolved against that directory and fail.

## Relationship to `scripts/bench_*.jl`

`scripts/bench_convolve.jl` and `scripts/bench_discretise.jl` stay where they are.
They are one-off diagnostics that time a superseded implementation against the current one to record why the current one was kept, so their comparison arms are deliberately dead code with a fixed answer.
Folding them in would mean carrying those reimplementations in the benchmark environment to report a ratio that never moves.
They also measure pure helpers below the component level this suite reports.

## Files

| File | What it does |
|---|---|
| `benchmarks.jl` | Builds `SUITE` from the shared fixtures; the entry point `benchpkg` discovers |
| `src/log_density.jl` | One evaluation per component |
| `src/ad_gradients.jl` | One gradient per component per backend |
| `run.jl` | Times one revision once, for local profiling |

## CI

`.github/workflows/benchmark.yaml` runs on pull requests that touch `src/`, `ext/`, `benchmark/` or `test/ad_fixtures.jl`.
It does not run on pushes to `main` and records no history.
There is no companion `benchmark-history.yaml`: the docs workflow already spends hours fitting and the runner queue has no room for a timeline.

The workflow is one job holding one `MilesCranmer/AirspeedVelocity.jl@action-v1` step, which installs Julia, caches the depot, benchmarks both revisions and posts the comment.
Fork pull requests are skipped: the workflow uses `pull_request` rather than `pull_request_target`, so a fork never runs with a write token, and it could not post the comment anyway.

The action has no input for `julia-actions/cache`, so the depot cache uses the action's own key rather than a pinned snapshot.
A run lands near 70 minutes, most of it compiling gradients twice.
