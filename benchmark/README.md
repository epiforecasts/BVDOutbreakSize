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
It benchmarks both revisions in one process, so there is no second runner to divide by.

## Why a driver rather than `benchpkg`

AirspeedVelocity's own entry points, the `action-v1` action and the `benchpkg` CLI, obtain a revision by handing it to `Pkg.add` as a git tree.
Pkg cannot check out a tree that declares a submodule:

```
GitError(Code:ERROR, Class:Submodule, cannot get submodules without a working tree)
```

This repository declares `external/bdbv-linelist-analysis`, so every route through `Pkg.add` fails, whether the package is named by url or by path.

`rev = "dirty"` is AirspeedVelocity's supported local mode.
It calls `Pkg.develop` on a path, which never clones and so never reaches the submodule.
`ci/run_pair.jl` uses that mode against two git worktrees, which `git worktree add` materialises with the submodule left uninitialised.
Nothing under `src/` reads the submodule, so an uninitialised one costs the benchmarks nothing.

Each worktree runs its own `benchmarks.jl`, its own `test/ad_fixtures.jl` and its own `Project.toml`.
So a pull request's changes to the suite are exercised by that pull request rather than the one after it, both arms resolve the versions their own revision pins, and a component that exists on one side only is reported as added or removed rather than failing the other arm.

## Why the comment is built here rather than by AirspeedVelocity

AirspeedVelocity runs the suite; `ci/comment.jl` reports it.

Two things the comment needs are not in AirspeedVelocity's own table.

It has no neutral band and no way to set one.
A fixed band is what let the old harness call noise a regression, so the band here is measured from the run.
It is the 90th percentile of the per-benchmark sample spread, floored at 2% and capped at 20%, and the comment states the number it measured.

That table also cannot show whether the benchmarks moved together.
Unrelated components share no cause, so one factor applied to all of them is an environment difference rather than the diff, and that is the signature that diagnosed the two-runner bias in the first place.
The comment reports the range of the ratios across benchmarks and warns when they all move as one.

Both are recoverable because AirspeedVelocity writes the raw per-sample times into its results JSON, not only a summary.

The ratio stays `PR / main`, so below 1 means the pull request is faster.

## What it still cannot resolve

The band is a lower bound.
Each revision is run once, one after the other, with no interleaving and no repeated rounds.
So the spread it can measure is dispersion within one revision's own samples, not drift between the two revisions, and those are separated by the twenty minutes it takes to compile the second one's gradients.

The large, systematic, cross-machine error is gone.
A smaller within-machine one is bounded from below rather than measured.
Treat a ratio inside the stated band as unresolved.

## Benchmark parameters

`benchmarks.jl` sets `BenchmarkTools.DEFAULT_PARAMETERS` before it builds the suite, because AirspeedVelocity calls `run(SUITE)` with no arguments and a benchmark's own parameters are the only place left to set a budget.
The budget is one second per benchmark rather than the default five.
`gctrial` and `gcsample` are off: both exist to keep a mean or a median honest, and this suite reports a minimum.
`run.jl` takes the same parameters, so a local run and a CI arm sample the same way.

## Running

```bash
task benchmark                        # results.json in the working directory
task benchmark -- out.json            # somewhere else
BVD_BENCH_JOINT=true task benchmark   # plus the full joint
BVD_BENCH_ENZYME=true task benchmark  # plus the Enzyme backend


# Reproduce the CI comparison: two revisions, one process, plus the comment
task benchmark-pair                   # main vs HEAD
task benchmark-pair -- v2.0.0 HEAD    # any two revisions
```

`task benchmark` times one revision, which is what you want while profiling a change in place.
`task benchmark-pair` is the CI comparison, and costs roughly twice as much.
It checks the two revisions out under `.benchmark-worktrees/` and writes both arms' results and the rendered comment to `benchmark-results/`.
Read it only from a quiet machine.

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

The suite is built once per revision, so the smoke test runs against each revision's own `src/`.
Because a pair is registered only when its smoke test passes, a pair that disappears between two revisions is a component that stopped differentiating.
`ci/comment.jl` calls that out separately from the timing tables.

## Shared fixtures

The component list lives in `test/ad_fixtures.jl`, not here, and this suite includes that file.
`test/test_ad_gradients.jl` asserts every component differentiates under Mooncake from the same list, so the benchmarked surface and the tested surface cannot drift apart.
Add a model there and it is both timed and asserted.

## Relationship to `scripts/bench_*.jl`

`scripts/bench_convolve.jl` and `scripts/bench_discretise.jl` stay where they are.
They are one-off diagnostics that time a superseded implementation against the current one to record why the current one was kept, so their comparison arms are deliberately dead code with a fixed answer.
Folding them in would mean carrying those reimplementations in the benchmark environment to report a ratio that never moves.
They also measure pure helpers below the component level this suite reports.

## Files

| File | What it does |
|---|---|
| `benchmarks.jl` | Builds `SUITE` from the shared fixtures |
| `src/log_density.jl` | One evaluation per component |
| `src/ad_gradients.jl` | One gradient per component per backend |
| `run.jl` | Times one revision once, for local profiling |
| `ci/run_pair.jl` | Times two worktrees in one process under AirspeedVelocity |
| `ci/comment.jl` | Turns the two results files into the PR comment |

`ci/Project.toml` carries the harness only, with no model dependency.
`Project.toml` is the environment the suite itself runs in, and is taken from each arm's own worktree.

## CI

`.github/workflows/benchmark.yaml` runs on pull requests that touch `src/`, `ext/`, `benchmark/` or `test/ad_fixtures.jl`.
It does not run on pushes to `main` and records no history.
There is no companion `benchmark-history.yaml`: the docs workflow already spends hours fitting and the runner queue has no room for a timeline.

The workflow is one job, and `scripts/run_benchmark_pair.sh` runs the same steps locally.
It checks the repository out with full history, materialises both revisions as worktrees, runs `ci/run_pair.jl` over them and posts the comment `ci/comment.jl` renders, which also goes to the job summary.
Fork pull requests are skipped: the workflow uses `pull_request` rather than `pull_request_target`, so a fork never runs with a write token, and it could not post the comment anyway.

A run lands near 70 minutes, most of it compiling gradients twice.
