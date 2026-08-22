# Refitting on line-list case counts

This model reads national counts from the INSP situation reports.
These scripts refit it on case counts derived from a DHIS2 case line list instead, holding the model fixed and changing only the data.
A difference against a fit of the situation reports at the same cut-off is then a difference between data sources rather than between models.

The comparison is carried by the `onsets` fit.
`confirmed` was tried and does not sample on line-list data; the measurements are under Scope.

## What is here

| File | |
|---|---|
| `manifest.jl` | Resolves the inputs and builds the observation manifest each fit reads: the released one with the line-list case streams substituted in, or the released one with only its cut-off moved. Shared by every fit so they cannot drift. |
| `delays.jl` | Turns the delay estimates fitted in `bvd-internal-cmmid` into priors. Reads them at run time; holds no values of its own. |
| `fit_single.jl` | One stream at a time. `onsets` (half an hour to an hour and a quarter) carries the comparison; `confirmed` (about half an hour) and `cases` still run but are out of it, for the reasons under Scope. Writes the daily reproduction-number trajectory. |
| `run_grid.sh` | Runs a whole comparison grid, detached and serially. Use it rather than calling Julia in a loop, since a closed terminal or a sleeping laptop is the likeliest way a multi-hour grid ends. |
| `plot_rt.jl` | Draws the trajectories the grid produced against each other. |
| `fit_joint.jl` | The headline joint fit on line-list cases. Five to seven hours. Not part of the R_t comparison; see Scope. |
| `run_joint_fit.sh` | Detached wrapper for `fit_joint.jl`. |
| `CONFIRMED_FIT_ISSUE.md` | Why `confirmed` is out of the comparison: the measurements, what was ruled out, and what is still open. |
| `REFIT_PLAN.md` | The refit the pinned breakpoint requires. Delete once done. |

## Inputs

Two files, produced by whoever holds the line list, in a directory named by `LINELIST_INPUT_DIR`.
There is no default: these are counts derived from individual patient records, so the operator says where they are rather than a default finding a directory that happens to exist.

`linelist_streams_known.csv` carries the case streams that replace the situation-report ones.

| Column | |
|---|---|
| `stream` | One of `confirmed_case_history`, `reported_case_history`, `suspected_daily_history`. Any other value is ignored. |
| `date` | Vintage date, ISO format. One row per stream per date, ordered arbitrarily. |
| `value` | National count. Cumulative for the two histories, daily new for the suspected series, matching what each block means in the released manifest. |

`onset_curve_scanned.csv` carries the onset-by-vintage reporting triangle, one block of onset dates per vintage.

| Column | |
|---|---|
| `sitrep` | Vintage identifier, one per reporting snapshot. |
| `report_date` | Date that snapshot was reported. |
| `onset_date` | Symptom-onset date the counts on this row belong to. |
| `confirmed_alive`, `confirmed_dead`, `confirmed_total` | Confirmed cases with that onset date, as that snapshot recorded them. |

The filename is not free.
`load_observations` reads the triangle from exactly `onset_curve_scanned.csv` beside the manifest it is handed, and `load_onset_curve` returns an empty no-op rather than an error when that file is absent.
A fit whose manifest and triangle sit in different directories therefore drops the onset stream and says nothing about it, so `place_onset_curve` copies the triangle next to the manifest and refuses an empty or missing one.

This repository also ships its own `data/onset_curve_scanned.csv`, digitised from the situation-report figures.
That is a different construction of the same quantity under the same name.
Which one a fit reads is named by `--data`, not inferred from whichever file is to hand: a situation-report fit handed the line-list triangle would be a mixture, and nothing downstream would show it.

The released observation manifest is the third input.
The scripts fetch it from this repository's own tagged results releases; `LINELIST_RELEASE` pins a tag and the default is the newest.
`LINELIST_RELEASED_MANIFEST` overrides with a local path instead, which is how the synthetic fixture runs without network access, and how a run pins itself to the working tree's own `data/observations.toml`.

`LINELIST_OUT_DIR` sets where the manifest and every result are written.
It defaults to `ignore/linelist/`, and any other path inside the repository is refused: `ignore/` is the only git-ignored output root, and this repository is public.

## How the case streams are indexed, and why it dominates

`--data` selects one of two observation sets.

| `--data` | Streams | Triangle |
|---|---|---|
| `sitrep` | The released manifest unchanged, with only its cut-off moved | `data/onset_curve_scanned.csv` |
| `linelist_known` | `linelist_streams_known.csv`, each case counted at the snapshot that first held it | the line list's |

Both are dated the way the situation reports date theirs: a case enters a series when the reporting system first records it as that kind of case.
A situation report published on a date prints the cumulative total known on that date, so `linelist_known` is indexed the same way the released manifest is, and the two are on the same footing.

An event-date construction, counting each case by the date its alert was notified, was tried and dropped.
For the confirmed stream it is right-truncated by both the notification-to-confirmation interval and the lag to the record reaching the export, so its last retained day holds about a quarter of the cases that will eventually carry that date, and the model has no truncation correction to tell that from a fall in transmission.
On the August data it gave a cut-off reproduction number around 0.55 against about 1.12 for the known-by construction of nearly the same cases: an artefact of the indexing, larger than any plausible difference between the two data sources.

## The breakpoint

The intervention breakpoint is the grid day before the first `reported_case_history` vintage: `default_breakpoint(obs)` is `obs.n - obs.who_first_sitrep_days`, and `who_first_sitrep_days` is `n - reported_history.days[1] + 1` (`src/data.jl`).

`reported_case_history` is one of the three streams the line list replaces, so taking the breakpoint from each run's own manifest moved it with the data source, even for `onsets`, which never conditions on that stream.
On the August data the released manifest's first reported vintage is 2026-05-18 and the line list's is 2026-05-26, which put the two arms eight days apart at breakpoints 94 and 102.

The breakpoint is a fact about the intervention timeline rather than about which series is being fitted, and the released manifest's first reported vintage is the first WHO joint situation report where the line list's is only the date its export begins.
So it is taken from the released manifest for every data mode, and both arms sit at 94.
For `--data=sitrep` that is the value the fit would have used anyway, so nothing about that arm changes.
`--breakpoint=N` overrides it.

It is part of the chain cache key, because it changes the model and a chain fitted at one breakpoint must not be reused under another.

The difference in triangle size between the two sources is deliberately not corrected.
The line-list triangle carries 904 cells over 35 loaded vintages against 419 over 19: that is part of what the comparison measures, not a confound to remove.

`LINELIST_AS_OF` pins the cut-off, and every run in a comparison must pass the same one.
The two constructions end on different days, so left to their own last day they sit on different grids and their trajectories are not on one axis.
It is required for `--data=sitrep`, which has no replacement streams to take a cut-off from.

## Delay configurations

`--delays` selects which delay distributions the fit is given.

| | Generation interval | Onset-to-report |
|---|---|---|
| `repo` | package default | package default |
| `cmmid_gi_any`, `cmmid_gi_case`, `cmmid_gi_diag` | cmmid, by transmission-pair definition | package default |
| `cmmid_rep` | package default | cmmid |
| `cmmid_any`, `cmmid_case`, `cmmid_diag` | cmmid, by transmission-pair definition | cmmid |

`repo` passes no overrides at all, so the package defaults stand by construction rather than by a restatement here that could drift from them.

The `onsets` fit reaches the generation interval and the incubation period, and then estimates its reporting delay from the triangle rather than taking a prior.
That is why `repo` is the configuration the comparison rests on: at `repo` the fit takes no cmmid estimate at all, so a difference between the two sources owes nothing to a generation interval fitted from line-list transmission pairs.
The `cmmid_gi_*` runs are the sensitivity around that, not the result.
The incubation period is not swapped: cmmid's own fits are unidentified on 17 to 64 contacts and it recommends the MacNeil et al. (2010) estimate, which is already this repository's prior.
Onset-to-death, the laboratory receipt delay and the treatment stays are not reached by either fit.

Running `cmmid_gi_*` and `cmmid_*` on the same data is what separates the two effects.
The cmmid onset-to-report delay is about nine days longer than the package's, which moves the trajectory against calendar time as well as changing its level, so a configuration that changes it and the generation interval together cannot say which did the work.

`cmmid_rep` is the other single-axis configuration: the report delay alone, the generation interval left at the package default.
It has no cell in the grid as it stands.
A report-delay override only means something to a fit that takes one as a prior, which is `confirmed`, and `confirmed` is out of the comparison; `onsets` refuses it.
It is kept, and tested, because it is what the grid would use if `confirmed` were ever usable on this data.

A `cmmid_*` configuration handed to the `onsets` fit is refused rather than silently ignored.

Two things to carry into any write-up.
The cmmid standard errors come from large samples, so these priors are near-fixed and carry almost no delay uncertainty into R_t; `BVD_DELAY_PRIOR_INFLATE` widens both spreads by a constant factor for a run that asks what the answer owes to the delay being exactly this.
And the cmmid onset-to-report kernel's 98% support is about 53 days against 16 for the package's, so every convolution using it is three times longer and the fit is correspondingly slower.

### Where the numbers come from

`BVD_DELAY_DIR` points at the `results/` directory of `bvd-internal-cmmid`, which is private.
Its disclosure rules permit fitted distribution parameters to be committed there and not shared onward, and this repository is public, so `delays.jl` reads them at run time and holds none of them.
`test/fixtures/linelist/delays/` holds invented numbers in the same shape, which is what the tests run against.

Each run writes `delay_provenance_<config>.csv` beside its outputs, recording every derived parameter and the SHA-256 of each file it came from, so a chain on disk can be traced back to the estimates that produced it.

## Running

```bash
export LINELIST_INPUT_DIR=$PWD/ignore/linelist/inputs
export LINELIST_RELEASED_MANIFEST=$PWD/data/observations.toml
export BVD_DELAY_DIR=<bvd-internal-cmmid>/results
export LINELIST_AS_OF=2026-08-10

scripts/linelist/run_grid.sh --stage=1 --dry-run   # what it would run
scripts/linelist/run_grid.sh --stage=1             # data source moves
scripts/linelist/run_grid.sh --stage=2             # delay priors move
julia --project=docs scripts/linelist/plot_rt.jl
```

Or one fit at a time:

```bash
julia -t 2 --project=docs scripts/linelist/fit_single.jl confirmed \
  --data=linelist_known --delays=repo
```

`-t 2` matters.
`nuts_sample` runs its chains with `MCMCThreads`, so a single thread runs them one after the other and doubles the wall clock.
Two is also the ceiling worth giving it, since the fit registry sets two chains.

Sampler settings come from this repository's own fit registry rather than being set here, so a refit stays comparable with the release it is compared against.
The only thing a run changes about the model is its delay configuration.

Each fit writes its chain to `chains/<fit>_<data>_<delays>_<as_of>.jls` and reuses it on rerun, so a correction to the trajectory export or the plots costs nothing rather than another fit.
`--refit` forces past it, and the key carries the cut-off so a rerun at new data refits rather than silently reusing an old grid.

Outputs are tagged `<fit>_<data>_<delays>`, so nothing in a grid overwrites anything else:

| | |
|---|---|
| `linelist_<tag>_rt.csv` | the daily R_t trajectory, median and 30/60/90% bounds, one row per established grid day |
| `linelist_<tag>_rt_draws.csv` | thinned per-draw daily R_t, so a downstream summary comes from draws rather than from these intervals |
| `linelist_<tag>_stream_estimates.csv` | the cut-off and basic reproduction numbers and the final size, in the shape the release's own `stream_estimates.csv` uses |
| `linelist_<tag>_delay.csv` | the fitted onset-to-report delay, for the fits that estimate one |
| `linelist_<tag>_diagnostics.csv` | worst R-hat, smallest bulk ESS, divergent-transition count |

## Reading the comparison

`plot_rt.jl` collects those per-run files and reduces them to the contrast.
It refits nothing, so a number can only change by rerunning the grid.

| | |
|---|---|
| `comparison.csv` | the file to read first: one row per fit, delay configuration and quantity, with the two data sources side by side, their difference and their ratio |
| `rt_confirmed.png`, `rt_onsets.png` | the trajectories overlaid, one series per data source and delay configuration |
| `delay_comparison.csv` | the onset-to-report delay each `onsets` fit estimated, against cmmid's independent fit of the same interval |
| `diagnostics_all.csv` | sampler quality per run, with a warning on stderr naming any run past the thresholds |
| `inputs_at_cutoff.csv` | what each observation set actually held at the cut-off, read back out of the manifest the fit was given |
| `rt_all.csv`, `rt_summary.csv`, `stream_estimates_all.csv` | the collected per-run detail the above is built from |

`comparison.csv` carries `Rt_cutoff` and the same quantity 7, 14 and 28 days earlier, then `R_T`, `R0` and `C_T`.
Only the medians are differenced.
The two fits are the same model on two constructions of one outbreak rather than independent samples of a population, so the difference describes how far apart the answers sit and is not a test; the per-source 90% bounds are carried so the gap can be read against the width of either fit.

Two guards run before anything is written.
Every series within a fit must end on the same date, since a series on a different grid is a run at a different cut-off and pairing it would contrast two questions.
And a configuration run on one source but not the other is left out rather than half-reported, so a row in `comparison.csv` always has both sides.

Read `diagnostics_all.csv` before `comparison.csv`.
A difference between two data sources is only attributable to the data if both sides sampled, and they have not sampled equally well: on every pair run so far the `linelist_known` arm records more divergent transitions than its `sitrep` pair, and one `confirmed` and one `onsets` run are past 1% of their draws.
The thresholds the warning fires on are R-hat above 1.01, bulk ESS below 400, and divergences above 1% of the run's own draws.

`cases` is excluded from the figures and from `comparison.csv`, and its rows are kept in `rt_all.csv` and `rt_summary.csv`; the script says so on stderr when it finds them.

To exercise the whole path without the real data, point the inputs at the synthetic fixture:

```bash
LINELIST_INPUT_DIR=test/fixtures/linelist \
LINELIST_RELEASED_MANIFEST=data/observations.toml \
LINELIST_OUT_DIR=$(mktemp -d) \
BVD_DELAY_DIR=test/fixtures/linelist/delays \
julia -t 2 --project=docs scripts/linelist/fit_single.jl confirmed \
  --pilot --data=linelist_known --delays=cmmid_any
```

The fixture's counts are invented, so the estimates mean nothing.
`test/test_linelist_manifest.jl` runs the manifest substitution, the triangle placement and the delay arithmetic on it in seconds rather than minutes.

## What these fits can and cannot say

`onsets` is single-stream, so ascertainment is not identified.
Level differences between the data sources land in `C_T` and `p_drc`; only shape differences reach R_t.
So the trajectory comparison is the one to read, and `C_T` from these fits is not a size estimate: the `confirmed` fit on situation-report data puts it at 52,141 against the release's joint 9,415.

The `onsets` fit estimates the onset-to-report delay rather than assuming it, and writes it.
On line-list data it is the same interval `bvd-internal-cmmid` fits independently, so the two are worth reading against each other; they have agreed at the median so far.

## Disclosure

This repository is public.
A case line list is individual patient records and must not reach it.

Nothing here reads a line list.
These scripts read files of national aggregate counts, produced elsewhere, from a directory outside the repository, and delay parameters from a private repository's results directory.
Neither those inputs nor any fit output derived from them at case level should be committed, and no delay parameter value belongs in a tracked file here.

The boundary as it stands: national counts by date, and by onset date and vintage, are what crosses into this repository as inputs.
Fitted parameters and posterior summaries are publishable.
Anything stratified more finely is not, and moving that boundary outward is a decision to take deliberately rather than by adding a file.

## Scope

`fit_joint.jl` is a partial refit and is not part of the R_t comparison.
Only three of the manifest's streams have a line-list source: cumulative confirmed cases, cumulative reported cases, and daily new suspected cases.
Deaths, laboratory volumes, isolation and treatment beds, recoveries, Uganda exports and travel have no line-list counterpart, since a case line list records no outcome.
Those streams stay as the situation reports gave them, so a joint refit mixes the two sources and cannot attribute a difference to either.

The `onsets` fit has no such mixture: it conditions on the reporting triangle alone, which the line list supplies in full.

### Why `confirmed` is not in the comparison

It does not sample on line-list data: max R-hat 2.63 and a bulk ESS of 2.4 out of 1000 draws, with zero divergences.
Both chains freeze, one pinned at the `inv_sqrt_k = 0` truncation boundary of `surveillance_dispersion_model` where `k` is `1/eps` exactly.
The reading is that the line-list confirmed history, 35 sparse vintages whose increments are multi-day sums, carries less dispersion than the negative-binomial observation model can express.

`scripts/linelist/CONFIRMED_FIT_ISSUE.md` has the measurements, what was ruled out and how, what is still open, and how to reproduce it.
Read that before proposing a fix.

`confirmed` is a partial exception and should be read as one.
It conditions on the confirmed-case history, which the line list supplies, and on the laboratory volumes, which it does not: `tests_analysed_history` and `tests_analysed_daily_history` are outside `LINELIST_BLOCKS` and so are byte-identical in both manifests.
The fit therefore pairs line-list confirmed counts with situation-report test denominators in the positivity likelihood.
At the cut-off the two sources differ by about 13% on the numerator (3,862 against 4,449) with the denominator held fixed.
No day implies a positivity above 100% under either source, so this is not on its own a broken likelihood, but a difference from the `confirmed` fit is a difference between confirmed-case series measured against one shared laboratory denominator, not between two independent constructions of the outbreak.
