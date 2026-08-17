# Refitting on line-list case counts

This model reads national counts from the INSP situation reports.
These scripts refit it on case counts derived from a DHIS2 case line list instead, holding the model fixed and changing only the data.
A difference against the released fit of the same id is then a difference between data sources rather than between models.

## What is here

| File | |
|---|---|
| `manifest.jl` | Resolves the inputs and substitutes the line-list case streams into a released observation manifest. Shared by both fits so they cannot drift. |
| `fit_joint.jl` | The headline joint fit, exactly as the released report runs it. Five to seven hours. |
| `fit_single.jl` | One stream at a time: `cases` (about 14 minutes) or `onsets` (about an hour and a quarter). |
| `run_joint_fit.sh` | Detached wrapper for `fit_joint.jl`. Use it rather than calling Julia directly, since a closed terminal or a sleeping laptop is the likeliest way a seven-hour run ends. |

## Inputs

Two files, produced by whoever holds the line list, in a directory named by `LINELIST_INPUT_DIR`.
There is no default: these are counts derived from individual patient records, so the operator says where they are rather than a default finding a directory that happens to exist.

`linelist_streams.csv` carries the case streams that replace the situation-report ones.

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
A fit whose manifest and triangle sit in different directories therefore drops the onset stream and says nothing about it, so `manifest.jl` copies the triangle next to the manifest and both fits refuse to start without it.

This repository ships its own `data/onset_curve_scanned.csv`, digitised from the situation-report figures.
That is a different construction of the same quantity under the same name, and it is not what these scripts fit.

The third input is the released observation manifest, which the scripts fetch themselves from this repository's own tagged results releases.
`LINELIST_RELEASE` pins a tag; the default is the newest.
`LINELIST_RELEASED_MANIFEST` overrides with a local path instead, which is how the synthetic fixture runs without network access.

`LINELIST_OUT_DIR` sets where the manifest and every result are written.
It defaults to `ignore/linelist/`, and any other path inside the repository is refused: `ignore/` is the only git-ignored output root, and this repository is public.

## Running

```bash
# joint fit, detached
LINELIST_INPUT_DIR=<dir> scripts/linelist/run_joint_fit.sh

# one stream, foreground
LINELIST_INPUT_DIR=<dir> julia -t 2 --project=docs scripts/linelist/fit_single.jl cases
LINELIST_INPUT_DIR=<dir> julia -t 2 --project=docs scripts/linelist/fit_single.jl onsets
```

`-t 2` matters.
`nuts_sample` runs its chains with `MCMCThreads`, so a single thread runs them one after the other and doubles the wall clock.
Two is also the ceiling worth giving it, since the fit registry sets two chains.

Sampler settings come from this repository's own fit registry rather than being set here, so a refit stays comparable with the release it is compared against.

To exercise the whole path without the real data, point both variables at the synthetic fixture:

```bash
LINELIST_INPUT_DIR=test/fixtures/linelist \
LINELIST_RELEASED_MANIFEST=data/observations.toml \
LINELIST_OUT_DIR=$(mktemp -d) \
julia -t 2 --project=docs scripts/linelist/fit_single.jl cases --pilot
```

The fixture's counts are invented, so the estimates mean nothing.
`test/test_linelist_manifest.jl` runs the same substitution on it and checks the result, and takes seconds rather than minutes.

## Disclosure

This repository is public.
A case line list is individual patient records and must not reach it.

Nothing here reads a line list.
These scripts read two files of national aggregate counts, produced elsewhere, from a directory outside the repository.
Neither those inputs nor any fit output derived from them at case level should be committed.

The boundary as it stands: national counts by date, and by onset date and vintage, are what crosses into this repository as inputs.
Fitted parameters and posterior summaries are publishable.
Anything stratified more finely is not, and moving that boundary outward is a decision to take deliberately rather than by adding a file.

## Scope

The joint refit is partial.
Only three of the manifest's streams have a line-list source: cumulative confirmed cases, cumulative reported cases, and daily new suspected cases.
Deaths, laboratory volumes, isolation and treatment beds, recoveries, Uganda exports and travel have no line-list counterpart, since a case line list records no outcome.
Those streams stay as the situation reports gave them.

The single-stream fits are the comparison worth trusting.
Each conditions on a stream the line list supplies in full, so a difference against the release's own fit of the same id is a difference between data sources.
The joint refit cannot say that, since most of its streams have no line-list source.

Counting each case by notification date leaves the last fortnight short by however long cases take to reach the export, and the model reads that ragged edge as a fall in transmission.
A construction that counts each case at the snapshot which first held it has no such edge.
`fit_single.jl` takes one via `--streams`, which tags the manifest and every output with the construction's name so two of them never overwrite each other.
