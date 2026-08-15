# Refitting on line-list case counts

Refits this model on case counts derived from the DHIS2 line list rather than
from the INSP situation reports it normally reads. The model is held fixed and
the data changed, so a difference against the released fit of the same id is a
difference between data sources, not between models.

## What is here

| File | |
|---|---|
| `manifest.jl` | Substitutes line-list case streams into a released observation manifest. Shared by both fits so they cannot drift. |
| `fit_joint.jl` | The headline joint fit, exactly as the released report runs it. Five to seven hours. |
| `fit_single.jl` | One stream at a time: `cases` (~14 min) or `onsets` (~1 h 15). |
| `run_joint_fit.sh` | Detached wrapper for `fit_joint.jl`. Use this rather than calling Julia directly — a closed terminal or a sleeping laptop is the likeliest way a seven-hour run ends. |

## Inputs

Produced outside this repo, by whatever has access to the line list. By default
they are read from a sibling `bvd-analysis` clone:

```
../bvd-analysis/ignore/report/
    bos_observations.toml         released manifest, from fetch_bos.sh
    bos_linelist_streams.csv      case streams built from the line list
    onset_curve_scanned.csv       onset-by-vintage reporting triangle
```

Point `LINELIST_INPUT_DIR` elsewhere to read them from anywhere else, and
`LINELIST_OUT_DIR` to write results somewhere other than `ignore/linelist/`.

`onset_curve_scanned.csv` is named for this repo's `load_observations()`, which
looks for exactly that filename beside the manifest it is handed. It is built
upstream as a plain reporting triangle and renamed on the way in. Note that
`load_onset_curve()` returns an empty no-op rather than an error when the file
is missing, so the stream disappears silently — `fit_single.jl` checks for it up
front rather than discovering it seven hours later.

## Running

```bash
# joint fit, detached
scripts/linelist/run_joint_fit.sh

# one stream, foreground
julia -t 2 --project=docs scripts/linelist/fit_single.jl cases
julia -t 2 --project=docs scripts/linelist/fit_single.jl onsets
```

`-t 2` matters: `nuts_sample` runs its chains with `MCMCThreads`, so a single
thread runs them one after the other and doubles the wall clock. Two is also the
ceiling worth giving it, since the fit registry sets two chains.

Sampler settings come from this repo's own fit registry rather than being set
here, so a refit stays comparable with the release it is compared against.

## Scope

The joint refit is partial. Only three of the manifest's streams have a
line-list source — cumulative confirmed cases, cumulative reported cases, and
daily new suspected cases. Deaths, laboratory volumes, isolation and treatment
beds, recoveries, Uganda exports and travel have no line-list counterpart; the
line list records no outcome, so deaths in particular cannot come from it. Those
streams stay as the situation reports gave them.

The single-stream fits are the comparison worth trusting: each conditions on a
stream the line list supplies in full, so a difference against the release's own
fit of the same id is a difference between data sources. The joint refit cannot
say that, since most of its streams have no line-list source.
