# Scripts

Helper scripts for the analysis.
Each runs against a specific Julia project.
Run it as `julia --project=<project> scripts/<name>.jl` from the repository root, unless noted otherwise below.

## Data pipeline (`--project=scripts`)

These are the scripts to reach for when a new SitRep lands, run in this order.
See `data/README.md` for the full data-update procedure, including the manual transcription steps these scripts do not cover.

| Script | What it does |
| --- | --- |
| `check_new_sitreps.jl` | Lists INSP SitReps not yet in `data/insp_sitrep_scanned.csv`. Exits non-zero if any are missing. |
| `download_sitreps.jl` | Downloads the INSP SitRep PDFs into `data/sitrep_pdfs/` (git-ignored). Also `task download-sitreps`. |
| `confirm_insp_data.jl` | Cross-checks the scanned confirmed-case and confirmed-death totals against the INRB-UMIE mirror, and reports the dates each source carries alone. Also `task confirm-data`. |
| `refresh_releases.jl` | Pulls each tagged results release's headline estimate into `data/released_estimates.csv`. Also `task refresh-releases`. |

`check_new_sitreps.jl`, `download_sitreps.jl` and `confirm_insp_data.jl` need no Julia packages beyond `Downloads`.
`refresh_releases.jl` also needs the `gh` CLI, authenticated against the repo.

### Onset-curve digitiser

The analytique SitReps carry a raster bar chart of confirmed cases by symptom-onset date and no data table.
`digitize_onset_curve.jl` reads the daily counts off that figure's pixels and writes one block per vintage to `data/onset_curve_scanned.csv`.
It needs no Julia packages beyond stdlib, plus poppler's `pdfimages`, `pdftotext` and `pdfinfo` on `PATH` (`apt install poppler-utils`, or `brew install poppler`).
`digitize_onset_curve.py` is its Python port for the automated data-updater, which has Python but not Julia access, run with `uv run scripts/digitize_onset_curve.py` (Pillow and numpy come from its PEP 723 metadata).
`audit_onset_curve.jl` checks the digitised file against the figures and `onset_check_panels.py` writes the crops for a read by eye.
`extract_dashboard_onsets.py` reads the INRB-UMIE dashboard's own onset charts into `data/onset_dashboard_history.csv`, the independent series the scan is checked against.
The Taskfile wraps the three regular runs as `task onset-digitise`, `task onset-port-check` and `task onset-audit`.

`audit_onset_curve.jl` writes `data/onset_curve_figures.csv`, one row per vintage with the figure's page, image size and md5, the pixel scales the digitiser calibrated, the n printed in the figure title, the digitised total, the gap between the two and the earlier vintage it reprints.
It also writes `output/onset_curve_audit.md` with the gap table, the settled-bar checks between consecutive snapshots (bars that fell, net change and the L1 distance at day shifts of up to two) and the worst vintages by gap and by settled-bar distance.
The printed n is part of the raster, not the PDF text layer, so it is read with tesseract from the title strip of the embedded image and checked against the source line.
Pass `--crops=DIR` to keep the upscaled strips each n was read from.
Each row also records where the figure's axis starts, where the digitiser's day loop starts and the first onset date in the block, with the cases the block leaves unread between them.
Run it with `task onset-audit` or `julia scripts/audit_onset_curve.jl`.

#### Adding a vintage

Run `task download-sitreps` so the new PDF is in `data/sitrep_pdfs/`.
Find the figure page with `pdftotext -layout` (the caption reads "date de début des symptômes") and extract the embedded image with `pdfimages -f P -l P <pdf> <prefix>`.
Open the extracted image and read the chart's own title, not the page caption.
The title must say "par date de début des symptômes" and print the `n`; SitRep 110's chart said "par date de notification" under the same caption and is left out.
Write down the printed `n`.
Read the y-axis labels: a 0/25/50/75 grid needs a `Y_AXIS_STEP` entry of 25 in both scripts, a 0/20/40/60 grid needs nothing.
Crop the right end of the x-axis at 8x (Pillow's `crop` and `resize` with nearest-neighbour sampling) and have two readers each write down the date under the rightmost weekly tick without seeing the other's answer or the previous vintage's `CONFIG` row.
Accept the date only when the two readings agree.
When they disagree, crop again at 12x and read again; if they still disagree, stop and do not add the vintage.
Append the `CONFIG` row (SitRep number, rapportage date, last-tick date) to both scripts, keeping the two tables identical, and add the y-axis step if it is 25.
Run `task onset-digitise`.
It is incremental, so it opens only the new PDF, and it prints the number of onset days and the total it read.
Run `task onset-port-check` to confirm the port reproduces the file.
Run `task onset-audit` and read the new vintage's rows in `output/onset_curve_audit.md`.

The block is accepted when all of the following hold.
In the gap table, `gap %` is within 2.1% of the printed `n` (the reader sits within 2.1% on every one of the 60 figures that print one, within 0.5% on 43 of them).
When the OCR left `printed n` empty, compare the digitised total against the `n` written down above by hand.
In the settled-bars table, the new pair's `best shift` is 0 and its `L1 0` column is the smallest of the five L1 columns.
Its `fell` count is a handful at most, each fall a single case: settled bars can only rise, so falls are scan noise and the noise floor is 0.40 cases per bar-day with a fall on 14% of days.
In the axis-coverage table, `cases before` is 0, so the block starts where the axis starts.
`reprint of` names an earlier vintage only when the figure is pixel-identical to that vintage's (same `image_md5` in `data/onset_curve_figures.csv`).

When a check fails, the row says which input is wrong.
A `best shift` of +1 or -1 means the new last-tick date or the previous vintage's is off by a day: re-read both ticks with two fresh readers before anything else.
A gap beyond 2.1% with shift 0 means the count scale: check `Y_AXIS_STEP`, then the printed `n` itself (the OCR misreads a digit now and then; the `note` column says when the title and source strips disagreed), then look at the check panels.
Many falls with shift 0 and a gap inside 2.1% mean the day grid inside the block has moved: compare `pixels_per_day` and `pixels_per_day_fit` for the vintage in `data/onset_curve_figures.csv` and run the vision check.
`cases before` above 0 means the reader's day loop started after the axis: the loop runs from a week before the first chain tick, so a lost tick at the left end is the usual cause.
A block that fails after the ticks and the scale have been re-read is not committed; remove its `CONFIG` row and open an issue with the audit rows.

#### Vision check

`uv run scripts/onset_check_panels.py SR` writes `output/onset_panels/blind_SR_k.png` and `check_SR_k.png`, one pair per twelve days, each an 8x crop with the day of month under every bar and a count ruler every five counts drawn from the digitiser's own calibration.
`--days`, `--scale`, `--out` change the panel width, the magnification and the directory.
Read the printed `n` from the title strip and five bars from the blind panels against the ruler, choosing bars across the range of heights and at least one in the pink incomplete-data band, before opening the check panels or the CSV.
Then compare with the check panels, where the green line is the digitised total and the magenta line the alive count.
A bar read more than one count away from the digitised value, or an `n` that disagrees with the one the audit read, refuses the vintage.
Run this check on every new vintage whose render size or layout differs from the previous one (`width` and `height` in `data/onset_curve_figures.csv`), and on every vintage after a change to the digitiser.

#### Rebuilding everything

Run `task onset-digitise -- --rebuild` to re-read every vintage, then `task onset-port-check` and `task onset-audit`.
A rebuild is required after any change to either digitiser, since an incremental run reuses the committed rows and would hide the change.
Keep the previous `output/onset_curve_audit.md` and compare the gap and settled-bars tables before and after: a change is an improvement only when the gaps tighten and the falls drop across vintages, not on the vintage that motivated it alone.
Commit the CSV, `data/onset_curve_figures.csv` and the two scripts together.

#### Keeping the port byte-identical

The Julia script is the reference and the committed CSV is its output.
The port mirrors it function by function, under the same names, with the same pixel classes, thresholds, tie-breaks and rounding; a disagreement is a fault in the port.
Change the reference first, then carry the same change into the port, then run `task onset-port-check`, which rebuilds the file with both and diffs them against each other and against the committed CSV.
`test/test_onset_digitiser.jl` runs the same comparison whenever `data/sitrep_pdfs` is present.
Three conventions carry the exact match.
The reference indexes pixels from 1 and the port does its bar-window arithmetic in that frame, converting to 0-based only when indexing, because round-half-to-even is not translation-invariant and a window edge on exactly .5 would otherwise land one column off.
Sums that feed the pixels-per-day slope are of integer-valued floats, so their order does not matter; a change that puts non-integer terms into a sum must fix the summation order in both scripts.
Ties are broken the same way on both sides: the first row with the longest run, the leftmost y-axis strip, the first-seen value in the modal height, and the floored mean in tick clustering.
When a cell differs, print both scripts' per-column run heights for that vintage and column window and chase the first differing pixel class or rounding; do not change the reference to match the port.

#### Failure modes to recognise in the audit

Each of these was found on the September renders and each shows in the audit tables before it shows in the fit.
A loop bound fixed at 105 days: `cases before` is above 0 and grows with each vintage as the axis extends back, `first onset` sits weeks after `axis start`, and the gap turns negative and widens vintage by vintage.
Washed JPEG chroma losing the dead segment: `confirmed_dead` is near zero on bars whose crimson segment is visible in the panels, the gap drops to -10% to -23% from the small renders at SitRep 106 on, and the alive count rises while the dead count falls between vintages.
A flood capped at the top y-axis tick: every bar taller than the top gridline reads the same height, the tallest bars fall between vintages whenever the render's axis range changes, and the gap is negative on vintages whose peak stands above the top tick.
Day-grid drift from an integer tick spacing: `best shift` leaves 0 on pairs whose ticks were verified, the preferred shift alternates direction between consecutive pairs, and `drift_days` in `data/onset_curve_figures.csv` is non-zero.

#### Dashboard cross-check

`data/onset_dashboard_history.csv` holds the INRB-UMIE dashboard's national and province onset curves, one block per dashboard build, read from the inline SVG charts by `extract_dashboard_onsets.py` (see `data/README.md` for the columns and the clone it reads from).
It is the one independent series on the same basis.
For a new vintage, take the dashboard block whose `snapshot_date` is nearest the report date, join the national `observed` column on onset date to the scan block, and check the correlation and the mean absolute difference per day against the values recorded in `data/README.md` (r above 0.96, mean absolute difference under 4 cases a day).
A vintage that agrees with the printed `n` but not with the dashboard has its days shifted or its scale wrong, and the vision check says which.

## Publishing the results (`--project=docs`)

| Script | What it does |
| --- | --- |
| `score_releases.jl` | Scores every past release's saved forecasts against the now-observed data and refreshes the scoring and per-release R_T/C_T/R0 overlay CSVs. |

It also runs under `--project=.`.
CI uses `--project=docs` because that environment is already instantiated at that point in the build.

## Entry points and reproduction

| Script | Project | What it does |
| --- | --- | --- |
| `run.jl` | `.` | Regenerates the published results by running the analysis and sensitivity pages. |
| `reproduce.jl` | none | Bootstraps a full reproduction from a fresh clone; run with `curl -fsSL https://raw.githubusercontent.com/epiforecasts/BVDOutbreakSize/main/scripts/reproduce.jl \| julia`. |
| `backfill_forecasts.jl` | none | Reconstructs the one-week-ahead forecast each past release made but never saved, for `score_releases.jl` to score. Needs only `Dates` itself: it checks out each release tag into its own worktree and runs there under that tag's own project. |
| `backfill_drivers/driver_v1.0.0.jl`, `driver_v1.1.0.jl` | that tag's own `docs` | Standalone drivers for the two release tags whose model predates the fit registry. Run by `backfill_forecasts.jl`, not directly. |

## Formatting

`run_formatter.sh` runs Runic over `src/`, `test/`, `docs/`, `scripts/`, `benchmark/` and `ext/` from its own isolated `test/formatter/` sub-environment.
`task format` and `task lint` both call it directly.
It takes no project flag.
The pre-commit hook does not: pre-commit builds the hook its own environment from the Runic version in `.pre-commit-config.yaml`, and `test/package/CodeFormatting.jl` checks that version against the pin in `test/formatter/Project.toml`.

## Developer diagnostics (`--project=.`)

These are not wired into CI, the Taskfile or any workflow.
Reach for them when a fit misbehaves or a hot path changes.
Each documents its own invocation in its header comment.

| Script | What it does |
| --- | --- |
| `bench_convolve.jl` | Times the Mooncake gradient of the renewal convolution, scalar loop against a vectorised rewrite. Its finding is cited from a docstring in `src/renewal.jl`, which justifies keeping the simpler scalar loop. Run it if you touch that code path. |
| `bench_discretise.jl` | Times the Mooncake gradient of the censored delay discretisation. Documents the reasoning behind the CDF-difference form in `src/renewal.jl`. |
| `summarise_chain.jl` | Prints posterior summaries and fit diagnostics for a saved chain. |
| `prior_vs_posterior.jl` | Samples the prior and sets it beside a saved posterior, so a parameter the data does not inform is visible. |
| `fit_joint_stream.jl` | Fits the full joint model outside the fit cache and CI, streaming live progress to `logs/joint_fit.log` so a long fit can be watched with `tail -f`. Writes the chain to `logs/joint_chain.jls`. The name refers to the streamed progress, not to fitting a single data stream. |
