# Paper figures

Every figure is drawn by a script under `paper/scripts/` from the assets of one `results-*` GitHub release of epiforecasts/BVDOutbreakSize, so a figure can be regenerated from the release tag alone.
The scripts run in the `paper/` environment (`paper/Project.toml`); instantiate it once with `julia --project=paper -e 'using Pkg; Pkg.instantiate()'` from the repository root.
The release assets are downloaded by `gh` into `paper/data/release/` (git-ignored) on the first run; the shared helpers, the default tag and the asset list live in `paper/scripts/release_assets.jl`.

## fig-current (Figure 2)

`paper/scripts/fig_current.jl` draws the current model's outputs at the pinned release into `fig-current.pdf` and `fig-current.png` (180 mm wide, 600 dpi).
Inputs: `onsets_over_time.csv` (panel A bands), `observations.toml` (the confirmed-case and confirmed-death histories in panels A and D and the cut-off date), `posterior_summary.csv` and `posterior_draws.csv` (national reproduction number and ascertainment intervals and medians in B and C), `stream_estimates.csv` (the national reproduction-number median in B), `forecast_validation.csv` and `forecast.csv` (panel D), and from `site.zip` the province tables on `estimates/province.html` (panel B) and the onset-report table on `estimates/national.html` (panel C).
The prior medians in panel C are read from `src/models/priors.jl` at the release's source commit through `git show`.
Regenerate with `julia --project=paper paper/scripts/fig_current.jl [tag]` from the repository root; the tag defaults to the one in `release_assets.jl`.

## numbers.yml

`paper/scripts/paper_numbers.jl` is not a figure script but shares the inputs above.
It writes `paper/generated/numbers.yml`, the quoted numbers the manuscript reads through Quarto meta shortcodes, each with a comment naming its source column.
Regenerate with `julia --project=paper paper/scripts/paper_numbers.jl [tag]`.

## fig-first (Figure 1)

`paper/scripts/fig_first.jl` draws the first model's outputs at the `results-v1.0.0` release (data cut-off 18 May 2026) into `fig-first.pdf` and `fig-first.png` (180 mm wide, 600 dpi).
It needs only CairoMakie, Dates and Statistics, so it runs in the root environment (`julia --project=. paper/scripts/fig_first.jl`) as well as in `paper/`; it downloads its assets into `paper/data/release/results-v1.0.0/` itself.
Inputs: `imperial_comparison.csv` (panel A: the two McCabe et al. headline scenarios with their reported 90% intervals, the Method 2 reproduction, and the joint fit to the report's 16 May data and to the 18 May data, each a median with a 90% interval), `scenario_coverage.csv` (panel A: the fifteen McCabe et al. scenario point estimates), `cumulative_cases_by_stream.csv` (panel B: the 30, 60 and 90% intervals of cumulative infections from each single-stream fit and the joint; the release carries no medians for these), `posterior_draws.csv` (the 400 thinned joint draws: the 50% interval of the headline row in A, and the growth rate and time since seeding per draw in C) and `observations.toml` (the four fitted counts named in B and marked in C, and the cut-off).
Panel C evaluates the released model's own latent curve, cumulative infections `exp(r s)` from a single seed `T` days before the cut-off, per draw on a daily grid and shows pointwise quantiles; the seeding-date marginal is the density of cut-off minus `T` over the same draws.
The release's posterior-predictive check for the four streams exists only as an embedded figure in its rendered report, and the thinned draws omit the detection window, traveller volume, delay and dispersion parameters, so that check is not redrawn here; panel B shows the single-stream fits instead.

## fig-journey (Figure 3)

`paper/scripts/fig_journey.jl` draws the development journey into `fig-journey.pdf` and `fig-journey.png` (180 mm wide, 600 dpi) from the CSVs in `paper/data/` described in `paper/data/README.md`, not from a release.
Unlike the other figure scripts it runs in the root environment, `julia --project=. paper/scripts/fig_journey.jl`, and needs only CairoMakie and the standard library.
Four panels share one date axis from 18 May 2026 to the end of the last week in `commits_weekly.csv`, with every release tag as a thin vertical line, labelled in the strip above panel A, and the three model versions as background bands.
Panel A steps `loc_src` and `loc_test` from `code_size.csv` between tags, with `n_streams_fitted` as a dashed step on the right axis.
Panel B stacks the weekly human, agent and other commits from `commits_weekly.csv` and draws the weekly total of merged pull requests as a line on the same count axis.
Panel C places each release's median and 90% interval from `release_estimates.csv` at its data cut-off, on a log axis; releases sharing a cut-off are drawn 0.7 days either side of it, and the three tags without a results release do not appear.
The strip under panel C marks each row of `data_events.csv` at its date, coloured by `human_decision`, with a hand-shortened label from the `event` column.
Panel D draws every `joint_gradient_ms` value in `fit_cost.csv` whose `model` is a production joint or the CI benchmark suite, open markers for the before arm and filled for the after arm of the same PR, joined by a line; the two `synthetic joint` rows of #810 are not drawn.
The one `joint_fit_minutes` value, 320 minutes at #716, is a labelled diamond whose height on the ms axis carries no value.
Colours are Makie's Wong palette, as in `src/plots.jl`.

## fig-evaluation (Figure 4)

`paper/scripts/fig_evaluation.jl` draws the evaluation into `fig-evaluation.pdf` and `fig-evaluation.png` (180 mm wide, 600 dpi).
It runs in the root environment, `julia --project=. paper/scripts/fig_evaluation.jl`, and needs only CairoMakie, Dates and TOML.
Panel A reads `data/forecast_scores.csv`, the cross-release scores `scripts/score_releases.jl` writes (last refreshed 4 September 2026 in PR #637, so the newest scored release is cut at 24 August).
It keeps the seven-day horizon only and divides the joint fit's CRPS by the persistence baseline's for every release and stream that carries both, so a value below one beats the baseline.
The still-reported streams are those with a baseline row: confirmed cases, confirmed deaths, isolation beds, onset reports and recovered; reported cases and suspected deaths have stopped being updated, carry no baseline and are absent.
Releases whose name carries the `(backfill)` label (v1.4.0, v1.5.0, v1.7.0, v1.8.0 and v1.9.0) are reconstructions of tagged model versions rather than forecasts published at the time, so they are drawn as hollow markers and left off the connecting lines; the v1.1.0 backfill scores only the two stopped streams and contributes nothing.
The strip under panel A marks whether the observation fell inside the joint's 90% interval, from the `coverage_90` column of the joint row.
Panels B to D read the assets of the pinned release (results-2416, data cut-off 22 September 2026) in `paper/data/release/` and the scenario constants in `src/constants.jl`.
Panel B is the `C_T` row of every fit in `stream_estimates.csv`, a median with nested 30, 60 and 90% intervals on a log axis.
Panel C is `REPORT_SCENARIOS_CI`, the McCabe et al. scenario means with their 95% confidence intervals at the 18 May, 20 May and 27 May cut-offs, against the joint fit frozen at 20 and 27 May from `frozen_matched_cutoffs.csv`, which publishes the 30, 60 and 90% bounds without a median.
The 18 May report has no frozen fit of its own, and the asset's 23 May and 8 June fits match no McCabe vintage, so they are not drawn.
The frozen joint estimates cumulative infections to the cut-off, including infections not yet symptomatic, whereas McCabe et al. estimate cumulative cases; the report's own comparison reads cumulative onsets off the current fit instead, a quantity the release does not publish.
Panel D is `CHAMLA_CONFIRMED_CENTRAL` and `CHAMLA_CONFIRMED_W12`, the confirmed-case projection of the joint fit frozen at 8 June read from the matched-date table on `sensitivity.html` inside the release's `site.zip` (it is not published as a CSV), and the observed confirmed cases from the release's `observations.toml`.
The script prints the numbers the text quotes: the stream-release pairs where the joint beat persistence at one week, the 90% coverage count, and the ratio of the largest to the smallest single-stream median.
