# Quoted numbers for paper/main.qmd, written to paper/generated/numbers.yml
# from the assets of one `results-*` release of epiforecasts/BVDOutbreakSize.
#
#     julia --project=paper paper/scripts/paper_numbers.jl [tag]
#
# `tag` defaults to `DEFAULT_TAG` in release_assets.jl. The assets are
# downloaded into paper/data/release/ (git-ignored) on the first run.
#
# Every key is a string already formatted for prose: two significant
# figures, thousands separators, intervals as "a to b". Each key's comment
# names the asset and column it comes from. A quantity the assets do not
# carry is omitted rather than estimated. Medians come from
# stream_estimates.csv where it carries the quantity and otherwise from the
# 200 thinned draws in posterior_draws.csv; the 90% bounds come from the
# full-chain summaries in posterior_summary.csv.

include(joinpath(@__DIR__, "release_assets.jl"))

using Distributions: Beta

const TAG = isempty(ARGS) ? DEFAULT_TAG : ARGS[1]

ensure_assets(TAG)
meta = release_meta(TAG)

summary = read_asset("posterior_summary.csv")
draws = read_asset("posterior_draws.csv")
streams = read_asset("stream_estimates.csv")
forecast = read_asset("forecast.csv")
validation = read_asset("forecast_validation.csv")
province_forecast = read_asset("province_forecast.csv")
obs = TOML.parsefile(asset_path("observations.toml"))
cutoff = Date(obs["as_of_date"])

## One entry per key: (key, value, source comment).
entries = Tuple{String, String, String}[]
add!(key, value, source) = push!(entries, (key, value, source))

## Lower/upper 90% of a posterior_summary.csv quantity.
function bounds90(q)
    r = summary_row(summary, q)
    return (r["Lower 90%"], r["Upper 90%"])
end

fmt_int(x) = fmt_count(round(Int, x))
fmt_day(d) = Dates.format(d, "d U Y")

## --- Release ---------------------------------------------------------------

add!("release_tag", meta.tag, "gh release view: tagName")
## Only the short commit: detect-secrets flags a full 40-character SHA.
add!(
    "release_commit_short", meta.commit[1:7],
    "gh release view: source commit, first seven characters"
)
add!("release_date", fmt_day(meta.date), "gh release view: publishedAt")
add!("data_cutoff", fmt_day(cutoff), "observations.toml: as_of_date")
add!("data_cutoff_iso", string(cutoff), "observations.toml: as_of_date")

## --- Headline size, growth and severity ------------------------------------

ct = stream_row(streams, "joint", "C_T")
ct90 = bounds90("C_T")
add!(
    "cumulative_infections_median", fmt2(ct.median),
    "stream_estimates.csv: fit joint, quantity C_T, median"
)
add!(
    "cumulative_infections_lower90", fmt2(ct90[1]),
    "posterior_summary.csv: C_T, Lower 90%"
)
add!(
    "cumulative_infections_upper90", fmt2(ct90[2]),
    "posterior_summary.csv: C_T, Upper 90%"
)
add!(
    "cumulative_infections_interval90", fmt_interval(ct90...),
    "posterior_summary.csv: C_T, Lower 90% to Upper 90%"
)

cdates, cvals = history(obs, "confirmed_case_history")
ddates, dvals = history(obs, "confirmed_death_history")
cdates[end] == cutoff ||
    error("confirmed_case_history does not end at the cut-off")
add!(
    "confirmed_cases_at_cutoff", fmt_count(cvals[end]),
    "observations.toml: confirmed_case_history, last value (dated as_of_date)"
)
add!(
    "confirmed_deaths_at_cutoff", fmt_count(dvals[end]),
    "observations.toml: confirmed_death_history, last value (dated as_of_date)"
)
add!(
    "ratio_infections_to_confirmed", fmt2(ct.median / cvals[end]),
    "cumulative_infections_median (unrounded) divided by " *
        "confirmed_cases_at_cutoff"
)
add!(
    "ratio_infections_to_confirmed_interval90",
    fmt_interval(ct90[1] / cvals[end], ct90[2] / cvals[end]),
    "posterior_summary.csv C_T 90% bounds divided by " *
        "confirmed_cases_at_cutoff"
)

rt = stream_row(streams, "joint", "R_T")
rt90 = bounds90("R_T")
add!(
    "rt_cutoff_median", fmt2(rt.median),
    "stream_estimates.csv: fit joint, quantity R_T, median"
)
add!(
    "rt_cutoff_lower90", fmt2(rt90[1]),
    "posterior_summary.csv: R_T, Lower 90%"
)
add!(
    "rt_cutoff_upper90", fmt2(rt90[2]),
    "posterior_summary.csv: R_T, Upper 90%"
)
add!(
    "rt_cutoff_interval90", fmt_interval(rt90...),
    "posterior_summary.csv: R_T, Lower 90% to Upper 90%"
)
add!(
    "rt_cutoff_prob_above_one", fmt2(100 * mean(draws.R_T .> 1)) * "%",
    "posterior_draws.csv: share of the thinned R_T draws above one"
)

r0 = stream_row(streams, "joint", "R0")
add!(
    "r0_median", fmt2(r0.median),
    "stream_estimates.csv: fit joint, quantity R0 (renewal walk start), " *
        "median"
)
add!(
    "r0_interval90", fmt_interval(r0.lo90, r0.hi90),
    "stream_estimates.csv: fit joint, quantity R0, lo90 to hi90"
)

growth90 = bounds90("r")
add!(
    "growth_rate_interval90", fmt_interval(growth90...),
    "posterior_summary.csv: r (latest daily growth rate), " *
        "Lower 90% to Upper 90%"
)
add!(
    "doubling_time_median", fmt2(median(draws.doubling_time)),
    "posterior_draws.csv: median of the thinned doubling_time draws (days; " *
        "the 90% interval spans zero growth so no bounds are quoted)"
)

cfr90 = bounds90("CFR")
add!(
    "cfr_median", fmt_pct(median(draws.CFR)),
    "posterior_draws.csv: median of the thinned CFR draws"
)
add!("cfr_lower90", fmt_pct(cfr90[1]), "posterior_summary.csv: CFR, Lower 90%")
add!("cfr_upper90", fmt_pct(cfr90[2]), "posterior_summary.csv: CFR, Upper 90%")
add!(
    "cfr_interval90", fmt_interval(cfr90...; f = fmt_pct),
    "posterior_summary.csv: CFR, Lower 90% to Upper 90%"
)

age90 = bounds90("T")
add!(
    "outbreak_age_days_interval90", fmt_interval(age90...; f = fmt_int),
    "posterior_summary.csv: T (days from the first infection to the " *
        "cut-off), Lower 90% to Upper 90%"
)
start90 = (
    cutoff - Day(round(Int, age90[2])), cutoff - Day(round(Int, age90[1])),
)
add!(
    "outbreak_start_interval90",
    fmt_interval(start90...; f = d -> Dates.format(d, "d U")),
    "as_of_date minus the posterior_summary.csv T bounds"
)

## --- Provinces --------------------------------------------------------------
## The province tables are published only on the site snapshot's province
## page (site.zip: estimates/province.html), one "Detail by province" table
## per patch under an h4 heading, and the infection shares and the
## probability that R is above one only in the summary table
## (site.zip: summary_assets/provinces.html).

const PROVINCES = [
    ("Ituri", "ituri"), ("Nord-Kivu", "nord_kivu"),
    ("Haut-Uele", "haut_uele"), ("Other provinces", "other"),
]
province_page = "site.zip estimates/province.html"
shares_page = "site.zip summary_assets/provinces.html"
detail = Dict(headed_tables(site_page("estimates/province.html"), "h4"))
shares_table = html_tables(site_page("summary_assets/provinces.html"))
for (label, key) in PROVINCES
    t = detail[label]
    ## The numeric columns `Lower 90%` .. `Upper 90%` follow the label.
    inf = table_row([t], "Cumulative infections")
    r = table_row([t], "Reproduction number")
    lo, hi = parse(Float64, inf[2]), parse(Float64, inf[7])
    add!(
        "$(key)_infections_interval90", fmt_interval(lo, hi),
        "$province_page: $label table, Cumulative infections, " *
            "Lower 90% to Upper 90%"
    )
    lo, hi = parse(Float64, r[2]), parse(Float64, r[7])
    add!(
        "$(key)_rt_lower90", fmt2(lo),
        "$province_page: $label table, Reproduction number, Lower 90%"
    )
    add!(
        "$(key)_rt_upper90", fmt2(hi),
        "$province_page: $label table, Reproduction number, Upper 90%"
    )
    add!(
        "$(key)_rt_interval90", fmt_interval(lo, hi),
        "$province_page: $label table, Reproduction number, " *
            "Lower 90% to Upper 90%"
    )
    ## Columns after the label: share of infections (%), R at the
    ## cut-off, P(R > 1), CFR (%), relative ascertainment.
    row = table_row(shares_table, label)
    m = match(r"^(\d+)–(\d+)$", row[2])
    m === nothing && error("no infection share for $label in $shares_page")
    add!(
        "$(key)_share_interval90",
        m.captures[1] * " to " * m.captures[2] * "%",
        "$shares_page: table, $label, Share of infections (%)"
    )
    add!(
        "$(key)_rt_prob_above_one", row[4],
        "$shares_page: table, $label, P(R > 1)"
    )
end
add!(
    "n_provinces", string(length(unique(province_forecast.province))),
    "province_forecast.csv: number of distinct values of province"
)

## --- Ascertainment by stream ---------------------------------------------
## Priors are read from src/models/priors.jl at the release commit: the
## pooled ascertainment centre `mu_prior = Normal(logit(p), 1.0)` puts the
## prior median of both p_drc and p_uganda at p, and the fraction tested
## has `fraction_tested_prior = Beta(a, b)`.

priors_src = read(
    `git -C $(repo_dir()) show $(meta.commit):src/models/priors.jl`, String
)
m = match(r"mu_prior = Normal\(logit\(([0-9.]+)\), ([0-9.]+)\)", priors_src)
m === nothing &&
    error("pooled ascertainment prior not found in priors.jl at $(meta.commit)")
asc_prior_median = parse(Float64, m.captures[1])
m = match(r"fraction_tested_prior = Beta\(([0-9.]+), ([0-9.]+)\)", priors_src)
m === nothing &&
    error("fraction tested prior not found in priors.jl at $(meta.commit)")
tau_prior_median = median(
    Beta(parse(Float64, m.captures[1]), parse(Float64, m.captures[2]))
)
asc_prior_source = "src/models/priors.jl at the release commit: " *
    "pooled_ascertainment_model mu_prior centre"

p90 = bounds90("p_drc")
add!(
    "suspected_ascertainment_median", fmt2(median(draws.p_drc)),
    "posterior_draws.csv: median of the thinned p_drc draws " *
        "(suspected-case ascertainment)"
)
add!(
    "suspected_ascertainment_interval90", fmt_interval(p90...),
    "posterior_summary.csv: p_drc, Lower 90% to Upper 90%"
)
add!(
    "suspected_ascertainment_prior_median", fmt2(asc_prior_median),
    asc_prior_source
)
p90 = bounds90("p_uganda")
add!(
    "export_ascertainment_median", fmt2(median(draws.p_uganda)),
    "posterior_draws.csv: median of the thinned p_uganda draws " *
        "(export detection in Uganda)"
)
add!(
    "export_ascertainment_interval90", fmt_interval(p90...),
    "posterior_summary.csv: p_uganda, Lower 90% to Upper 90%"
)
add!(
    "export_ascertainment_prior_median", fmt2(asc_prior_median),
    asc_prior_source
)
p90 = bounds90("tau_test")
add!(
    "tested_fraction_interval90", fmt_interval(p90...),
    "posterior_summary.csv: tau_test (fraction of suspects sampled), " *
        "Lower 90% to Upper 90%"
)
add!(
    "tested_fraction_prior_median", fmt2(tau_prior_median),
    "src/models/priors.jl at the release commit: fraction_tested_prior " *
        "Beta median"
)
national_tables = html_tables(site_page("estimates/national.html"))
onset = table_row(national_tables, "median modelled ascertainment")
add!(
    "onset_ascertainment_interval90",
    fmt_interval(parse(Float64, onset[2]), parse(Float64, onset[7])),
    "site.zip estimates/national.html: onset-report table, median " *
        "modelled ascertainment, Lower 90% to Upper 90%"
)

## --- Sampler and diagnostics ---------------------------------------------
## The sampler budget is read from docs/fits/registry.jl at the release
## commit (`joint_sampler_args` and the `chains` default of
## `build_fit_specs`); the diagnostics from the site snapshot's diagnostics
## table (site.zip: summary_assets/diagnostics.html, row `joint`).

registry_src = read(
    `git -C $(repo_dir()) show $(meta.commit):docs/fits/registry.jl`, String
)
function registry_value(pattern, what)
    m = match(pattern, registry_src)
    m === nothing &&
        error("$what not found in registry.jl at $(meta.commit)")
    return m.captures[1]
end
registry = "docs/fits/registry.jl at the release commit"
add!(
    "sampler_chains", registry_value(r"chains::Integer = (\d+)", "chains"),
    "$registry: build_fit_specs chains default"
)
add!(
    "sampler_draws",
    fmt_count(
        parse(Int, registry_value(r"joint_samples\((\d+)\)", "samples"))
    ),
    "$registry: joint_sampler_args samples (per chain, after warm-up)"
)
add!(
    "sampler_adaptation",
    fmt_count(
        parse(Int, registry_value(r"joint_warmup\((\d+)\)", "warm-up"))
    ),
    "$registry: joint_sampler_args n_adapts"
)
add!(
    "sampler_target_accept",
    registry_value(
        r"BVD_JOINT_TARGET_ACCEPT\", \"([0-9.]+)\"", "target accept"
    ),
    "$registry: joint_target_accept default"
)
add!(
    "sampler_max_depth",
    registry_value(r"BVD_JOINT_MAX_DEPTH\", \"(\d+)\"", "max depth"),
    "$registry: joint_max_depth default"
)

## --- Convergence gate and engineering figures ---------------------------
## The gate thresholds are read from docs/fits/convergence.jl at the release
## commit. The gradient and build-time figures are the ones recorded in the
## pull requests that made the changes, quoted here with their source; they
## are not measured from the release assets.

convergence_src = read(
    `git -C $(repo_dir()) show $(meta.commit):docs/fits/convergence.jl`,
    String
)
function convergence_value(field)
    m = match(
        Regex("const CONVERGENCE_FAIL = \\([^)]*?" * field * " = ([0-9.]+)"),
        convergence_src
    )
    m === nothing && error("CONVERGENCE_FAIL.$field not found")
    return m.captures[1]
end
gate = "docs/fits/convergence.jl at the release commit: CONVERGENCE_FAIL"
add!("convergence_fail_rhat", convergence_value("rhat"), "$gate.rhat")
add!(
    "convergence_fail_ess_bulk",
    fmt_count(round(Int, parse(Float64, convergence_value("ess_bulk")))),
    "$gate.ess_bulk"
)
add!(
    "convergence_fail_divergent_fraction",
    string(round(Int, 100 * parse(Float64, convergence_value("divergent_fraction")))) * "%",
    "$gate.divergent_fraction, as a percentage of draws"
)
add!(
    "gradient_speedup_joint_percent", "20",
    "PR #810 body and news v2.2.0: joint gradient under Mooncake about 20% faster from the hand-written rules"
)
add!(
    "gradient_speedup_kernels_min", "2.7",
    "PR #837 body: smallest per-kernel gradient speed-up (onset-reporting CDF table, 262 to 97 us)"
)
add!(
    "gradient_speedup_kernels_max", "17",
    "PR #837 body: largest per-kernel gradient speed-up (abscond thinning, 5 to 17x)"
)
add!(
    "precompile_cold_build_factor", "3.8",
    "PR #791 and news v2.2.0: headline joint fit cold build 1095 s to 292 s"
)

diag_page = "site.zip summary_assets/diagnostics.html"
diag = table_row(
    html_tables(site_page("summary_assets/diagnostics.html")), "joint"
)
add!("max_rhat", diag[2], "$diag_page: row joint, max_rhat")
add!("min_ess_bulk", diag[3], "$diag_page: row joint, min_ess_bulk")
add!("n_divergent", diag[4], "$diag_page: row joint, divergences")

## --- One-week-ahead forecasts --------------------------------------------
## forecast_validation.csv holds the one-week-back validation forecast (made
## a week before the cut-off, scored against the cut-off report);
## forecast.csv the forecast made at the cut-off, not yet observable.
## Values are new confirmed cases or deaths over the seven days.

vmade = only(unique(validation.made_date))
vtarget = vmade + Day(7)
add!(
    "forecast_validation_made", fmt_day(vmade),
    "forecast_validation.csv: made_date"
)
add!(
    "forecast_validation_target", fmt_day(vtarget),
    "forecast_validation.csv: made_date plus seven days"
)
for (stream, key, hist, dates, values) in (
        (
            "confirmed cases", "confirmed_cases", "confirmed_case_history",
            cdates, cvals,
        ),
        (
            "confirmed deaths", "confirmed_deaths", "confirmed_death_history",
            ddates, dvals,
        ),
    )
    v = forecast_draws(validation, stream; made = vmade, horizon = 7)
    add!(
        "forecast_validation_$(key)_median", fmt_int(median(v)),
        "forecast_validation.csv: $stream, horizon 7, median of draws"
    )
    add!(
        "forecast_validation_$(key)_interval90",
        fmt_interval(q90(v)...; f = fmt_int),
        "forecast_validation.csv: $stream, horizon 7, 5th to 95th " *
            "percentile of draws"
    )
    o = increment_over(dates, values, vmade, vtarget)
    o === missing || add!(
        "observed_$(key)_validation_week", fmt_count(o),
        "observations.toml: $hist, value at the validation target minus " *
            "value at the validation origin"
    )
    f = forecast_draws(forecast, stream; made = cutoff, horizon = 7)
    add!(
        "forecast_$(key)_next_week_median", fmt_int(median(f)),
        "forecast.csv: $stream, horizon 7, median of draws"
    )
    add!(
        "forecast_$(key)_next_week_interval90",
        fmt_interval(q90(f)...; f = fmt_int),
        "forecast.csv: $stream, horizon 7, 5th to 95th percentile of draws"
    )
end

## --- Released estimates across the journey --------------------------------
## paper/data/release_estimates.csv (written by release_estimates.jl) holds
## the headline median and 90% bounds of every results-<tag> release. The
## McCabe et al. headline scenarios are the two rows of the first release's
## imperial_comparison.csv that fig_first.jl draws filled.

released = CSV.read(
    joinpath(paper_dir(), "data", "release_estimates.csv"), DataFrame
)
released_row(tag) = released[only(findall(==(tag), released.tag)), :]
first_rel = released_row("v1.0.0")
add!(
    "first_release_median", fmt_int(first_rel.median),
    "release_estimates.csv: row v1.0.0, median"
)
add!(
    "first_release_interval90",
    fmt_interval(first_rel.lower90, first_rel.upper90; f = fmt_int),
    "release_estimates.csv: row v1.0.0, lower90 to upper90"
)
last_rel = released[findlast(!ismissing, released.median), :]
add!(
    "last_release_tag", last_rel.tag,
    "release_estimates.csv: last row with a median"
)
add!(
    "last_release_median", fmt_int(last_rel.median),
    "release_estimates.csv: row $(last_rel.tag) (the last with a median), " *
        "median"
)
## The largest step between consecutive released medians.
with_median = released[.!ismissing.(released.median), :]
step_ratios = with_median.median[2:end] ./ with_median.median[1:(end - 1)]
istep = argmax(step_ratios)
step_from, step_to = with_median.tag[istep], with_median.tag[istep + 1]
all(==("closed-form"), with_median.model_version[istep:(istep + 1)]) ||
    error("the largest release step is no longer within the closed-form model")
add!(
    "closed_form_step_from", step_from,
    "release_estimates.csv: tag before the largest ratio of consecutive medians"
)
add!(
    "closed_form_step_to", step_to,
    "release_estimates.csv: tag after the largest ratio of consecutive medians"
)
add!(
    "closed_form_step_ratio", fmt2(step_ratios[istep]),
    "release_estimates.csv: median of $step_to divided by median of $step_from"
)

const FIRST_TAG = "results-v1.0.0"
ensure_tag_assets(FIRST_TAG, ["imperial_comparison.csv"])
comparison = CSV.read(
    tag_asset_path(FIRST_TAG, "imperial_comparison.csv"), DataFrame
)
function comparison_row(source)
    i = findfirst(==(source), comparison.Source)
    i === nothing && error("no `$source` row in imperial_comparison.csv")
    return comparison[i, :]
end
comparison_src = "$FIRST_TAG imperial_comparison.csv"
add!(
    "mccabe_headline_exports",
    fmt_int(
        comparison_row("McCabe Method 1 (Ituri, w=15 d)")["Central estimate"]
    ),
    "$comparison_src: McCabe Method 1 (Ituri, w=15 d), Central estimate"
)
add!(
    "mccabe_headline_deaths",
    fmt_int(
        comparison_row("McCabe Method 2 (τ=14 d, CFR 30%)")["Central estimate"]
    ),
    "$comparison_src: McCabe Method 2 (τ=14 d, CFR 30%), Central estimate"
)

## --- One-week-ahead skill across releases --------------------------------
## data/forecast_scores.csv in the repository holds the cross-release
## scores; forecast_skill (release_assets.jl) pairs each joint score at
## horizon seven with its persistence baseline, as fig_evaluation.jl does.

scores = CSV.read(
    joinpath(repo_dir(), "data", "forecast_scores.csv"), DataFrame
)
skill = forecast_skill(scores)
beat(rs) = count(<(1), rs.skill)
skill_src = "data/forecast_scores.csv: horizon 7, joint CRPS below the " *
    "baseline's"
add!(
    "skill_pairs_beat", string(beat(skill)),
    "$skill_src, all release-stream pairs"
)
add!(
    "skill_pairs_total", string(nrow(skill)),
    "data/forecast_scores.csv: release-stream pairs with a joint and a " *
        "baseline score at horizon 7"
)
live = skill[.!skill.backfill, :]
add!(
    "skill_pairs_beat_live", string(beat(live)),
    "$skill_src, releases not marked (backfill)"
)
add!(
    "skill_pairs_total_live", string(nrow(live)),
    "data/forecast_scores.csv: scored pairs at releases not marked (backfill)"
)
for (stream, key) in (
        ("confirmed cases", "confirmed_cases"),
        ("confirmed deaths", "confirmed_deaths"),
        ("isolation beds", "isolation_beds"),
        ("onset reports", "onset_reports"), ("recovered", "recovered"),
    )
    rs = skill[skill.stream .== stream, :]
    add!(
        "skill_$key", string(beat(rs), " of ", nrow(rs)),
        "$skill_src, stream $stream, over its scored releases"
    )
end
add!(
    "coverage90_covered", string(count(skill.covered)),
    "data/forecast_scores.csv: horizon 7 joint coverage_90 true"
)
add!(
    "coverage90_total", string(nrow(skill)),
    "data/forecast_scores.csv: scored pairs, as skill_pairs_total"
)

## --- Single-stream fits at the cut-off -----------------------------------
## The single-stream fit ids of stream_estimates.csv, as the paper names
## them (the list fig_evaluation.jl draws in its panel B).

const FIT_LABELS = [
    ("confirmed", "confirmed cases"), ("confirmed_deaths", "confirmed deaths"),
    ("cases", "reported cases"), ("deaths", "reported deaths"),
    ("treatment", "isolation beds"), ("onsets", "onset reports"),
    ("exports", "exports"),
]
singles = [stream_row(streams, fit, "C_T").median for (fit, _) in FIT_LABELS]
imin, imax = argmin(singles), argmax(singles)
single_src = "stream_estimates.csv: quantity C_T, median"
add!(
    "single_stream_min_fit", FIT_LABELS[imin][2],
    "$single_src, single-stream fit with the smallest median"
)
add!(
    "single_stream_min_median", fmt_int(singles[imin]),
    "$single_src, fit $(FIT_LABELS[imin][1]) (the smallest)"
)
add!(
    "single_stream_max_fit", FIT_LABELS[imax][2],
    "$single_src, single-stream fit with the largest median"
)
add!(
    "single_stream_max_median", fmt_int(singles[imax]),
    "$single_src, fit $(FIT_LABELS[imax][1]) (the largest)"
)
add!(
    "single_stream_ratio",
    string(round(singles[imax] / singles[imin]; digits = 2)),
    "$single_src, largest single-stream median over the smallest"
)
add!(
    "joint_median_unrounded", fmt_int(ct.median),
    "stream_estimates.csv: fit joint, quantity C_T, median to the nearest " *
        "infection"
)

## --- Chamla et al. comparison -----------------------------------------------
## Their week-12 projections are constants in src/constants.jl at the
## release commit; the joint fit frozen at their calibration date is
## published only as the comparison table on the site snapshot's
## sensitivity page (site.zip: sensitivity.html).

constants_src = read(
    `git -C $(repo_dir()) show $(meta.commit):src/constants.jl`, String
)
function source_const(name)
    m = match(Regex("const $name = (\\[.*?\\n\\])", "s"), constants_src)
    m === nothing &&
        error("no `const $name` in src/constants.jl at $(meta.commit)")
    return eval(Meta.parse(m.captures[1]))
end
chamla_w12 = source_const("CHAMLA_CONFIRMED_W12")
chamla_central = source_const("CHAMLA_CONFIRMED_CENTRAL")
w12 = only(filter(c -> occursin("central", c[1]), chamla_w12))
w12_date = Date(
    only(filter(c -> c[2:4] == w12[2:4], chamla_central))[1]
)
constants = "src/constants.jl at the release commit"
add!(
    "chamla_week12_date", fmt_day(w12_date),
    "$constants: CHAMLA_CONFIRMED_CENTRAL row matching the " *
        "CHAMLA_CONFIRMED_W12 central scenario"
)
add!(
    "chamla_central_median", fmt_int(w12[2]),
    "$constants: CHAMLA_CONFIRMED_W12 central scenario, median"
)
add!(
    "chamla_central_interval90", fmt_interval(w12[3], w12[4]; f = fmt_int),
    "$constants: CHAMLA_CONFIRMED_W12 central scenario, lower_90 to upper_90"
)
sensitivity_page = "site.zip sensitivity.html"
sensitivity_tables = html_tables(site_page("sensitivity.html"))
chamla_table = findfirst(
    t -> any(r -> "Chamla central (90% PI)" in r, t), sensitivity_tables
)
chamla_table === nothing &&
    error("no Chamla comparison table in $sensitivity_page")
chamla_header = sensitivity_tables[chamla_table][1]
iours = findfirst(==("Our projection (90% CrI)"), chamla_header)
frozen_row = table_row(
    [sensitivity_tables[chamla_table]], Dates.format(w12_date, "d U")
)
m = match(r"(\d+) \((\d+)–(\d+)\)", frozen_row[iours])
m === nothing && error("no projection at $w12_date in $sensitivity_page")
add!(
    "chamla_frozen_projection_median", fmt_count(parse(Int, m.captures[1])),
    "$sensitivity_page: Chamla comparison table, row " *
        "$(Dates.format(w12_date, "d U")), Our projection median"
)
add!(
    "chamla_frozen_projection_interval90",
    fmt_interval(
        parse(Int, m.captures[2]), parse(Int, m.captures[3]); f = fmt_count
    ),
    "$sensitivity_page: Chamla comparison table, row " *
        "$(Dates.format(w12_date, "d U")), Our projection 90% CrI"
)
iw12 = findfirst(==(w12_date), cdates)
iw12 === nothing && error("confirmed_case_history has no $w12_date point")
add!(
    "observed_confirmed_week12", fmt_count(cvals[iw12]),
    "observations.toml: confirmed_case_history, value at chamla_week12_date"
)

## --- Write ---------------------------------------------------------------

out = joinpath(paper_dir(), "generated", "numbers.yml")
mkpath(dirname(out))
open(out, "w") do io
    println(
        io,
        "# Quoted numbers for paper/main.qmd, read with Quarto meta shortcodes."
    )
    println(
        io, "# Written by paper/scripts/paper_numbers.jl from release ",
        meta.tag, " (commit ", meta.commit[1:7], ", published ", meta.date,
        ");"
    )
    println(io, "# never edited by hand. Each key's comment names its source.")
    for (key, value, source) in entries
        println(io, "# ", source)
        println(io, key, ": \"", replace(value, "\"" => "\\\""), "\"")
    end
end
println("wrote ", length(entries), " keys to ", out)
