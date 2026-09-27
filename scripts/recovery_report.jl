# Summarise the parameter-recovery runs (scripts/recovery.jl) across seeds as
# Markdown: one verdict over every seed, a table per level with each
# quantity's error relative to the truth across seeds, and each seed's fit
# in a dropdown. Writes the body to `<dir>/report.md` and the overall status
# to `<dir>/status`, which scripts/recovery_report.sh posts.
#
# Usage:
#   julia --project=docs scripts/recovery_report.jl <dir with recovery_*.csv>

using BVDOutbreakSize, CSV, DataFrames
using Printf: @sprintf
using Statistics: median

dir = get(ARGS, 1, joinpath("output", "recovery"))
read_all(pattern) = begin
    files = sort(filter(f -> occursin(pattern, f), readdir(dir)))
    isempty(files) ? DataFrame() :
        reduce(vcat, [CSV.read(joinpath(dir, f), DataFrame) for f in files])
end
params = read_all(r"^recovery_\d+\.csv$")
forecasts = read_all(r"^forecast_recovery_\d+\.csv$")
isempty(params) && error("no recovery results in $dir")

seeds = recovery_seed_verdicts(params)
status = recovery_overall(seeds.status)
summary = recovery_summary(params)

function fmt(x)
    ismissing(x) && return "–"
    a = abs(x)
    a >= 1000 && return BVDOutbreakSize._recovery_tick(x)
    a >= 10 && return string(round(x; digits = 1))
    a >= 0.1 && return @sprintf("%.2f", x)
    return @sprintf("%.3g", x)
end
f2(x) = ismissing(x) ? "–" : @sprintf("%.2f", x)
spread(m, lo, hi) = "$(f2(m)) ($(f2(lo)) to $(f2(hi)))"
function table(rows)
    lines = [
        "| Quantity | Truth | Relative error | z | Truth in 90% | Largest z |",
        "|---|---|---|---|---|---|",
    ]
    for r in eachrow(rows)
        truth = r.truth_min == r.truth_max ? fmt(r.truth_min) :
            "$(fmt(r.truth_min)) to $(fmt(r.truth_max))"
        flag = r.outside > 0 ? " **outside 99% in $(r.outside)**" : ""
        push!(
            lines,
            "| `$(r.quantity)` | $truth | " *
                spread(r.rel_error_median, r.rel_error_min, r.rel_error_max) *
                " | $(spread(r.z_median, r.z_min, r.z_max)) | " *
                "$(r.covered_90)/$(r.n_seeds)$flag | " *
                "seed $(r.worst_seed) ($(f2(r.worst_z))) |"
        )
    end
    return join(lines, "\n")
end
national = summary[.!occursin.("[", summary.quantity), :]
province = summary[occursin.("[", summary.quantity), :]

n_pairs = nrow(params)
n_in = count(params.covered_90)
n_out = sum(summary.outside)
n_conv = count(!=(:unconverged), seeds.status)
first_seed = params[params.seed .== first(params.seed), :]
n_allowed = recovery_verdict(first_seed).outside_allowed
skill = isempty(forecasts) ? "No forecast was scored." :
    "Forecasts: $(nrow(forecasts)) scored, median CRPS relative to " *
    "persistence $(fmt(median(forecasts.relative_crps))) (below one beats it)."

seed_rows = [
    "| $(r.seed) | $(r.status) | $(fmt(r.coverage_90)) | " *
        "$(isempty(r.outside) ? "none" : r.outside) | $(r.fit_minutes) | " *
        "$(@sprintf("%.3f", r.max_rhat)) | $(round(Int, r.min_ess_bulk)) | " *
        "$(r.divergences) |"
        for r in eachrow(seeds)
]

body = """
### Parameter recovery: $status ($(nrow(seeds)) seeds)

$n_conv of $(nrow(seeds)) fits converged (R-hat at most 1.1, bulk ESS at least 30).
An unconverged seed's recovery is shown but not judged.
The truth is inside the 90% interval in $n_in of $n_pairs quantity-seed pairs ($(round(Int, 100 * n_in / n_pairs))%) and outside the 99% interval in $(n_out == 0 ? "none" : n_out).
$skill

**National**

$(table(national))

<details><summary>Province</summary>

$(table(province))

</details>

<details><summary>Fit per seed</summary>

| Seed | Verdict | 90% coverage | Outside 99% | Fit (min) | Max R-hat | Min bulk ESS | Divergences |
|---|---|---|---|---|---|---|---|
$(join(seed_rows, "\n"))

</details>

A seed fails when more quantities lie outside their 99% interval than the 99th percentile of a Binomial with a 1% chance each ($n_allowed of $(nrow(first_seed))), or fewer than 60% lie inside their 90% interval.
Relative error is (posterior median − truth) / |truth| and z is (posterior mean − truth) / posterior SD, each the median across seeds with the range in brackets.
"""
write(joinpath(dir, "report.md"), body)
write(joinpath(dir, "status"), string(status))
println(status)
