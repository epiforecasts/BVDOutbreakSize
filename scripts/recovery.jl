# Parameter recovery and simulated forecasts for the headline joint.
#
# One run simulates a dataset from the model itself at a fixed seed, fits the
# model to it with short sampler settings, and writes how well the fit
# recovers the values that generated the data, and how its forecasts score
# against the simulated future and a persistence baseline. See
# `src/recovery.jl` for the method.
#
# Usage:
#   julia --project=docs scripts/recovery.jl <seed> [out_dir]
#
# Settings, from the environment:
#   BVD_RECOVERY_SAMPLES (200), BVD_RECOVERY_WARMUP (200),
#   BVD_RECOVERY_CHAINS (2), BVD_RECOVERY_MAX_DEPTH (8),
#   BVD_RECOVERY_HORIZON (14)
#   BVD_RECOVERY_DRY_RUN (false): simulate and check the density, but do not
#   fit

using BVDOutbreakSize, CSV, DataFrames
using Statistics: median
using BVDOutbreakSize: _draws, _draw_vectors

seed = parse(Int, get(ARGS, 1, "1"))
out_dir = get(ARGS, 2, joinpath("output", "recovery"))
env_int(k, d) = parse(Int, get(ENV, k, string(d)))
samples = env_int("BVD_RECOVERY_SAMPLES", 200)
warmup = env_int("BVD_RECOVERY_WARMUP", 200)
chains = env_int("BVD_RECOVERY_CHAINS", 2)
max_depth = env_int("BVD_RECOVERY_MAX_DEPTH", 8)
horizon = env_int("BVD_RECOVERY_HORIZON", 14)
mkpath(out_dir)

## National and province quantities checked, as the report reads them.
const SCALARS = ["C_T", "T", "R_T", "r", "CFR", "p_drc", "tau_test", "lambda_bg"]
const PER_PROVINCE = ["C_T_patch", "R_T_patch", "CFR_patch", "province_ascertainment"]
## Forecast columns scored, with the in-window stream each one's persistence
## baseline is read from: the data variable names and the observation days.
const FORECASTS = [
    (;
        column = "confirmed_new",
        prefixes = [
            "confirmed_state.early_increments.increments",
            "confirmed_state.late_increments.increments",
        ],
        history = :confirmed_history,
    ),
    (;
        column = "confirmed_deaths_new",
        prefixes = ["confirmed_deaths_state.cdeath_increments.increments"],
        history = :confirmed_deaths_history,
    ),
]

obs = load_observations()
breakpoint = default_breakpoint(obs)
generator = generator_joint(obs; breakpoint)
observed = recovery_observed_varnames(
    generator, production_joint(obs; breakpoint)
)

## Keep a draw whose outbreak is within a factor of five of the one observed,
## so the check runs on an epidemic like this one rather than on the prior's
## tails.
function plausible(truth)
    c = only(_draws(truth, :C_T))
    return isfinite(c) && obs.confirmed_cases / 5 <= c <= 5 * obs.confirmed_cases
end

t0 = time()
sim = simulate_recovery(
    generator, observed; seed, horizon, accept = plausible
)
model = generator_joint(obs; breakpoint, simulated = recovery_data(sim))
density = recovery_density_check(generator, model, sim.truth)
println(
    "seed $seed: simulated in ", round(time() - t0; digits = 1),
    " s; log joint ", density.from_generator, " (generator) and ",
    density.from_data, " (rebuilt with the simulated data)"
)
get(ENV, "BVD_RECOVERY_DRY_RUN", "false") == "true" && exit(0)
fit = recovery_fit(
    model; samples, n_adapts = warmup, chains, max_depth, seed = seed + 1
)
fit_minutes = round((time() - t0) / 60; digits = 1)

## Parameter recovery.
truth = Dict{String, Float64}()
draws = Dict{String, Vector{Float64}}()
for k in SCALARS
    s = Symbol(k)
    truth[k] = only(_draws(sim.truth, s))
    draws[k] = vec(_draws(fit.chain, s))
end
for k in PER_PROVINCE
    s = Symbol(k)
    tv = only(_draw_vectors(sim.truth, s))
    dv = _draw_vectors(fit.chain, s)
    for p in eachindex(tv)
        name = "$(k)[$(PROVINCE_LABELS[p])]"
        truth[name] = tv[p]
        draws[name] = [v[p] for v in dv]
    end
end
params = recovery_table(truth, draws)
params.seed .= seed
params.fit_minutes .= fit_minutes
CSV.write(joinpath(out_dir, "recovery_$(seed).csv"), params)

## Forecasts from the fit, scored against the simulated future.
## The data values of a stream in index order.
function stream_values(data, prefix)
    hits = [
        (parse(Int, m[1]), v) for (vn, v) in data
            for m in (
                match(
                    Regex(
                        "^" * replace(prefix, "." => "\\.") *
                        "\\[(\\d+)\\]\$"
                    ), string(vn)
                ),
            ) if m !== nothing
    ]
    return [v for (_, v) in sort(hits; by = first)]
end
## The stream's last `h` days of simulated increments carried forward.
function persistence(spec, h)
    values = reduce(vcat, [stream_values(sim.data, p) for p in spec.prefixes])
    days = getproperty(obs, spec.history).days
    length(values) == length(days) || return nothing
    return sum(
        v for (v, d) in zip(values, days) if obs.n - h < d <= obs.n;
        init = 0
    )
end
pp = forecast_draws(fit.model, fit.chain; horizon, seed = seed + 2)
zero_kw = (;
    obs_cases = 0, obs_deaths = 0, obs_confirmed = 0,
    obs_confirmed_deaths = 0,
)
forecast_rows = DataFrame[]
for h in unique([7, horizon])
    h <= horizon || continue
    model_fc = forecast_reported(pp; horizon = h, zero_kw...)
    truth_fc = forecast_reported(sim.truth; horizon = h, zero_kw...)
    ftruth = Dict{String, Float64}()
    fdraws = Dict{String, Vector{Float64}}()
    fbase = Dict{String, Float64}()
    for spec in FORECASTS
        col = Symbol(spec.column)
        (col in propertynames(model_fc) && col in propertynames(truth_fc)) ||
            continue
        base = persistence(spec, h)
        base === nothing && continue
        ftruth[spec.column] = only(truth_fc[!, col])
        fdraws[spec.column] = Float64.(model_fc[!, col])
        fbase[spec.column] = base
    end
    tab = forecast_recovery_table(ftruth, fdraws, fbase)
    tab.horizon .= h
    push!(forecast_rows, tab)
end
forecasts = reduce(vcat, forecast_rows; init = DataFrame())
forecasts.seed .= seed
CSV.write(joinpath(out_dir, "forecast_recovery_$(seed).csv"), forecasts)

verdict = recovery_verdict(params)
## One line per seed for the report job: the verdict, then the summary.
skill = isempty(forecasts) ? "no forecast scored" :
    "median relative CRPS " *
    string(round(median(forecasts.relative_crps); digits = 2))
write(
    joinpath(out_dir, "verdict_$(seed).txt"),
    (verdict.pass ? "pass" : "fail") * "\t" *
        "seed $seed: 90% coverage $(round(verdict.coverage_90; digits = 2)), " *
        "outside the 99% interval: " *
        (isempty(verdict.outside) ? "none" : join(verdict.outside, ", ")) *
        "; forecasts $skill; fit $fit_minutes min\n"
)
println(
    "seed $seed: recovery ", verdict.pass ? "passes" : "FAILS",
    " (90% coverage ", round(verdict.coverage_90; digits = 2),
    ", outside the 99% interval: ",
    isempty(verdict.outside) ? "none" : join(verdict.outside, ", "),
    "), fit took $fit_minutes min"
)
