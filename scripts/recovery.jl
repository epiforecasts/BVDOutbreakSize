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
#   BVD_RECOVERY_SAMPLES (1000), BVD_RECOVERY_WARMUP (500),
#   BVD_RECOVERY_CHAINS (2), BVD_RECOVERY_MAX_DEPTH (10): the headline
#   joint's own settings, since a shorter run does not converge on this
#   model and would test the sampler rather than the model,
#   BVD_RECOVERY_HORIZON (14)
#   BVD_RECOVERY_DRY_RUN (false): simulate and check the density, but do not
#   fit

using BVDOutbreakSize, CSV, DataFrames, Turing
using Random: MersenneTwister
using BVDOutbreakSize: _draws, _draw_vectors

seed = parse(Int, get(ARGS, 1, "1"))
out_dir = get(ARGS, 2, joinpath("output", "recovery"))
env_int(k, d) = parse(Int, get(ENV, k, string(d)))
samples = env_int("BVD_RECOVERY_SAMPLES", 1000)
warmup = env_int("BVD_RECOVERY_WARMUP", 500)
chains = env_int("BVD_RECOVERY_CHAINS", 2)
max_depth = env_int("BVD_RECOVERY_MAX_DEPTH", 10)
horizon = env_int("BVD_RECOVERY_HORIZON", 14)
mkpath(out_dir)

## National and province quantities checked, as the report reads them.
const SCALARS = ["C_T", "T", "R_T", "r", "CFR", "p_drc", "tau_test", "lambda_bg"]
const PER_PROVINCE = ["C_T_patch", "R_T_patch", "CFR_patch", "province_ascertainment"]

obs = load_observations()
breakpoint = default_breakpoint(obs)
generator = generator_joint(obs; breakpoint)
observed = recovery_observed_varnames(
    generator, production_joint(obs; breakpoint)
)

## Forecast columns scored, with the in-window streams each one's persistence
## baseline is read from: the data variable names and their observation days.
## The confirmed cases are split into windows at the first and last
## laboratory dates, each drawn as its own stream.
windows = BVDOutbreakSize.confirmed_positivity_windows(
    obs.confirmed_history, obs.lab_history, obs.lab_daily_history,
    obs.confirmed_break_days
)
FORECASTS = [
    (;
        column = "confirmed_new",
        parts = [
            "confirmed_state.early_increments.increments" =>
                windows.early_days,
            "confirmed_state.confirmed_positives.positives" =>
                windows.obs_days,
            "confirmed_state.late_increments.increments" => windows.late_days,
        ],
    ),
    (;
        column = "confirmed_deaths_new",
        parts = [
            "confirmed_deaths_state.cdeath_increments.increments" =>
                obs.confirmed_deaths_history.days,
        ],
    ),
]

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
flush(stdout)
get(ENV, "BVD_RECOVERY_DRY_RUN", "false") == "true" && exit(0)
fit = recovery_fit(
    model; samples, n_adapts = warmup, chains, max_depth, seed = seed + 1
)
fit_minutes = round((time() - t0) / 60; digits = 1)

## Parameter recovery. The draws of every checked quantity by name.
function tracked(chain)
    out = Dict{String, Vector{Float64}}()
    for k in SCALARS
        out[k] = vec(_draws(chain, Symbol(k)))
    end
    for k in PER_PROVINCE
        dv = _draw_vectors(chain, Symbol(k))
        for p in eachindex(PROVINCE_LABELS)
            out["$(k)[$(PROVINCE_LABELS[p])]"] = [v[p] for v in dv]
        end
    end
    return out
end
truth = Dict(k => only(v) for (k, v) in tracked(sim.truth))
draws = tracked(fit.chain)
params = recovery_table(truth, draws)
diagnostics = fit_diagnostics(fit.chain)
params.seed .= seed
params.fit_minutes .= fit_minutes
params.max_rhat .= diagnostics.max_rhat
params.min_ess_bulk .= diagnostics.min_ess_bulk
params.divergences .= diagnostics.n_divergent
CSV.write(joinpath(out_dir, "recovery_$(seed).csv"), params)
## Every fifth posterior draw and a prior sample of the fitted model, one
## column per quantity, for the report's figures. The prior is the same for
## every seed, so each seed's sample adds to one pooled prior.
thinned(d) = DataFrame(
    [k => v[1:5:end] for (k, v) in sort(collect(d); by = first)]
)
CSV.write(joinpath(out_dir, "recovery_draws_$(seed).csv"), thinned(draws))
prior = sample(MersenneTwister(seed), model, Prior(), 1000; progress = false)
CSV.write(
    joinpath(out_dir, "recovery_prior_$(seed).csv"), thinned(tracked(prior))
)

## Forecasts from the fit, scored against the simulated future.
## The data values of a stream in index order, drawn element by element or
## as one vector.
function stream_values(data, prefix)
    for (vn, v) in data
        string(vn) == prefix && return collect(v)
    end
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
    total = 0
    for (prefix, days) in spec.parts
        values = stream_values(sim.data, prefix)
        isempty(days) && continue
        length(values) == length(days) || return nothing
        total += sum(
            v for (v, d) in zip(values, days) if obs.n - h < d <= obs.n;
            init = 0
        )
    end
    return total
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

verdict = recovery_verdict(params; diagnostics)
println(
    "seed $seed: recovery ", verdict.status,
    " (90% coverage ", round(verdict.coverage_90; digits = 2),
    ", outside the 99% interval: ",
    isempty(verdict.outside) ? "none" : join(verdict.outside, ", "),
    "), ", nrow(forecasts), " forecasts scored, fit took $fit_minutes min"
)
flush(stdout)
