# Parameter recovery and simulated forecasts for the headline joint. A
# dataset is drawn from the model itself: fixed-seed parameters from the
# prior, then one draw of the model run past its cut-off with those
# parameters fixed, which gives the in-window observations and the future
# they lead to from one latent path. The model conditioned on the in-window
# observations is then fitted, and its estimates and forecasts are scored
# against the values that generated the data.

"""
$(TYPEDSIGNATURES)

The headline patch [`bvd_joint`](@ref) for a `load_observations()` result
with every stream's counts dropped and its observation grid kept, so the
model draws the counts itself (the predictive-generator path). Its sampled
parameters match [`production_joint`](@ref)'s one for one, so a chain from
either replays through the other.

The confirmed-case, confirmed-death and laboratory histories keep their
counts: they define the confirmed windows, the positivity random effect and
the published break discrepancies, while the `missing` cut-off totals gate
their generator paths. The onset triangle and the province compositions
keep their cell grids with `missing` increments.
"""
function generator_joint(obs; breakpoint)
    days_only(h) = (; days = h.days, counts = Int[])
    prov_cases = province_increment_matrix(
        obs.province_confirmed_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES)
    )
    prov_deaths = province_increment_matrix(
        obs.province_death_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES)
    )
    return bvd_joint(
        obs.n, missing, missing, missing, missing, missing, missing;
        confirmed_deaths = missing,
        recovered_cases = missing,
        deaths_history = days_only(obs.deaths_history),
        reported_history = days_only(obs.reported_history),
        suspected_daily_history = days_only(obs.suspected_daily_history),
        suspected_daily_deaths_history =
            days_only(obs.suspected_daily_deaths_history),
        isolation_history = days_only(obs.isolation_history),
        bed_capacity_history = days_only(obs.bed_capacity_history),
        occupancy_break_days = obs.occupancy_break_days,
        recovered_history = days_only(obs.recovered_history),
        treatment_admissions_history =
            days_only(obs.treatment_admissions_history),
        treatment_deaths_history = days_only(obs.treatment_deaths_history),
        treatment_ruleout_history = days_only(obs.treatment_ruleout_history),
        treatment_absconded_history =
            days_only(obs.treatment_absconded_history),
        treatment_confirmed_incare_history =
            days_only(obs.treatment_confirmed_incare_history),
        treatment_suspect_incare_history =
            days_only(obs.treatment_suspect_incare_history),
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        confirmed_break_days = obs.confirmed_break_days,
        confirmed_break_gross_cases = obs.confirmed_break_gross_cases,
        confirmed_break_gross_deaths = obs.confirmed_break_gross_deaths,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        onset_curve_history = (;
            onset_days = obs.onset_curve_history.onset_days,
            report_days = obs.onset_curve_history.report_days,
            prev_report_days = obs.onset_curve_history.prev_report_days,
            increments = missing,
        ),
        breakpoint = breakpoint,
        background_pooling = background_pooling_model,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days,
        n_patches = length(PROVINCE_NAMES),
        province_increments = missing,
        province_days = prov_cases.days,
        province_testing_covariate =
            province_testing_covariate(obs.province_lab_daily_history),
        province_death_increments = missing,
        province_death_days = prov_deaths.days
    )
end

"""
$(TYPEDSIGNATURES)

How well posterior draws recover known values. `truth` maps each quantity's
name to the value that generated the data, and `draws` maps the same names
to that quantity's posterior draws. Returns one row per quantity with its
true value, the posterior median and 50% and 90% equal-tailed intervals,
whether each interval covers the truth, the truth's quantile in the
posterior (`truth_quantile`, 0.5 when it sits at the median) and its
z-score against the posterior mean and standard deviation.
"""
function recovery_table(
        truth::AbstractDict{<:AbstractString, <:Real},
        draws::AbstractDict{<:AbstractString, <:AbstractVector{<:Real}}
    )
    names = sort([k for k in keys(truth) if haskey(draws, k)])
    rows = map(names) do k
        x = collect(Float64, draws[k])
        t = Float64(truth[k])
        q = quantile(x, [0.05, 0.25, 0.5, 0.75, 0.95])
        s = std(x)
        (;
            quantity = k, truth = t, median = q[3],
            lower_90 = q[1], upper_90 = q[5],
            lower_50 = q[2], upper_50 = q[4],
            covered_50 = q[2] <= t <= q[4],
            covered_90 = q[1] <= t <= q[5],
            truth_quantile = mean(x .<= t),
            z = s > 0 ? (mean(x) - t) / s : 0.0,
        )
    end
    return DataFrame(rows)
end

"""
$(TYPEDSIGNATURES)

Whether a [`recovery_table`](@ref) shows the model failing badly to recover
the values that generated its data. It fails when any quantity's truth lies
outside the posterior's central `outer` interval (the 99% interval by
default), or when fewer than `min_coverage_90` of the quantities have the
truth inside their 90% interval. Returns `(; pass, outside, coverage_90)`:
the verdict, the quantities outside the outer interval, and the share of
quantities covered at 90%.
"""
function recovery_verdict(
        tab::DataFrame; outer::Real = 0.99,
        min_coverage_90::Real = 0.5
    )
    tail = (1 - outer) / 2
    outside = String[
        r.quantity for r in eachrow(tab)
            if r.truth_quantile < tail || r.truth_quantile > 1 - tail
    ]
    coverage_90 = isempty(tab.covered_90) ? 1.0 : mean(tab.covered_90)
    pass = isempty(outside) && coverage_90 >= min_coverage_90
    return (; pass, outside, coverage_90)
end

"""
$(TYPEDSIGNATURES)

Forecasts from a fit to simulated data, scored against the simulated future
and against a persistence baseline. `truth` maps each forecast quantity to
its simulated value, `draws` maps it to the fitted forecast's draws and
`baseline` maps it to the persistence forecast, the last observed stretch
of that stream carried forward. Returns one row per quantity with the
forecast's CRPS ([`score_draws`](@ref)), the baseline's absolute error (a
point forecast's CRPS), their ratio (`relative_crps`, below one when the
model beats the baseline) and whether the forecast's 90% interval covers
the truth.
"""
function forecast_recovery_table(
        truth::AbstractDict{<:AbstractString, <:Real},
        draws::AbstractDict{<:AbstractString, <:AbstractVector{<:Real}},
        baseline::AbstractDict{<:AbstractString, <:Real}
    )
    names = sort(
        [k for k in keys(truth) if haskey(draws, k) && haskey(baseline, k)]
    )
    rows = map(names) do k
        t = Float64(truth[k])
        s = score_draws(t, collect(Float64, draws[k]))
        base = abs(Float64(baseline[k]) - t)
        (;
            quantity = k, truth = t, baseline = Float64(baseline[k]),
            crps = s.crps, baseline_crps = base,
            relative_crps = base > 0 ? s.crps / base :
                (s.crps == 0 ? 1.0 : Inf),
            covered_90 = s.coverage_90,
        )
    end
    return DataFrame(rows)
end

"""
$(TYPEDSIGNATURES)

The variables `generator` samples that `fitted` observes: the names in a
draw of the generator that a draw of the fitted model does not carry. These
are the simulated data a recovery fit is conditioned on.
"""
function recovery_observed_varnames(generator::Model, fitted::Model)
    rng = MersenneTwister(1)
    gen = keys(last(init!!(rng, generator, VarInfo(), InitFromPrior())))
    fit = Set(keys(last(init!!(rng, fitted, VarInfo(), InitFromPrior()))))
    return [vn for vn in gen if !(vn in fit)]
end

"""
$(TYPEDSIGNATURES)

One simulated dataset from `generator` (see [`generator_joint`](@ref)): a
single prior draw of the model run `horizon` days past its cut-off, with a
fixed `seed`, redrawn up to `attempts` times until `accept(truth)` holds.
One draw carries the parameters, every deterministic, the in-window
observations and the future counts, all from one latent path.

Returns `(; truth, data)`: the draw as a one-draw chain, which the chain
readers and [`forecast_reported`](@ref) read like a fitted chain, and its
in-window observations (the names in `observed`) keyed by `VarName`, ready
for `condition`.
"""
function simulate_recovery(
        generator::Model, observed::AbstractVector;
        seed::Integer, horizon::Integer = 14, accept = _ -> true,
        attempts::Integer = 50
    )
    rng = MersenneTwister(seed)
    extended = with_horizon(generator, horizon)
    for _ in 1:attempts
        truth = sample(rng, extended, Prior(), 1; progress = false)
        accept(truth) || continue
        data = Dict(vn => only(truth[vn]) for vn in observed)
        return (; truth, data)
    end
    throw(
        ErrorException(
            "no prior draw passed the acceptance check in $attempts attempts."
        )
    )
end

"""
$(TYPEDSIGNATURES)

Fit `generator` conditioned on a [`simulate_recovery`](@ref) dataset's
in-window observations. Returns `(; model, chain)`: the conditioned model,
which [`forecast_draws`](@ref) runs past the cut-off, and its NUTS chain.
The sampler settings default to a short run (two chains of 200 draws after
200 warmup steps, tree depth at most 8), enough to check recovery without
the cost of the headline fit. Other keywords pass to
[`nuts_sample`](@ref).
"""
function recovery_fit(
        generator::Model, sim; samples::Integer = 200,
        n_adapts::Integer = 200, chains::Integer = 2,
        max_depth::Integer = 8, seed::Integer = 20260518, kwargs...
    )
    model = condition(generator, sim.data)
    chain = nuts_sample(
        model; samples, n_adapts, chains, max_depth, seed, kwargs...
    )
    return (; model, chain)
end
