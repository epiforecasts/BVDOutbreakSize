# Parameter recovery and simulated forecasts for the headline joint. A
# dataset is one fixed-seed prior draw of the model run past its cut-off,
# which gives the parameters, the in-window observations and the future
# they lead to from one latent path. The same model built with those
# in-window observations as its data is then fitted, and its estimates and
# forecasts are scored against the values that generated the data.

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

With `simulated` (a [`recovery_data`](@ref) result) the same model is built
with that simulated dataset in place of each `missing` observation: each
stream's counts through `simulated_data`, and the onset triangle and the
province compositions as their own arguments. Its structure is unchanged,
so its density at the simulated data is the generator's.
"""
function generator_joint(obs; breakpoint, simulated = nothing)
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
            increments = simulated === nothing ? missing : simulated.onsets,
        ),
        breakpoint = breakpoint,
        background_pooling = background_pooling_model,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days,
        n_patches = length(PROVINCE_NAMES),
        province_increments =
            simulated === nothing ? missing : simulated.province_cases,
        province_days = prov_cases.days,
        province_death_increments =
            simulated === nothing ? missing : simulated.province_deaths,
        province_death_days = prov_deaths.days,
        province_care_args(obs; simulated)...,
        simulated_data = simulated === nothing ? nothing : simulated.streams
    )
end

## The province laboratory, occupancy and bed rows at the thinning
## `patch_fit_args` fits them with, their counts `missing` so the generator
## draws them. The occupancy and bed rows keep each day's printed sum, which
## the splits are conditional on. With `simulated` the laboratory counts are
## the simulated ones; the occupancy and bed counts reach the treatment
## stream through its simulated observations.
function province_care_args(obs; simulated = nothing)
    lab = province_lab_increment_matrix(
        obs.province_lab_daily_history, PROVINCE_NAMES,
        length(PROVINCE_NAMES); every = 7
    )
    function generated(rows)
        totals = [sum(rows.counts[rows.days .== d]) for d in rows.days]
        return (; rows.days, rows.patches, counts = missing, totals)
    end
    return (;
        province_lab_increments =
            simulated === nothing ? missing : simulated.province_lab,
        province_lab_days = lab.days, province_lab_bins = lab.bins,
        province_isolation = generated(
            province_care_observations(
                obs.province_isolation_history, PROVINCE_NAMES; every = 7
            )
        ),
        province_capacity = generated(
            province_care_observations(
                obs.province_bed_capacity_history, PROVINCE_NAMES;
                changes_only = true
            )
        ),
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

The verdict on a [`recovery_table`](@ref), and on the fit behind it when
`diagnostics` (a [`fit_diagnostics`](@ref) result) is given. `status` is
one of:

- `:unconverged`: the fit's worst R-hat is above `max_rhat` or its smallest
  bulk ESS below `min_ess`, so its intervals say nothing about the model
  and the recovery is not judged.
- `:fail`: a quantity's true value lies outside its posterior's central
  `outer` interval (99% by default), or fewer than `coverage_fail` of the
  quantities have the truth inside their 90% interval.
- `:warn`: fewer than `coverage_warn` do. A correct model misses a 90%
  interval about one time in ten by chance, so the bar sits below 0.9.
- `:pass` otherwise.

Returns `(; status, pass, outside, coverage_90, converged)`, where `pass`
is `status` being `:pass` or `:warn` and `outside` lists the quantities
outside the outer interval.
"""
function recovery_verdict(
        tab::DataFrame; diagnostics = nothing, outer::Real = 0.99,
        coverage_fail::Real = 0.6, coverage_warn::Real = 0.8,
        max_rhat::Real = 1.05, min_ess::Real = 100
    )
    tail = (1 - outer) / 2
    outside = String[
        r.quantity for r in eachrow(tab)
            if r.truth_quantile < tail || r.truth_quantile > 1 - tail
    ]
    coverage_90 = isempty(tab.covered_90) ? 1.0 : mean(tab.covered_90)
    converged = diagnostics === nothing || (
        diagnostics.max_rhat <= max_rhat &&
            diagnostics.min_ess_bulk >= min_ess
    )
    status = !converged ? :unconverged :
        (!isempty(outside) || coverage_90 < coverage_fail) ? :fail :
        coverage_90 < coverage_warn ? :warn : :pass
    return (;
        status, pass = status in (:pass, :warn), outside, coverage_90,
        converged,
    )
end

"""
$(TYPEDSIGNATURES)

The [`recovery_verdict`](@ref) of each seed in `params`, the
[`recovery_table`](@ref) rows of every seed stacked with the `seed`,
`fit_minutes`, `max_rhat`, `min_ess_bulk` and `divergences` columns that
`scripts/recovery.jl` adds. Returns one row per seed with its `status`,
90% coverage, the quantities outside the 99% interval and its
convergence diagnostics. Other keywords pass to [`recovery_verdict`](@ref).
"""
function recovery_seed_verdicts(params::DataFrame; kwargs...)
    rows = map(sort(unique(params.seed))) do s
        tab = params[params.seed .== s, :]
        d = (;
            max_rhat = first(tab.max_rhat),
            min_ess_bulk = first(tab.min_ess_bulk),
        )
        v = recovery_verdict(tab; diagnostics = d, kwargs...)
        (;
            seed = s, status = v.status, coverage_90 = v.coverage_90,
            outside = join(v.outside, ", "),
            fit_minutes = first(tab.fit_minutes), d.max_rhat,
            d.min_ess_bulk, divergences = first(tab.divergences),
        )
    end
    return DataFrame(rows)
end

"""
$(TYPEDSIGNATURES)

One verdict over several seeds' statuses (from
[`recovery_seed_verdicts`](@ref)): the worst of them, in the order `:fail`,
`:unconverged`, `:warn`, `:pass`.
"""
function recovery_overall(statuses::AbstractVector{Symbol})
    for s in (:fail, :unconverged, :warn)
        s in statuses && return s
    end
    return :pass
end

"""
$(TYPEDSIGNATURES)

Recovery across seeds, one row per quantity, from the stacked
[`recovery_table`](@ref) rows of every seed (with a `seed` column). Each row
has the number of seeds, the range of the true values, the posterior
median's error relative to the truth (`(median - truth) / |truth|`, missing
for a true value of zero) and the z-score as their median and range across
seeds, how many seeds have the truth inside the 90% interval
(`covered_90`) and outside the central `outer` interval (`outside`), and
the seed with the largest absolute z (`worst_seed`, `worst_z`).
"""
function recovery_summary(params::DataFrame; outer::Real = 0.99)
    tail = (1 - outer) / 2
    rows = map(unique(params.quantity)) do q
        tab = params[params.quantity .== q, :]
        rel = [
            t == 0 ? missing : (m - t) / abs(t)
                for (m, t) in zip(tab.median, tab.truth)
        ]
        known = collect(skipmissing(rel))
        worst = argmax(abs.(tab.z))
        (;
            quantity = q, n_seeds = nrow(tab),
            truth_min = minimum(tab.truth), truth_max = maximum(tab.truth),
            rel_error_median = isempty(known) ? missing : median(known),
            rel_error_min = isempty(known) ? missing : minimum(known),
            rel_error_max = isempty(known) ? missing : maximum(known),
            z_median = median(tab.z), z_min = minimum(tab.z),
            z_max = maximum(tab.z), covered_90 = count(tab.covered_90),
            outside = count(
                q -> q < tail || q > 1 - tail, tab.truth_quantile
            ),
            worst_seed = tab.seed[worst], worst_z = tab.z[worst],
        )
    end
    return DataFrame(rows)
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
are the simulated data a recovery fit is fitted to.
"""
function recovery_observed_varnames(generator::Model, fitted::Model)
    rng = MersenneTwister(1)
    gen = keys(last(init!!(rng, generator, VarInfo(), InitFromPrior())))
    fit = Set(keys(last(init!!(rng, fitted, VarInfo(), InitFromPrior()))))
    return [vn for vn in gen if !(vn in fit)]
end

## One draw's value of `vn` from a one-draw chain. A whole-row variable,
## such as a composition's `obs_increments[p, :]`, is stored cell by cell
## (`obs_increments[p, i]`), so its row is gathered from those cells.
function _draw_value(chain, vn)
    name = string(vn)
    m = match(r"^(.*)\[(\d+), :\]$", name)
    m === nothing && return only(chain[vn])
    prefix, row = m[1], parse(Int, m[2])
    ## Chain keys print wrapped, as `Parameter(<name>)`.
    pat = Regex(
        "(?:^|\\()" * replace(prefix, "." => "\\.") *
            "\\[$row, (\\d+)\\]\\)?\$"
    )
    cells = [
        (parse(Int, c[1]), only(chain[k]))
            for k in keys(chain) for c in (match(pat, string(k)),)
            if c !== nothing
    ]
    return [v for (_, v) in sort(cells; by = first)]
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
in-window observations (the names in `observed`) keyed by `VarName`,
which [`recovery_data`](@ref) arranges for [`generator_joint`](@ref).
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
        data = Dict(vn => _draw_value(truth, vn) for vn in observed)
        return (; truth, data)
    end
    throw(
        ErrorException(
            "no prior draw passed the acceptance check in $attempts attempts."
        )
    )
end

## The simulated values of the draw by name: `(stream, observation) =>`
## the values in index order, and the whole-row composition draws by row.
function _grouped_observations(data::AbstractDict)
    streams = Dict{Symbol, Dict{Symbol, Vector{Tuple{Int, Any}}}}()
    onsets = Tuple{Int, Any}[]
    rows = Dict{Symbol, Vector{Tuple{Int, Any}}}()
    for (vn, v) in data
        name = string(vn)
        m = match(r"^(\w+)\.obs_increments\[(\d+), :\]$", name)
        if m !== nothing
            push!(
                get!(rows, Symbol(m[1]), Tuple{Int, Any}[]),
                (parse(Int, m[2]), v)
            )
            continue
        end
        m = match(r"^onset_report_state\.increments(?:\[(\d+)\])?$", name)
        if m !== nothing
            ## Element by element, or the whole triangle as one vector.
            m[1] === nothing ? append!(onsets, collect(enumerate(v))) :
                push!(onsets, (parse(Int, m[1]), v))
            continue
        end
        m = match(r"^(\w+)\.(\w+)\.\w+\[(\d+)\]$", name)
        if m !== nothing
            obs = get!(
                streams, Symbol(m[1]), Dict{Symbol, Vector{Tuple{Int, Any}}}()
            )
            push!(
                get!(obs, Symbol(m[2]), Tuple{Int, Any}[]),
                (parse(Int, m[3]), v)
            )
            continue
        end
        ## A stream drawn as one whole vector rather than element by element.
        m = match(r"^(\w+)\.(\w+)\.\w+$", name)
        m === nothing && throw(
            ArgumentError("no rule for the simulated observation `$name`.")
        )
        obs = get!(streams, Symbol(m[1]), Dict{Symbol, Vector{Tuple{Int, Any}}}())
        append!(get!(obs, Symbol(m[2]), Tuple{Int, Any}[]), collect(enumerate(v)))
    end
    in_order(xs) = [v for (_, v) in sort(xs; by = first)]
    return (; streams, onsets = in_order(onsets), rows, in_order)
end

## A composition's simulated counts as the full province-by-vintage matrix:
## the drawn rows, then the last province as the remainder of the recorded
## totals, which is how the predictive path fills it in.
function _composition_matrix(rows, totals)
    drawn = reduce(vcat, [reshape(Int.(round.(collect(r))), 1, :) for r in rows])
    last_row = reshape(
        max.(Int.(round.(totals)) .- vec(sum(drawn; dims = 1)), 0), 1, :
    )
    return vcat(drawn, last_row)
end

"""
$(TYPEDSIGNATURES)

The observations of a [`simulate_recovery`](@ref) dataset in the form
[`generator_joint`](@ref) takes them: each stream's counts by observation
name (`streams`), the onset triangle increments (`onsets`) and the province
case and death compositions as province-by-vintage count matrices
(`province_cases`, `province_deaths`). Each composition's last province is
the remainder of its recorded totals, as the predictive path fills it in.
"""
function recovery_data(sim)
    g = _grouped_observations(sim.data)
    ## The occupancy and bed splits are drawn a position at a time, so
    ## their counts are read whole from the draw's recorded `split_counts`.
    splits = (:occupancy_split, :capacity_split)
    streams = Dict{Symbol, Any}(
        state => NamedTuple(
            name => (
                name in splits ?
                    Int.(
                        only(
                            _draw_vectors(
                                sim.truth,
                                Symbol("$(state).$(name).split_counts")
                            )
                        )
                    ) :
                    Int.(round.(g.in_order(xs)))
            ) for (name, xs) in obs
        )
            for (state, obs) in g.streams
    )
    composition(state) = haskey(g.rows, state) ? _composition_matrix(
            g.in_order(g.rows[state]),
            only(_draw_vectors(sim.truth, Symbol("$(state).composition_totals")))
        ) : missing
    return (;
        streams, onsets = isempty(g.onsets) ? missing : Float64.(g.onsets),
        province_cases = composition(:composition_state),
        province_deaths = composition(:death_composition_state),
        province_lab = composition(:lab_composition_state),
    )
end

"""
$(TYPEDSIGNATURES)

Check that `model`, the generator rebuilt with a simulated dataset
([`generator_joint`](@ref) with `simulated`), scores the true draw `truth`
exactly as `generator` does. The generator's log joint at the true draw
counts the simulated observations as draws; the rebuilt model counts them
as data. The two agree only when every simulated observation reached the
stream it came from, so a mismatch throws rather than letting a recovery
fit run on data the model did not simulate. Returns the two log joints.
"""
function recovery_density_check(
        generator::Model, model::Model, truth; rtol::Real = 1.0e-8
    )
    from_generator = only(logjoint(generator, truth))
    from_data = only(logjoint(model, truth))
    isapprox(from_generator, from_data; rtol) || throw(
        ErrorException(
            "the rebuilt model scores the true draw at $from_data, not the " *
                "generator's $from_generator: a simulated observation is " *
                "not reaching its stream."
        )
    )
    return (; from_generator, from_data)
end

"""
$(TYPEDSIGNATURES)

Fit the headline joint to a [`simulate_recovery`](@ref) dataset with NUTS:
`model` is the generator rebuilt with the simulated observations
([`generator_joint`](@ref) with `simulated = recovery_data(sim)`), after
[`recovery_density_check`](@ref). Returns `(; model, chain)`;
[`forecast_draws`](@ref) runs the model past its cut-off. The sampler
settings default to the headline joint's (two chains of 1000 draws after
500 warmup steps, tree depth at most 10): a shorter run does not converge
on this model, so its recovery would test the sampler rather than the
model. Other keywords pass to [`nuts_sample`](@ref).
"""
function recovery_fit(
        model::Model; samples::Integer = 1000, n_adapts::Integer = 500,
        chains::Integer = 2, max_depth::Integer = 10,
        seed::Integer = 20260518, kwargs...
    )
    chain = nuts_sample(
        model; samples, n_adapts, chains, max_depth, seed, kwargs...
    )
    return (; model, chain)
end
