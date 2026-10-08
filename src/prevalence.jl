# Point prevalence of infection: the people infected at the cut-off who are
# still in the community, by stage of infection. An infection's stage follows
# from its age and the model's delays (incubation, onset to detection, onset
# to death or recovery). Prevalence is therefore a weighted sum of the recent
# daily infections, computed for each posterior draw of the joint fit and of
# the health-zone stage.

# P(delay > s) for s = 0, …, horizon − 1 from a PMF indexed from lag 0.
function _survival(pmf::AbstractVector, horizon::Integer)
    F = cumsum(pmf)
    return [
        max(1 - (s < length(F) ? F[s + 1] : F[end]), 0.0)
            for s in 0:(horizon - 1)
    ]
end

"""
Probability that a person infected `s` days before a reference day is, on
that day, in each infection state, at index `s + 1` for
`s = 0, …, horizon − 1`.

`incubation`, `detection`, `death` and `recovery` are delay PMFs indexed
from lag 0: infection to onset, onset to detection as a suspected case,
onset to death and onset to recovery. A symptomatic infection is detected
with probability `ascertainment` and leaves the community at detection.
One that is not detected stays until it dies, with probability `cfr`, or
recovers.

Returns `(; incubating, symptomatic)`. `incubating` is the probability
that onset has not yet happened. `symptomatic` is the probability of being
past onset and still undetected, alive and unrecovered.
"""
function infection_state_weights(
        incubation::AbstractVector, detection::AbstractVector,
        death::AbstractVector, recovery::AbstractVector;
        ascertainment::Real, cfr::Real, horizon::Integer
    )
    horizon >= 1 || throw(ArgumentError("horizon must be at least 1"))
    (0 <= ascertainment <= 1 && 0 <= cfr <= 1) || throw(
        ArgumentError("ascertainment and cfr must lie in [0, 1]")
    )
    incubating = _survival(incubation, horizon)
    since_onset = ascertainment .* _survival(detection, horizon) .+
        (1 - ascertainment) .* (
        cfr .* _survival(death, horizon) .+
            (1 - cfr) .* _survival(recovery, horizon)
    )
    symptomatic = convolve_delay(since_onset, incubation)
    return (; incubating, symptomatic)
end

"""
People in a state on the last day of `infections`, the daily infections
over the day grid, given the probability `weights[s + 1]` that an infection
`s` days old is in it ([`infection_state_weights`](@ref)).
"""
function point_prevalence(
        infections::AbstractVector, weights::AbstractVector
    )
    n = length(infections)
    return sum(
        infections[n - s] * weights[s + 1]
            for s in 0:(min(n, length(weights)) - 1);
        init = 0.0
    )
end

const PREVALENCE_STATES = (:incubating, :symptomatic, :community)

# Per-draw delay parameters and probabilities from a joint chain, one named
# tuple per draw.
function _prevalence_parameters(chn)
    d(vn) = Float64.(vec(collect(chn[vn])))
    asc = d(:p_drc)
    ## A chain without the province death split holds only the national CFR.
    cfr = _has_key(chn, :CFR_patch) ? _draw_vectors(chn, :CFR_patch) :
        [[c] for c in d(:CFR)]
    cols = (;
        inc_mean = d(@varname(inc_state.delay_mean)),
        inc_sd = d(@varname(inc_state.delay_sd)),
        report_alpha = d(@varname(cases_state.report_state.α)),
        report_theta = d(@varname(cases_state.report_state.θ)),
        admission_alpha = d(@varname(deaths_state.od_state.oa.α)),
        admission_theta = d(@varname(deaths_state.od_state.oa.θ)),
        death_alpha = d(@varname(deaths_state.od_state.ad.α)),
        death_theta = d(@varname(deaths_state.od_state.ad.θ)),
        recovery_mean = d(
            @varname(treatment_state.recovery_los_state.delay_mean)
        ),
        recovery_sd = d(@varname(treatment_state.recovery_los_state.delay_sd)),
    )
    return [
        (; map(c -> c[i], cols)..., ascertainment = asc[i], cfr = cfr[i])
            for i in eachindex(asc)
    ]
end

# The delay PMFs of one draw's parameters, discretised as the model
# discretises them. Death and recovery run from onset through admission.
function _prevalence_delays(p, horizon::Integer)
    pmf(dist) = discretise_censored(dist, horizon)
    admission = pmf(Distributions.Gamma(p.admission_alpha, p.admission_theta))
    return (;
        incubation = pmf(lognormal_meansd(p.inc_mean, p.inc_sd)),
        detection = pmf(Distributions.Gamma(p.report_alpha, p.report_theta)),
        death = convolve_pmf(
            admission, pmf(Distributions.Gamma(p.death_alpha, p.death_theta))
        ),
        recovery = convolve_pmf(
            admission, pmf(lognormal_meansd(p.recovery_mean, p.recovery_sd))
        ),
    )
end

# Prevalence in each state, one vector of draws per area. `series[a]` holds
# area `a`'s daily infections, one row per draw, and `params[j]` the
# parameters of draw `j`.
function _prevalence_draws(series, params, patch_of, horizon)
    na = length(series)
    nd = size(first(series), 1)
    out = NamedTuple{PREVALENCE_STATES}(
        ntuple(_ -> [zeros(nd) for _ in 1:na], length(PREVALENCE_STATES))
    )
    delays = [_prevalence_delays(p, horizon) for p in params]
    for a in 1:na, j in 1:nd
        p, dl = params[j], delays[j]
        w = infection_state_weights(
            dl.incubation, dl.detection, dl.death, dl.recovery;
            ascertainment = p.ascertainment,
            cfr = p.cfr[min(patch_of[a], end)], horizon
        )
        x = view(series[a], j, :)
        out.incubating[a][j] = point_prevalence(x, w.incubating)
        out.symptomatic[a][j] = point_prevalence(x, w.symptomatic)
        out.community[a][j] = out.incubating[a][j] + out.symptomatic[a][j]
    end
    return out
end

"""
Point prevalence at the cut-off by province, draw by draw, from a joint
chain of [`bvd_joint`](@ref) with `n_patches` patches. Each state is one
vector of draws per patch: `incubating`, `symptomatic` (past onset and
undetected) and `community` (incubating or symptomatic and undetected), as in
[`infection_state_weights`](@ref). Detection is the national suspected-case
ascertainment and the case-fatality ratio is the patch's own, or the
national one when the chain has no province death split.
"""
function patch_prevalence(
        chn, n_patches::Integer = length(PROVINCE_NAMES);
        horizon::Integer = 120
    )
    n = length(first(_draw_vectors(chn, :infections_patch))) ÷ n_patches
    daily = _patch_daily(chn, :infections_patch, n_patches, n)
    series = [reduce(vcat, permutedims.(d)) for d in daily]
    return _prevalence_draws(
        series, _prevalence_parameters(chn), 1:n_patches, horizon
    )
end

"""
Point prevalence at the cut-off by health zone, draw by draw, on the daily
zone infections of [`zone_infections`](@ref), in the order of
`inputs.zone_keys`. Each zone draw takes its incubation, detection,
admission and death delays from its own draw of the joint's delays. The
zone stage does not estimate the recovery stay, the detection probability
or the case-fatality ratio, so zone draw `j` takes them from draw `j` of
`joint_chn`, the joint the zone stage was fitted from, cycling when the zone
chain is longer. Returns the states of [`patch_prevalence`](@ref).
"""
function zone_prevalence(
        zone_chn, inputs, joint_chn;
        horizon::Integer = 120
    )
    series = zone_infections(zone_chn, inputs)
    joint = _prevalence_parameters(joint_chn)
    θ = _draw_vectors(zone_chn, :delay_parameters_zone)
    at(k) = findfirst(==(k), _ZONE_DELAY_KEYS)
    own = (
        :inc_mean, :inc_sd, :report_alpha, :report_theta, :admission_alpha,
        :admission_theta, :death_alpha, :death_theta,
    )
    params = [
        merge(
            joint[mod1(j, length(joint))],
            NamedTuple{own}(map(k -> θ[j][at(k)], own))
        ) for j in eachindex(θ)
    ]
    return _prevalence_draws(
        series, params, inputs.patch_of_zone, horizon
    )
end

"""
One row per area of a [`patch_prevalence`](@ref) or
[`zone_prevalence`](@ref) result: `label`, `patch`, `population`, and for
each state the median and 90% credible interval of the number of people
(`<state>_median`, `_lo90`, `_hi90`) and of the rate per 100,000
(`<state>_rate_median`, and so on). Rows are sorted by the median
community rate, highest first.
"""
function prevalence_overview(
        prev, populations::AbstractVector;
        labels::AbstractVector, patch::AbstractVector{<:Integer}
    )
    na = length(labels)
    sizes = (length(populations), length(patch), length(first(prev)))
    all(==(na), sizes) || throw(
        DimensionMismatch(
            "labels, patch, populations and the prevalence draws must " *
                "cover the same areas"
        )
    )
    df = DataFrame(; label = labels, patch, population = populations)
    for s in PREVALENCE_STATES, (suffix, scale) in (
                ("", ones(na)), ("_rate", 1.0e5 ./ populations),
            )
        draws = [prev[s][a] .* scale[a] for a in 1:na]
        df[!, "$(s)$(suffix)_median"] = median.(draws)
        df[!, "$(s)$(suffix)_lo90"] = quantile.(draws, 0.05)
        df[!, "$(s)$(suffix)_hi90"] = quantile.(draws, 0.95)
    end
    return sort!(df, :community_rate_median; rev = true)
end

"""
The [`prevalence_overview`](@ref) rows as a display table, each cell a
median with its 90% credible interval: the rate per 100,000 in each state
and the number of people in the community. `max_rows` keeps the first rows.
"""
function prevalence_table(
        overview::DataFrame; max_rows = nothing,
        area::AbstractString = "Area"
    )
    rows = max_rows === nothing ? nrow(overview) :
        min(max_rows, nrow(overview))
    cell(r, s; digits) = string(
        _fmt_prevalence(r["$(s)_median"], digits), " (",
        _fmt_prevalence(r["$(s)_lo90"], digits), "–",
        _fmt_prevalence(r["$(s)_hi90"], digits), ")"
    )
    df = DataFrame(
        area => String[], "Population" => String[],
        "Incubating" => String[],
        "Symptomatic, undetected" => String[],
        "In the community" => String[],
        "People in the community" => String[]
    )
    for r in eachrow(overview[1:rows, :])
        push!(
            df, [
                r.label, string(r.population),
                cell(r, "incubating_rate"; digits = 1),
                cell(r, "symptomatic_rate"; digits = 1),
                cell(r, "community_rate"; digits = 1),
                cell(r, "community"; digits = 0),
            ]
        )
    end
    return df
end

_fmt_prevalence(x, digits) = digits <= 0 ? string(round(Int, x)) :
    string(round(x; digits))

"""
Dot plot of point prevalence per 100,000 by area, from a
[`prevalence_overview`](@ref): the left panel is the rate incubating, the
right the rate in the community, each a median with its 90% credible
interval, one row per area coloured by patch and sorted by the community
rate. `max_rows` keeps the highest.
"""
function plot_prevalence_ranking(
        overview::DataFrame;
        patch_labels::AbstractVector = PROVINCE_LABELS,
        patch_colours = _ZONE_PATCH_COLOURS, max_rows = nothing,
        title::AbstractString = "Point prevalence at the cut-off"
    )
    rows = max_rows === nothing ? nrow(overview) :
        min(max_rows, nrow(overview))
    ov = overview[1:rows, :]
    ys = Float64.(rows:-1:1)
    fig = Figure(; size = (760, max(26 * rows, 180) + 170))
    ax1 = Axis(
        fig[1, 1]; xlabel = "Incubating per 100,000",
        yticks = (ys, String.(ov.label))
    )
    ax2 = Axis(fig[1, 2]; xlabel = "In the community per 100,000")
    CairoMakie.linkyaxes!(ax1, ax2)
    CairoMakie.hideydecorations!(ax2; grid = false)
    for (ax, s) in ((ax1, "incubating_rate"), (ax2, "community_rate"))
        for k in 1:rows
            c = patch_colours[mod1(ov.patch[k], length(patch_colours))]
            lines!(
                ax, [ov[k, "$(s)_lo90"], ov[k, "$(s)_hi90"]], [ys[k], ys[k]];
                color = (c, 0.6), linewidth = 2
            )
            scatter!(
                ax, [ov[k, "$(s)_median"]], [ys[k]];
                color = c, markersize = 9
            )
        end
        CairoMakie.xlims!(ax, 0, nothing)
    end
    CairoMakie.ylims!(ax1, 0.4, rows + 0.6)
    _patch_legend!(fig, (2, 1:2), ov.patch, patch_labels, patch_colours)
    CairoMakie.Label(
        fig[3, 1:2],
        "Median and 90% credible interval per 100,000 residents. " *
            "Incubating: infected, before symptom onset. In the community: " *
            "incubating, or symptomatic and not yet detected.";
        fontsize = 12, word_wrap = true, padding = (0, 0, 0, 6)
    )
    CairoMakie.Label(fig[0, 1:2], title; fontsize = 16, font = :bold)
    return fig
end
