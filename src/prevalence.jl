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

# Suspected-case ascertainment by patch: the national `p_drc` times the
# relative ascertainment `rel`, rescaled so that the patches weighted by
# `w` keep the national rate, as the joint splits the national suspect
# volume by ascertainment-weighted incidence. The weights are the patches'
# expected confirmations on the cut-off day
# ([`_cutoff_confirmations`](@ref)), the incidence the joint splits on that
# day. The joint uses the relative ascertainment only to split national
# counts, so nothing in it bounds the product, and it is capped at one.
function _patch_ascertainment(p_drc::Real, rel::AbstractVector, w::AbstractVector)
    weighted = sum(rel .* w)
    scale = weighted > 0 ? sum(w) / weighted : 1.0
    return [min(p_drc * r * scale, 1.0) for r in rel]
end

# Each area's expected confirmations on the last day: the infections
# `x[a]` through the infection-to-confirmation PMF `f`
# (incubation, onset to report and report to receipt), the onset-to-
# confirmation kernel the joint applies to each patch's onsets.
_cutoff_confirmations(x, f::AbstractVector) = [point_prevalence(xa, f) for xa in x]

# The infection-to-confirmation PMF of one draw of the joint, discretised as
# the delays of [`_prevalence_delays`](@ref).
function _confirmation_pmf(p, receipt_mean::Real, receipt_sd::Real, horizon::Integer)
    pmf(dist) = discretise_censored(dist, horizon)
    return convolve_pmf(
        convolve_pmf(
            pmf(lognormal_meansd(p.inc_mean, p.inc_sd)),
            pmf(Distributions.Gamma(p.report_alpha, p.report_theta))
        ),
        pmf(lognormal_meansd(receipt_mean, receipt_sd))
    )
end

# Per-draw delay parameters and probabilities from a joint chain, one named
# tuple per draw. The ascertainment is one value per patch, weighted by the
# patches' expected confirmations on the cut-off day, or the
# national value alone when the chain has no province ascertainment.
function _prevalence_parameters(chn, horizon::Integer = 120)
    d(vn) = Float64.(vec(collect(chn[vn])))
    p_drc = d(:p_drc)
    ## A chain without the province death split holds only the national CFR.
    ## The joint builds `CFR_patch` on the logit scale, so it lies in (0, 1).
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
    delays = [(; map(c -> c[i], cols)...) for i in eachindex(p_drc)]
    asc = if _has_key(chn, :province_ascertainment) &&
            _has_key(chn, :infections_patch)
        rel = _draw_vectors(chn, :province_ascertainment)
        inf = _draw_vectors(chn, :infections_patch)
        receipt_mean = d(_ZONE_PARENT_KEYS.receipt_mean)
        receipt_sd = d(_ZONE_PARENT_KEYS.receipt_sd)
        map(eachindex(p_drc)) do i
            np = length(rel[i])
            f = _confirmation_pmf(
                delays[i], receipt_mean[i], receipt_sd[i], horizon
            )
            ## `infections_patch` is day-major, so patch `q` is every `np`th
            ## entry from `q`.
            w = _cutoff_confirmations([inf[i][q:np:end] for q in 1:np], f)
            _patch_ascertainment(p_drc[i], rel[i], w)
        end
    else
        [[a] for a in p_drc]
    end
    return [
        (; delays[i]..., ascertainment = asc[i], cfr = cfr[i])
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
            ascertainment = p.ascertainment[min(patch_of[a], end)],
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
[`infection_state_weights`](@ref). The suspected-case ascertainment and the
case-fatality ratio are the patch's own, or the national ones when the
chain has no province composition or death split. The ascertainment is
capped at one.
"""
function patch_prevalence(
        chn, n_patches::Integer = length(PROVINCE_NAMES);
        horizon::Integer = 120
    )
    n = length(first(_draw_vectors(chn, :infections_patch))) ÷ n_patches
    daily = _patch_daily(chn, :infections_patch, n_patches, n)
    series = [reduce(vcat, permutedims.(d)) for d in daily]
    return _prevalence_draws(
        series, _prevalence_parameters(chn, horizon), 1:n_patches, horizon
    )
end

# Per-draw parameters for the zones: zone draw `j` takes the joint's draw
# `j` (cycling), its own delays from `θ[j]`, and one ascertainment and one
# case-fatality ratio per zone. The ascertainment is the province's times the
# zone's relative ascertainment `asc[j]`, rescaled within the province by the
# weights `w[j]`, the zones' expected confirmations on the cut-off day,
# as [`_patch_ascertainment`](@ref) rescales the provinces. The case-fatality
# ratio is the province's times the zone's relative fatality `sev[j]`, the
# multiplier the zone stage fits the within-province deaths with. That
# multiplier is not bounded, so the zone ratio is capped at one, as the
# ascertainment is.
function _zone_prevalence_parameters(joint, θ, asc, sev, w, patch_of)
    at(k) = findfirst(==(k), _ZONE_DELAY_KEYS)
    own = (
        :inc_mean, :inc_sd, :report_alpha, :report_theta, :admission_alpha,
        :admission_theta, :death_alpha, :death_theta,
    )
    return map(eachindex(θ)) do j
        p = joint[mod1(j, length(joint))]
        ascertainment = zeros(length(patch_of))
        for pr in unique(patch_of)
            zs = findall(==(pr), patch_of)
            ascertainment[zs] = _patch_ascertainment(
                p.ascertainment[min(pr, end)], asc[j][zs], w[j][zs]
            )
        end
        cfr = [
            min(p.cfr[min(patch_of[a], end)] * sev[j][a], 1.0)
                for a in eachindex(patch_of)
        ]
        merge(p, NamedTuple{own}(map(k -> θ[j][at(k)], own)), (; ascertainment, cfr))
    end
end

"""
Point prevalence at the cut-off by health zone, draw by draw, on the daily
zone infections of [`zone_infections`](@ref), in the order of
`inputs.zone_keys`. Each zone draw takes its incubation, detection,
admission and death delays from its own draw of the joint's delays. The
zone stage does not estimate the recovery stay, so zone draw `j` takes it
from draw `j` of `joint_chn`, the joint the zone stage was fitted from,
cycling when the zone chain is longer. It estimates ascertainment and
fatality only relative to the zone's province
(`zone_ascertainment_relative`, `zone_severity_relative`). A zone's
case-fatality ratio is its province's from that joint draw times the zone's
contrast. Its suspected-case ascertainment is the province's times the
zone's contrast, rescaled so that the province's zones, weighted by their
expected confirmations on the cut-off day, keep the province's rate. Both are capped at
one. The ascertainment contrast is fitted to confirmed cases, so this assumes
suspected-case detection varies between the zones of a province as
confirmed-case detection does. Returns the states of
[`patch_prevalence`](@ref).
"""
function zone_prevalence(
        zone_chn, inputs, joint_chn;
        horizon::Integer = 120
    )
    series = zone_infections(zone_chn, inputs)
    θ = _draw_vectors(zone_chn, :delay_parameters_zone)
    w = [
        _cutoff_confirmations(
            [view(x, j, :) for x in series], zone_delay_pmfs(θ[j]).f
        ) for j in eachindex(θ)
    ]
    params = _zone_prevalence_parameters(
        _prevalence_parameters(joint_chn, horizon), θ,
        _draw_vectors(zone_chn, :zone_ascertainment_relative),
        _draw_vectors(zone_chn, :zone_severity_relative), w,
        inputs.patch_of_zone
    )
    return _prevalence_draws(
        series, params, eachindex(inputs.patch_of_zone), horizon
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
        "Infected but not yet symptomatic per 100,000" => String[],
        "Symptomatic, undetected per 100,000" => String[],
        "In the community per 100,000" => String[],
        "Number in the community" => String[]
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
Dot plot of the rate in the community per 100,000 by area, from a
[`prevalence_overview`](@ref): each a median with its 90% credible interval,
one row per area coloured by patch and sorted by the rate. `max_rows` keeps
the highest.
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
    fig = Figure(; size = (680, max(26 * rows, 180) + 170))
    ax = Axis(
        fig[1, 1]; xlabel = "In the community per 100,000",
        yticks = (ys, String.(ov.label))
    )
    for k in 1:rows
        c = patch_colours[mod1(ov.patch[k], length(patch_colours))]
        lines!(
            ax, [ov[k, "community_rate_lo90"], ov[k, "community_rate_hi90"]],
            [ys[k], ys[k]]; color = (c, 0.6), linewidth = 2
        )
        scatter!(
            ax, [ov[k, "community_rate_median"]], [ys[k]];
            color = c, markersize = 9
        )
    end
    CairoMakie.xlims!(ax, 0, nothing)
    CairoMakie.ylims!(ax, 0.4, rows + 0.6)
    _patch_legend!(fig, (2, 1), ov.patch, patch_labels, patch_colours)
    CairoMakie.Label(
        fig[3, 1],
        "Median and 90% credible interval per 100,000 residents. In the " *
            "community: infected but not yet symptomatic, or symptomatic and " *
            "not yet detected.";
        fontsize = 12, word_wrap = true, padding = (0, 0, 0, 6)
    )
    CairoMakie.Label(fig[0, 1], title; fontsize = 16, font = :bold)
    return fig
end
