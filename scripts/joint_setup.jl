# Shared setup for the deaths-background fitting harnesses
# (`fit_deaths_background.jl`, `cascade_priors.jl`): builds the joint
# fit-args from the loaded observations, the cut-off-aware growth submodel,
# moment-matching helpers for cascaded priors, and the streaming-callback
# watchdog. Kept dependency-light so both harnesses can `include` it.

using BVDOutbreakSize
import Turing.AbstractMCMC as AbstractMCMC
using Distributions: truncated, Normal, Beta, Gamma
using Statistics: mean, std, var
using Printf: @sprintf

## ------------------------------------------------------------------
## Fit-arg construction (mirrors docs/examples/analysis.jl `joint_obs`):
## per-vintage increment streams for deaths, reported and confirmed cases,
## the lab received/analysed denominators, the single confirmed-death total
## and the travel-gated export series.
## ------------------------------------------------------------------
function _inc(v)
    d = similar(v, Int)
    prev = 0
    for i in eachindex(v)
        d[i] = v[i] - prev
        prev = v[i]
    end
    return d
end

function build_fit_args(o)
    rh = o.reported_case_history
    dh = o.death_history
    ch = o.confirmed_case_history
    sa = o.tests_analysed_history
    sr = o.tests_received_history

    rep, rep_off = _inc(rh.values), rh.offsets
    dth, dth_off = _inc(dh.values), dh.offsets

    ## FULL confirmed series for the missing-lab-data queue: every vintage
    ## 18 May-3 June enters as a between-vintage increment. The 23-28 May
    ## windows carry the published analysed/received denominators (with the
    ## 25 May analysis stall merged into the next window so every observed
    ## denominator is positive); the early 18-22 May AND late 29 May-3 June
    ## windows have no national analysed total, so they are DARK
    ## (`missing` analysed/received) and the queue owns them via the
    ## Poisson-thinned marginal. Offsets come straight from the vintage
    ## histories, so they track the as_of_date (3 June) automatically.
    coff = collect(ch.offsets)
    conf = Union{Missing, Int}[_inc(ch.values)...]
    ## Observed-denominator windows, stall-merged: keep the first analysed
    ## vintage and any with strictly more cumulative analysed than the
    ## previous kept one (drops the 24-25 May flat 295 stall), then
    ## difference to between-vintage increments.
    keep = [i == 1 || sa.values[i] > sa.values[i - 1]
            for i in eachindex(sa.values)]
    aoff = collect(sa.offsets)[keep]
    a_inc = _inc(sa.values[keep])
    ridx = [findfirst(==(off), sr.offsets) for off in aoff]
    r_inc = _inc([sr.values[i] for i in ridx])
    ## Map observed increments onto the kept offsets; every other confirmed
    ## vintage (early + late dark) gets `missing`.
    aobs = Dict(off => a_inc[k] for (k, off) in enumerate(aoff))
    robs = Dict(off => r_inc[k] for (k, off) in enumerate(aoff))
    analysed = Union{Missing, Int}[get(aobs, off, missing) for off in coff]
    received = Union{Missing, Int}[get(robs, off, missing) for off in coff]

    ## Travel-gated exports truncated at the most recent reported import.
    ec_full = o.exported_cases_daily
    last_import = isempty(ec_full) ? nothing : findlast(!=(0), ec_full)
    export_last_offset = last_import === nothing ? 0 :
                         length(ec_full) - last_import
    _trunc(v) = v[1:max(length(v) - export_last_offset, 0)]
    edaily = _trunc(o.export_deaths_daily)
    ecases = isempty(ec_full) ? ec_full : _trunc(ec_full)

    ## Confirmed deaths span their full dark window too (26 May-3 June): the
    ## stream is a lab/positivity process on the death specimens, not gated
    ## on the analysed denominator, so every vintage enters as an increment
    ## at its own offset.
    cdh = o.confirmed_death_history
    cdeath = Union{Missing, Int}[_inc(cdh.values)...]
    cdeath_off = collect(cdh.offsets)

    return (deaths = dth, reported = rep, export_deaths = edaily,
        kw = (; reported_offsets = rep_off, death_offsets = dth_off,
            confirmed_cases = conf, confirmed_offsets = coff,
            samples_analysed = analysed, samples_received = received,
            confirmed_deaths = cdeath, confirmed_death_offsets = cdeath_off,
            exported_cases_daily = ecases,
            export_last_offset = export_last_offset,
            confirmed_queue = true, confirmed_epi_exclusion = nothing,
            tests_analysed = o.cumulative_tests_analysed, tests_offset = 0))
end

## Cut-off-aware growth submodel (size prior centre advances with the
## cut-off date), matching the report.
function growth_for(o)
    exponential_growth_model(
        m_prior = truncated(Normal(m_prior_centre(o.as_of_date), 3.0);
        lower = 0))
end

## ------------------------------------------------------------------
## Moment-matching helpers for cascaded priors: fit a prior in the same
## family the submodel uses to a posterior draw vector.
## ------------------------------------------------------------------
"""Moment-match a `Beta` to draws on (0, 1) (mean/variance)."""
function fit_beta(x)
    m = clamp(mean(x), 1e-3, 1 - 1e-3)
    v = max(var(x), 1e-6)
    c = max(m * (1 - m) / v - 1, 1e-2)
    return Beta(m * c, (1 - m) * c)
end

"""Moment-match a lower-truncated `Normal` to draws (mean/sd)."""
function fit_trunc_normal(x; lower = 0.0)
    truncated(Normal(mean(x), max(std(x), 1e-3)); lower = lower)
end

## ------------------------------------------------------------------
## Watchdog: abort a fit once divergences run away early.
## ------------------------------------------------------------------
struct EarlyKill <: Exception
    msg::String
end

"""
Watchdog callback: throws `EarlyKill` once the divergent fraction exceeds
`max_div_frac` after at least `min_iter` kept draws, killing a fit stuck in
a bad basin in seconds rather than running to completion. Reads the
divergence flag through the sampler-agnostic `ParamsWithStats` interface.
"""
function watchdog_callback(; min_iter::Integer = 60, max_div_frac::Real = 0.4)
    lk = ReentrantLock()
    ndiv = Ref(0)
    niter = Ref(0)
    return function (rng, model, sampler, transition, state, iteration;
            kwargs...)
        stats = try
            AbstractMCMC.ParamsWithStats(
                model, sampler, transition, state; stats = true).stats
        catch
            (;)
        end
        Base.@lock lk begin
            niter[] += 1
            if get(stats, :numerical_error, false) === true
                ndiv[] += 1
            end
            if niter[] >= min_iter && ndiv[] / niter[] > max_div_frac
                throw(EarlyKill(@sprintf("watchdog: %d/%d (%.0f%%) divergent — killing fit",
                    ndiv[], niter[], 100 * ndiv[] / niter[])))
            end
        end
        return nothing
    end
end

"""Compose several AbstractMCMC callbacks into one (run in order)."""
compose_callbacks(cbs...) = (args...; kwargs...) -> begin
    for cb in cbs
        cb(args...; kwargs...)
    end
    return nothing
end
