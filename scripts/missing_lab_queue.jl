# Fit the joint model over ALL confirmed vintages with the condensed
# laboratory-throughput queue (Poisson-thinned denominator for the dark
# windows that lack a published analysed total). Run under the test env:
#
#   julia --project=test scripts/missing_lab_queue.jl
#
# Streams to TensorBoard (logs/queue) and appends a summary block to
# scripts/out/queue_result.txt after each fit (headline e = 0 and the
# e ~ Beta(2, 12) sensitivity).

using BVDOutbreakSize
using BVDOutbreakSize: bvd_joint, load_observations, nuts_sample,
                       epi_exclusion_model, confirmed_q_re_model,
                       fit_diagnostics, posterior_summary,
                       tensorboard_callback
using TensorBoardLogger
import FlexiChains
using Statistics: median

const OUT = joinpath(@__DIR__, "out", "queue_result.txt")

## Sample / chain counts are overridable from the command line
## (`julia scripts/missing_lab_queue.jl 400 4`) so the fit can be run lean
## under heavy CPU/memory contention. Defaults are 400 samples, 4 chains.
const SAMPLES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 400
const CHAINS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

_inc(values) = begin
    out = similar(collect(values), Int)
    prev = 0
    for i in eachindex(out)
        out[i] = values[i] - prev
        prev = values[i]
    end
    out
end

## Build the FULL confirmed series straight from the vintage histories:
## every confirmed vintage (18 May-3 June) enters as a between-vintage
## increment at its own offset. The 23-28 May windows carry the published
## analysed/received denominators; every other vintage (the early
## 18-22 May and late 29 May onward windows that lack a national analysed
## total) is DARK, carrying `missing` analysed and received counts.
## Offsets come from the data, so the series tracks the as_of_date.
function full_extension(obs)
    ch = obs.confirmed_case_history
    sa = obs.tests_analysed_history
    sr = obs.tests_received_history
    coff = collect(ch.offsets)
    confirmed = _inc(collect(ch.values))
    ## Keep only observed analysed windows whose cumulative strictly
    ## increased: the 24-25 May vintage is flat (295 -> 295), so its
    ## analysed increment is 0 and a Binomial(0, p) cannot observe the
    ## nonzero confirmed increment. Collapsing it onto the next window
    ## keeps every observed denominator positive.
    keep = [i == 1 || sa.values[i] > sa.values[i - 1]
            for i in eachindex(sa.values)]
    aoff = collect(sa.offsets)[keep]
    a_inc = _inc(collect(sa.values)[keep])
    sroff = collect(sr.offsets)
    srval = collect(sr.values)
    ridx = [findfirst(==(off), sroff) for off in aoff]
    r_inc = _inc([srval[i] for i in ridx])
    aobs = Dict(off => a_inc[k] for (k, off) in enumerate(aoff))
    robs = Dict(off => r_inc[k] for (k, off) in enumerate(aoff))
    a_vals = Union{Missing, Int}[get(aobs, off, missing) for off in coff]
    r_vals = Union{Missing, Int}[get(robs, off, missing) for off in coff]
    return (offsets = coff, confirmed = confirmed,
        analysed = a_vals, received = r_vals)
end

function build_model(obs, c; epi_exclusion = nothing)
    rep = obs.reported_case_history
    dh = obs.death_history
    return bvd_joint(obs.exported_cases, _inc(dh.values), _inc(rep.values),
        obs.export_deaths_daily;
        reported_offsets = rep.offsets,
        death_offsets = dh.offsets,
        confirmed_cases = c.confirmed,
        confirmed_offsets = c.offsets,
        samples_analysed = c.analysed,
        samples_received = c.received,
        tests_analysed = obs.cumulative_tests_analysed,
        tests_offset = 0,
        first_export_detection_delta = obs.first_export_detection_delta,
        confirmed_q_random_effect = confirmed_q_re_model,
        confirmed_queue = true,
        confirmed_epi_exclusion = epi_exclusion)
end

ival(s) = "$(round(s.lo90; digits = 0)) to $(round(s.hi90; digits = 0))"

function summarise(io, label, chn, c, obs)
    d = fit_diagnostics(chn)
    ess_tail = FlexiChains.ess(chn; kind = :tail)
    et = filter(isfinite, collect(Iterators.flatten(values(ess_tail))))
    min_ess_tail = isempty(et) ? NaN : minimum(et)
    C = vec(Array(chn[:cumulative_cases]))
    Cs = posterior_summary(C)
    dark = vec(Array(chn[:dark_analysed_total]))
    nw = length(c.offsets)
    ## Per-window predicted analysed mean and received mean (median over
    ## draws), aligned with the confirmed offsets. The deterministics are
    ## stored as one column of per-draw length-`nw` vectors.
    μA_draws = vec(Array(chn[:μ_A_pred]))
    recv_draws = vec(Array(chn[:recv_pred]))
    μA_med = [median(getindex.(μA_draws, w)) for w in 1:nw]
    recv_med = [median(getindex.(recv_draws, w)) for w in 1:nw]
    println(io, "=== ", label, " ===")
    println(io, "n confirmed vintages: ", length(c.offsets))
    println(io, "max R-hat: ", round(d.max_rhat; digits = 3))
    println(io, "min ESS bulk: ", round(d.min_ess_bulk; digits = 1))
    println(io, "min ESS tail: ", round(min_ess_tail; digits = 1))
    println(io, "divergences: ", d.n_divergent)
    println(io, "C_T median: ", round(median(C); digits = 0),
        "  90% interval: ", ival(Cs))
    println(io, "dark analysed total (median): ",
        round(median(dark); digits = 0))
    println(io, "offsets:        ", c.offsets)
    println(io, "analysed (obs): ", c.analysed)
    println(io, "pred μ_A (med):  ",
        round.(μA_med; digits = 1))
    println(io, "received (obs): ", c.received)
    println(io, "pred recv (med): ", round.(recv_med; digits = 1))
    println(io, "")
    flush(io)
end

function main()
    obs = load_observations()
    c = full_extension(obs)
    open(OUT, "w") do io
        println(io, "Condensed lab-throughput queue: ALL confirmed ",
            "vintages with Poisson-thinned dark denominators")
        println(io, "confirmed offsets: ", c.offsets)
        println(io, "confirmed increments: ", c.confirmed)
        println(io, "observed analysed anchors (23-28 May): ",
            filter(!ismissing, c.analysed))
        println(io, "observed received anchors (23-28 May): ",
            filter(!ismissing, c.received))
        println(io, "")
        flush(io)

        ## Headline: e = 0 (forward fraction 1).
        m0 = build_model(obs, c; epi_exclusion = nothing)
        chn0 = nuts_sample(m0;
            callback = tensorboard_callback("logs/queue/e0"),
            samples = SAMPLES, chains = CHAINS, seed = 1)
        summarise(io, "headline e=0", chn0, c, obs)

        ## Sensitivity: e ~ Beta(2, 12).
        me = build_model(obs, c; epi_exclusion = epi_exclusion_model())
        chne = nuts_sample(me;
            callback = tensorboard_callback("logs/queue/ebeta"),
            samples = SAMPLES, chains = CHAINS, seed = 1)
        summarise(io, "sensitivity e~Beta(2,12)", chne, c, obs)
    end
    println("Wrote ", OUT)
end

main()
