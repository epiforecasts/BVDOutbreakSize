# Capacity-driven analysed-denominator imputation for the missing lab
# vintages. Fits the joint model over ALL confirmed-case vintages
# (18-28 May), with the early 18-22 May windows — which carry no
# published national analysed total — given a MECHANISTIC denominator
# derived from the lab capacity-times-backlog throughput the model
# already runs for the observed 23-28 May windows
# (`confirmed_capacity_impute = true`). See
# `confirmed_cases_model` / `bvd_joint` for the mechanism.
#
# Run under the test project so TensorBoardLogger (a weak dep) resolves:
#   julia --project=test scripts/missing_lab_capacity.jl
#
# Streams the live log-density / step-size / divergence traces to
# `logs/capacity` for TensorBoard, and writes a small diagnostics report
# to `scripts/out/capacity_result.txt`.

using BVDOutbreakSize
using BVDOutbreakSize: bvd_joint, load_observations, nuts_sample,
                       confirmed_q_re_model, fit_diagnostics,
                       posterior_summary
using TensorBoardLogger
using Statistics: median
import FlexiChains

const SAMPLES = get(ENV, "CAP_SAMPLES", "400")
const CHAINS = get(ENV, "CAP_CHAINS", "4")
const SEED = 20260518

_inc(v) = begin
    out = similar(v, Int)
    p = 0
    for i in eachindex(v)
        out[i] = v[i] - p
        p = v[i]
    end
    out
end

"""
Build full-window confirmed vectors from `load_observations()`: every
confirmed vintage (18-28 May) enters, the analysed / received
denominators present only where the national lab series publishes them
(23-28 May), `missing` otherwise. Mirrors `early_extension` in
`test/test_confirmed_impute.jl` but keeps the whole confirmed series.
"""
function full_window_confirmed(obs)
    ch = obs.confirmed_case_history
    sa = obs.tests_analysed_history
    sr = obs.tests_received_history
    offs = collect(ch.offsets)
    ccum = ch.values
    sa_inc = _inc(sa.values)
    sr_inc = _inc(sr.values)
    analysed = Union{Missing, Int}[]
    received = Union{Missing, Int}[]
    for off in offs
        k = findfirst(==(off), sa.offsets)
        if k === nothing
            push!(analysed, missing)
            push!(received, missing)
        else
            push!(analysed, sa_inc[k])
            push!(received, sr_inc[k])
        end
    end
    return (offsets = offs, confirmed = _inc(ccum),
        analysed = analysed, received = received,
        confirmed_cum = ccum)
end

function build_model(obs, c)
    rep = obs.reported_case_history
    dh = obs.death_history
    return bvd_joint(obs.exported_cases, _inc(dh.values),
        _inc(rep.values), obs.export_deaths_daily;
        reported_offsets = rep.offsets,
        death_offsets = dh.offsets,
        confirmed_cases = c.confirmed,
        confirmed_offsets = c.offsets,
        samples_analysed = c.analysed,
        samples_received = c.received,
        tests_analysed = obs.cumulative_tests_analysed,
        tests_offset = 0,
        first_export_detection_delta = obs.first_export_detection_delta,
        report_onset_offset = report_onset_offset(obs.as_of_date),
        confirmed_q_random_effect = confirmed_q_re_model,
        confirmed_capacity_impute = true)
end

function main()
    obs = load_observations()
    c = full_window_confirmed(obs)
    n_missing = count(ismissing, c.analysed)
    n_total = length(c.confirmed)
    n_obs = n_total - n_missing
    @info "confirmed vintages" total=n_total missing=n_missing obs=n_obs

    model = build_model(obs, c)

    samples = parse(Int, SAMPLES)
    chains = parse(Int, CHAINS)
    isdir("logs") || mkpath("logs")
    cb = tensorboard_callback("logs/capacity"; every = 10)
    @info "fitting" samples chains seed = SEED
    chn = nuts_sample(model; samples = samples, chains = chains,
        seed = SEED, warmup = true, callback = cb, progress = false)

    diag = fit_diagnostics(chn)

    ## Posterior C(T) (cumulative cases) and imputed denominators.
    CT = vec(Array(chn[:cumulative_cases]))
    ct_s = posterior_summary(CT)

    ## The model exposes the per-edge analysed denominator only through
    ## the throughput; reconstruct the posterior imputed ΔA for the
    ## missing windows via the generated quantity if present, else report
    ## the capacity at the cut-off as a proxy.
    report = IOBuffer()
    println(report, "Capacity-driven analysed-denominator imputation")
    println(report, "=" ^ 50)
    println(report, "samples=$samples chains=$chains seed=$SEED")
    println(report, "confirmed vintages fitted = ",
        length(c.confirmed),
        " (missing-denominator = ", n_missing, ")")
    println(report)
    println(report, "Diagnostics")
    println(report, "  max R-hat      = ", round(diag.max_rhat; digits = 4))
    println(report, "  min ESS (bulk) = ",
        round(diag.min_ess_bulk; digits = 1))
    println(report, "  divergences    = ", diag.n_divergent)
    println(report)
    println(report, "cumulative_cases C(T)")
    println(report, "  median   = ", round(median(CT); digits = 1))
    println(report, "  90% CrI  = [", round(ct_s.lo90; digits = 1), ", ",
        round(ct_s.hi90; digits = 1), "]")
    println(report)
    println(report, "Observed analysed increments (23-28 May):")
    sa_inc = _inc(obs.tests_analysed_history.values)
    println(report, "  offsets ", collect(obs.tests_analysed_history.offsets),
        " ΔA = ", sa_inc)
    println(report)
    miss_idx = [i for i in eachindex(c.offsets) if ismissing(c.analysed[i])]
    println(report, "Missing-denominator vintages (offsets): ",
        [c.offsets[i] for i in miss_idx])
    println(report)

    ## Derived `:=` quantities are not bare Symbols in
    ## `FlexiChains.parameters`, so fetch them by direct indexing and skip
    ## gracefully if absent.
    _try_get(sym) =
        try
            Array(chn[sym])
        catch
            nothing
        end

    ## Per-edge imputed/observed analysed denominators (median over draws).
    ## `analysed_denominator` is a length-n vector quantity: hcat over draws
    ## gives an n_edges × n_draws matrix.
    AD = _try_get(:analysed_denominator)
    if AD !== nothing
        admat = reduce(hcat, vec(AD))
        println(report, "Imputed/observed analysed denominator ΔA ",
            "(median, 90% CrI) per vintage:")
        for i in eachindex(c.offsets)
            col = admat[i, :]
            s = posterior_summary(col)
            tag = ismissing(c.analysed[i]) ? "imputed " : "observed="
            ov = ismissing(c.analysed[i]) ? "" : string(c.analysed[i])
            println(report, "  offset ", c.offsets[i], "  ", tag, ov,
                "  median ", round(median(col); digits = 1), "  [",
                round(s.lo90; digits = 1), ", ",
                round(s.hi90; digits = 1), "]")
        end
        println(report)
    end

    ## Received posterior-predictive vs observed (the 431→662 surge, #182).
    RM_raw = _try_get(:received_mean)
    if RM_raw !== nothing
        RM = reduce(hcat, vec(RM_raw))  # n × draws
        sr_inc = _inc(obs.tests_received_history.values)
        println(report, "Received PP mean ΔR vs observed (offsets with ",
            "data):")
        for (k, off) in enumerate(obs.tests_received_history.offsets)
            j = findfirst(==(off), c.offsets)
            j === nothing && continue
            col = RM[j, :]
            s = posterior_summary(col)
            println(report, "  offset ", off, "  observed ΔR=", sr_inc[k],
                "  PP median ", round(median(col); digits = 1), "  [",
                round(s.lo90; digits = 1), ", ",
                round(s.hi90; digits = 1), "]")
        end
        println(report)
    end

    cap_raw = _try_get(:capacity_cutoff)
    if cap_raw !== nothing
        cap = vec(cap_raw)
        cs = posterior_summary(cap)
        println(report, "capacity at cut-off (samples/day): median ",
            round(median(cap); digits = 1), " 90% [",
            round(cs.lo90; digits = 1), ", ",
            round(cs.hi90; digits = 1), "]")
    end
    tf_raw = _try_get(Symbol("τ_forward_out"))
    if tf_raw !== nothing
        tf = vec(tf_raw)
        ts = posterior_summary(tf)
        println(report, "τ_forward: median ", round(median(tf); digits = 3),
            " 90% [", round(ts.lo90; digits = 3), ", ",
            round(ts.hi90; digits = 3), "]")
    end

    txt = String(take!(report))
    print(txt)
    open("scripts/out/capacity_result.txt", "w") do io
        print(io, txt)
    end
    @info "wrote scripts/out/capacity_result.txt"
    return chn
end

main()
