# Prior-vs-posterior sense check for ALL sampled parameters of the joint.
# Loads the saved posterior chain, samples the SAME model under the prior,
# and tabulates prior vs posterior median + 90% interval per parameter,
# flagging which parameters the data actually update (identified) vs those
# that stay at the prior (unidentified / prior-dominated).
# Run: julia --project=. scripts/prior_vs_posterior.jl

using BVDOutbreakSize
using Turing: sample, Prior
using Random: MersenneTwister
using Statistics: median, quantile
using Serialization: deserialize

post = deserialize("logs/joint_chain.jls")
obs = load_observations()
const BP = obs.n - obs.who_first_sitrep_days

model = bvd_joint(
    obs.n, obs.exported_cases, obs.total_deaths,
    obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
    obs.tests_analysed;
    confirmed_deaths = obs.confirmed_deaths,
    deaths_history = obs.deaths_history,
    reported_history = obs.reported_history,
    confirmed_history = obs.confirmed_history,
    confirmed_deaths_history = obs.confirmed_deaths_history,
    lab_history = obs.lab_history,
    lab_daily_history = obs.lab_daily_history,
    tests_received_history = obs.tests_received_history,
    export_case_days = obs.export_case_days,
    export_death_days = obs.export_death_days,
    breakpoint = BP, background_re = true,
    genetic = genetic_seeding_model, tmrca_days = obs.tmrca_days)

prior = sample(MersenneTwister(20260605), model, Prior(), 2000; progress = false)

stats(v) = (median(v), quantile(v, 0.05), quantile(v, 0.95))
draws(chn, k) = vec(collect(chn[k]))

## Common scalar parameter keys (intersection), excluding vector/derived
## entries that don't reduce to a scalar.
postkeys = collect(keys(post))
priorkeys = Set(string.(keys(prior)))
rows = []
for k in postkeys
    string(k) in priorkeys || continue
    pv = draws(post, k)
    qv = draws(prior, k)
    (eltype(pv) <: Real && length(pv) > 0) || continue
    pm, plo, phi = stats(pv)
    qm, qlo, qhi = stats(qv)
    priorw = qhi - qlo
    postw = phi - plo
    shrink = priorw > 0 ? postw / priorw : NaN
    push!(rows,
        (name = replace(string(k), "Parameter(" => "", ")" => ""),
            prior = (qm, qlo, qhi), post = (pm, plo, phi), shrink = shrink))
end
sort!(rows, by = r -> r.shrink)

function fmt(t)
    string(round(t[1]; sigdigits = 3), " [", round(t[2]; sigdigits = 3),
        ", ", round(t[3]; sigdigits = 3), "]")
end
println(rpad("parameter", 34), rpad("prior med[90%]", 26),
    rpad("posterior med[90%]", 26), "shrink(post/prior width)")
println("-"^110)
for r in rows
    flag = isnan(r.shrink) ? "" :
           (r.shrink > 0.85 ? "  <- prior-dominated" :
            (r.shrink < 0.4 ? "  <- well identified" : ""))
    println(rpad(r.name, 34), rpad(fmt(r.prior), 26),
        rpad(fmt(r.post), 26), round(r.shrink; digits = 2), flag)
end
