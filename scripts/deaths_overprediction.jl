# Investigation: confirmed-deaths over-prediction in the joint, and whether
# a wider death-background prior (λ_bg_death) fixes the level without breaking
# mixing. Relaxed settings: 600 samples, 2 chains, target_accept = 0.9.
#
# ONE scenario per invocation (results buffer otherwise and are lost on a
# long crawl). Pass the scenario key as ARGS[1]: baseline | wide | centred.
# Each run appends a result block to scripts/deaths_results.txt and flushes.
using BVDOutbreakSize
using BVDOutbreakSize: bvd_joint, load_observations,
    exponential_growth_model, genetic_seeding_model,
    report_onset_offset, m_prior_centre, confirmed_q_re_model,
    death_background_model, nuts_sample, fit_diagnostics
using Statistics, Printf
using Distributions: truncated, Normal
using Turing: predict, @varname

const NS = 600
const NC = 2
const TA = 0.9
const RESULTS = joinpath(@__DIR__, "deaths_results.txt")

scenario = isempty(ARGS) ? "baseline" : ARGS[1]
death_priors = Dict(
    "baseline" => ("BASELINE Normal(0,0.25)",
        truncated(Normal(0.0, 0.25); lower = 0)),
    "wide" => ("WIDER half-normal Normal(0,1.0)",
        truncated(Normal(0.0, 1.0); lower = 0)),
    "centred" => ("CENTRED half-normal Normal(1.24,0.6)",
        truncated(Normal(1.24, 0.6); lower = 0)))
label, lambda_prior = death_priors[scenario]

logf(io, s) = (println(io, s); flush(io); println(s); flush(stdout))

sample_fit(model) = nuts_sample(model; samples = NS, chains = NC,
    target_accept = TA, check_model = false)

q(xs) = (med = median(xs), lo = quantile(xs, 0.05), hi = quantile(xs, 0.95))

function getparam(chn, name)
    try
        return vec(Array(chn[Symbol(name)]))
    catch
        return nothing
    end
end

obs = load_observations()
growth_now = exponential_growth_model(
    m_prior = truncated(Normal(m_prior_centre(obs.as_of_date), 3.0);
        lower = 0))
genetic_seeding = T -> genetic_seeding_model(T, obs.genetic_tmrca_days;
    tmrca_days_sd = obs.genetic_tmrca_days_sd)

include(joinpath(@__DIR__, "_joint_obs.jl"))
fit_args = joint_obs(obs)
conf_deaths_obs = collect(skipmissing(fit_args.kw.confirmed_deaths))

death_bg = death_background_model(; lambda_prior = lambda_prior)

build_joint(dbg) = bvd_joint(obs.exported_cases, fit_args.deaths,
    fit_args.reported, fit_args.export_deaths; fit_args.kw...,
    growth = growth_now,
    first_export_detection_delta = obs.first_export_detection_delta,
    report_onset_offset = report_onset_offset(obs.as_of_date),
    confirmed_q_random_effect = confirmed_q_re_model,
    genetic = genetic_seeding,
    death_background = dbg)

function pp_confirmed_deaths(dbg, chn)
    kw = merge(fit_args.kw,
        (; confirmed_deaths = fill(missing,
            length(fit_args.kw.confirmed_deaths))))
    genmodel = bvd_joint(obs.exported_cases, fit_args.deaths,
        fit_args.reported, fit_args.export_deaths; kw...,
        growth = growth_now,
        first_export_detection_delta = obs.first_export_detection_delta,
        report_onset_offset = report_onset_offset(obs.as_of_date),
        confirmed_q_random_effect = confirmed_q_re_model,
        genetic = genetic_seeding,
        death_background = dbg)
    pred = predict(genmodel, chn)
    totals = vec(sum.(pred[@varname(confirmed_deaths)]))
    return q(totals)
end

open(RESULTS, "a") do io
    logf(io, "")
    logf(io, "================ $label ================")
    logf(io, "scenario_key=$scenario  as_of=$(obs.as_of_date)")
    logf(io, "suspected_deaths_T=$(obs.total_deaths)  " *
            "confirmed_deaths_T=$(obs.confirmed_death_history.values[end])  " *
            "conf_rate=" *
            "$(round(obs.confirmed_death_history.values[end] / obs.total_deaths; digits=3))")
    logf(io, "confirmed_death_incrs=$(fit_args.kw.confirmed_deaths) " *
            "sum=$(sum(conf_deaths_obs))")

    chn = sample_fit(build_joint(death_bg))
    d = fit_diagnostics(chn)
    logf(io, @sprintf("DIAG max_rhat=%.4f min_ess_bulk=%.0f div=%d",
        d.max_rhat, d.min_ess_bulk, d.n_divergent))
    for nm in ["CFR", "p_deaths", "λ_bg_death", "q_death_cutoff",
        "p_pos_death_cutoff", "τ_death", "expected_confirmed_deaths_total",
        "expected_deaths_T", "cumulative_infections"]
        xs = getparam(chn, nm)
        if xs === nothing
            logf(io, @sprintf("PARAM %-32s NA", nm))
        else
            p = q(xs)
            logf(io, @sprintf("PARAM %-32s med=%.4f [%.4f, %.4f]",
                nm, p.med, p.lo, p.hi))
        end
    end
    pp = pp_confirmed_deaths(death_bg, chn)
    logf(io, @sprintf("PP_CONF_DEATHS med=%.1f [%.1f, %.1f] (obs=%d)",
        pp.med, pp.lo, pp.hi, sum(conf_deaths_obs)))
    logf(io, "DONE $scenario")
end
