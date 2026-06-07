# Testing-selection prototype: does a PERSISTENT selection floor δ∞ in the
# composition link let λ_bg float high and bring the joint size down toward
# the single-stream band (~2280), with acceptable mixing? δ∞=0 (tight floor
# prior) is the current composition link (baseline). One scenario per run via
# ARGS[1]: link (floor~0) | select (persistent floor). Appends to results.
using BVDOutbreakSize
using BVDOutbreakSize: bvd_joint, load_observations, exponential_growth_model,
    genetic_seeding_model, report_onset_offset, m_prior_centre,
    confirmed_q_re_model, severity_enrichment_model, nuts_sample,
    fit_diagnostics
using Statistics, Printf
using Distributions: truncated, Normal

const RESULTS = joinpath(@__DIR__, "selection_results.txt")
scenario = isempty(ARGS) ? "link" : ARGS[1]
floors = Dict(
    "link"   => ("COMPOSITION LINK floor~0 N(0,0.05)",
        truncated(Normal(0.0, 0.05); lower = 0)),
    "select" => ("PERSISTENT SELECTION floor N(1.0,0.5)",
        truncated(Normal(1.0, 0.5); lower = 0)),
    "wide"   => ("WIDE SELECTION floor N(0,1.0)",
        truncated(Normal(0.0, 1.0); lower = 0)))
label, floor_prior = floors[scenario]

logf(io, s) = (println(io, s); flush(io); println(s); flush(stdout))
q(xs) = (med = median(xs), lo = quantile(xs, 0.05), hi = quantile(xs, 0.95))
getp(chn, nm) = try; vec(Array(chn[Symbol(nm)])); catch; nothing; end

obs = load_observations()
growth_now = exponential_growth_model(
    m_prior = truncated(Normal(m_prior_centre(obs.as_of_date), 3.0); lower = 0))
genetic = T -> genetic_seeding_model(T, obs.genetic_tmrca_days;
    tmrca_days_sd = obs.genetic_tmrca_days_sd)
include(joinpath(@__DIR__, "_joint_obs.jl"))
fa = joint_obs(obs)
enrich = severity_enrichment_model(floor_prior = floor_prior)

chn = nuts_sample(bvd_joint(obs.exported_cases, fa.deaths, fa.reported,
        fa.export_deaths; fa.kw..., growth = growth_now,
        first_export_detection_delta = obs.first_export_detection_delta,
        report_onset_offset = report_onset_offset(obs.as_of_date),
        confirmed_q_random_effect = confirmed_q_re_model, genetic = genetic,
        confirmed_severity_enrichment = enrich);
    samples = 600, chains = 2, target_accept = 0.9, check_model = false)

d = fit_diagnostics(chn)
open(RESULTS, "a") do io
    logf(io, "")
    logf(io, "============ $label ============")
    logf(io, @sprintf("DIAG max_rhat=%.4f min_ess_bulk=%.0f div=%d",
        d.max_rhat, d.min_ess_bulk, d.n_divergent))
    for nm in ["cumulative_infections", "m", "T", "λ_bg", "δ∞", "δ0",
        "q_cutoff", "positivity"]
        xs = getp(chn, nm)
        if xs === nothing
            logf(io, @sprintf("PARAM %-24s NA", nm))
        else
            p = q(xs)
            logf(io, @sprintf("PARAM %-24s med=%.4f [%.4f, %.4f]",
                nm, p.med, p.lo, p.hi))
        end
    end
    logf(io, "DONE $scenario")
end
