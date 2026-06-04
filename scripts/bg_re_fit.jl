# Short NUTS sanity fit comparing the bg-random-effect joint to the
# renewal baseline. Run with `--baseline` to fit the scalar-λ_bg joint
# (the renewal default), or no argument to fit with the per-vintage
# background random effect enabled. Prints max R-hat, min ESS,
# divergences, and the C_T / λ_bg summaries.

using BVDOutbreakSize
using BVDOutbreakSize: bvd_joint, background_re_model,
    test_positivity_model, deaths_model, reported_cases_model,
    confirmed_deaths_model
using Turing: summarize
using FlexiChains: FlexiChain
using Statistics: mean, median, quantile
using Printf: @printf

const USE_RE = !("--baseline" in ARGS)
const SAMPLES = 300
const CHAINS = 2

obs = load_observations()
bp = obs.n - obs.who_first_sitrep_days

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
    tests_received_history = obs.tests_received_history,
    breakpoint = bp,
    genetic = nothing,
    tmrca_days = obs.tmrca_days,
    background_re = USE_RE)

@info "Fitting" use_re = USE_RE samples = SAMPLES chains = CHAINS

chn = nuts_sample(model; samples = SAMPLES, chains = CHAINS,
    progress = false)

## Diagnostics across all sampled parameters.
diag = fit_diagnostics(chn)
@info "diagnostics" max_rhat = diag.max_rhat min_ess = diag.min_ess_bulk
ndiv = diag.n_divergent

function q(sym)
    v = vec(Array(chn[sym]))
    (median(v), quantile(v, 0.05), quantile(v, 0.95))
end

println("=================================================")
println(USE_RE ? "RANDOM-EFFECT BACKGROUND" : "SCALAR BASELINE")
println("=================================================")
@printf "max R-hat : %.4f\n" diag.max_rhat
@printf "min ESS   : %.1f\n" diag.min_ess_bulk
println("divergences: ", ndiv)
for s in (:C_T, :lambda_bg, :T, :r, :CFR, :p_drc)
    try
        m, lo, hi = q(s)
        @printf "%-12s median %.4g  (90%% %.4g - %.4g)\n" String(s) m lo hi
    catch e
        println(s, " : n/a (", e, ")")
    end
end
if USE_RE
    for s in (:bg_sigma, :bg_death_sigma, :lambda_bg_death)
        try
            m, lo, hi = q(s)
            @printf "%-16s median %.4g  (90%% %.4g - %.4g)\n" String(s) m lo hi
        catch
        end
    end
end
