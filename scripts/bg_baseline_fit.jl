# Renewal-baseline short NUTS fit of the full joint (scalar λ_bg, no
# background RE). Captures the comparison numbers before the random
# effect is wired in. Prints max R-hat, min ESS, divergences and the
# C_T / λ_bg summaries.

using BVDOutbreakSize
using BVDOutbreakSize: bvd_joint
using Statistics: median, quantile
using Printf: @printf

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
    tmrca_days = obs.tmrca_days)

@info "Fitting baseline" samples=SAMPLES chains=CHAINS
chn = nuts_sample(model; samples = SAMPLES, chains = CHAINS,
    progress = false)

diag = fit_diagnostics(chn)

function q(sym)
    v = vec(Array(chn[sym]))
    (median(v), quantile(v, 0.05), quantile(v, 0.95))
end

println("=================================================")
println("SCALAR BASELINE (renewal)")
println("=================================================")
@printf "max R-hat : %.4f\n" diag.max_rhat
@printf "min ESS   : %.1f\n" diag.min_ess_bulk
println("divergences: ", diag.n_divergent)
for s in (:C_T, :lambda_bg, :T, :r, :CFR, :p_drc)
    try
        m, lo, hi = q(s)
        @printf "%-12s median %.4g  (90%% %.4g - %.4g)\n" String(s) m lo hi
    catch e
        println(s, " : n/a")
    end
end
