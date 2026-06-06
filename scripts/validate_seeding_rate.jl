# Short joint validation for the R0-derived cryptic growth rate cleanup.
# Composition-linked confirmed positivity, dated exports, SCALAR
# background, genetic seeding bound. 300 samples x 4 chains.
# Run: julia --project=test scripts/validate_seeding_rate.jl

using BVDOutbreakSize
using BVDOutbreakSize: bvd_joint, nuts_sample, genetic_seeding_model,
                       fit_diagnostics
using Statistics: median, quantile, std
using Serialization: serialize

obs = load_observations()
_BREAKPOINT = obs.n - obs.who_first_sitrep_days

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
    export_case_days = obs.export_case_days,
    export_death_days = obs.export_death_days,
    breakpoint = _BREAKPOINT,
    background_re = false,                  # SCALAR background
    confirmed_positivity_link = :composition,
    genetic = genetic_seeding_model,
    tmrca_days = obs.tmrca_days)

chn = nuts_sample(model; samples = 300, chains = 4, progress = false)
serialize("logs/validate_seeding_rate.jls", chn)

d = fit_diagnostics(chn)
draws(sym) = vec(collect(chn[Symbol(sym)]))
function q(sym)
    v = draws(sym)
    (median(v), quantile(v, 0.05), quantile(v, 0.95), std(v))
end

println("=== Convergence ===")
println("  max R-hat    : ", round(d.max_rhat; digits = 4))
println("  min bulk ESS : ", round(d.min_ess_bulk; digits = 1))
try
    println("  min tail ESS : ", round(d.min_ess_tail; digits = 1))
catch
end
println("  divergences  : ", d.n_divergent)

println("\n=== Headline (median [5%, 95%] SD) ===")
for sym in (:T, :C_T, :R0, :R_T, :r0, :r, :doubling_time_initial,
    :doubling_time, :expected_exports_T, :CFR, :lambda_bg)
    try
        m, lo, hi, s = q(sym)
        println("  ", rpad(string(sym), 22),
            rpad(round(m; digits = 3), 12),
            "[", round(lo; digits = 3), ", ", round(hi; digits = 3),
            "]  SD=", round(s; digits = 3))
    catch e
        println("  ", rpad(string(sym), 22), "(absent)")
    end
end
println("\n  observed exports = 3 (Uganda)")
