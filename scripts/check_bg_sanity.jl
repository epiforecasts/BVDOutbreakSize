# Scratch: does the widened SD-5 background over-predict? Short joint fit,
# then print C_T, the total background over the grid vs the observed
# suspected cases, λ_bg, suspected positivity and divergences. The
# over-prediction guard is: background_total should stay a minority of the
# observed suspected total (reported_cases), with the composition link
# identifying λ_bg rather than the prior.

using BVDOutbreakSize
using Statistics: quantile, median
using FlexiChains: FlexiChains

obs = load_observations()
const BP = obs.n - obs.who_first_sitrep_days

model = bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
    obs.reported_cases, obs.exports_deaths, obs.confirmed_cases,
    obs.tests_analysed; confirmed_deaths = obs.confirmed_deaths,
    deaths_history = obs.deaths_history,
    reported_history = obs.reported_history,
    confirmed_history = obs.confirmed_history,
    confirmed_deaths_history = obs.confirmed_deaths_history,
    lab_history = obs.lab_history,
    lab_daily_history = obs.lab_daily_history,
    tests_received_history = obs.tests_received_history,
    breakpoint = BP, background_re = true,
    confirmed_positivity_link = :composition,
    export_case_days = obs.export_case_days,
    export_death_days = obs.export_death_days,
    genetic = genetic_seeding_model, tmrca_days = obs.tmrca_days)

chn = nuts_sample(model; samples = 400, chains = 2, progress = false)

function q(v)
    (round(median(v); digits = 1), round(quantile(v, 0.05); digits = 1),
        round(quantile(v, 0.95); digits = 1))
end
draws(sym) = vec(Array(chn[sym]))

CT = draws(:C_T)
bg = draws(:background_total)
lbg = draws(:lambda_bg)
pos = draws(:suspected_positivity)

println("observed suspected (reported_cases) = ", obs.reported_cases)
println("T (outbreak age)  median/90% = ", q(draws(:T)),
    "  (genetic floor tmrca_days = ", obs.tmrca_days, ")")
println("C_T               median/90% = ", q(CT))
println("background_total  median/90% = ", q(bg))
println("  background_total / observed (median) = ",
    round(median(bg) / obs.reported_cases; digits = 3))
println("  P(background_total > observed)       = ",
    round(sum(bg .> obs.reported_cases) / length(bg); digits = 3))
println("lambda_bg (per day) median/90% = ", q(lbg))
println("suspected_positivity median/90% = ", q(pos))
println("test_positivity      median/90% = ", q(draws(:test_positivity)))
println("expected_exports_T   median/90% = ", q(draws(:expected_exports_T)),
    "  (observed ", obs.exported_cases, ")")
println("expected_exports_deaths_T med/90% = ",
    q(draws(:expected_exports_deaths_T)), "  (observed ", obs.exports_deaths, ")")
println("R_T                  median/90% = ", q(draws(:R_T)))
let d = BVDOutbreakSize.fit_diagnostics(chn)
    println("max R-hat = ", round(d.max_rhat; digits = 3),
        "   min bulk ESS = ", round(d.min_ess_bulk; digits = 1))
end
