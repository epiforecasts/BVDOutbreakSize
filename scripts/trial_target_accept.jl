# Trial the lowered NUTS target_accept on the real-data renewal joint.
#
# Fits the full `bvd_joint` (the analysis configuration) at target_accept
# 0.90 (the new default) and 0.95 (the previous value) with a reduced
# sample budget, and compares divergences, convergence (max R-hat, min bulk
# ESS) and the headline C_T. A lower target acceptance should cut leapfrog
# steps without reintroducing divergences if the geometry is benign.
#
# Run: julia --project=. scripts/trial_target_accept.jl

using BVDOutbreakSize
using Turing: summarystats
using Statistics: median, quantile

const obs = load_observations()
const BREAKPOINT = obs.n - obs.who_first_sitrep_days

function build_joint()
    bvd_joint(
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
        breakpoint = BREAKPOINT,
        genetic = genetic_seeding_model,
        tmrca_days = obs.tmrca_days)
end

ndiv(chn) =
    try
        Int(sum(Array(chn[:numerical_error])))
    catch
        -1
    end

function convergence(chn)
    nt = summarystats(chn).nt
    col(name) = hasproperty(nt, name) ?
                filter(isfinite, collect(getproperty(nt, name))) : Float64[]
    rhat = col(:rhat)
    ess = isempty(col(:ess_bulk)) ? col(:ess) : col(:ess_bulk)
    return (maxrhat = isempty(rhat) ? NaN : maximum(rhat),
        miness = isempty(ess) ? NaN : minimum(ess))
end

function run_one(ta; samples, chains)
    chn = nuts_sample(build_joint(); samples, chains, target_accept = ta,
        progress = false)
    c = convergence(chn)
    CT = vec(Array(chn[:C_T]))
    println("--- target_accept = $ta ---")
    println("  divergences (of $(samples * chains)) : ", ndiv(chn))
    println("  max R-hat                  : ", round(c.maxrhat; digits = 4))
    println("  min ESS                    : ", round(c.miness; digits = 1))
    println("  C_T median [5%, 95%]       : ",
        round(median(CT); digits = 0), " [",
        round(quantile(CT, 0.05); digits = 0), ", ",
        round(quantile(CT, 0.95); digits = 0), "]")
    return nothing
end

const SAMPLES = 500
const CHAINS = 2
println("Joint target_accept trial: $(SAMPLES)x$(CHAINS), n=$(obs.n)\n")
run_one(0.90; samples = SAMPLES, chains = CHAINS)
run_one(0.95; samples = SAMPLES, chains = CHAINS)
