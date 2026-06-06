# Streamed full-joint fit for sanity-checking the renewal model after the
# shared-σ background + 75% ascertainment changes. Writes live progress
# (iteration, log-density, divergences) to logs/joint_fit.log via
# `progress_callback` — tail it with `tail -f logs/joint_fit.log` — and
# prints a posterior summary of the headline quantities, with C_T (the
# cumulative-infection outbreak size) front and centre.
#
# Run: JULIA_NUM_THREADS=4 julia --project=. scripts/fit_joint_stream.jl [samples] [chains]

using BVDOutbreakSize
using BVDOutbreakSize: fit_diagnostics
using Statistics: median, quantile
using Serialization: serialize

const SAMPLES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 500
const CHAINS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
const AD = length(ARGS) >= 3 ? ARGS[3] : "mooncake"
const PLINK = length(ARGS) >= 4 ? Symbol(ARGS[4]) : :free
const BG_RE = length(ARGS) >= 5 ? parse(Bool, ARGS[5]) : true

## Enzyme is opt-in (~3x faster than the Mooncake default on the joint).
## Its weak-dep extension registers only once `Enzyme` is loaded.
adtype = if AD == "enzyme"
    @eval using Enzyme
    BVDOutbreakSize.enzyme_adtype()
else
    BVDOutbreakSize.default_adtype()
end

const obs = load_observations()
const BREAKPOINT = obs.n - obs.who_first_sitrep_days

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
    breakpoint = BREAKPOINT,
    export_last_offset = obs.export_last_offset,
    background_re = BG_RE,
    confirmed_positivity_link = PLINK,
    genetic = genetic_seeding_model,
    tmrca_days = obs.tmrca_days)

isdir("logs") || mkdir("logs")
cb = progress_callback(; path = "logs/joint_fit.log", every = 25)

println("Fitting joint ($(SAMPLES)x$(CHAINS), n=$(obs.n), background_re=true, ",
    "ascertainment≈75%, AD=$(AD)). Tail logs/joint_fit.log for live progress.\n")
chn = nuts_sample(model; samples = SAMPLES, chains = CHAINS,
    adtype = adtype, callback = cb, progress = false)

serialize("logs/joint_chain.jls", chn)

d = fit_diagnostics(chn)

function q(sym)
    v = vec(collect(chn[sym]))
    return (median(v), quantile(v, 0.05), quantile(v, 0.95))
end

println("\n=== Convergence ===")
println("  divergences (of $(SAMPLES * CHAINS)) : ", d.n_divergent)
println("  max R-hat                  : ", round(d.max_rhat; digits = 4))
println("  min bulk ESS               : ", round(d.min_ess_bulk; digits = 1))

println("\n=== Headline posteriors  (median [5%, 95%]) ===")
for sym in (:C_T, :CFR, :r, :R_T, :T, :p_drc, :p_uganda,
    :lambda_bg, :bg_sigma, :lambda_bg_death,
    :suspected_positivity, :test_positivity,
    :expected_reports_T, :expected_deaths_T, :expected_exports_T,
    :expected_confirmed_T)
    try
        m, lo, hi = q(sym)
        println("  ", rpad(string(sym), 20),
            rpad(round(m; digits = 3), 12),
            "[", round(lo; digits = 3), ", ", round(hi; digits = 3), "]")
    catch
        println("  ", rpad(string(sym), 20), "(absent)")
    end
end

## Observed anchors for the sanity check.
println("\n=== Observed (cut-off) for comparison ===")
println("  suspected cases  : ", obs.reported_cases)
println("  suspected deaths : ", obs.total_deaths)
println("  confirmed cases  : ", obs.confirmed_cases)
println("  exports          : ", obs.exported_cases)
