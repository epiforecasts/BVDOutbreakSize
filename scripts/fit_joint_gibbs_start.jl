# Gibbs joint fit with an EXPLICIT, sampled outbreak start time `T`.
# Samples the outbreak age `T` in its own Gibbs block (HMC by default)
# and NUTS on every other continuous parameter, replacing the renewal
# model's post-hoc `seeding_age` crossing (which produced a two-stage
# spike and a tiny-outbreak multimodal collapse). Writes live progress to
# logs/gibbs_start_fit.log and prints a posterior summary with the `T`
# posterior front and centre.
#
# Run: JULIA_NUM_THREADS=4 julia --project=. \
#        scripts/fit_joint_gibbs_start.jl [samples] [chains] [plink]

using BVDOutbreakSize
using BVDOutbreakSize: fit_diagnostics, genetic_seeding_model
using Statistics: median, quantile, std
using Serialization: serialize

const SAMPLES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 400
const CHAINS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
const PLINK = length(ARGS) >= 3 ? Symbol(ARGS[3]) : :composition
const BG_RE = length(ARGS) >= 4 ? parse(Bool, ARGS[4]) : true
const TAG = length(ARGS) >= 5 ? ARGS[5] : "start"

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
    background_re = BG_RE,
    confirmed_positivity_link = PLINK,
    fit_start = true,
    genetic = genetic_seeding_model,
    tmrca_days = obs.tmrca_days)

isdir("logs") || mkdir("logs")
const LOGFILE = "logs/gibbs_$(TAG)_fit.log"
const CHAINFILE = "logs/gibbs_$(TAG)_chain.jls"
cb = progress_callback(; path = LOGFILE, every = 10)

println("Gibbs joint w/ explicit T ($(SAMPLES)x$(CHAINS), n=$(obs.n), ",
    "tmrca_days=$(obs.tmrca_days), plink=$(PLINK)). ",
    "Tail logs/gibbs_start_fit.log for progress.\n")

chn = gibbs_sample(model; samples = SAMPLES, chains = CHAINS,
    callback = cb, progress = false)

serialize(CHAINFILE, chn)

d = fit_diagnostics(chn)

function q(sym)
    v = vec(collect(chn[sym]))
    return (median(v), quantile(v, 0.05), quantile(v, 0.95), std(v))
end

println("\n=== Convergence ===")
println("  divergences (of $(SAMPLES * CHAINS)) : ", d.n_divergent)
println("  max R-hat                  : ", round(d.max_rhat; digits = 4))
println("  min bulk ESS               : ", round(d.min_ess_bulk; digits = 1))

println("\n=== T (outbreak age) posterior ===")
try
    m, lo, hi, sd = q(:T)
    println("  T median ", round(m; digits = 1),
        " d  [5%, 95%] = [", round(lo; digits = 1), ", ",
        round(hi; digits = 1), "]  SD ", round(sd; digits = 1))
    println("  (tmrca_days = ", obs.tmrca_days, ", grid n = ", obs.n, ")")
catch e
    println("  (T absent: ", e, ")")
end

println("\n=== Headline posteriors  (median [5%, 95%]) ===")
for sym in (:C_T, :CFR, :r, :R_T, :T, :p_drc, :p_uganda,
    :lambda_bg, :lambda_bg_death,
    :suspected_positivity, :test_positivity,
    :expected_reports_T, :expected_deaths_T, :expected_exports_T,
    :expected_confirmed_T)
    try
        m, lo, hi, _ = q(sym)
        println("  ", rpad(string(sym), 20),
            rpad(round(m; digits = 3), 12),
            "[", round(lo; digits = 3), ", ", round(hi; digits = 3), "]")
    catch
        println("  ", rpad(string(sym), 20), "(absent)")
    end
end

println("\n=== Observed (cut-off) for comparison ===")
println("  suspected cases  : ", obs.reported_cases)
println("  suspected deaths : ", obs.total_deaths)
println("  confirmed cases  : ", obs.confirmed_cases)
println("  exports          : ", obs.exported_cases)
