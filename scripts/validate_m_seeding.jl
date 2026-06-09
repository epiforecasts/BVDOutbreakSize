# Short-joint validation for the m-induced two-phase renewal seeding.
# Composition positivity link (default) + dated exports (already in base),
# scalar background (background_re = false), genetic bound on. Reports the
# T width, per-chain C_T basins, divergences, R-hat and ESS, reading the
# numbers directly so they can be relayed without re-derivation.
#
# Run: JULIA_NUM_THREADS=4 julia --project=. scripts/validate_m_seeding.jl \
#        [samples] [chains]

using BVDOutbreakSize
using BVDOutbreakSize: fit_diagnostics
using Statistics: median, quantile, mean, std
using Serialization: serialize

const SAMPLES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 300
const CHAINS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

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
    lab_daily_history = obs.lab_daily_history,
    tests_received_history = obs.tests_received_history,
    export_case_days = obs.export_case_days,
    export_death_days = obs.export_death_days,
    breakpoint = BREAKPOINT,
    background_re = false,
    confirmed_positivity_link = :composition,
    genetic = genetic_seeding_model,
    tmrca_days = obs.tmrca_days)

isdir("logs") || mkdir("logs")
cb = progress_callback(; path = "logs/validate_m_seeding.log", every = 25)

println("Fitting joint ($(SAMPLES)x$(CHAINS), n=$(obs.n), ",
    "composition link, scalar background, genetic bound on, ",
    "tmrca_days=$(obs.tmrca_days)).\n")
chn = nuts_sample(model; samples = SAMPLES, chains = CHAINS,
    callback = cb, progress = false)

serialize("logs/validate_m_seeding_chain.jls", chn)

d = fit_diagnostics(chn)

function q(sym)
    v = vec(collect(chn[sym]))
    return (median(v), quantile(v, 0.05), quantile(v, 0.95))
end

println("\n=== Convergence ===")
println("  divergences (of $(SAMPLES * CHAINS)) : ", d.n_divergent)
println("  max R-hat                  : ", round(d.max_rhat; digits = 4))
println("  min bulk ESS               : ", round(d.min_ess_bulk; digits = 1))

## R0 = the SINGLE established reproduction number (= the first R_t), the
## one growth source for both the cryptic phase (via Euler–Lotka) and the
## renewal. THE headline failure signature of the earlier R0-rate attempt
## was R0 collapsing toward 1.0 and going bimodal across chains, because the
## seed back-scaling put R0 in both the seed and the renewal growth. With
## the 2^m seed magnitude r-free, R0 should stay near its prior (≈1.6),
## single-basin across chains.
R0v = vec(collect(chn[:R0]))
println("\n=== R0 (single established reproduction number) ===")
println("  median ", round(median(R0v); digits = 4),
    "  90% [", round(quantile(R0v, 0.05); digits = 4), ", ",
    round(quantile(R0v, 0.95); digits = 4), "]",
    "  SD ", round(std(R0v); digits = 4))
println("  quantiles 1/10/25/50/75/90/99 % : ",
    join(
        round.(quantile(R0v,
                [0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99]); digits = 3), "  "))
println("  fraction of draws below 1.1 (collapse check) : ",
    round(mean(R0v .< 1.1); digits = 4))
println("  per-chain R0 medians (bimodality check):")
R0all = collect(chn[:R0])
for c in 1:CHAINS
    rows = ((c - 1) * SAMPLES + 1):(c * SAMPLES)
    R0c = vec(R0all)[rows]
    println("    chain ", c, " : median ", round(median(R0c); digits = 4),
        "  mean ", round(mean(R0c); digits = 4),
        "  [", round(quantile(R0c, 0.05); digits = 3), ", ",
        round(quantile(R0c, 0.95); digits = 3), "]")
end

## T width: median / 90% interval / SD. The headline question is whether
## the wide m prior leaves T genuinely uncertain (vs the SD-3 tightness of
## the fixed-anchor placement in PR #216).
Tv = vec(collect(chn[:T]))
println("\n=== T (outbreak age, days) ===")
println("  median ", round(median(Tv); digits = 2),
    "  90% [", round(quantile(Tv, 0.05); digits = 2), ", ",
    round(quantile(Tv, 0.95); digits = 2), "]",
    "  SD ", round(std(Tv); digits = 2))
println("  tmrca_days (genetic floor) = ", obs.tmrca_days)

## Per-chain basins for T and C_T: one outbreak-size basin or several?
nper = SAMPLES
println("\n=== Per-chain means (basin check) ===")
println("  chain :  T mean   C_T mean   C_T median")
Tall = collect(chn[:T])
Call = collect(chn[:C_T])
for c in 1:CHAINS
    rows = ((c - 1) * nper + 1):(c * nper)
    Tc = vec(Tall)[rows]
    Cc = vec(Call)[rows]
    println("    ", c, "   : ", rpad(round(mean(Tc); digits = 2), 8), " ",
        rpad(round(mean(Cc); digits = 1), 10), " ",
        round(median(Cc); digits = 1))
end

println("\n=== Headline posteriors  (median [5%, 95%]) ===")
for sym in (:C_T, :R0, :r0, :CFR, :r, :R_T, :T, :p_drc, :p_uganda,
    :lambda_bg, :lambda_bg_death,
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

println("\n=== Observed (cut-off) for comparison ===")
println("  suspected cases  : ", obs.reported_cases)
println("  suspected deaths : ", obs.total_deaths)
println("  confirmed cases  : ", obs.confirmed_cases)
println("  exports          : ", obs.exported_cases)
