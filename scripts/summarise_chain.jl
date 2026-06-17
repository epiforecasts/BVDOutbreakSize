# Summarise a saved joint chain (logs/joint_chain.jls) without re-fitting,
# using the package's FlexiChains-aware helpers.
# Run: julia --project=. scripts/summarise_chain.jl [path]

using BVDOutbreakSize
using BVDOutbreakSize: fit_diagnostics, posterior_summary
using Serialization: deserialize
using Statistics: median, quantile

path = length(ARGS) >= 1 ? ARGS[1] : "logs/joint_chain.jls"
chn = deserialize(path)
obs = load_observations()

d = fit_diagnostics(chn)
println("=== Convergence ===")
println("  max R-hat      : ", round(d.max_rhat; digits = 4))
println("  min bulk ESS   : ", round(d.min_ess_bulk; digits = 1))
println("  divergences    : ", d.n_divergent)

draws(sym) = vec(collect(chn[Symbol(sym)]))

println("\n=== Headline posteriors (median [5%, 95%]) ===")
for sym in (:C_T, :CFR, :r, :R_T, :T, :p_drc, :p_uganda, :lambda_bg,
    :bg_sigma, :death_ascertainment, :background_cfr, :lambda_bg_death,
    :tau_death, :death_composition, :death_confirmation,
    :background_total, :background_death_total,
    :expected_reports_T, :expected_deaths_T, :expected_exports_T,
    :expected_confirmed_T, :expected_confirmed_deaths_T)
    try
        v = draws(sym)
        println("  ", rpad(string(sym), 22),
            rpad(round(median(v); digits = 3), 12),
            "[", round(quantile(v, 0.05); digits = 3), ", ",
            round(quantile(v, 0.95); digits = 3), "]")
    catch e
        println("  ", rpad(string(sym), 22), "(absent: ", e, ")")
    end
end

println("\n=== Observed (cut-off) ===")
println("  suspected cases  : ", obs.reported_cases)
println("  suspected deaths : ", obs.total_deaths)
println("  confirmed cases  : ", obs.confirmed_cases)
println("  exports          : ", obs.exported_cases)
