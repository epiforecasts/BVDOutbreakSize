# Post-fit report for the explicit-T Gibbs joint: loads the serialised
# chain, prints the T posterior + convergence diagnostics, and writes a
# start-date pair plot. Run after scripts/fit_joint_gibbs_start.jl.
#
# Run: julia --project=. scripts/report_gibbs_start.jl

using BVDOutbreakSize
using BVDOutbreakSize: fit_diagnostics, plot_start_date_pair
using Statistics: median, quantile, std, mean
using Serialization: deserialize
import CairoMakie

const chn = deserialize("logs/gibbs_start_chain.jls")
const obs = load_observations()

d = fit_diagnostics(chn)
T = vec(collect(chn[:T]))

println("=== Gibbs explicit-T fit ===")
println("  draws            : ", length(T))
println("  divergences      : ", d.n_divergent)
println("  max R-hat        : ", round(d.max_rhat; digits = 4))
println("  min bulk ESS     : ", round(d.min_ess_bulk; digits = 1))
println()
println("=== T (outbreak age, days) ===")
println("  mean   ", round(mean(T); digits = 1))
println("  median ", round(median(T); digits = 1))
println("  SD     ", round(std(T); digits = 1))
println("  90% CI [", round(quantile(T, 0.05); digits = 1), ", ",
    round(quantile(T, 0.95); digits = 1), "]")
println("  range  [", round(minimum(T); digits = 1), ", ",
    round(maximum(T); digits = 1), "]")
println("  (tmrca_days = ", obs.tmrca_days, ", grid n = ", obs.n, ")")

## Per-chain T means to gauge cross-chain agreement (no second mode).
try
    arr = Array(chn[:T])           # iterations × chains
    if ndims(arr) == 2
        println("\n  per-chain T means: ",
            [round(mean(arr[:, c]); digits = 1) for c in 1:size(arr, 2)])
    end
catch
end

isdir("figures") || mkdir("figures")
fig = plot_start_date_pair(chn; as_of_date = string(obs.cutoff))
CairoMakie.save("figures/gibbs_start_T.png", fig)
println("\nWrote figures/gibbs_start_T.png")
