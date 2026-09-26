# Smoke fit of the headline joint: the screened prior starts each chain
# uses, then a short multi-chain NUTS run at the production sampler
# settings, reporting the sampler and mixing diagnostics that decide a
# CI fit. Run before pushing a model change:
#
#   SMOKE_CHAINS=4 SMOKE_WARMUP=150 SMOKE_DRAWS=150 \
#     julia --project=docs -t 4 scripts/smoke_joint.jl
using BVDOutbreakSize
using BVDOutbreakSize: parameter_diagnostics, fit_diagnostics
using DataFrames: DataFrame, first, sort!
using Printf: @sprintf
using Random: MersenneTwister
using Statistics: mean, median, std
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "fits", "cache.jl"))

say(args...) = (println(args...); flush(stdout))
warm = parse(Int, get(ENV, "SMOKE_WARMUP", "150"))
draws = parse(Int, get(ENV, "SMOKE_DRAWS", "150"))
nchains = parse(Int, get(ENV, "SMOKE_CHAINS", "4"))
max_depth = parse(Int, get(ENV, "SMOKE_MAX_DEPTH", "10"))

obs = load_observations()
say("onset cells: ", length(obs.onset_curve_history.increments))
model = BVDOutbreakSize.production_joint(
    obs; breakpoint = default_breakpoint(obs)
)

## The screened prior starts. A start tens of thousands of log units
## below the posterior is the signature of a chain that never leaves it.
rng = MersenneTwister(20260518)
ldf = BVDOutbreakSize.LogDensityFunction(model)
for c in 1:nchains
    init = BVDOutbreakSize.viable_prior_init(rng, model; ldf)
    vi = BVDOutbreakSize.VarInfo(rng, model, init)
    say(@sprintf("chain %d start logjoint %.1f", c, BVDOutbreakSize.getlogjoint(vi)))
end

t0 = time()
chn = nuts_sample(
    model; samples = draws, chains = nchains, n_adapts = warm,
    target_accept = 0.8, max_depth = max_depth, callback = nothing
)
say(@sprintf("sampling wall clock: %.1f min", (time() - t0) / 60))
chn = repair_chain_keys(chn)
FC = parentmodule(typeof(chn))
for c in 1:nchains
    ss = vec(chn[FC.Extra(:step_size)][:, c])
    td = vec(chn[FC.Extra(:tree_depth)][:, c])
    ns = vec(chn[FC.Extra(:n_steps)][:, c])
    ar = vec(chn[FC.Extra(:acceptance_rate)][:, c])
    ne = vec(chn[FC.Extra(:numerical_error)][:, c])
    lp = vec(chn[FC.Extra(:logjoint)][:, c])
    ct = vec(chn[:C_T][:, c])
    line = @sprintf(
        "chain %d: step %.5f  depth>=%d %.2f  steps %.0f  accept %.2f  div %d  lp %.0f±%.0f  C_T %.0f",
        c, median(ss), max_depth, mean(td .>= max_depth), mean(ns), mean(ar),
        sum(ne), median(lp), std(lp), median(ct)
    )
    say(line)
end
say("headline: ", fit_diagnostics(chn))
df = sort!(parameter_diagnostics(chn), :ess_bulk)
say("worst 12 by ESS bulk:")
show(stdout, first(df, 12); allrows = true, allcols = true)
println()
## Parameters that have split between chains in past fits.
for n in [
        "cases_state.report_state.α", "onset_to_sample_mean",
        "onset_report_state.inv_sqrt_k", "rt_state.sigma_rw",
        "patch_state.σ_drift",
    ]
    try
        v = chn[Symbol(n)]
        parts = [
            @sprintf(
                "c%d %.3g", c,
                median([e isa AbstractVector ? e[1] : e for e in vec(v[:, c])])
            ) for c in 1:nchains
        ]
        say(rpad(n, 30), join(parts, "  "))
    catch e
        say(n, ": ", sprint(showerror, e)[1:min(end, 60)])
    end
end
