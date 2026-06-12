# Diagnose the systematic first-bin overestimate across the per-vintage
# streams. Fits a small joint, then for each cumulative stream (suspected
# cases, suspected deaths, confirmed cases) extracts the posterior expected
# FIRST-bin increment (the modelled cumulative from grid day 1 to the first
# vintage day) and compares it with the observed first cumulative value.
# Decomposes the first suspected-case bin into the BVD signal vs the
# accumulated non-BVD background to test the background-accumulation
# hypothesis. Also reports the same modelled-vs-observed ratio for the LATER
# increments so the first-bin bias can be contrasted with the rest.
#
# Run: JULIA_NUM_THREADS=6 julia --project=. scripts/diagnose_first_bin.jl \
#        [samples] [chains]

using BVDOutbreakSize
using BVDOutbreakSize: bin_increments
using Turing: predict
using Statistics: median, quantile, mean
using Serialization: serialize, deserialize

const SAMPLES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 200
const CHAINS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2
const CHAIN_PATH = length(ARGS) >= 3 ? ARGS[3] : "logs/diag_chain.jls"

obs = load_observations()
const BP = obs.n - obs.who_first_sitrep_days

println("n (grid length)        = ", obs.n)
println("seeding day            = ", obs.seeding, " (grid day 1)")
println("cut-off                = ", obs.cutoff, " (grid day ", obs.n, ")")
println("first suspected vintage day = ", obs.reported_history.days[1],
    " (", obs.seeding, " + ", obs.reported_history.days[1] - 1, "d)")
println("first death vintage day     = ", obs.deaths_history.days[1])
println("first confirmed vintage day = ", obs.confirmed_history.days[1])
println()

_days_only(h) = (; days = h.days, counts = Int[])

function build_model(; fit = true)
    ## In predict mode (`fit = false`) strip the cumulative histories'
    ## counts so the per-vintage increments become latent and `predict`
    ## resamples them; with the counts present `vintage_obs` keeps them
    ## conditioned and predict cannot draw them.
    bvd_joint(obs.n, obs.exported_cases, obs.total_deaths,
        fit ? obs.reported_cases : missing, obs.exports_deaths,
        fit ? obs.confirmed_cases : missing,
        fit ? obs.tests_analysed : missing;
        confirmed_deaths = fit ? obs.confirmed_deaths : missing,
        deaths_history = fit ? obs.deaths_history :
                         _days_only(obs.deaths_history),
        reported_history = fit ? obs.reported_history :
                           _days_only(obs.reported_history),
        confirmed_history = obs.confirmed_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        suspected_daily_history = obs.suspected_daily_history,
        export_case_days = obs.export_case_days,
        export_death_days = obs.export_death_days,
        breakpoint = BP, background_re = true,
        confirmed_positivity_link = :composition,
        genetic = genetic_seeding_model, tmrca_days = obs.tmrca_days)
end

if isfile(CHAIN_PATH)
    println("Loading cached chain from ", CHAIN_PATH)
    global chn = deserialize(CHAIN_PATH)
else
    println("Fitting joint ($(SAMPLES)x$(CHAINS)) ...")
    global chn = nuts_sample(build_model(; fit = true);
        samples = SAMPLES, chains = CHAINS, progress = false)
    serialize(CHAIN_PATH, chn)
end

d = BVDOutbreakSize.fit_diagnostics(chn)
println("\n=== Convergence ===")
println("  divergences : ", d.n_divergent)
println("  max R-hat   : ", round(d.max_rhat; digits = 4))
println("  min ESS     : ", round(d.min_ess_bulk; digits = 1))

draws(sym) = vec(Array(chn[sym]))
q3(v) = (median(v), quantile(v, 0.05), quantile(v, 0.95))
function showq(label, v)
    m, lo, hi = q3(v)
    println("  ", rpad(label, 34), rpad(round(m; digits = 2), 12),
        "[", round(lo; digits = 2), ", ", round(hi; digits = 2), "]")
end

println("\n=== Headline ===")
for s in (:C_T, :CFR, :T, :p_drc, :lambda_bg, :background_total)
    showq(string(s), draws(s))
end

## Reconstruct per-vintage modelled increments per draw, using the
## composer deterministics already on the chain. The suspected-case daily
## series is `reports_daily = p_drc .* bvd_reports_daily .+ bg_daily`; the
## composer exposes `cumulative_onsets` and `lambda_bg` but not the daily
## reports series directly, so rebuild the increments via predict() means.

println("\n=== Per-vintage modelled vs observed (predictive expected) ===")
pp = predict(build_model(; fit = false), chn)

function modelled_increment_means(prefix)
    ks = [k for k in keys(pp) if occursin("$prefix.increments", string(k))]
    isempty(ks) && return Vector{Vector{Float64}}()
    ## Two possible layouts: one key per vintage (trailing [i]) each a flat
    ## draw vector, OR a single key whose draws are per-vintage vectors.
    has_idx = all(occursin(r"\[(\d+)\]", string(k)) for k in ks)
    if has_idx && length(ks) > 1
        idx(k) = parse(Int, match(r"\[(\d+)\]", string(k)).captures[1])
        ks = sort(ks; by = idx)
        return [Float64.(vec(collect(pp[k]))) for k in ks]
    else
        draws_of_vecs = vec(collect(pp[ks[1]]))
        nv = length(first(draws_of_vecs))
        return [Float64[d[i] for d in draws_of_vecs] for i in 1:nv]
    end
end

function report_stream(name, prefix, obs_inc)
    reps = modelled_increment_means(prefix)
    isempty(reps) && (println("  [", name, "] no predicted increments");
        return)
    nrep = min(length(reps), length(obs_inc))
    println("\n  -- ", name, " (", length(reps), " vintages) --")
    println("  ", rpad("vintage", 9), rpad("obs_inc", 10),
        rpad("model_mean", 12), rpad("model_med", 11),
        rpad("model_90%", 22), "ratio(mean/obs)")
    for i in 1:nrep
        mu = mean(reps[i])
        m, lo, hi = q3(reps[i])
        r = mu / max(obs_inc[i], 1)
        tag = i == 1 ? " <== FIRST" : ""
        println("  ", rpad(i, 9), rpad(obs_inc[i], 10),
            rpad(round(mu; digits = 1), 12),
            rpad(round(m; digits = 1), 11),
            rpad(string("[", round(lo; digits = 1), ", ",
                    round(hi; digits = 1), "]"), 22),
            round(r; digits = 2), tag)
    end
    return reps, obs_inc
end

## Observed between-vintage increments.
rc_obs = diff(vcat(0, collect(Int.(obs.reported_history.counts))))
dd_obs = diff(vcat(0, collect(Int.(obs.deaths_history.counts))))
## Confirmed early-window observed increments (vintages up to the first lab
## day), from the same windowing the model uses.
wins = BVDOutbreakSize.confirmed_positivity_windows(
    obs.confirmed_history, obs.lab_history, obs.lab_daily_history)
ca_obs = wins.early_increments

rc = report_stream("Suspected cases", "cases_state.reported_increments",
    rc_obs)
dd = report_stream("Suspected deaths", "deaths_state.death_increments",
    dd_obs)
ca = report_stream("Confirmed (early windows)",
    "confirmed_state.early_increments", ca_obs)

## Decompose the first suspected-case bin into BVD signal vs background.
## reports_daily = p_drc * bvd_reports_daily + bg_daily, summed over
## days 1:days[1]. The composer exposes lambda_bg (the daily background
## baseline) and the first vintage day; background contribution to bin 1 is
## approximately lambda_bg * days[1] (constant background) or the integral
## of the per-vintage RE over that span. Use background_total scaled by the
## fraction of the grid before the first vintage as an approximation, but
## better: rebuild from bg_daily is not on the chain. Use lambda_bg * d1.
println("\n=== First-bin background decomposition (suspected cases) ===")
d1 = obs.reported_history.days[1]
lbg = draws(:lambda_bg)
bg_first = lbg .* d1
obs_first = diff(vcat(0, collect(Int.(obs.reported_history.counts))))[1]
println("  first vintage day d1            = ", d1)
println("  observed first increment        = ", obs_first)
showq("lambda_bg (per day)", lbg)
showq("background in bin1 (lambda_bg*d1)", bg_first)
if rc !== nothing
    model_first = rc[1][1]
    over = model_first .- obs_first
    println()
    showq("modelled first bin", model_first)
    showq("over-prediction (model-obs)", over)
    ## fraction of the over-prediction explained by accumulated background
    frac = bg_first ./ max.(over, 1e-6)
    showq("bg_in_bin1 / over-prediction", clamp.(frac, -5, 5))
    ## background as a share of the modelled first bin
    showq("bg_in_bin1 / modelled_first_bin", bg_first ./ model_first)
end

println("\nDone.")
