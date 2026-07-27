# Four-cell attribution of the confirmed harmonisation break, on the cheap
# single-stream probes.
#
# Two mechanisms landed together on the confirmed streams: the DE-ANCHOR
# (dropping the break day's published 24h analysed denominator, which removes a
# BetaBinomial term and loosens `λ_bg` and the ascertainment) and the STEP
# fitted into that day's modelled mean. The mechanism is local to the confirmed
# streams, so `confirmed_only_model` and `confirmed_deaths_only_model` carry it
# without the joint's multi-hour cost.
#
# Cells:
#   off        no break day declared: the pre-feature baseline.
#   deanchor   break day declared with the step pinned to zero (`gross` set to
#              the whole observed increment, so the centre is 0, and sd = 0 so
#              nothing is sampled). Isolates the de-anchor.
#   determ     break day declared with the step FIXED at the published
#              discrepancy (real `gross`, sd = 0): the deterministic
#              correction, two parameters fewer than `sampled`.
#   sampled    break day declared with the step sampled around that discrepancy
#              (real `gross`, sd = 25): the shipped behaviour.
#
# `off -> deanchor` isolates the de-anchor, `deanchor -> determ` the correction
# itself, and `determ -> sampled` the cost of making it a parameter. Reports
# divergences, min bulk ESS, worst R-hat, the step marginal and the step's
# posterior correlation with `λ_bg`, the ascertainment and `C_T`.
#
# Run (from the repo root):
#   BREAK_SAMPLES=500 julia --project=. -t 4 \
#     scripts/diag_break_attribution.jl 2>&1 | tee logs/break_attr.log

using BVDOutbreakSize
using BVDOutbreakSize: FlexiChains
using Turing: @varname
using Statistics: median, quantile, cor
using Printf: @printf

const M = BVDOutbreakSize
const SAMPLES = parse(Int, get(ENV, "BREAK_SAMPLES", "500"))
const CHAINS = parse(Int, get(ENV, "BREAK_CHAINS", "2"))

mkpath("logs")
obs = M.load_observations()
const BP = obs.n - obs.who_first_sitrep_days

# Observed between-vintage increment on each break day, per stream: the whole
# step the likelihood sees. Passing it AS the gross count centres the step on
# zero, which is how the `deanchor` cell pins the step out of the way.
function net_increments(history, break_days)
    inc = diff(vcat(0, collect(history.counts)))
    days = collect(history.days)
    return Int[inc[findfirst(==(d), days)] for d in break_days]
end

const BRK = obs.confirmed_break_days
const NET_CASES = net_increments(obs.confirmed_history, BRK)
const NET_DEATHS = net_increments(obs.confirmed_deaths_history, BRK)

@info "Break days" BRK NET_CASES NET_DEATHS obs.confirmed_break_gross_cases obs.confirmed_break_gross_deaths

# (break_days, gross_cases, gross_deaths, sd) per cell.
function cells(net_c, net_d, gross_c, gross_d)
    [
        ("off", Int[], Int[], Int[], 25.0),
        ("deanchor", BRK, net_c, net_d, 0.0),
        ("determ", BRK, gross_c, gross_d, 0.0),
        ("sampled", BRK, gross_c, gross_d, 25.0)
    ]
end

# Pooled finite draws of a scalar parameter, or empty when absent.
function getdraws(chn, sym)
    v = try
        vec(Array(chn[sym]))
    catch
        return Float64[]
    end
    return Float64[x for x in v if x isa Real && isfinite(x)]
end

# First element of a (submodel-prefixed) vector-valued variable, per draw.
function getvecdraws(chn, vn)
    raw = try
        collect(chn[FlexiChains.Prefixed(vn)])
    catch
        return Float64[]
    end
    out = Float64[]
    for v in vec(raw)
        v === missing && continue
        vv = collect(v)
        isempty(vv) && continue
        isfinite(vv[1]) && push!(out, Float64(vv[1]))
    end
    return out
end

function summarise(d)
    isempty(d) && return (med = NaN, lo = NaN, hi = NaN)
    (med = median(d), lo = quantile(d, 0.05), hi = quantile(d, 0.95))
end

function corr(x, y)
    (isempty(x) || isempty(y) || length(x) != length(y)) && return NaN
    keep = isfinite.(x) .& isfinite.(y)
    sum(keep) < 10 && return NaN
    return cor(x[keep], y[keep])
end

# One probe: a label, a model builder taking the cell's break settings, and the
# name of the step variable that probe samples.
function run_probe(label, build, stepvn)
    results = Vector{NamedTuple}()
    for (cell, bd, gc, gd, sd) in cells(NET_CASES, NET_DEATHS,
        obs.confirmed_break_gross_cases, obs.confirmed_break_gross_deaths)
        println("\n", "="^72)
        println("$label / $cell   ($SAMPLES draws x $CHAINS chains, sd=$sd)")
        println("="^72)
        flush(stdout)
        t0 = time()
        chn = M.nuts_sample(build(bd, gc, gd, sd);
            samples = SAMPLES, chains = CHAINS, check_model = false)
        mins = (time() - t0) / 60
        dg = M.fit_diagnostics(chn)
        step = getvecdraws(chn, stepvn)
        ct = getdraws(chn, :C_T)
        pdrc = getdraws(chn, :p_drc)
        lbg = getdraws(chn, :λ_bg)
        rec = (; cell, mins,
            ndiv = dg.n_divergent, miness = dg.min_ess_bulk,
            maxrhat = dg.max_rhat,
            CT = summarise(ct), step = summarise(step),
            nstep = length(step),
            cor_step_lbg = corr(step, lbg),
            cor_step_pdrc = corr(step, pdrc),
            cor_step_CT = corr(step, ct),
            cor_pdrc_CT = corr(pdrc, ct))
        push!(results, rec)
        @printf("DONE %-9s %5.1f min  div=%4d  minESS=%6.0f  rhat=%.3f  C_T=%.0f  step=%.2f  cor(step,λ_bg)=%.2f cor(step,p_drc)=%.2f cor(step,C_T)=%.2f\n",
            cell, mins, rec.ndiv, rec.miness, rec.maxrhat, rec.CT.med,
            rec.step.med, rec.cor_step_lbg, rec.cor_step_pdrc,
            rec.cor_step_CT)
        flush(stdout)
    end
    return results
end

function report(label, results)
    println("\n\n", "#"^96)
    println("BREAK ATTRIBUTION — $label  ($SAMPLES draws x $CHAINS chains)")
    println("#"^96)
    @printf("%-9s %6s %6s %8s %6s %9s %20s %8s %8s %8s %8s\n",
        "cell", "mins", "div", "minESS", "rhat", "C_T",
        "step 90% CI", "nstep", "cor_lbg", "cor_pdrc", "cor_CT")
    for r in results
        @printf("%-9s %6.1f %6d %8.0f %6.3f %9.0f  %6.1f [%5.1f,%5.1f] %8d %8.2f %8.2f %8.2f\n",
            r.cell, r.mins, r.ndiv, r.miness, r.maxrhat, r.CT.med,
            r.step.med, r.step.lo, r.step.hi, r.nstep,
            r.cor_step_lbg, r.cor_step_pdrc, r.cor_step_CT)
    end
    println("#"^96)
    flush(stdout)
end

function cases_probe(bd, gc, gd, sd)
    M.confirmed_only_model(obs.n, obs.confirmed_cases;
        confirmed_history = obs.confirmed_history,
        lab_history = obs.lab_history,
        lab_daily_history = obs.lab_daily_history,
        tests_analysed = obs.tests_analysed,
        breakpoint = BP,
        confirmed_positivity_link = :composition,
        confirmed_break_days = bd,
        confirmed_break_gross_cases = gc,
        confirmed_break_sd = sd)
end

function deaths_probe(bd, gc, gd, sd)
    M.confirmed_deaths_only_model(
        obs.n, obs.confirmed_deaths, obs.total_deaths;
        deaths_history = obs.deaths_history,
        confirmed_deaths_history = obs.confirmed_deaths_history,
        breakpoint = BP,
        confirmed_break_days = bd,
        confirmed_break_gross_deaths = gd,
        confirmed_break_sd = sd)
end

cases_results = run_probe("confirmed cases", cases_probe,
    @varname(confirmed_step))
report("confirmed cases (confirmed_only_model)", cases_results)

deaths_results = run_probe("confirmed deaths", deaths_probe,
    @varname(cdeath_step))
report("confirmed deaths (confirmed_deaths_only_model)", deaths_results)
