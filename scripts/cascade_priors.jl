# Cascade single-stream posteriors into the joint to speed up fitting.
#
# Two composable speedups, both driven by cheap single-stream fits:
#
#   1. Cascaded PRIORS for the transferable, ~outbreak-size-independent
#      physical nuisances — the onset-to-death delay (from deaths_only) and
#      the onset-to-report delay (from cases_only). These are stream-local
#      Gamma-delay shapes, identified in isolation, and ~orthogonal to C(T),
#      so moment-matching their sub-fit marginals into tighter priors shrinks
#      the joint's nuisance space without biasing the outbreak size. This is
#      a deliberate cut/empirical-Bayes approximation (the sub-fit data also
#      appears in the joint), justified only for parameters the joint barely
#      updates beyond the single stream.
#
#   2. InitFromParams INITIALISATION of the joint chains at the single-stream
#      posterior medians (r, m, CFR, the delays, dispersion, λ_bg, τ_forward),
#      with the rest from the prior. This is inference-neutral — it only moves
#      the warmup starting point off the prior into the data-informed region,
#      cutting warmup. C(T)/r are used ONLY as a starting point here, never as
#      a prior (that would double-count the data).
#
# Compares baseline (prior init, default priors) against cascaded+init on
# matched short fits: min bulk ESS per wall-second and the C(T) posterior.
# Run with: julia --project=test -t 4 scripts/cascade_priors.jl

using BVDOutbreakSize
using Turing: Turing
import Turing.DynamicPPL as DynamicPPL
using Statistics: median
using Printf: @printf

const REPO = pkgdir(BVDOutbreakSize)
include(joinpath(@__DIR__, "joint_setup.jl"))

obs = load_observations()
fit_args = build_fit_args(obs)
growth_now = growth_for(obs)
genetic_seeding = T -> genetic_seeding_model(T, obs.genetic_tmrca_days;
    tmrca_days_sd = obs.genetic_tmrca_days_sd)

med(chn, sym) = median(vec(Array(chn[sym])))

## ------------------------------------------------------------------
## Stage A — cheap single-stream sub-fits.
## ------------------------------------------------------------------
const SUB_SAMPLES = 500
const SUB_CHAINS = 4
println("=== Stage A: single-stream sub-fits ===")
t_sub = @elapsed begin
    chn_deaths = nuts_sample(deaths_only_model(obs.total_deaths;
            growth = growth_now);
        samples = SUB_SAMPLES, chains = SUB_CHAINS, seed = 1, progress = false)
    chn_cases = nuts_sample(cases_only_model(obs.reported_cases;
            growth = growth_now);
        samples = SUB_SAMPLES, chains = SUB_CHAINS, seed = 1, progress = false)
end
@printf("sub-fits done in %.1fs\n", t_sub)
@printf("  deaths-only:  r=%.3f m=%.2f CFR=%.3f  (α=%.2f θ=%.2f)\n",
    med(chn_deaths, :r), med(chn_deaths, :m), med(chn_deaths, :CFR),
    med(chn_deaths, :α), med(chn_deaths, :θ))
@printf("  cases-only:   r=%.3f m=%.2f λ_bg=%.3f τ_fwd=%.3f (α_rep=%.2f θ_rep=%.2f)\n",
    med(chn_cases, :r), med(chn_cases, :m), med(chn_cases, :λ_bg),
    med(chn_cases, :τ_forward), med(chn_cases, :α_rep),
    med(chn_cases, :θ_rep))

## ------------------------------------------------------------------
## Stage B — build cascaded priors + InitFromParams.
## ------------------------------------------------------------------
## Cascaded delay priors (transferable physical nuisances). Moment-match
## the sub-fit marginals to the same truncated-Normal families the delay
## submodels use. The onset-to-death delay lives inside `deaths_model`, so
## it is injected through a `deaths` submodel wrapper (the report delay is a
## direct `bvd_joint` kwarg).
delay_cascade = delay_model(
    alpha_prior = fit_trunc_normal(vec(Array(chn_deaths[:α]))),
    theta_prior = fit_trunc_normal(vec(Array(chn_deaths[:θ]))))
report_cascade = report_delay_model(
    alpha_prior = fit_trunc_normal(vec(Array(chn_cases[:α_rep]))),
    theta_prior = fit_trunc_normal(vec(Array(chn_cases[:θ_rep]))))
deaths_cascade = (td, gs, k, edges;
    kwargs...) -> deaths_model(td, gs, k, edges; delay = delay_cascade, kwargs...)

## Initialisation NamedTuple from sub-fit medians (rest from the prior).
## Deaths anchor outbreak size, cases anchor the reporting nuisances.
init_nt = (; r = med(chn_deaths, :r), m = med(chn_deaths, :m),
    CFR = med(chn_deaths, :CFR), α = med(chn_deaths, :α),
    θ = med(chn_deaths, :θ), inv_sqrt_k = med(chn_deaths, :inv_sqrt_k),
    α_rep = med(chn_cases, :α_rep), θ_rep = med(chn_cases, :θ_rep),
    λ_bg = med(chn_cases, :λ_bg), τ_forward = med(chn_cases, :τ_forward))
cascade_init = DynamicPPL.InitFromParams(init_nt, DynamicPPL.InitFromPrior())

## ------------------------------------------------------------------
## Stage C — baseline vs cascaded+init joint fits (matched, short).
## ------------------------------------------------------------------
const SAMPLES = 500
const CHAINS = 4

function make_joint(; kw...)
    bvd_joint(obs.exported_cases, fit_args.deaths,
        fit_args.reported, fit_args.export_deaths; fit_args.kw...,
        growth = growth_now,
        first_export_detection_delta = obs.first_export_detection_delta,
        report_onset_offset = report_onset_offset(obs.as_of_date),
        genetic = genetic_seeding, kw...)
end

function summarise(label, chn, secs)
    d = fit_diagnostics(chn)
    C = vec(Array(chn[:cumulative_cases]))
    ess_per_s = d.min_ess_bulk / secs
    @printf("%-14s %6.1fs  minESS=%-7.1f maxR=%-6.3f div=%-4d ESS/s=%-6.2f  C(T)=%.0f [%.0f, %.0f]\n",
        label, secs, d.min_ess_bulk, d.max_rhat, d.n_divergent, ess_per_s,
        median(C), posterior_summary(C).lo90, posterior_summary(C).hi90)
    ## Suspected-case background: is the posterior pushing against the
    ## prior's upper tail (evidence the λ_bg prior is too tight)? The
    ## default prior is truncated(Normal(0, 1); lower = 0), 97.5% ≈ 2.24/day.
    λ = vec(Array(chn[:λ_bg]))
    s = posterior_summary(λ)
    ## Background cases ≈ λ_bg · T (T the latent seeding-to-cut-off time).
    Tmed = median(vec(Array(chn[:T])))
    bg_cases = median(λ) * Tmed
    @printf("    lambda_bg median=%.3f 90%% CI [%.3f, %.3f]  (prior 97.5%% ~ 2.24/day; ~%.0f bg cases, %.0f%% of %d suspected)\n",
        median(λ), s.lo90, s.hi90, bg_cases,
        100 * bg_cases / obs.reported_cases, obs.reported_cases)
    return ess_per_s
end

println("\n=== Stage C: baseline vs cascaded+init joint ===")
wd() = watchdog_callback(; min_iter = 120, max_div_frac = 0.5)

t_base = @elapsed chn_base = try
    nuts_sample(make_joint(); samples = SAMPLES, chains = CHAINS, seed = 1,
        callback = wd(), progress = false)
catch e
    e isa EarlyKill ? (println("BASELINE KILLED: ", e.msg); rethrow(e)) :
    rethrow(e)
end
base_eps = summarise("baseline", chn_base, t_base)

t_casc = @elapsed chn_casc = try
    nuts_sample(
        make_joint(deaths = deaths_cascade, report_delay = report_cascade);
        samples = SAMPLES, chains = CHAINS, seed = 1, init = cascade_init,
        callback = wd(), progress = false)
catch e
    e isa EarlyKill ? (println("CASCADE KILLED: ", e.msg); rethrow(e)) :
    rethrow(e)
end
casc_eps = summarise("cascade+init", chn_casc, t_casc)

@printf("\nESS/s speedup (cascade+init / baseline): %.2fx\n",
    casc_eps / base_eps)
@printf("wall-clock incl. sub-fits: cascade %.1fs vs baseline %.1fs\n",
    t_sub + t_casc, t_base)
println("done.")
