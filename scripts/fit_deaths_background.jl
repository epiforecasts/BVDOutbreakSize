# Fitting harness for the deaths-background change (issue #193).
#
# Tracks the key inputs and outputs of the new death streams while fitting
# the joint model, and kills fits that are not working quickly:
#
#   * a TensorBoard callback streams every kept draw (params + diagnostics)
#     so the run can be watched live with `tensorboard --logdir <logdir>`;
#   * a watchdog callback aborts the fit early (throws) once the divergent
#     fraction runs away past a threshold, so a broken geometry does not
#     burn a full multi-hour fit before it is spotted.
#
# Two stages, cheapest first (the "fast MCMC iteration" workflow):
#   1. deaths submodel group in isolation (confirmed_deaths_only_model) as a
#      canary — if the death lab/positivity geometry is broken it shows here
#      in seconds;
#   2. the full joint with the confirmed-deaths stream on, streamed to
#      TensorBoard with the watchdog armed.
#
# Self-contained (no analysis.jl include, which needs the docs/plot env).
# TensorBoardLogger lives in the test environment, so run with:
#   julia --project=test -t 4 scripts/fit_deaths_background.jl

using BVDOutbreakSize
using Turing: Turing, @varname
using TensorBoardLogger
using Statistics: median
using Printf: @printf

const REPO = pkgdir(BVDOutbreakSize)
include(joinpath(@__DIR__, "joint_setup.jl"))

obs = load_observations()
fit_args = build_fit_args(obs)
growth_now = growth_for(obs)
genetic_seeding = T -> genetic_seeding_model(T, obs.genetic_tmrca_days;
    tmrca_days_sd = obs.genetic_tmrca_days_sd)

logroot = get(ENV, "TB_LOGDIR", joinpath(REPO, "logs", "deaths_background"))
mkpath(logroot)
println("TensorBoard log root: ", logroot)
println("  view with: tensorboard --logdir ", logroot)

## ------------------------------------------------------------------
## Stage 1 — deaths submodel-group canary (confirmed_deaths_only).
## ------------------------------------------------------------------
susp_deaths = obs.total_deaths
conf_deaths = obs.confirmed_death_history.values[end]
@printf("\n=== Stage 1: deaths canary (susp=%d, confirmed=%d) ===\n",
    susp_deaths, conf_deaths)

stage1_cb = compose_callbacks(
    tensorboard_callback(joinpath(logroot, "canary"); every = 25),
    watchdog_callback(; min_iter = 80, max_div_frac = 0.5))

chn_canary = try
    nuts_sample(confirmed_deaths_only_model(susp_deaths, conf_deaths);
        samples = 500, chains = 4, seed = 1, callback = stage1_cb)
catch e
    e isa EarlyKill ? (println("STAGE 1 KILLED: ", e.msg); rethrow(e)) :
    rethrow(e)
end

let div = sum(vec(Array(chn_canary[:numerical_error]))),
    τd = vec(Array(chn_canary[:τ_death])), C = vec(Array(chn_canary[:cumulative_cases]))

    @printf("canary: divergences=%d  τ_death median=%.3f  C(T) median=%.0f\n",
        Int(div), median(τd), median(C))
end

## ------------------------------------------------------------------
## Stage 2 — full joint with confirmed-deaths stream on.
## ------------------------------------------------------------------
println("\n=== Stage 2: joint fit (confirmed-deaths stream on) ===")
joint_cb = compose_callbacks(
    tensorboard_callback(joinpath(logroot, "joint"); every = 20),
    watchdog_callback(; min_iter = 100, max_div_frac = 0.4))

joint_model = bvd_joint(obs.exported_cases, fit_args.deaths,
    fit_args.reported, fit_args.export_deaths; fit_args.kw...,
    growth = growth_now,
    first_export_detection_delta = obs.first_export_detection_delta,
    report_onset_offset = report_onset_offset(obs.as_of_date),
    genetic = genetic_seeding)

chn_joint = try
    nuts_sample(joint_model; samples = 1_000, chains = 4, seed = 1,
        callback = joint_cb)
catch e
    e isa EarlyKill ? (println("STAGE 2 KILLED: ", e.msg); rethrow(e)) :
    rethrow(e)
end

## ------------------------------------------------------------------
## Report the key inputs/outputs of the death streams.
## ------------------------------------------------------------------
function row(chn, sym)
    d = vec(Array(chn[sym]))
    s = posterior_summary(d)
    @printf("  %-22s median=%-9.3f 90%% CI [%.3f, %.3f]\n",
        string(sym), median(d), s.lo90, s.hi90)
end

div = sum(vec(Array(chn_joint[:numerical_error])))
@printf("\njoint divergences: %d\n", Int(div))
println("key outputs:")
for sym in (:cumulative_cases, :cumulative_infections, :r,
    :λ_bg_death, :τ_death, :CFR, :s_test, :spec_test)
    try
        row(chn_joint, sym)
    catch
        @printf("  %-22s (unavailable)\n", string(sym))
    end
end

## Posterior-predictive confirmed deaths must recentre on the observed 17.
pp_model = bvd_joint(obs.exported_cases, fit_args.deaths, fit_args.reported,
    fit_args.export_deaths; fit_args.kw...,
    growth = growth_now,
    first_export_detection_delta = obs.first_export_detection_delta,
    report_onset_offset = report_onset_offset(obs.as_of_date),
    genetic = genetic_seeding,
    confirmed_deaths = Union{Missing, Int}[missing],
    confirmed_death_offsets = [0])
pp = Turing.predict(pp_model, chn_joint)
cd = reduce(vcat, vec(Array(pp[@varname(confirmed_deaths)])))
@printf("\nconfirmed-death PP: median=%.1f (observed=%d)\n",
    median(cd), conf_deaths)
println("done.")
