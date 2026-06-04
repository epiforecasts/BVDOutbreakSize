# Severity-selection (q0) prior sweep.
#
# The confirmed-positivity dynamics can be explained two ways that the
# model confounds: (i) severe-first TESTING SELECTION — the first-tested
# batch is enriched for BVD, positivity then falls as testing widens — or
# (ii) a large non-BVD BACKGROUND among suspected cases. The default
# `q0 ~ Beta(20, 1.5)` (mean ≈ 0.93) is a radical selection prior that
# pins explanation (i) and starves the background.
#
# This sweep weakens q0 (less radical selection) to see whether the
# explanation shifts onto the suspected-case background — λ_bg should grow,
# C(T) should fall toward the deaths, and the confirmed-death PP should
# move toward the observed 17. Even with severity-induced sampling the
# enrichment cannot be extreme (other haemorrhagic/febrile illness is also
# severe), so a moderate q0 is the more defensible prior.
#
# Run with: julia --project=test -t 4 scripts/severity_sweep.jl

using BVDOutbreakSize
using Turing: Turing, @varname
using Distributions: Beta
using Statistics: median
using Printf: @printf

const REPO = pkgdir(BVDOutbreakSize)
include(joinpath(@__DIR__, "joint_setup.jl"))

obs = load_observations()
fit_args = build_fit_args(obs)
growth_now = growth_for(obs)
genetic_seeding = T -> genetic_seeding_model(T, obs.genetic_tmrca_days;
    tmrca_days_sd = obs.genetic_tmrca_days_sd)
conf_deaths = obs.confirmed_death_history.values[end]

med(chn, s) = median(vec(Array(chn[s])))

## q0 prior variants: baseline (radical, near-pure first batch) → moderate
## → weak (the first-tested batch is only modestly BVD-enriched).
variants = [
    ("baseline Beta(20,1.5) m=.93", Beta(20.0, 1.5)),
    ("moderate Beta(4,2)   m=.67", Beta(4.0, 2.0)),
    ("weak     Beta(2,2)   m=.50", Beta(2.0, 2.0))
]

const SAMPLES = 500
const CHAINS = 4

function make_joint(q0p; observe_cdeath = true)
    sel = test_selection_model(q0_prior = q0p)
    cdeath = observe_cdeath ? fit_args.kw.confirmed_deaths :
             Union{Missing, Int}[missing]
    bvd_joint(obs.exported_cases, fit_args.deaths, fit_args.reported,
        fit_args.export_deaths; fit_args.kw...,
        confirmed_deaths = cdeath,
        growth = growth_now,
        first_export_detection_delta = obs.first_export_detection_delta,
        report_onset_offset = report_onset_offset(obs.as_of_date),
        genetic = genetic_seeding,
        confirmed_q_random_effect = confirmed_q_re_model,
        test_selection = sel)
end

println("q0 severity-selection sweep (4x500, watchdog armed)\n")
@printf("%-30s %8s %8s %8s %7s %8s %7s %6s %5s\n",
    "q0 prior", "C(T)", "λ_bg", "%susp", "q0post", "posit", "cdPP",
    "div", "R̂")

for (label, q0p) in variants
    chn = try
        nuts_sample(make_joint(q0p); samples = SAMPLES, chains = CHAINS,
            seed = 1, callback = watchdog_callback(; min_iter = 120,
                max_div_frac = 0.6),
            progress = false)
    catch e
        e isa EarlyKill ? (println(label, "  KILLED: ", e.msg); continue) :
        rethrow(e)
    end
    d = fit_diagnostics(chn)
    C = med(chn, :cumulative_cases)
    λ = med(chn, :λ_bg)
    Tm = med(chn, :T)
    pct = 100 * λ * Tm / obs.reported_cases
    q0post = med(chn, :q0)
    posit = med(chn, :p_positive)
    ## Confirmed-death posterior predictive (recentre check vs 17).
    pp = Turing.predict(make_joint(q0p; observe_cdeath = false), chn)
    cdpp = median(reduce(vcat, vec(Array(pp[@varname(confirmed_deaths)]))))
    @printf("%-30s %8.0f %8.3f %7.0f%% %7.3f %8.3f %7.1f %6d %5.3f\n",
        label, C, λ, pct, q0post, posit, cdpp, d.n_divergent, d.max_rhat)
end
println("\ndone.")
