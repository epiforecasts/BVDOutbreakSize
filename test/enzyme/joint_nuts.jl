# Exploratory harness: can Enzyme drive a SHORT NUTS fit of the full
# renewal `bvd_joint` to completion, or does warmup hit the
# `EnzymeNonScalarReturn` blocker seen on the earlier 29-dim joint?
#
# The single-gradient path is already known to work and match Mooncake
# (see joint_explore.jl). This exercises trajectory exploration into the
# extreme warmup regions where the model's isfinite-guard branches return
# a constant, which is where Enzyme previously failed.
#
# Run from the isolated Enzyme sub-environment:
#   julia --project=test/enzyme test/enzyme/joint_nuts.jl
#
# Strategy: tiny fit first (50 samples, 1 chain) to see if warmup
# completes at all; only if that works, a longer 300x2 fit. Each Enzyme
# fit is paired with the same Mooncake fit for a wall-clock comparison and
# a sanity check that the chains are not all-divergent.
#
# STATUS as of the renewal merge (`background_re` switch): the Enzyme
# tiny fit no longer reaches NUTS — it fails at the first gradient with
# the boxed-closure `TypeError` documented in joint_explore.jl, so the
# 300x2 fit is skipped. This harness is kept for when the model-side
# closure box is resolved; until then it records the failure.

using Enzyme    # loads BVDOutbreakSizeEnzymeExt so enzyme_adtype() exists
using Mooncake  # default backend
using BVDOutbreakSize: default_adtype, enzyme_adtype, bvd_joint,
                       load_observations, genetic_seeding_model, nuts_sample
using Turing: DynamicPPL

obs = load_observations()
breakpoint = obs.n - obs.who_first_sitrep_days

build_joint() = bvd_joint(
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
    breakpoint = breakpoint,
    background_re = true,
    genetic = genetic_seeding_model,
    tmrca_days = obs.tmrca_days)

# Fraction of divergent transitions across the post-warmup chain, a quick
# all-divergent sanity flag. Returns NaN if the diagnostic is absent.
function divergent_fraction(chn)
    try
        d = chn[:numerical_error]
        return sum(d) / length(d)
    catch
        return NaN
    end
end

# Run one short fit, returning (:ok, chain, seconds) or (:error, msg, t).
function try_fit(label, adtype; samples, chains, target_accept = 0.9)
    println("--- $label: NUTS samples=$samples chains=$chains ---")
    t0 = time()
    try
        chn = nuts_sample(build_joint();
            adtype = adtype, samples = samples, chains = chains,
            target_accept = target_accept, check_model = false)
        dt = time() - t0
        df = divergent_fraction(chn)
        println("  OK in $(round(dt; digits = 1))s, divergent fraction = ",
            round(df; digits = 3))
        return (:ok, chn, dt)
    catch err
        dt = time() - t0
        msg = sprint(showerror, err)
        first_line = first(split(msg, '\n'))
        println("  ERROR after $(round(dt; digits = 1))s: ", first_line)
        println("  --- full error head ---")
        println(first(msg, 1500))
        println("  --- backtrace head ---")
        for (i, fr) in enumerate(stacktrace(catch_backtrace()))
            i > 20 && break
            println("    ", fr)
        end
        return (:error, first_line, dt)
    end
end

println("base joint dimension check")
let m = build_joint()
    vi = DynamicPPL.link(DynamicPPL.VarInfo(m), m)
    println("  dim = ", length(collect(vi[:])))
end

# 1. Tiny Enzyme fit: does warmup complete at all?
status_tiny, res_tiny, t_en_tiny = try_fit(
    "Enzyme tiny", enzyme_adtype(); samples = 50, chains = 1)

# Paired tiny Mooncake fit for a baseline wall-clock.
_, _, t_mc_tiny = try_fit(
    "Mooncake tiny", default_adtype(); samples = 50, chains = 1)

println()
println("=== tiny fit summary ===")
println("  Enzyme warmup+sample (50x1): ",
    status_tiny == :ok ? "OK $(round(t_en_tiny; digits = 1))s" :
    "FAILED ($res_tiny)")
println("  Mooncake (50x1): $(round(t_mc_tiny; digits = 1))s")

# 2. Longer fit only if the tiny Enzyme fit completed.
if status_tiny == :ok
    status_big, _, t_en_big = try_fit(
        "Enzyme 300x2", enzyme_adtype(); samples = 300, chains = 2)
    _, _, t_mc_big = try_fit(
        "Mooncake 300x2", default_adtype(); samples = 300, chains = 2)
    println()
    println("=== 300x2 fit summary ===")
    println("  Enzyme: ",
        status_big == :ok ? "OK $(round(t_en_big; digits = 1))s" :
        "FAILED")
    println("  Mooncake: $(round(t_mc_big; digits = 1))s")
else
    println()
    println("Skipping 300x2 fit: tiny Enzyme fit did not complete.")
end
