# Mooncake gradient-timing diagnostic for the province composition's split.
#
# `province_composition_model` (src/models/observations.jl) scores the
# per-province split with sequential `BetaBinomial`s at a shared `ρ`
# (`stick_breaking_loglik`, the "stick-breaking" of the model's docstring).
# `zone_composition_model` scores the same kind of split a different way,
# in closed form as a single Dirichlet-multinomial `@addlogprob!`
# (`dirichlet_multinomial_logpdf`, src/models/zone.jl).
#
# The two agree only at two patches: with `np = 2` there is one free split
# per vintage, and `BetaBinomial(N, π₁, ρ)` already *is* the
# Dirichlet-multinomial margin. At `np ≥ 3`, as the fitted model has (Ituri,
# Nord-Kivu, Sud-Kivu), they diverge: the exact Dirichlet-multinomial
# decomposition would shrink each later step's concentration by its
# remaining tail, while the current form keeps `ρ` fixed at every step. So
# the closed form is not a drop-in at the province model's patch count;
# swapping it in would change what the likelihood says about a province's
# overdispersion, not just how it is computed, and would need a refit, not
# a parameter translation.
#
# This checks that divergence, then times both forms' Mooncake compile and
# steady-state gradient at `np = 3`. Each form's compile runs in its own
# fresh process, invoked below with `--form`, for the reason
# `benchmark/compile_one.jl` gives: the first Mooncake build in a session
# pays a fixed floor that every later build in the same session avoids, so
# timing both forms one after another in one process would charge that
# floor to whichever ran first and hide it from whichever ran second.
#
# Run: julia --project=. scripts/bench_province_composition.jl

using BVDOutbreakSize
using BVDOutbreakSize: stick_breaking_loglik, dirichlet_multinomial_logpdf,
    composition_shares
using Mooncake: Mooncake
using Statistics: median
using Printf: @printf

## Current form: one `stick_breaking_loglik` call over all patches and
## vintages, the vectorised form `composition_split_model` scores with on
## the fitting path.
function stick_breaking_objective(
        obs::AbstractMatrix{<:Integer}, modelled::AbstractMatrix,
        weight::AbstractVector, ρ
    )
    shares = composition_shares(weight, modelled)
    np, nv = size(shares)
    groups = repeat(1:nv; inner = np)
    return stick_breaking_loglik(groups, vec(obs), vec(shares), ρ)
end

## Candidate form: one closed-form Dirichlet-multinomial per vintage, at
## concentration vector `κ .* shares[:, v]`, with `κ = (1 - ρ) / ρ` the same
## ρ → concentration map `safe_betabinomial` uses.
function dirichlet_multinomial_objective(
        obs::AbstractMatrix{<:Integer}, modelled::AbstractMatrix,
        weight::AbstractVector, ρ
    )
    shares = composition_shares(weight, modelled)
    np, nv = size(shares)
    κ = (1 - ρ) / ρ
    lp = zero(promote_type(eltype(shares), typeof(κ)))
    for v in 1:nv
        lp += dirichlet_multinomial_logpdf(obs[:, v], κ .* shares[:, v])
    end
    return lp
end

## Synthetic (patches × vintages) data at the order of magnitude of the
## province confirmed table: a dominant patch 1, a smaller patch 2 and,
## for `np = 3`, a thin patch 3, growing over the vintages.
function synthetic_data(np::Integer, nv::Integer)
    base_weights = np == 2 ? [0.82, 0.18] : [0.78, 0.19, 0.03]
    modelled = [
        base_weights[p] * (20.0 + 3.0 * v + 2.0 * sin(v / 3))
            for p in 1:np, v in 1:nv
    ]
    shares = composition_shares(ones(np), modelled)
    totals = [round(Int, sum(@view modelled[:, v])) for v in 1:nv]
    obs = Matrix{Int}(undef, np, nv)
    for v in 1:nv
        remaining = totals[v]
        for p in 1:(np - 1)
            c = round(Int, shares[p, v] * totals[v])
            obs[p, v] = c
            remaining -= c
        end
        obs[np, v] = remaining
    end
    return obs, modelled
end

## At `np` patches and `nv` vintages, the two forms' value at a shared ρ.
function compare_values(np::Integer, nv::Integer, ρ)
    obs, modelled = synthetic_data(np, nv)
    weight = ones(np)
    v_old = stick_breaking_objective(obs, modelled, weight, ρ)
    v_new = dirichlet_multinomial_objective(obs, modelled, weight, ρ)
    return (; v_old, v_new)
end

## One form's compile (cache construction plus the Mooncake rule build) and
## the median of `reps` further gradients once that is paid.
function time_gradient(f, args...; reps = 50)
    local cache
    compile = @elapsed begin
        cache = Mooncake.prepare_gradient_cache(f, args...)
        Mooncake.value_and_gradient!!(cache, f, args...)
    end
    times = Float64[]
    for _ in 1:reps
        push!(times, @elapsed Mooncake.value_and_gradient!!(cache, f, args...))
    end
    return (; compile, steady = median(times))
end

## One form's gradient timing at the model's own patch count, printed as one
## parseable line. Called in a fresh subprocess by `main` below.
function run_form(form::AbstractString)
    obs, modelled = synthetic_data(3, 100)
    weight = ones(3)
    ρ = 0.05
    f = form == "stick-breaking" ?
        ((w, r) -> stick_breaking_objective(obs, modelled, w, r)) :
        ((w, r) -> dirichlet_multinomial_objective(obs, modelled, w, r))
    t = time_gradient(f, weight, ρ)
    @printf("RESULT %.6f %.9f\n", t.compile, t.steady)
    return nothing
end

function main()
    length(ARGS) == 1 && return run_form(ARGS[1])

    println(
        "Province composition: stick-breaking vs closed-form " *
            "Dirichlet-multinomial\n"
    )

    println("1. Value at a shared ρ, by patch count (20 vintages):")
    for np in (2, 3)
        v = compare_values(np, 20, 0.05)
        @printf(
            "   np=%d  stick-breaking %10.3f   Dirichlet-multinomial %10.3f   Δ %8.3f\n",
            np, v.v_old, v.v_new, v.v_old - v.v_new
        )
    end
    println(
        "   → equal at np=2 (one free split is already the DM margin), " *
            "diverge at np=3.\n"
    )

    println("2. Gradient cost at the model's own patch count (np=3, 100 vintages):")
    println("   (each form in its own fresh process, see header)")
    timings = Dict{String, NamedTuple}()
    for form in ("stick-breaking", "dirichlet-multinomial")
        cmd = `$(Base.julia_cmd()) --project=$(dirname(@__DIR__)) $(@__FILE__) $form`
        out = read(cmd, String)
        m = match(r"RESULT ([\d.eE+-]+) ([\d.eE+-]+)", out)
        m === nothing && error("no RESULT line from $form:\n$out")
        timings[form] = (
            compile = parse(Float64, m.captures[1]),
            steady = parse(Float64, m.captures[2]),
        )
    end
    t_old = timings["stick-breaking"]
    t_new = timings["dirichlet-multinomial"]
    @printf(
        "   stick-breaking         compile %7.1f ms   gradient %7.2f µs\n",
        t_old.compile * 1.0e3, t_old.steady * 1.0e6
    )
    @printf(
        "   Dirichlet-multinomial  compile %7.1f ms   gradient %7.2f µs  (%.2fx)\n",
        t_new.compile * 1.0e3, t_new.steady * 1.0e6, t_old.steady / t_new.steady
    )
    println(
        "   → both forms already emit one `@addlogprob!`/`~` term rather " *
            "than one per cell, so neither compile is the hour the model's " *
            "docstring describes; a value divergence from np=3 on, not " *
            "compile time, is what would need settling before any swap."
    )
    return nothing
end

main()
