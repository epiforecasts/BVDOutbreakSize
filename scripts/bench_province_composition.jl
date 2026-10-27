# Mooncake gradient-timing diagnostic for the province composition's split.
#
# `province_composition_model` (src/models/observations.jl) scores the
# per-province split with sequential `BetaBinomial`s at a shared `ρ`
# (`stick_breaking_loglik`). `zone_composition_model` (src/models/zone.jl)
# scores the same kind of split in closed form, as a Dirichlet-multinomial
# added once with `@addlogprob!` (`dirichlet_multinomial_logpdf`).
#
# The two agree only at two patches. With two patches there is one free
# split per vintage, and `BetaBinomial(N, π₁, ρ)` is the Dirichlet-multinomial
# margin. From three patches on they are different distributions: the
# Dirichlet-multinomial's sequential decomposition shrinks each later step's
# concentration by its remaining tail, while the stick-breaking keeps `ρ`
# fixed at every step. The fitted province model has eight patches, so the
# closed form is not a drop-in. Swapping it in changes the overdispersion the
# likelihood describes and needs a refit.
#
# This checks that divergence on synthetic data, then times both forms'
# Mooncake compile and steady-state gradient at the fitted model's patch
# count. Each form's timing runs in its own fresh process, for the reason
# `benchmark/compile_one.jl` gives: the first Mooncake build in a session
# pays a fixed cost that later builds in the same session avoid.
#
# Run: julia --project=. scripts/bench_province_composition.jl

using BVDOutbreakSize
using BVDOutbreakSize: stick_breaking_loglik, dirichlet_multinomial_logpdf,
    composition_shares, PROVINCE_NAMES
using Mooncake: Mooncake
using Statistics: median
using Printf: @printf

## Patch count of the fitted province composition.
const N_PATCHES = length(PROVINCE_NAMES)

## Stick-breaking: one `stick_breaking_loglik` call over every patch and
## vintage, as `composition_split_model` scores it.
function stick_breaking_objective(
        obs::AbstractMatrix{<:Integer}, modelled::AbstractMatrix,
        weight::AbstractVector, ρ
    )
    shares = composition_shares(weight, modelled)
    np, nv = size(shares)
    groups = repeat(1:nv; inner = np)
    return stick_breaking_loglik(groups, vec(obs), vec(shares), ρ)
end

## Closed form: one Dirichlet-multinomial per vintage at concentration
## `κ .* shares[:, v]`, with `κ = (1 - ρ) / ρ`, the map from `ρ` to a
## concentration that `safe_betabinomial` uses.
function dirichlet_multinomial_objective(
        obs::AbstractMatrix{<:Integer}, modelled::AbstractMatrix,
        weight::AbstractVector, ρ
    )
    shares = composition_shares(weight, modelled)
    nv = size(shares, 2)
    κ = (1 - ρ) / ρ
    lp = zero(promote_type(eltype(shares), typeof(κ)))
    for v in 1:nv
        lp += dirichlet_multinomial_logpdf(obs[:, v], κ .* shares[:, v])
    end
    return lp
end

## Synthetic (patches × vintages) data: a dominant first patch and
## geometrically thinner later ones, growing over the vintages.
function synthetic_data(np::Integer, nv::Integer)
    raw = [0.25^(p - 1) for p in 1:np]
    base_weights = raw ./ sum(raw)
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
            c = min(round(Int, shares[p, v] * totals[v]), remaining)
            obs[p, v] = c
            remaining -= c
        end
        obs[np, v] = remaining
    end
    return obs, modelled
end

## Both forms' value at a shared ρ, at `np` patches and `nv` vintages.
function compare_values(np::Integer, nv::Integer, ρ)
    obs, modelled = synthetic_data(np, nv)
    weight = ones(np)
    stick = stick_breaking_objective(obs, modelled, weight, ρ)
    closed = dirichlet_multinomial_objective(obs, modelled, weight, ρ)
    return (; stick, closed)
end

## One form's compile (cache construction plus the first gradient) and the
## median of `reps` further gradients.
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

## One form's gradient timing at the fitted patch count, printed as one
## parseable line. `main` runs this in a fresh subprocess per form.
function run_form(form::AbstractString)
    obs, modelled = synthetic_data(N_PATCHES, 100)
    weight = ones(N_PATCHES)
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
        "Province composition: stick-breaking against closed-form " *
            "Dirichlet-multinomial\n"
    )

    println("1. Log likelihood at a shared ρ = 0.05, by patch count (20 vintages):")
    for np in unique((2, 3, N_PATCHES))
        v = compare_values(np, 20, 0.05)
        @printf(
            "   np=%d  stick-breaking %10.3f   Dirichlet-multinomial %10.3f   difference %8.3f\n",
            np, v.stick, v.closed, v.stick - v.closed
        )
    end
    println()

    println(
        "2. Gradient cost at the fitted patch count " *
            "(np=$(N_PATCHES), 100 vintages), each form in a fresh process:"
    )
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
    t_stick = timings["stick-breaking"]
    t_closed = timings["dirichlet-multinomial"]
    @printf(
        "   stick-breaking         compile %8.1f ms   gradient %8.2f µs\n",
        t_stick.compile * 1.0e3, t_stick.steady * 1.0e6
    )
    @printf(
        "   Dirichlet-multinomial  compile %8.1f ms   gradient %8.2f µs\n",
        t_closed.compile * 1.0e3, t_closed.steady * 1.0e6
    )
    @printf(
        "   stick-breaking gradient time over closed form: %.2f\n",
        t_stick.steady / t_closed.steady
    )
    return nothing
end

main()
