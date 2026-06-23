# A/B: current `_pmf_from_dic` (two comprehensions, `c` then `raw`) vs a
# fused single-loop variant that keeps only the previous CDF value, removing
# the intermediate `c` array from the Mooncake tape. Same nmax+2 censored-CDF
# evaluations (preserving #327), same values; the question is whether dropping
# one allocated array speeds the reverse pass.
#
# Run: julia --project=. scripts/bench_pmf_fused.jl

using BVDOutbreakSize: lognormal_meansd
using CensoredDistributions: double_interval_censored
using Distributions: Gamma, cdf, pdf
using Mooncake: Mooncake
using Statistics: median
using Printf: @printf

## Current form (mirrors src/renewal.jl `_pmf_from_dic`).
function pmf_current(dist, nmax::Integer)
    dic = double_interval_censored(dist; interval = 1.0, upper = float(nmax))
    c = [cdf(dic, float(b)) for b in 0:(nmax + 1)]
    z0 = zero(eltype(c))
    raw = [max(c[i + 1] - c[i], z0) for i in 1:(nmax + 1)]
    s = sum(raw)
    if !isfinite(s) || s <= zero(s)
        z = zero(pdf(dist, oneunit(float(nmax))))
        return fill(one(z) / (nmax + 1), nmax + 1)
    end
    return raw ./ s
end

## Fused form: one loop, scalar `prev` carry, single output array.
function pmf_fused(dist, nmax::Integer)
    dic = double_interval_censored(dist; interval = 1.0, upper = float(nmax))
    prev = cdf(dic, 0.0)
    z0 = zero(prev)
    raw = Vector{typeof(prev)}(undef, nmax + 1)
    @inbounds for i in 1:(nmax + 1)
        cur = cdf(dic, float(i))
        raw[i] = max(cur - prev, z0)
        prev = cur
    end
    s = sum(raw)
    if !isfinite(s) || s <= zero(s)
        z = zero(pdf(dist, oneunit(float(nmax))))
        return fill(one(z) / (nmax + 1), nmax + 1)
    end
    return raw ./ s
end

function pmf_mean(pmf)
    s = zero(eltype(pmf))
    @inbounds for i in eachindex(pmf)
        s += (i - 1) * pmf[i]
    end
    return s
end

function time_gradient(f, args...; reps = 80)
    cache = Mooncake.prepare_gradient_cache(f, args...)
    Mooncake.value_and_gradient!!(cache, f, args...)
    times = Float64[]
    for _ in 1:reps
        push!(times, @elapsed Mooncake.value_and_gradient!!(cache, f, args...))
    end
    return median(times)
end

function run_case(label, build, a, b, nmaxes)
    println(label)
    for nmax in nmaxes
        f_cur = (m, s) -> pmf_mean(pmf_current(build(m, s), nmax))
        f_fus = (m, s) -> pmf_mean(pmf_fused(build(m, s), nmax))
        v_cur = f_cur(a, b)
        v_fus = f_fus(a, b)
        t_cur = time_gradient(f_cur, a, b)
        t_fus = time_gradient(f_fus, a, b)
        @printf("  nmax=%-3d  current %8.2f µs  fused %8.2f µs  (%.2fx)  Δval=%.1e\n",
            nmax, t_cur * 1e6, t_fus * 1e6, t_cur / t_fus, abs(v_cur - v_fus))
    end
    println()
    return nothing
end

run_case("LogNormal incubation (6.3, 3.5):", lognormal_meansd, 6.3, 3.5,
    (16, 40))
run_case("Gamma onset->report (1.178, 3.694):", Gamma, 1.178, 3.694, (16, 31))
run_case("Gamma generation interval (2.71, 5.65):", Gamma, 2.71, 5.65, (40,))
run_case("LogNormal recovery (14, 8):", lognormal_meansd, 14.0, 8.0, (42,))
