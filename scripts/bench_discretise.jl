# Mooncake gradient-timing diagnostic for the censored-delay discretisation.
#
# `discretise_censored` → `_pmf_from_dic` (src/renewal.jl) is the documented
# per-gradient hot path: every delay submodel discretises a
# `double_interval_censored` distribution once per draw, and the reverse pass
# walks the censored-CDF evaluations. The PMF was built as
# `[pdf(dic, d) for d in 0:nmax]`, which evaluates every interior integer
# boundary CDF TWICE (lag `d` reads `cdf(d)` and `cdf(d+1)`, lag `d+1` reads
# `cdf(d+1)` and `cdf(d+2)`, …). Differencing one CDF path over `0:nmax+1`
# evaluates each boundary once and halves those evaluations, numerically
# identical to the `pdf` differences it replaces.
#
# This times the OLD pdf-loop against the CURRENT `discretise_censored`
# (cdf-difference) under Mooncake, for the delays the model uses, confirming
# the values match and reporting the gradient-cost delta.
#
# Run: julia --project=. scripts/bench_discretise.jl

using BVDOutbreakSize: discretise_censored, lognormal_meansd
using CensoredDistributions: double_interval_censored
using Distributions: Gamma, cdf, pdf
using Mooncake: Mooncake
using Statistics: median
using Printf: @printf

## Old implementation: a `pdf` per lag, overlapping boundary CDFs.
function pmf_pdf_loop(dist, nmax::Integer)
    dic = double_interval_censored(dist; interval = 1.0, upper = float(nmax))
    raw = [pdf(dic, float(d)) for d in 0:nmax]
    return raw ./ sum(raw)
end

## Objective: PMF mean (a non-trivial functional of the censored CDF path) so
## the gradient w.r.t. the delay parameters is exercised in full.
function pmf_mean(pmf)
    s = zero(eltype(pmf))
    @inbounds for i in eachindex(pmf)
        s += (i - 1) * pmf[i]
    end
    return s
end

function time_gradient(f, args...; reps = 50)
    cache = Mooncake.prepare_gradient_cache(f, args...)
    Mooncake.value_and_gradient!!(cache, f, args...)        # warmup / compile
    times = Float64[]
    for _ in 1:reps
        push!(times, @elapsed Mooncake.value_and_gradient!!(cache, f, args...))
    end
    return median(times)
end

function run_case(label, build, a, b, nmaxes)
    println(label)
    for nmax in nmaxes
        f_old = (m, s) -> pmf_mean(pmf_pdf_loop(build(m, s), nmax))
        f_new = (m, s) -> pmf_mean(discretise_censored(build(m, s), nmax))
        v_old = f_old(a, b)
        v_new = f_new(a, b)
        t_old = time_gradient(f_old, a, b)
        t_new = time_gradient(f_new, a, b)
        @printf("  nmax=%-3d  pdf-loop %8.2f µs  cdf-diff %8.2f µs  (%.2fx)  Δval=%.1e\n",
            nmax, t_old * 1e6, t_new * 1e6, t_old / t_new, abs(v_old - v_new))
    end
    println()
    return nothing
end

run_case("LogNormal(mean=6.3, sd=3.5) [incubation]:",
    lognormal_meansd, 6.3, 3.5, (20, 30))
run_case("LogNormal(mean=4.5, sd=4.0) [lab receipt / ruleout]:",
    lognormal_meansd, 4.5, 4.0, (20, 30))
run_case("Gamma(1.178, 3.694) [onset->report / detection]:",
    Gamma, 1.178, 3.694, (20, 30))
run_case("Gamma(3.33, 3.83) [onset->death atomic]:",
    Gamma, 3.33, 3.83, (40, 60))
