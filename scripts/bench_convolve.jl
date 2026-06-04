# Mooncake gradient-timing diagnostic for the renewal per-draw hot path.
#
# Two questions, each isolated to a pure function so the gradient is the
# real model cost without the Turing machinery:
#
#  1. Is `convolve_delay` worth restructuring? A vectorised lag-AXPY form
#     was tried against the scalar double loop. Mooncake differentiates the
#     scalar loop about as fast (≈parity), so the convolution is NOT the
#     bottleneck and the simple loop stays.
#  2. Where does the per-gradient cost actually sit? The delay discretisation
#     (`discretise_censored`, censored-distribution CDF evaluations over
#     0:nmax lags) dominates, and its cost scales with `nmax` — which is why
#     trimming the oversized onset-to-death `nmax` (60→40) is the real win.
#
# Run: julia --project=. scripts/bench_convolve.jl

using BVDOutbreakSize
using BVDOutbreakSize: discretise_censored, lognormal_meansd
using Mooncake: Mooncake
using Statistics: median
using Printf: @printf

## --- Convolution: scalar loop vs vectorised lag-AXPY -------------------
function conv_old(x::AbstractVector, delay::AbstractVector)
    n = length(x)
    Tp = promote_type(eltype(x), eltype(delay))
    y = zeros(Tp, n)
    @inbounds for t in 1:n
        acc = zero(Tp)
        dmax = min(t - 1, length(delay) - 1)
        for d in 0:dmax
            acc += x[t - d] * delay[d + 1]
        end
        y[t] = acc
    end
    return y
end

function conv_new(x::AbstractVector, delay::AbstractVector)
    n = length(x)
    L = length(delay)
    Tp = promote_type(eltype(x), eltype(delay))
    y = zeros(Tp, n)
    @inbounds for d in 0:(L - 1)
        w = delay[d + 1]
        @views y[(d + 1):n] .+= w .* x[1:(n - d)]
    end
    return y
end

function make_conv_objective(conv, kernels)
    return function (x)
        s = x
        acc = zero(eltype(x))
        for k in kernels
            s = conv(s, k)
            acc += sum(s)
        end
        return acc
    end
end

## Single Mooncake gradient wall-clock, median of `reps` after a warmup.
function time_gradient(f, args...; reps = 50)
    cache = Mooncake.prepare_gradient_cache(f, args...)
    Mooncake.value_and_gradient!!(cache, f, args...)        # warmup / compile
    times = Float64[]
    for _ in 1:reps
        push!(times, @elapsed Mooncake.value_and_gradient!!(cache, f, args...))
    end
    return median(times)
end

function bench_convolution()
    n = 93
    klens = (30, 40, 30, 30, 30)
    kernels = [(w = abs.(sin.(1.0:l)) .+ 0.1; w ./ sum(w))
               for l in klens]
    x = abs.(cos.(1.0:n)) .* 10 .+ 1
    f_old = make_conv_objective(conv_old, kernels)
    f_new = make_conv_objective(conv_new, kernels)
    @assert f_old(x) ≈ f_new(x)
    t_old = time_gradient(f_old, x)
    t_new = time_gradient(f_new, x)
    println("1. convolve_delay gradient, 5-convolution chain (n=$n):")
    @printf("   scalar double loop  : %8.2f µs\n", t_old * 1e6)
    @printf("   vectorised lag-AXPY : %8.2f µs  (%.2fx)\n",
        t_new * 1e6, t_old / t_new)
    println("   → no speedup; the scalar loop is kept.\n")
    return nothing
end

## --- Discretisation: gradient cost as a function of nmax ----------------
## Objective is the PMF mean (a non-trivial functional of the censored CDF
## path) so the gradient w.r.t. the delay mean/sd is exercised in full.
function pmf_mean(mean, sd, nmax)
    pmf = discretise_censored(lognormal_meansd(mean, sd), nmax)
    s = zero(eltype(pmf))
    @inbounds for i in eachindex(pmf)
        s += (i - 1) * pmf[i]
    end
    return s
end

function bench_discretisation()
    println("2. discretise_censored gradient (∂ PMF-mean / ∂(mean, sd)):")
    base = nothing
    for nmax in (40, 60)
        f = (m, s) -> pmf_mean(m, s, nmax)
        t = time_gradient(f, 11.2, 5.4)
        base === nothing && (base = t)
        tag = nmax == 40 ? "(onset→death, trimmed)" :
              nmax == 60 ? "(onset→death, original)" : ""
        @printf("   nmax=%-3d : %8.2f µs  %s\n", nmax, t * 1e6, tag)
    end
    println("   → the 60→40 trim removes ~1/3 of the per-delay CDF evals.")
    return nothing
end

bench_convolution()
bench_discretisation()
