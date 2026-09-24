## Textbook loops for the kernels whose package bodies are written for speed:
## BLAS calls, a shared onset column walk, a shared occupancy forward and
## summed likelihoods. Each is the formula as a plain loop, with no guards
## beyond what the formula needs, so a faster body in `src/` cannot change
## the model's values unnoticed. Included by `test/test_reference_kernels.jl`
## and never by the package.

using Distributions: BetaBinomial, NegativeBinomial, TDist, logpdf
using StatsFuns: logistic

## `y[t] = Σ_{d = 0}^{t − 1} w[d + 1] x[t − d]`.
function ref_convolve_delay(x, w)
    n = length(x)
    y = zeros(n)
    for t in 1:n, d in 0:min(t - 1, length(w) - 1)
        y[t] += w[d + 1] * x[t - d]
    end
    return y
end

## `(a ⊕ b)[k] = Σ_{i + j = k + 1} a[i] b[j]`.
function ref_convolve_pmf(a, b)
    (isempty(a) || isempty(b)) && return Float64[]
    y = zeros(length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        y[i + j - 1] += a[i] * b[j]
    end
    return y
end

## `I_t` is the seed for `t ≤ L`, then `R_t f_t` with the force
## `f_t = Σ_{s = 1}^{min(t − 1, G)} I_{t − s} g_s`.
function ref_renewal_infections_with_force(Rt, g, seed)
    n = length(Rt)
    I = zeros(n)
    force = zeros(n)
    for t in 1:n
        if t <= length(seed)
            I[t] = seed[t]
        else
            for s in 1:min(t - 1, length(g))
                force[t] += I[t - s] * g[s]
            end
            I[t] = Rt[t] * force[t]
        end
    end
    return I, force
end

## Patch `p` generates `G[p, t] = R[p, t] Σ_{s = 1}^{min(t − 1, G)} I[p, t − s] g[s]`
## after its seed. It keeps `(1 − ε[p, t] Σ_{r ≠ p} K[r, p]) G[p, t]` and
## receives the arrivals `A[p, t] = Σ_{q ≠ p} ε[q, t] K[p, q] G[q, t]`.
function ref_patch_infections(R, g, seeds, K, ε)
    np, n = size(R)
    I = zeros(np, n)
    A = zeros(np, n)
    I[:, 1:min(size(seeds, 2), n)] = seeds[:, 1:min(size(seeds, 2), n)]
    for t in (size(seeds, 2) + 1):n
        G = [
            R[p, t] * sum(I[p, t - s] * g[s] for s in 1:min(t - 1, length(g)))
                for p in 1:np
        ]
        for p in 1:np
            sent = sum(K[r, p] for r in 1:np if r != p; init = 0.0)
            A[p, t] = sum(
                ε[q, t] * K[p, q] * G[q] for q in 1:np if q != p; init = 0.0
            )
            I[p, t] = (1 - ε[p, t] * sent) * G[p] + A[p, t]
        end
    end
    return (; infections = I, importation = A)
end

## The occupancy balance of `accumulate_occupancy`'s docstring, day by day.
function ref_accumulate_occupancy(A_bvd, A_bg, deaths, recover, ruleout, κ, h)
    n = length(A_bvd)
    out = (;
        demand = zeros(n), O_bvd = zeros(n), O_conf = zeros(n),
        O_susp = zeros(n), abscond = zeros(n),
    )
    bvd, bg, conf, susp = 0.0, 0.0, 0.0, 0.0
    for t in 1:n
        discharged = deaths[t] + recover[t]
        ab = κ * susp
        unconf = max(bvd - conf, 0.0)
        denom = max(susp, eps())
        share = bvd > 0 ? conf / bvd : 0.0
        bvd_t = max(bvd + A_bvd[t] - discharged - ab * unconf / denom, 0.0)
        bg_t = max(bg + A_bg[t] - ruleout[t] - ab * bg / denom, 0.0)
        conf_t = clamp(conf + h[t] * unconf - discharged * share, 0.0, bvd_t)
        susp_t = max(bvd_t + bg_t - conf_t, 0.0)
        out.demand[t] = bvd_t + bg_t
        out.O_bvd[t] = bvd_t
        out.O_conf[t] = conf_t
        out.O_susp[t] = susp_t
        out.abscond[t] = ab
        bvd, bg, conf, susp = bvd_t, bg_t, conf_t, susp_t
    end
    return out
end

## `c = min(max(x, 0), O_bvd)`, `O_susp = max(D − c, 0)`, absconds
## `κ O_susp(t − 1)` from day 2, total `D + Δ` and suspect census
## `max(D + Δ − c, 0)`.
function ref_incare_census(D, O_bvd, x, κ, Δ)
    c = min.(max.(x, 0.0), O_bvd)
    O_susp = max.(D .- c, 0.0)
    return (;
        confirmed = c, suspect = max.(D .+ Δ .- c, 0.0),
        abscond = [t == 1 ? 0.0 : κ * O_susp[t - 1] for t in eachindex(D)],
        total = D .+ Δ,
    )
end

## Reporting hazard at delay `j` for onset date `u`, with the calendar index
## held at the nearest edge of `γ`.
ref_onset_hazard(lh, γ, gs, u, j) = logistic(
    lh[j + 1] + γ[clamp(u + j - gs + 1, 1, length(γ))]
)

## `table[d + 1, k] = 1 − Π_{j = 0}^{d} (1 − h(u, j))` for onset date
## `u = u_lo + k − 1`.
function ref_onset_report_cdf_table(lh, γ, gs, u_lo, u_hi)
    D = length(lh)
    table = zeros(D, max(u_hi - u_lo + 1, 0))
    for k in axes(table, 2), d in 0:(D - 1)
        u = u_lo + k - 1
        surv = prod(1 - ref_onset_hazard(lh, γ, gs, u, j) for j in 0:d)
        table[d + 1, k] = 1 - surv
    end
    return table
end

## `Σ_{u = 1}^{min(as_of, n)} onsets[u] α(u) G(u, as_of − u)` with the
## normalised delay CDF `G(u, δ) = F(u, min(δ, D − 1)) / F(u, D − 1)`.
function ref_onset_report_expected_total(onsets, lh, γ, gs, alpha, as_of)
    D = length(lh)
    total = 0.0
    for u in 1:min(as_of, length(onsets))
        δ = as_of - u
        F = ref_onset_report_cdf_table(lh, γ, gs, u, u)
        α = alpha[clamp(u - gs + 1, 1, length(alpha))]
        total += onsets[u] * α * F[min(δ, D - 1) + 1] / F[D]
    end
    return total
end

## Cell `i` reads scan `s_i` now and scan `p_i` before, each at that scan's
## level, or one for an index naming no scan. Its mean is
## `ℓ_cur c_s − ℓ_prev c_p` and its scale `sqrt(max(mean, 0) + pixel_sd² r)`
## over `r = 2` reads, or one read when there is no previous report.
function ref_onset_scanned_cells(lc, lp, c, s, p, prev_report, pixel_sd)
    level(k) = 1 <= k <= length(c) ? c[k] : 1.0
    means = [lc[i] * level(s[i]) - lp[i] * level(p[i]) for i in eachindex(lc)]
    reads = [prev_report[i] > 0 ? 2 : 1 for i in eachindex(lc)]
    scales = sqrt.(max.(means, 0) .+ pixel_sd^2 .* reads)
    return (; means, scales)
end

## One `logpdf` per count, with in-range parameters.
function ref_nbinomial_loglik(k, μ, x)
    return sum(
        logpdf(NegativeBinomial(k, k / (k + μ[i])), x[i])
            for i in eachindex(x);
        init = 0.0
    )
end

function ref_studentt_loglik(μ, σ, x, ν)
    return sum(
        logpdf(TDist(ν), (x[i] - μ[i]) / σ[i]) - log(σ[i])
            for i in eachindex(x);
        init = 0.0
    )
end

## Concentration `c = (1 − ρ) / ρ`, so `α = c p` and `β = c (1 − p)`.
function ref_betabinomial_loglik(n, p, ρ, x)
    c = (1 - ρ) / ρ
    return sum(
        logpdf(BetaBinomial(n[i], c * p[i], c * (1 - p[i])), x[i])
            for i in eachindex(x);
        init = 0.0
    )
end
