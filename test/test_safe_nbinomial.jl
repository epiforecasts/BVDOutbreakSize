## Tests for `safe_nbinomial`, the guarded NegativeBinomial used by the
## surveillance likelihoods (`src/models/observations.jl`). A fit can drive
## the dispersion `k` to zero (or a non-finite value); the distribution must
## stay valid rather than throwing `DomainError: r > 0` and aborting the fit.

@testitem "safe_nbinomial: zero/negative/non-finite dispersion is valid" begin
    using BVDOutbreakSize: safe_nbinomial
    using Distributions: NegativeBinomial, params, succprob, mean

    ## An unguarded k = 0 would throw `DomainError: r > 0`.
    for k in (0.0, -1.0, -eps(), NaN, Inf)
        d = safe_nbinomial(k, 5.0)
        @test d isa NegativeBinomial
        @test params(d)[1] > 0          # r floored strictly positive
        @test 0 < succprob(d) < 1
        @test isfinite(mean(d))
    end
end

@testitem "safe_nbinomial: valid dispersion passes through" begin
    using BVDOutbreakSize: safe_nbinomial
    using Distributions: NegativeBinomial, params, mean

    ## A normal dispersion is used as-is, and the mean matches the requested
    ## μ (the NegativeBinomial mean is r(1-p)/p with p = r/(r+μ)).
    d = safe_nbinomial(2.0, 10.0)
    @test params(d)[1] == 2.0
    @test mean(d) ≈ 10.0 rtol = 1.0e-6
end

@testitem "safe_nbinomial: zero mean is handled" begin
    using BVDOutbreakSize: safe_nbinomial
    using Distributions: succprob

    d = safe_nbinomial(3.0, 0.0)
    @test 0 < succprob(d) < 1
end

@testitem "summed NegativeBinomial likelihoods match one term per count" begin
    using BVDOutbreakSize: safe_nbinomial, safe_rate,
        vintage_increments_model, censored_occupancy_model, NegBinomialVector
    using Distributions: logpdf, censored
    using Random: Xoshiro
    using Turing: DynamicPPL
    using Turing.DynamicPPL: @varname

    loglik(m) = DynamicPPL.loglikelihood(m, DynamicPPL.VarInfo(Xoshiro(1), m))
    μ = [12.0, 0.0, 40.5, 7.2, 90.0]
    x = [10, 0, 44, 3, 60]
    k = 6.5

    ## A fully observed vector is one summed term, equal to a `~` per count.
    per_term = sum(logpdf(safe_nbinomial(k, safe_rate(μ[i])), x[i]) for i in 1:5)
    @test loglik(vintage_increments_model(μ, x, k)) ≈ per_term
    @test loglik(vintage_increments_model(Float64[], Int[], k)) == 0

    ## Counts at the ceiling score the censored tail, the rest the pmf.
    ceilings = [1.0e6, 1.0e6, 44.0, 1.0e6, 60.0]
    per_term_c = sum(
        logpdf(
            censored(
                safe_nbinomial(k, safe_rate(μ[i]));
                upper = safe_rate(ceilings[i])
            ), x[i]
        ) for i in 1:5
    )
    @test loglik(censored_occupancy_model(μ, ceilings, x, k)) ≈ per_term_c

    ## A `missing` vector samples as one variable under the whole-vector
    ## key the predictive path reads, with each entry drawn from its own
    ## kernel. A draw above a ceiling returns the ceiling.
    keyset(m) = Set(keys(DynamicPPL.VarInfo(Xoshiro(1), m)))
    @test keyset(vintage_increments_model(μ, missing, k)) ==
        Set([@varname(increments)])
    @test keyset(censored_occupancy_model(μ, ceilings, missing, k)) ==
        Set([@varname(obs)])
    draw = vintage_increments_model(μ, missing, k)(Xoshiro(2)).increments
    @test draw isa Vector{Int} && length(draw) == 5 && all(>=(0), draw)
    low = [100.0, 100.0, 5.0]
    frac = [44.5, 60.0, 1.0e6]
    draw = censored_occupancy_model(low, frac, missing, k)(Xoshiro(2)).obs
    @test draw isa Vector{Float64} && all(draw .<= frac)
    ## The censored vector draws what a `censored` NegativeBinomial per entry
    ## draws from the same stream.
    entries = [
        censored(safe_nbinomial(k, safe_rate(low[i])); upper = frac[i])
            for i in 1:3
    ]
    rng = Xoshiro(3)
    @test rand(Xoshiro(3), censored(NegBinomialVector(k, low); upper = frac)) ==
        [rand(rng, e) for e in entries]
    ## No vintages: no variable.
    @test isempty(keyset(vintage_increments_model(Float64[], missing, k)))
    @test isempty(keyset(censored_occupancy_model(Float64[], Float64[], missing, k)))
end

@testitem "nbinomial_logtail matches the censored logpdf at the ceiling" begin
    using BVDOutbreakSize: nbinomial_logtail, NegBinomialVector
    using Distributions: NegativeBinomial, logpdf, censored

    ## The reference is Distributions' `censored` logpdf, which takes the tail
    ## from Rmath. The grid includes tails too small for a normal float.
    for r in (0.3, 2.0, 8.3, 200.0), μ in (5.0, 40.0, 400.0), u in (1, 40, 300)
        d = NegativeBinomial(r, r / (r + μ))
        @test nbinomial_logtail(d, u) ≈
            logpdf(censored(d; upper = float(u)), u) rtol = 1.0e-10
    end
    d = NegativeBinomial(8.3, 0.2)
    @test nbinomial_logtail(d, 0) == 0
    @test nbinomial_logtail(d, 44.5) ≈
        logpdf(censored(d; upper = 44.5), 44.5) rtol = 1.0e-10
    ## A count above its ceiling has no probability.
    @test logpdf(
        censored(NegBinomialVector(8.3, [30.0]); upper = [40.0]), [41]
    ) == -Inf
end
