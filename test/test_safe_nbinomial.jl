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
        vintage_increments_model, censored_occupancy_model
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
    ## No vintages: no variable.
    @test isempty(keyset(vintage_increments_model(Float64[], missing, k)))
    @test isempty(keyset(censored_occupancy_model(Float64[], Float64[], missing, k)))
end
