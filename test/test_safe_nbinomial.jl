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
    @test mean(d) ≈ 10.0 rtol = 1e-6
end

@testitem "safe_nbinomial: zero mean is handled" begin
    using BVDOutbreakSize: safe_nbinomial
    using Distributions: succprob

    d = safe_nbinomial(3.0, 0.0)
    @test 0 < succprob(d) < 1
end
