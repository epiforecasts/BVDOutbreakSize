## Count draws past `typemax(Int)`. A forecast from a loosely constrained fit
## can put a count's mean above what an `Int` holds, where the Distributions
## samplers throw `InexactError` on converting the draw.

@testitem "SafePoisson: an extreme mean draws a saturated count" begin
    using BVDOutbreakSize: SafePoisson
    using Distributions: Poisson
    using Random: Xoshiro
    ## The mean the exports forecast reached in #893.
    λ = 9.480115359493294e18
    @test_throws InexactError rand(Xoshiro(1), Poisson(λ))
    @test rand(Xoshiro(1), SafePoisson(λ)) == typemax(Int)
    x = rand(Xoshiro(1), SafePoisson(2.0^62))
    @test x isa Int
    @test isapprox(x, 2.0^62; rtol = 1.0e-6)
end

@testitem "SafePoisson: Poisson draws and density within range" begin
    using BVDOutbreakSize: SafePoisson
    using Distributions: Poisson, logpdf
    using Random: Xoshiro
    for λ in (0.0, 0.3, 5.9, 6.1, 250.0, 1.0e9)
        @test rand(Xoshiro(2), SafePoisson(λ)) == rand(Xoshiro(2), Poisson(λ))
        for x in (0, 3, 7, 300)
            @test logpdf(SafePoisson(λ), x) == logpdf(Poisson(λ), x)
        end
    end
end

@testitem "count vectors: an extreme mean draws a saturated count" begin
    using BVDOutbreakSize: NegBinomialVector, safe_nbinomial
    using Distributions: censored
    using Random: Xoshiro
    ## The success probability is floored at `eps`, so the mean reaches past
    ## `typemax(Int)` only with a large dispersion.
    k = 1.0e4
    μ = [3.0, 1.0e22, 40.0]
    @test_throws InexactError rand(Xoshiro(3), safe_nbinomial(k, 1.0e22))
    x = rand(Xoshiro(3), NegBinomialVector(k, μ))
    @test x isa Vector{Int}
    @test x[2] == typemax(Int)
    y = rand(Xoshiro(3), censored(NegBinomialVector(k, μ); upper = fill(1.0e6, 3)))
    @test y[2] == 1.0e6
end

@testitem "count vectors: draws within range match the scalar sampler" begin
    using BVDOutbreakSize: NegBinomialVector, safe_nbinomial
    using Distributions: censored
    using Random: Xoshiro
    μ = [0.0, 0.4, 8.0, 3000.0, 2.0e7]
    x = rand(Xoshiro(4), NegBinomialVector(2.5, μ))
    rng = Xoshiro(4)
    @test x == [rand(rng, safe_nbinomial(2.5, m)) for m in μ]
    up = [10.0, 10.0, 10.0, 10.0, 10.0]
    y = rand(Xoshiro(5), censored(NegBinomialVector(2.5, μ); upper = up))
    rng = Xoshiro(5)
    @test y == [
        rand(rng, censored(safe_nbinomial(2.5, m); upper = u))
            for (m, u) in zip(μ, up)
    ]
end

@testitem "dated exports forecast from an extreme mean" begin
    using BVDOutbreakSize: dated_poisson_model
    using Random: Xoshiro
    draw = rand(Xoshiro(6), dated_poisson_model([2.0, 9.480115359493294e18], missing))
    counts = [draw[k] for k in keys(draw)]
    @test counts[2] == typemax(Int)
    @test 0 <= counts[1] < 100
end
