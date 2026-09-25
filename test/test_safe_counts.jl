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
    y = rand(
        Xoshiro(3), censored(NegBinomialVector(k, μ); upper = fill(1.0e6, 3))
    )
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
    draw = rand(
        Xoshiro(6), dated_poisson_model([2.0, 9.480115359493294e18], missing)
    )
    counts = [draw[k] for k in keys(draw)]
    @test counts[2] == typemax(Int)
    @test 0 <= counts[1] < 100
end

@testitem "forecast sums saturate rather than wrap" begin
    using BVDOutbreakSize: _saturating_sum
    @test _saturating_sum([typemax(Int), 5, 7]) == typemax(Int)
    @test sum([typemax(Int), 5, 7]) < 0
    @test _saturating_sum([1, 2, 3]) == 6
    @test _saturating_sum([1.5, 2.5]) == 4.0
end

@testitem "every observation draw saturates at an extreme mean" begin
    using BVDOutbreakSize: dated_poisson_model, vintage_increments_model,
        censored_occupancy_model, late_confirmed_model, exports_model,
        exports_deaths_model, composition_split_model, SafePoisson,
        SafeNegBinomial
    using Random: Xoshiro
    ## Each submodel that can sample a missing or forecast count, at a mean
    ## past `typemax(Int)`. Every draw must lie in `[0, typemax(Int)]`. The
    ## censored occupancy takes its float ceiling where it binds.
    big = 1.0e22
    k = 1.0e4
    pmf = [0.2, 0.5, 0.3]
    infections = fill(big, 30)
    models = [
        "dated Poisson" => dated_poisson_model(fill(big, 3), missing),
        "vintage increments" => vintage_increments_model(
            fill(big, 3), missing, k
        ),
        "censored occupancy" => censored_occupancy_model(
            fill(big, 3), fill(1.0e6, 3), missing, k
        ),
        "late confirmed" => late_confirmed_model(
            missing, fill(big, 3), [0, 5, 0], fill(0.1, 3), k
        ),
        "cumulative exports" => exports_model(
            missing, infections, 0.5; incubation_pmf = pmf
        ),
        "dated exports" => exports_model(
            missing, infections, 0.5; incubation_pmf = pmf,
            export_case_days = [10, 20], pre_detection_exports = missing
        ),
        "cumulative export deaths" => exports_deaths_model(
            missing, infections, 0.5, pmf, pmf
        ),
        "dated export deaths" => exports_deaths_model(
            missing, infections, 0.5, pmf, pmf;
            export_death_days = [10, 20], pre_death_exports = missing
        ),
        "composition split" => composition_split_model(
            missing, fill(1 / 3, 3, 2), [typemax(Int), 7], 0.1
        ),
    ]
    for (name, m) in models
        @testset "$name" begin
            draw = rand(Xoshiro(7), m)
            counts = [
                x for (key, v) in pairs(draw)
                    if !occursin(r"state|daily_travellers", string(key))
                    for x in v
            ]
            @test !isempty(counts)
            @test all(x -> isfinite(x) && 0 <= x <= typemax(Int), counts)
        end
    end
    @test eltype(SafePoisson(1.0)) === Int
    @test eltype(SafeNegBinomial(k, 1.0)) === Int
end

@testitem "beta-binomial draws: huge trial counts stay in range" begin
    using BVDOutbreakSize: BetaBinomialVector, safe_betabinomial
    using Random: Xoshiro
    trials = [5, 2^52 - 1, typemax(Int)]
    x = rand(Xoshiro(8), BetaBinomialVector(trials, fill(0.4, 3), 0.1))
    @test all(0 .<= x .<= trials)
    ## Within range the draw is the Distributions one.
    y = rand(Xoshiro(9), BetaBinomialVector([5, 40], [0.3, 0.6], 0.1))
    rng = Xoshiro(9)
    @test y == [
        rand(rng, safe_betabinomial(5, 0.3, 0.1)),
        rand(rng, safe_betabinomial(40, 0.6, 0.1)),
    ]
end
