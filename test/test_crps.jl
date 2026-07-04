## Tests for the sample-based CRPS proper score (`crps_sample`), matching the
## `scoringutils::crps_sample` energy-form estimator.

@testitem "crps_sample matches known closed-form values" begin
    using BVDOutbreakSize: crps_sample

    ## A deterministic forecast (every draw equal) reduces to the absolute
    ## error between the point forecast and the observation.
    @test crps_sample(1.0, fill(5.0, 7)) ≈ 4.0
    @test crps_sample(5.0, fill(5.0, 7)) ≈ 0.0

    ## Two-point predictive {0, 2}, observation at 1: mean|x−y| = 1, and the
    ## mean pairwise |x_i − x_j| = 1, so CRPS = 1 − 0.5·1 = 0.5. Cross-checked
    ## against the definition E|X−y| − ½E|X−X'|.
    @test crps_sample(1.0, [0.0, 2.0]) ≈ 0.5
    ## The observation sitting on one of the two support points.
    @test crps_sample(0.0, [0.0, 2.0]) ≈ 1.0 - 0.5 * 1.0

    ## Empty predictive sample is undefined, mirroring `bias_sample`.
    @test isnan(crps_sample(1.0, Float64[]))
end

@testitem "crps_sample is non-negative and order-invariant" begin
    using BVDOutbreakSize: crps_sample
    using Random: MersenneTwister, randperm

    rng = MersenneTwister(11)
    for _ in 1:50
        x = randn(rng, 40)
        y = randn(rng)
        c = crps_sample(y, x)
        @test c >= -1e-9
        ## The sorted estimator must not depend on the input ordering.
        @test crps_sample(y, x) ≈ crps_sample(y, reverse(x))
        @test crps_sample(y, x) ≈ crps_sample(y, x[randperm(rng, length(x))])
    end
end

@testitem "crps_sample sorted identity equals the brute-force double sum" begin
    using BVDOutbreakSize: crps_sample
    using Random: MersenneTwister

    ## Independent O(m^2) reference implementation of the energy-form estimator
    ## CRPS = mean|x_i − y| − (1/2m^2) ΣΣ|x_i − x_j|.
    function crps_bruteforce(y, x)
        m = length(x)
        mae = sum(abs.(x .- y)) / m
        pair = 0.0
        for xi in x, xj in x
            pair += abs(xi - xj)
        end
        return mae - pair / (2 * m^2)
    end

    rng = MersenneTwister(23)
    for _ in 1:100
        m = rand(rng, 2:60)
        x = 100 .* rand(rng, m)
        y = 100 * rand(rng)
        @test crps_sample(y, x) ≈ crps_bruteforce(y, x) atol = 1e-8
    end
end

@testitem "crps_sample rewards a sharper, well-centred forecast" begin
    using BVDOutbreakSize: crps_sample
    using Random: MersenneTwister

    rng = MersenneTwister(7)
    y = 0.0
    sharp = 0.3 .* randn(rng, 5_000)      # tight, centred on the observation
    diffuse = 3.0 .* randn(rng, 5_000)    # wide, same centre
    biased = 5.0 .+ 0.3 .* randn(rng, 5_000)  # tight but off-centre
    ## A strictly proper score prefers the sharp, calibrated forecast over both
    ## a diffuse one and a sharp-but-biased one.
    @test crps_sample(y, sharp) < crps_sample(y, diffuse)
    @test crps_sample(y, sharp) < crps_sample(y, biased)
end
