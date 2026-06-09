## Tests for the genetic seeding submodel from `src/models/priors.jl`.
## The submodel encodes the molecular-clock TMRCA as an upper-censored,
## noisy reading of the seeding time `T`: a soft one-sided bound, flat
## above the TMRCA and decaying below it.

@testitem "genetic_seeding adds the soft lower-bound log density" begin
    using Distributions: Normal, logcdf
    using Turing: logjoint
    using BVDOutbreakSize: genetic_seeding_model

    g, sd = 80.0, 20.0
    at(T) = logjoint(genetic_seeding_model(T, g; tmrca_days_sd = sd), (;))

    ## Matches the analytic one-sided (censored) likelihood P(read ≥ g).
    @test at(112.0) ≈ logcdf(Normal(0.0, sd), 112.0 - g)
    @test at(56.0) ≈ logcdf(Normal(0.0, sd), 56.0 - g)

    ## Monotone increasing in T: earlier seeding is penalised more.
    @test at(40.0) < at(80.0) < at(160.0)

    ## Flat well above the bound, steep decay well below it.
    @test at(200.0) - at(160.0) < at(40.0) - at(20.0)
end

@testitem "genetic_seeding: missing tmrca_days is a no-op" begin
    using Turing: logjoint
    using BVDOutbreakSize: genetic_seeding_model

    ## When tmrca_days = missing the model adds zero log-probability.
    lp = logjoint(genetic_seeding_model(100.0, missing; tmrca_days_sd = 15.0),
        (;))
    @test lp == 0.0
end

@testitem "genetic_seeding composes into bvd_joint via `genetic` kwarg" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, genetic_seeding_model

    seed = genetic_seeding_model
    chn = sample(
        bvd_joint(40, 2, 18;
            tmrca_days = 30.0,
            tmrca_days_sd = 15.0,
            genetic = seed),
        Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false
    )
    T_draws = vec(Array(chn[:T]))
    @test length(T_draws) == 50
    @test all(isfinite, T_draws)
end
