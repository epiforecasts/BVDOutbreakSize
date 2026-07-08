## The headline confirmation-sensitivity prior credits the repeat-control
## confirmation process (Beta(38, 2), mean 0.95) rather than one analytical
## assay draw (Beta(10, 1.76), mean 0.85); see issue #374.

@testitem "headline process prior is higher and tighter than single-assay" tags=[:slow] begin
    using Turing: sample, Prior
    using Random: MersenneTwister
    using Statistics: mean, std
    using Distributions: Beta
    using BVDOutbreakSize: test_sensitivity_model

    headline = sample(MersenneTwister(1), test_sensitivity_model(),
        Prior(), 4_000; progress = false)
    single = sample(MersenneTwister(1),
        test_sensitivity_model(sensitivity_prior = Beta(10.0, 1.76)),
        Prior(), 4_000; progress = false)
    s_headline = vec(Array(headline[:s_test]))
    s_single = vec(Array(single[:s_test]))
    @test mean(s_headline) > mean(s_single)
    @test std(s_headline) < std(s_single)
    @test mean(s_headline) > 0.9
    @test all(0 .< s_headline .< 1)
end
