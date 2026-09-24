## Parameter recovery helpers, on a toy model with the same shape as the
## joint: a parameter, a deterministic, and a submodel observation that
## samples when its data are `missing` (the generator path).
@testsnippet RecoveryToy begin
    using BVDOutbreakSize, Turing, Distributions, DataFrames
    using Turing.DynamicPPL: VarInfo, getlogjoint, init!!, InitFromVector,
        LogDensityFunction, condition
    @model function toy_obs(mu, y)
        if ismissing(y)
            y = Vector{Union{Missing, Float64}}(missing, 3)
        end
        for i in eachindex(y)
            y[i] ~ Normal(mu, 1)
        end
        return y
    end
    @model function toy(y; forecast = nothing)
        mu ~ Normal(0, 1)
        total := 2 * mu
        obs_state ~ to_submodel(toy_obs(mu, y))
        return total
    end
end

@testitem "recovery_observed_varnames picks the generated data" setup = [RecoveryToy] begin
    observed = recovery_observed_varnames(toy(missing), toy([0.1, 0.2, 0.3]))
    @test sort(string.(observed)) ==
        ["obs_state.y[1]", "obs_state.y[2]", "obs_state.y[3]"]
end

@testitem "simulate_recovery keeps the truth and conditions like the fit" setup = [RecoveryToy] begin
    gen = toy(missing)
    observed = recovery_observed_varnames(gen, toy([0.1, 0.2, 0.3]))
    sim = simulate_recovery(gen, observed; seed = 3, horizon = 0)
    mu = only(BVDOutbreakSize._draws(sim.truth, :mu))
    @test only(BVDOutbreakSize._draws(sim.truth, :total)) ≈ 2 * mu
    y = [sim.data[vn] for vn in sort(collect(keys(sim.data)); by = string)]
    @test length(y) == 3
    ## Conditioning the generator on the simulated data scores the same log
    ## joint as the model built with those data, at any parameter value.
    cond = condition(gen, sim.data)
    direct = toy(Float64.(y))
    for m in (0.0, mu, -1.3)
        @test getlogjoint(last(init!!(cond, VarInfo(), InitFromVector([m], LogDensityFunction(cond))))) ≈
            getlogjoint(last(init!!(direct, VarInfo(), InitFromVector([m], LogDensityFunction(direct)))))
    end
    ## A fixed seed gives the same dataset.
    again = simulate_recovery(gen, observed; seed = 3, horizon = 0)
    @test only(BVDOutbreakSize._draws(again.truth, :mu)) == mu
    ## An acceptance check that never passes errors rather than looping.
    @test_throws ErrorException simulate_recovery(
        gen, observed; seed = 3, horizon = 0, accept = _ -> false,
        attempts = 2
    )
end

@testitem "recovery_table and recovery_verdict" begin
    using BVDOutbreakSize, DataFrames, Random
    rng = MersenneTwister(1)
    draws = Dict("a" => randn(rng, 4000), "b" => randn(rng, 4000))
    tab = recovery_table(Dict("a" => 0.0, "b" => 1.0, "c" => 5.0), draws)
    @test tab.quantity == ["a", "b"]
    @test tab.covered_50 == [true, false]
    @test tab.covered_90 == [true, true]
    @test isapprox(tab.truth_quantile[1], 0.5; atol = 0.03)
    @test isapprox(tab.z[2], -1.0; atol = 0.05)
    @test recovery_verdict(tab).pass
    bad = recovery_table(Dict("a" => 4.0, "b" => 0.0), draws)
    v = recovery_verdict(bad)
    @test !v.pass
    @test v.outside == ["a"]
    @test v.coverage_90 == 0.5
    ## Most quantities outside their 90% interval fails even with none
    ## outside the 99% interval.
    edge = recovery_table(Dict("a" => 1.9, "b" => -1.9), draws)
    @test isempty(recovery_verdict(edge).outside)
    @test !recovery_verdict(edge).pass
end

@testitem "forecast_recovery_table scores against persistence" begin
    using BVDOutbreakSize, DataFrames
    draws = Dict("cases_new" => fill(10.0, 200), "deaths_new" => fill(3.0, 200))
    tab = forecast_recovery_table(
        Dict("cases_new" => 10.0, "deaths_new" => 5.0),
        draws, Dict("cases_new" => 4.0, "deaths_new" => 5.0)
    )
    @test tab.quantity == ["cases_new", "deaths_new"]
    @test tab.crps[1] ≈ 0.0
    @test tab.baseline_crps == [6.0, 0.0]
    @test tab.relative_crps[1] ≈ 0.0
    @test tab.relative_crps[2] == Inf
    @test tab.covered_90 == [true, false]
end
