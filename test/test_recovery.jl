## Parameter recovery helpers, on a toy model with the same shape as the
## joint: a parameter, a deterministic, and a submodel observation that
## samples when its data are `missing` (the generator path).
@testsnippet RecoveryToy begin
    using BVDOutbreakSize, Turing, Distributions, DataFrames
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

@testitem "simulate_recovery keeps the truth with the data" setup = [RecoveryToy] begin
    gen = toy(missing)
    observed = recovery_observed_varnames(gen, toy([0.1, 0.2, 0.3]))
    sim = simulate_recovery(gen, observed; seed = 3, horizon = 0)
    mu = only(BVDOutbreakSize._draws(sim.truth, :mu))
    @test only(BVDOutbreakSize._draws(sim.truth, :total)) ≈ 2 * mu
    @test length(sim.data) == 3
    @test all(v isa Real for v in values(sim.data))
    ## A fixed seed gives the same dataset.
    again = simulate_recovery(gen, observed; seed = 3, horizon = 0)
    @test only(BVDOutbreakSize._draws(again.truth, :mu)) == mu
    ## An acceptance check that never passes errors rather than looping.
    @test_throws ErrorException simulate_recovery(
        gen, observed; seed = 3, horizon = 0, accept = _ -> false,
        attempts = 2
    )
end

@testitem "simulated observations group by stream" begin
    using BVDOutbreakSize
    g = BVDOutbreakSize._grouped_observations(
        Dict(
            "cases_state.reported_increments.increments[2]" => 5,
            "cases_state.reported_increments.increments[1]" => 3,
            "confirmed_state.confirmed_positives.positives[1]" => 7,
            "treatment_state.bed_capacity.increments" => [4, 6, 8],
            "onset_report_state.increments[2]" => 1.5,
            "onset_report_state.increments[1]" => 0.5,
            "composition_state.obs_increments[2, :]" => [1, 2],
            "composition_state.obs_increments[1, :]" => [3, 4],
        )
    )
    @test g.in_order(g.streams[:cases_state][:reported_increments]) == [3, 5]
    @test g.in_order(g.streams[:confirmed_state][:confirmed_positives]) == [7]
    ## A stream drawn as one vector keeps its order.
    @test g.in_order(g.streams[:treatment_state][:bed_capacity]) == [4, 6, 8]
    @test g.onsets == [0.5, 1.5]
    ## The onset triangle drawn as one vector.
    whole = BVDOutbreakSize._grouped_observations(
        Dict("onset_report_state.increments" => [2.0, 3.0])
    )
    @test whole.onsets == [2.0, 3.0]
    rows = g.in_order(g.rows[:composition_state])
    @test rows == [[3, 4], [1, 2]]
    ## The last province is the remainder of the recorded totals.
    m = BVDOutbreakSize._composition_matrix(rows, [10.0, 9.0])
    @test m == [3 4; 1 2; 6 3]
    @test vec(sum(m; dims = 1)) == [10, 9]
    ## An observation with no rule is an error, not a silent drop.
    @test_throws ArgumentError BVDOutbreakSize._grouped_observations(
        Dict("oddly_named" => 1)
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
    v = recovery_verdict(tab)
    @test v.pass && v.status == :pass && v.converged
    bad = recovery_table(Dict("a" => 4.0, "b" => 0.0), draws)
    v = recovery_verdict(bad)
    @test v.status == :fail && !v.pass
    @test v.outside == ["a"]
    @test v.coverage_90 == 0.5
    ## Most quantities outside their 90% interval fails even with none
    ## outside the 99% interval.
    edge = recovery_table(Dict("a" => 1.9, "b" => -1.9), draws)
    @test isempty(recovery_verdict(edge).outside)
    @test recovery_verdict(edge).status == :fail
    ## Coverage between the two bars warns but passes.
    five = Dict("q$i" => randn(rng, 4000) for i in 1:5)
    mid = recovery_table(
        Dict("q1" => 0.0, "q2" => 0.0, "q3" => 0.0, "q4" => 1.8, "q5" => -1.8),
        five
    )
    @test recovery_verdict(mid).status == :warn
    @test recovery_verdict(mid).pass
    ## A fit that did not converge is not judged, whatever its intervals.
    stuck = (; max_rhat = 1.3, min_ess_bulk = 8.0)
    v = recovery_verdict(tab; diagnostics = stuck)
    @test v.status == :unconverged && !v.pass && !v.converged
    fine = (; max_rhat = 1.01, min_ess_bulk = 400.0)
    @test recovery_verdict(tab; diagnostics = fine).status == :pass
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

@testitem "recovery summaries across seeds" begin
    using BVDOutbreakSize, DataFrames
    ## Two seeds of two quantities: a growth rate with a negative truth and a
    ## count, with the count's truth outside the 99% interval in seed 2.
    row(seed, q, truth, med, tq, z; covered = true) = (;
        quantity = q, truth, median = med, lower_90 = med - 1,
        upper_90 = med + 1, lower_50 = med - 0.5, upper_50 = med + 0.5,
        covered_50 = covered, covered_90 = covered, truth_quantile = tq, z,
        seed, fit_minutes = 10.0, max_rhat = seed == 1 ? 1.01 : 1.2,
        min_ess_bulk = 400.0, divergences = 0,
    )
    params = DataFrame(
        [
            row(1, "r", -0.02, -0.01, 0.4, 0.3),
            row(1, "C_T", 100.0, 120.0, 0.3, 0.5),
            row(2, "r", -0.04, -0.05, 0.6, -0.4),
            row(2, "C_T", 100.0, 40.0, 0.999, -3.0; covered = false),
        ]
    )
    s = recovery_summary(params)
    @test s.quantity == ["r", "C_T"]
    ## Relative error is taken against the truth's magnitude, so a
    ## negative truth keeps the sign of the error.
    @test s.rel_error_min[1] ≈ -0.25
    @test s.rel_error_max[1] ≈ 0.5
    @test s.rel_error_median[2] ≈ (0.2 - 0.6) / 2
    @test s.covered_90 == [2, 1]
    @test s.outside == [0, 1]
    @test s.worst_seed == [2, 2]
    @test s.worst_z[2] == -3.0
    @test (s.truth_min[1], s.truth_max[1]) == (-0.04, -0.02)
    ## A zero truth has no relative error.
    zero_truth = copy(params)
    zero_truth.truth[zero_truth.quantity .== "r"] .= 0.0
    z = recovery_summary(zero_truth)
    @test ismissing(z.rel_error_median[1])
    @test !ismissing(z.rel_error_median[2])

    v = recovery_seed_verdicts(params)
    @test v.seed == [1, 2]
    @test v.status == [:pass, :unconverged]
    @test v.outside[1] == ""
    @test recovery_overall(v.status) == :unconverged
    @test recovery_overall([:pass, :fail, :unconverged]) == :fail
    @test recovery_overall([:pass, :warn]) == :warn
    @test recovery_overall([:pass]) == :pass
end
