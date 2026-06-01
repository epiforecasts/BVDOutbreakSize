## Unit tests for the discrete-time renewal primitives in src/renewal.jl.

@testitem "lognormal_meansd: moments match the requested mean and sd" begin
    using Distributions: LogNormal, mean, std
    using BVDOutbreakSize: lognormal_meansd

    for (m, s) in ((5.0, 2.0), (10.0, 3.0), (1.0, 0.5))
        d = lognormal_meansd(m, s)
        @test d isa LogNormal
        @test isapprox(mean(d), m; rtol = 1e-8)
        @test isapprox(std(d), s; rtol = 1e-8)
    end
end

@testitem "discretise_censored: PMF sums to 1 and is non-negative" begin
    using Distributions: LogNormal
    using BVDOutbreakSize: discretise_censored, lognormal_meansd

    d = lognormal_meansd(5.0, 2.0)
    pmf = discretise_censored(d, 30)
    @test length(pmf) == 31
    @test all(pmf .>= 0)
    @test isapprox(sum(pmf), 1.0; atol = 1e-10)
end

@testitem "discretise_censored: PMF contract (non-negative, sums to 1)" begin
    using BVDOutbreakSize: discretise_censored, lognormal_meansd
    using Distributions: LogNormal

    ## The discretised delay is a proper PMF over lags 0…nmax: non-negative
    ## and summing to one, length nmax+1.
    d = lognormal_meansd(2.0, 1.0)
    pmf = discretise_censored(d, 60)
    @test length(pmf) == 61
    @test all(pmf .>= 0)
    @test isapprox(sum(pmf), 1.0; atol = 1e-10)
end

@testitem "euler_lotka_r: round-trips R → r → R" begin
    using BVDOutbreakSize: euler_lotka_r
    using Distributions: LogNormal
    using BVDOutbreakSize: lognormal_meansd, discretise_censored

    ## Build a PMF for the generation interval (drop lag-0 bin).
    gi_raw = discretise_censored(lognormal_meansd(15.3, 9.3), 40)
    g = gi_raw[2:end] ./ sum(gi_raw[2:end])

    for R in (0.8, 1.0, 1.5, 2.0, 3.0)
        r = euler_lotka_r(R, g; steps = 5)
        ## Verify Euler-Lotka identity: R · Σ g_s e^{-r s} = 1
        S = sum(g[s] * exp(-r * s) for s in eachindex(g))
        @test isapprox(R * S, 1.0; rtol = 1e-5)
    end
end

@testitem "euler_lotka_r: r > 0 when R > 1, r < 0 when R < 1" begin
    using BVDOutbreakSize: euler_lotka_r, lognormal_meansd, discretise_censored

    gi_raw = discretise_censored(lognormal_meansd(15.3, 9.3), 40)
    g = gi_raw[2:end] ./ sum(gi_raw[2:end])

    @test euler_lotka_r(1.5, g) > 0
    @test euler_lotka_r(0.8, g) < 0
    ## R = 1 → r ≈ 0
    r_one = euler_lotka_r(1.0, g; steps = 10)
    @test abs(r_one) < 1e-4
end

@testitem "doubling_time: log(2)/r" begin
    using BVDOutbreakSize: doubling_time

    r = log(2) / 14.0
    @test isapprox(doubling_time(r), 14.0; rtol = 1e-10)
    r2 = log(2) / 7.0
    @test isapprox(doubling_time(r2), 7.0; rtol = 1e-10)
end

@testitem "seed_infections: ends at I0, grows exponentially" begin
    using BVDOutbreakSize: seed_infections

    I0 = 5.0
    r = 0.1
    len = 20
    s = seed_infections(I0, r, len)

    @test length(s) == len
    @test isapprox(s[end], I0; rtol = 1e-10)
    ## Check the exponential shape: s[j] = I0 * exp(r*(j - len))
    for j in 1:len
        @test isapprox(s[j], I0 * exp(r * (j - len)); rtol = 1e-10)
    end
end

@testitem "renewal_infections: hand-calculation on a tiny example" begin
    using BVDOutbreakSize: renewal_infections

    ## g = [1.0] (all infectivity at lag 1), R_t = 2 for all t.
    ## seed = [1.0], so I[2] = 2 * I[1] * g[1] = 2,
    ## I[3] = 2 * I[2] * g[1] = 4, etc.
    g = [1.0]
    Rt = fill(2.0, 5)
    seed = [1.0]
    I = renewal_infections(Rt, g, seed)

    @test I[1] ≈ 1.0
    @test I[2] ≈ 2.0
    @test I[3] ≈ 4.0
    @test I[4] ≈ 8.0
    @test I[5] ≈ 16.0
end

@testitem "renewal_infections: grows under R > 1" begin
    using BVDOutbreakSize: renewal_infections, lognormal_meansd,
                           discretise_censored

    gi_raw = discretise_censored(lognormal_meansd(15.3, 9.3), 40)
    g = gi_raw[2:end] ./ sum(gi_raw[2:end])
    n = 60
    seed = ones(length(g))
    Rt = fill(2.0, n)
    I = renewal_infections(Rt, g, seed)
    ## Under R > 1 the trajectory must grow on average.
    @test I[n] > I[1]
    @test all(I .>= 0)
end

@testitem "renewal_infections: declines under R < 1" begin
    using BVDOutbreakSize: renewal_infections, lognormal_meansd,
                           discretise_censored

    gi_raw = discretise_censored(lognormal_meansd(15.3, 9.3), 40)
    g = gi_raw[2:end] ./ sum(gi_raw[2:end])
    n = 60
    seed = fill(10.0, length(g))
    Rt = fill(0.5, n)
    I = renewal_infections(Rt, g, seed)
    ## Under R < 1 the trajectory must eventually fall below seed level.
    @test I[n] < seed[end]
end

@testitem "convolve_delay: hand-calculation" begin
    using BVDOutbreakSize: convolve_delay

    ## x = [1, 2, 3], delay = [0.5, 0.5] (lag 0 and lag 1)
    ## y[1] = 1 * 0.5 = 0.5
    ## y[2] = 2 * 0.5 + 1 * 0.5 = 1.5
    ## y[3] = 3 * 0.5 + 2 * 0.5 = 2.5
    x = [1.0, 2.0, 3.0]
    delay = [0.5, 0.5]
    y = convolve_delay(x, delay)
    @test y ≈ [0.5, 1.5, 2.5]
end

@testitem "convolve_delay: all-lag-0 delay is identity" begin
    using BVDOutbreakSize: convolve_delay

    x = [3.0, 1.0, 4.0, 1.0, 5.0]
    delay = [1.0]
    @test convolve_delay(x, delay) ≈ x
end

@testitem "knot_days: first is 1, last is n, spacing ≤ week" begin
    using BVDOutbreakSize: knot_days

    for n in (7, 14, 20, 40, 60)
        days = knot_days(n; week = 7)
        @test days[1] == 1
        @test days[end] == n
        @test issorted(days)
        @test allunique(days)
        diffs = diff(days)
        @test all(diffs .<= 7)
    end
    ## Single-day grid
    @test knot_days(1) == [1]
end

@testitem "interpolate_knots: linear between knots" begin
    using BVDOutbreakSize: interpolate_knots, knot_days

    ## Two knots: value 0 at day 1, value 10 at day 11.
    ## Midpoint day 6 should be 5.
    vals = [0.0, 10.0]
    days = [1, 11]
    out = interpolate_knots(vals, days, 11)
    @test out[1] ≈ 0.0
    @test out[11] ≈ 10.0
    @test out[6] ≈ 5.0
end

@testitem "interpolate_knots: constant when all knots equal" begin
    using BVDOutbreakSize: interpolate_knots, knot_days

    n = 20
    days = knot_days(n)
    vals = fill(3.14, length(days))
    out = interpolate_knots(vals, days, n)
    @test all(out .≈ 3.14)
end

@testitem "sigmoid_ramp: logistic shape, missing gives zeros" begin
    using BVDOutbreakSize: sigmoid_ramp
    using StatsFuns: logistic

    n = 40
    day = 20
    ramp_val = 7.0
    r = sigmoid_ramp(n, day; ramp = ramp_val)
    @test length(r) == n
    ## Logistic at the breakpoint should be ≈ 0.5
    @test isapprox(r[day], logistic(0.0); atol = 1e-8)
    ## Well before breakpoint should be close to 0
    @test r[1] < 0.1
    ## Well after breakpoint should be close to 1
    @test r[n] > 0.9
    ## Monotone non-decreasing
    @test issorted(r)

    ## Missing day gives all zeros
    @test all(sigmoid_ramp(n, missing) .== 0.0)
end

@testitem "seeding_age: returns n when cumulative never reaches 1" begin
    using BVDOutbreakSize: seeding_age

    n = 50
    ## Cumulative that stays below 1
    cum = fill(0.1, n)
    @test seeding_age(cum, n) == float(n)
end

@testitem "seeding_age: round-trips on a monotone growing series" begin
    using BVDOutbreakSize: seeding_age

    ## Exponential cumulative: cumsum of exp(0.1 * t)
    n = 50
    infections = [exp(0.1 * t) for t in 1:n]
    ## First infection is already > 1, so it crosses at t ≈ 1
    cum = cumsum(infections)
    age = seeding_age(cum, n)
    ## Age should be positive and less than n
    @test age > 0
    @test age <= n
end
