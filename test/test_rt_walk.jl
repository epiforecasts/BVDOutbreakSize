## The national reproduction-number walk samples the total log-Rt at each
## knot after the first, intervention ramp included, and centres the walk on
## those knots net of the ramp. These items pin that structure and check it
## against the standard-normal innovation form it reparameterises.

@testsnippet RtWalkFixtures begin
    using BVDOutbreakSize
    using BVDOutbreakSize: knot_days, sigmoid_ramp, interpolate_knots,
        RT_INTERVENTION_RAMP
    using Turing: @model, to_submodel, returned, logjoint, @varname,
        DynamicPPL
    using Distributions: Normal, Poisson, product_distribution, truncated
    using Random: Xoshiro

    ## Cut-off day 60, the walk from day 12 and the ramp centred on day 30,
    ## so the ramp moves across the knots.
    const N = 60
    const WS = 12
    const BP = 30
    const OBS = collect(20:5:N)

    ## Reference: the walk as standard-normal innovations `z` scaled by
    ## `sigma_rw` and summed from `log_R0`, with the ramp added daily.
    @model function innovation_walk(
            n::Integer, log_R0_base::Real;
            breakpoint, rt_start::Integer = 1, week::Integer = 7,
            ramp::Real = RT_INTERVENTION_RAMP
        )
        days = knot_days(n; week, start = rt_start)
        nb = length(days)
        log_R0 := log_R0_base
        sigma_rw ~ truncated(Normal(0, 0.1); lower = 0)
        z ~ product_distribution(fill(Normal(0, 1), nb - 1))
        intervention_effect ~ truncated(Normal(0, 0.4); upper = 0)
        log_R = log_R0 .+ vcat(zero(log_R0), cumsum(sigma_rw .* z))
        log_Rt = interpolate_knots(log_R, days, n) .+
            intervention_effect .* sigmoid_ramp(n, breakpoint; ramp)
        return (; Rt = exp.(log_Rt), log_R)
    end

    ## The infection process scored against counts, so the log joint carries
    ## a likelihood that reads the whole Rt path through the renewal.
    @model function scored(rt, y)
        lat ~ to_submodel(
            infection_model(
                N; breakpoint = BP, rt_start = 1, rt_walk_start = WS, rt
            ), false
        )
        y ~ product_distribution(Poisson.(lat.infections[OBS] .+ 1.0e-8))
        return lat
    end

    const KNOTS = @varname(rt_state.log_Rt_knots)
    const Z = @varname(rt_state.z)

    ## Map a draw `θ` of the production walk onto the innovation
    ## coordinates, filling a draw `θold` of the reference model.
    function to_innovations(θ, θold, log_R0)
        days = knot_days(N; start = WS)
        σ = θ[@varname(rt_state.sigma_rw)]
        e = θ[@varname(rt_state.intervention_effect)]
        net = θ[KNOTS] .- e .* sigmoid_ramp(N, BP)[days[2:end]]
        z = diff(vcat(log_R0, net)) ./ σ
        for k in keys(θ)
            k == KNOTS && continue
            θold = DynamicPPL.setindex!!(θold, θ[k], k)
        end
        return DynamicPPL.setindex!!(θold, z, Z), σ, length(z)
    end
end

@testitem "RandomWalkVector: the density is a chain of Normal steps" begin
    using BVDOutbreakSize: RandomWalkVector
    using Distributions: Normal, logpdf
    using Random: Xoshiro
    using Statistics: mean, var

    offset = [0.0, -0.1, -0.3, -0.4]
    d = RandomWalkVector(0.3, 0.1, offset)
    @test length(d) == 4
    x = [0.35, 0.1, 0.0, 0.0]
    level = x .- offset
    expected = logpdf(Normal(0.3, 0.1), level[1]) +
        sum(logpdf(Normal(level[i - 1], 0.1), level[i]) for i in 2:4)
    @test logpdf(d, x) ≈ expected rtol = 1.0e-12

    ## Draws are the walk plus the offset, with variance growing by one step
    ## per knot.
    rng = Xoshiro(7)
    draws = [rand(rng, d) for _ in 1:20_000]
    @test mean(x[4] for x in draws) ≈ 0.3 + offset[4] atol = 0.01
    @test var(x[1] for x in draws) ≈ 0.1^2 rtol = 0.1
    @test var(x[4] for x in draws) ≈ 4 * 0.1^2 rtol = 0.1

    ## A zero SD still gives a finite density, peaked on the offset path.
    d0 = RandomWalkVector(0.0, 0.0, [0.1, 0.2])
    @test isfinite(logpdf(d0, [0.1, 0.2]))
    @test logpdf(d0, [0.1, 0.2]) > logpdf(d0, [0.1 + 1.0e-6, 0.2])
end

@testitem "rt_walk_model: the knots are the total log Rt on knot days" setup = [
    RtWalkFixtures,
] begin
    m = rt_walk_model(N, log(1.4); breakpoint = BP, rt_start = WS)
    days = knot_days(N; start = WS)
    for seed in 1:5
        θ = rand(Xoshiro(seed), m)
        @test haskey(θ, @varname(log_Rt_knots))
        @test !haskey(θ, @varname(z))
        s = returned(m, θ)
        ## Ramp included: the sampled knots are the daily total at each knot
        ## day after the first, and the returned walk is net of the ramp.
        @test log.(s.Rt[days[2:end]]) ≈ θ[@varname(log_Rt_knots)]
        e = θ[@varname(intervention_effect)]
        @test s.log_R ≈ log.(s.Rt[days]) .- e .* sigmoid_ramp(N, BP)[days]
        @test s.log_R[1] == log(1.4)
    end
end

@testitem "rt_walk_model: the same model as the innovation form" setup = [
    RtWalkFixtures,
] begin
    y = [3, 5, 8, 12, 20, 30, 45, 60, 80]
    new = scored(rt_walk_model, y)
    old = scored(innovation_walk, y)
    for seed in 1:5
        θ = rand(Xoshiro(seed), new)
        lat = returned(new, θ)
        θold, σ, m = to_innovations(θ, rand(Xoshiro(0), old), log(lat.R0))
        lat_old = returned(old, θold)
        ## The same Rt path and infections at mapped points.
        @test lat_old.Rt ≈ lat.Rt rtol = 1.0e-12
        @test lat_old.infections ≈ lat.infections rtol = 1.0e-10
        ## Knot levels are `log_R0` plus `σ` times the summed innovations, so
        ## the change of coordinates has Jacobian `σ^m`. The joint density
        ## agrees once it is accounted for.
        jac = m * log(σ)
        @test logjoint(new, θ) ≈ logjoint(old, θold) - jac rtol = 1.0e-10
    end
end

@testitem "rt_walk_model: the Mooncake gradient matches ForwardDiff" tags = [
    :ad,
] setup = [RtWalkFixtures] begin
    using ADTypes: AutoForwardDiff
    using LogDensityProblems: logdensity_and_gradient
    using BVDOutbreakSize: default_adtype
    import ForwardDiff

    y = [3, 5, 8, 12, 20, 30, 45, 60, 80]
    for model in (
            rt_walk_model(N, log(1.4); breakpoint = BP, rt_start = WS),
            scored(rt_walk_model, y),
        )
        vi = DynamicPPL.link(DynamicPPL.VarInfo(Xoshiro(3), model), model)
        x = collect(vi[:])
        grads = map((default_adtype(), AutoForwardDiff())) do adtype
            ldf = DynamicPPL.LogDensityFunction(
                model, DynamicPPL.getlogjoint, vi; adtype
            )
            logdensity_and_gradient(ldf, x)
        end
        (lp, g), (lp_fd, g_fd) = grads
        @test isfinite(lp)
        @test lp ≈ lp_fd
        @test all(isfinite, g)
        @test g ≈ g_fd rtol = 1.0e-8
    end
end
