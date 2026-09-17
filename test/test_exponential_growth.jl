## Tests for the molecular-clock growth-and-size prior and the two-phase
## renewal-start seeding it drives in `infection_model`. The growth
## submodel samples the cryptic rate `r` (the prior sits on the genetic
## doubling time) along with the generation count `m`, so the cryptic
## duration `T = m·G` is prior-dominated and does not depend on `r`; the
## established `R0` is derived forward from `r` in `infection_model`. The
## renewal-start seed is the daily incidence `C_T = exp(r·T)` reached at the
## end of the cryptic phase, and the total age `T + τ_obs` carries the
## genetic bound while the renewal sets the realized size.

@testitem "exponential_growth_model: deterministics T, C_T, τ" begin
    using Turing: @model, to_submodel, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: exponential_growth_model

    ## `r` is sampled (the growth-rate prior), along with `m`.
    ## `to_submodel(x, false)` re-exposes `r`, `m` and the `:=` (`τ`, `G`,
    ## `T`, `C_T`) names at the parent.
    gi = [0.02, 0.08, 0.16, 0.22, 0.20, 0.15, 0.10, 0.07]
    gi = gi ./ sum(gi)
    G_true = sum(i * gi[i] for i in eachindex(gi))

    @model function _wrap()
        st ~ to_submodel(exponential_growth_model(gi), false)
        return st
    end

    chn = sample(_wrap(), Prior(), 400;
        chain_type = FlexiChains.VNChain, progress = false)
    T = vec(Array(chn[:T]))
    C_T = vec(Array(chn[:C_T]))
    τ = vec(Array(chn[:τ]))
    G = vec(Array(chn[:G]))
    m = vec(Array(chn[:m]))
    r = vec(Array(chn[:r]))

    @test all(isfinite, T) && all(T .> 0)
    @test all(isfinite, C_T) && all(C_T .> 0)
    ## τ = log(2)/r, T = m·G and C_T = exp(r·T) hold draw-by-draw.
    @test all(isapprox.(τ, log(2) ./ r; rtol = 1e-8))
    @test all(isapprox.(G, G_true; rtol = 1e-8))
    @test all(isapprox.(T, m .* G_true; rtol = 1e-8))
    @test all(isapprox.(C_T, exp.(r .* T); rtol = 1e-8))
    ## The point of counting generations: elapsed time per generation is the
    ## generation interval, so it is the same for every draw regardless of the
    ## growth rate. Counting doublings would make it log(2)/r and vary.
    @test all(isapprox.(T ./ m, G_true; rtol = 1e-8))
end

@testitem "exponential_growth_model: r and m priors are wide" begin
    using Statistics: mean, std, quantile
    using Turing: @model, to_submodel, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: exponential_growth_model

    gi = [0.02, 0.08, 0.16, 0.22, 0.20, 0.15, 0.10, 0.07]
    gi = gi ./ sum(gi)
    G_true = sum(i * gi[i] for i in eachindex(gi))

    @model function _wrap()
        st ~ to_submodel(exponential_growth_model(gi), false)
        return st
    end

    chn = sample(_wrap(), Prior(), 4000;
        chain_type = FlexiChains.VNChain, progress = false)
    m = vec(Array(chn[:m]))
    T = vec(Array(chn[:T]))
    r = vec(Array(chn[:r]))

    ## Centre near 2.79, SD near 1.16 (truncated Normal(2.75, 1.2); lower 0).
    @test 2.7 < mean(m) < 2.9
    @test 1.06 < std(m) < 1.26
    ## The growth rate is centred on the BEAST X 11.7-day doubling
    ## (r ≈ 0.059).
    @test 0.05 < mean(r) < 0.08
    ## `r` is LogNormal, so `log r` pins both the centre and the spread
    ## directly. The centre is log(log 2 / 11.7) ≈ -2.826 and the log-SD is
    ## 0.40. At 4000 draws the standard error on each is under 0.007, so
    ## these bounds hold the prior to the documented values rather than
    ## merely to the right order of magnitude.
    @test -2.86 < mean(log.(r)) < -2.79
    @test 0.37 < std(log.(r)) < 0.43
    ## The induced doubling time τ = log 2 / r is LogNormal(log 11.7, 0.40),
    ## a 95% interval of 5.3-25.6 d. This is the interval the analysis text
    ## quotes, so it is guarded here.
    τ = log(2) ./ r
    @test 4.9 < quantile(τ, 0.025) < 5.8
    @test 23.6 < quantile(τ, 0.975) < 27.8
    ## The cryptic duration T = m·G is `m`'s own spread scaled by the fixed
    ## generation interval, and carries none of `r`'s.
    @test isapprox(std(T), std(m) * G_true; rtol = 1e-8)
    @test std(T) > 4.0
end

@testitem "infection_model: two-phase renewal-start seeding" tags=[:slow] begin
    using Statistics: mean, std
    using Turing: @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: infection_model, nuts_sample

    n = 101
    renewal_start = 38   # genetic-TMRCA day + renewal-start lead
    τ_obs = n - renewal_start

    @model function _wrap()
        st ~ to_submodel(infection_model(n; rt_start = renewal_start), false)
        T := st.T
        C_T := st.C_T
        m := st.m
        sa := st.seed_at_renewal_start
        r0 := st.r0
        G := sum(i * st.g[i] for i in eachindex(st.g))
        return st
    end

    chn = nuts_sample(_wrap(); samples = 200, chains = 2, progress = false)
    T = vec(Array(chn[:T]))
    C_T = vec(Array(chn[:C_T]))
    sa = vec(Array(chn[:sa]))
    m = vec(Array(chn[:m]))
    r0 = vec(Array(chn[:r0]))
    G = vec(Array(chn[:G]))

    @test all(isfinite, T) && all(isfinite, C_T) && all(C_T .> 0)
    ## Total age T = m·G + τ_obs is prior-dominated by the m prior and
    ## ≥ τ_obs by construction (the cryptic phase adds m·G ≥ 0).
    @test std(T) > 12.0
    @test all(T .>= τ_obs)
    @test mean(T) > τ_obs    # origin sits before the renewal start (cryptic)
    ## Elapsed cryptic time is `m` generation intervals, set by the sampled
    ## generation interval alone. Counting doublings instead would make it
    ## log(2)/r and tie the origin date to the growth rate.
    @test all(isapprox.(T .- τ_obs, m .* G; rtol = 1e-6))
    ## The renewal-start seed is the daily incidence the cryptic phase reaches,
    ## `exp(r·m·G)`, so the origin date and the seed share one growth rate.
    @test all(sa .> 0)
    @test all(isapprox.(sa, exp.(r0 .* (T .- τ_obs)); rtol = 1e-6))
end

@testitem "infection_model: rt_walk_start holds Rt flat before the walk" begin
    using Turing: returned
    using Random: default_rng, seed!
    using BVDOutbreakSize: infection_model

    ## The renewal/seeding starts at `rt_start` (the genetic-TMRCA renewal
    ## start), but the reproduction-number random walk can be held flat at the
    ## established `R0` until a later day `rt_walk_start` (the first situation
    ## report), so no unsupported Rt drift happens over the pre-surveillance
    ## window where the dynamics are unidentified. With `breakpoint = missing`
    ## there is no intervention ramp, so `Rt` must be exactly constant on every
    ## day from day 1 up to the first walk knot at `rt_walk_start`, while the
    ## renewal still seeds and grows from the earlier `rt_start`.
    n = 109
    rt_start = 45
    rt_walk_start = 81
    seed!(13)
    m = infection_model(n; rt_start, rt_walk_start)
    rng = default_rng()
    for _ in 1:50
        s = returned(m, rand(rng, m))
        ## Flat (= the established R0) from day 1 through the walk start.
        @test all(s.Rt[1:rt_walk_start] .≈ s.Rt[1])
        ## The renewal still ran over the full grid: a finite, positive
        ## trajectory of the right length (seeding unaffected by the walk
        ## start).
        @test length(s.infections) == n
        @test all(isfinite, s.infections)
        @test all(>(0), s.infections)
    end

    ## Default `rt_walk_start = rt_start`: the walk starts at the renewal
    ## start, so Rt is free to vary from `rt_start` onward whenever the walk
    ## start is not decoupled from it.
    seed!(13)
    m0 = infection_model(n; rt_start)
    s0 = returned(m0, rand(rng, m0))
    @test all(s0.Rt[1:rt_start] .≈ s0.Rt[1])
end

@testitem "infection_model: current r and R_T agree in sign per draw" begin
    using Turing: returned
    using Random: default_rng, seed!
    using BVDOutbreakSize: infection_model, doubling_time, euler_lotka_r

    ## The reported current growth rate `r` and the cut-off reproduction
    ## number `R_T = Rt[n]` describe the same instant, so they must never
    ## disagree in sign: `r < 0` must imply `R_T < 1` and `r >= 0` must imply
    ## `R_T >= 1`. The realised last-two-days slope `log I[n] - log I[n-1]`
    ## can disagree (the intervention ramp can depress the final renewal step
    ## so `I[n] < I[n-1]` while `Rt[n] >= 1`), so `r` is instead derived from
    ## `R_T` and the generation interval via Euler–Lotka, which keeps the
    ## relationship exact by construction.
    n = 109
    rt_start = 45       # genetic-TMRCA day + renewal-start lead
    breakpoint = 84     # intervention ramp near the cut-off

    seed!(11)
    m = infection_model(n; breakpoint, rt_start)
    rng = default_rng()
    ndraws = 1500
    disagree = 0
    for _ in 1:ndraws
        s = returned(m, rand(rng, m))
        r = s.r
        R_T = s.Rt[n]
        ## Sign agreement at the cut-off.
        ((r < 0) && (R_T >= 1)) && (disagree += 1)
        ((r >= 0) && (R_T < 1)) && (disagree += 1)
        ## `r` is the Euler–Lotka growth implied by `R_T` and `g`, and the
        ## reported doubling time is `log 2 / r`, both consistent with `R_T`.
        @test isapprox(r, euler_lotka_r(R_T, s.g); atol = 1e-8)
        @test isapprox(s.doubling_time, doubling_time(r); rtol = 1e-8)
    end
    @test disagree == 0
end
