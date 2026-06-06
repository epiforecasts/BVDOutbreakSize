## Tests for the molecular-clock growth-and-size prior and the two-phase
## anchor-day renewal seeding it drives in `infection_model`. The growth
## submodel takes the cryptic rate `r` from the single `R0` (Euler–Lotka)
## and samples only the doubling count `m`, so the cryptic duration `m·τ`
## is prior-dominated; the anchor seed magnitude is `2^m` DIRECTLY (no
## back-scaling), and the total age `T = m·τ + τ_obs` carries the genetic
## bound while the renewal sets the realized size.

@testitem "exponential_growth_model: deterministics T, C_T, τ" begin
    using Turing: @model, to_submodel, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: exponential_growth_model

    ## `r` is injected (the R0-implied cryptic rate), not sampled; the
    ## submodel samples only `m`. `to_submodel(x, false)` re-exposes the
    ## sampled `m` and the `:=` (`τ`, `T`, `C_T`) names at the parent.
    r_inj = log(2) / 20      # ≈ the molecular-clock rate (20-day doubling)
    @model function _wrap()
        st ~ to_submodel(exponential_growth_model(r_inj), false)
        return st
    end

    chn = sample(_wrap(), Prior(), 400;
        chain_type = FlexiChains.VNChain, progress = false)
    T = vec(Array(chn[:T]))
    C_T = vec(Array(chn[:C_T]))
    τ = vec(Array(chn[:τ]))
    m = vec(Array(chn[:m]))

    @test all(isfinite, T) && all(T .> 0)
    @test all(isfinite, C_T) && all(C_T .> 0)
    ## τ = log(2)/r is fixed by the injected rate; T = m·τ (cryptic
    ## duration), C_T = 2^m hold draw-by-draw.
    @test all(isapprox.(τ, log(2) / r_inj; rtol = 1e-8))
    @test all(isapprox.(T, m .* τ; rtol = 1e-8))
    @test all(isapprox.(C_T, 2.0 .^ m; rtol = 1e-8))
end

@testitem "exponential_growth_model: m prior is WIDE" begin
    using Statistics: mean, std
    using Turing: @model, to_submodel, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: exponential_growth_model

    r_inj = log(2) / 20
    @model function _wrap()
        st ~ to_submodel(exponential_growth_model(r_inj), false)
        return st
    end

    chn = sample(_wrap(), Prior(), 4000;
        chain_type = FlexiChains.VNChain, progress = false)
    m = vec(Array(chn[:m]))
    T = vec(Array(chn[:T]))

    ## Centre near 2, SD near 3 (truncated Normal(2, 3); lower 0): `m` now
    ## counts only the CRYPTIC doublings, so the centre is lower than the
    ## integral model's. The lower truncation lifts the mean above 2. The
    ## prior is deliberately wide so the cryptic duration stays uncertain.
    @test 2.0 < mean(m) < 4.0
    @test 1.9 < std(m) < 2.7
    ## The induced cryptic duration T = m·τ is correspondingly wide.
    @test std(T) > 30.0
end

@testitem "infection_model: two-phase anchor seeding" tags=[:slow] begin
    using Statistics: mean, std
    using Turing: @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: infection_model, nuts_sample

    n = 101
    anchor = 38          # genetic-TMRCA day + seeding-anchor lead
    τ_obs = n - anchor

    @model function _wrap()
        st ~ to_submodel(infection_model(n; rt_start = anchor), false)
        T := st.T
        C_T := st.C_T
        m := st.m
        sa := st.seed_at_anchor
        r0 := st.r0
        return st
    end

    chn = nuts_sample(_wrap(); samples = 200, chains = 2, progress = false)
    T = vec(Array(chn[:T]))
    C_T = vec(Array(chn[:C_T]))
    sa = vec(Array(chn[:sa]))
    m = vec(Array(chn[:m]))

    @test all(isfinite, T) && all(isfinite, C_T) && all(C_T .> 0)
    ## Total age T = m·τ + τ_obs is WIDE (prior-dominated by the m prior)
    ## and ≥ τ_obs by construction (the cryptic phase adds m·τ ≥ 0).
    @test std(T) > 20.0
    @test all(T .>= τ_obs)
    @test mean(T) > τ_obs    # origin sits before the anchor (cryptic phase)
    ## The anchor seed magnitude is `2^m` DIRECTLY, r-independent (no
    ## back-scaling): keeps the single R0 out of the seed magnitude.
    @test all(sa .> 0)
    @test all(isapprox.(sa, 2.0 .^ m; rtol = 1e-6))
end
