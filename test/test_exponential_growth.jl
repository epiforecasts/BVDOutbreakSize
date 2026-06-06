## Tests for the molecular-clock growth-and-size prior and the two-phase
## anchor-day renewal seeding it drives in `infection_model`. The growth
## submodel is ported from the integral model: a WIDE prior on the
## doubling count `m`, so the outbreak age `T = m·τ` is prior-dominated and
## carries the genetic bound, while the renewal sets the realized size.

@testitem "exponential_growth_model: deterministics T, C_T, τ" begin
    using Turing: @model, to_submodel, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: exponential_growth_model

    ## `to_submodel(x, false)` already re-exposes the submodel's sampled
    ## (`r`, `m`) and `:=` (`τ`, `T`, `C_T`) names at the parent, so they
    ## surface as bare chain keys without re-declaration here.
    @model function _wrap()
        st ~ to_submodel(exponential_growth_model(), false)
        return st
    end

    chn = sample(_wrap(), Prior(), 400;
        chain_type = FlexiChains.VNChain, progress = false)
    T = vec(Array(chn[:T]))
    C_T = vec(Array(chn[:C_T]))
    τ = vec(Array(chn[:τ]))
    m = vec(Array(chn[:m]))
    r = vec(Array(chn[:r]))

    @test all(isfinite, T) && all(T .> 0)
    @test all(isfinite, C_T) && all(C_T .> 0)
    ## τ = log(2)/r and T = m·τ, C_T = 2^m hold draw-by-draw.
    @test all(isapprox.(τ, log(2) ./ r; rtol = 1e-8))
    @test all(isapprox.(T, m .* τ; rtol = 1e-8))
    @test all(isapprox.(C_T, 2.0 .^ m; rtol = 1e-8))
end

@testitem "exponential_growth_model: m prior is WIDE" begin
    using Statistics: mean, std
    using Turing: @model, to_submodel, sample, Prior
    import FlexiChains
    using BVDOutbreakSize: exponential_growth_model

    @model function _wrap()
        st ~ to_submodel(exponential_growth_model(), false)
        return st
    end

    chn = sample(_wrap(), Prior(), 4000;
        chain_type = FlexiChains.VNChain, progress = false)
    m = vec(Array(chn[:m]))
    T = vec(Array(chn[:T]))

    ## Centre near 9, SD near 3 (truncated Normal(9, 3)): the prior is
    ## deliberately wide so T stays uncertain.
    @test 8.0 < mean(m) < 10.0
    @test 2.5 < std(m) < 3.3
    ## The induced T is correspondingly wide (SD well above the SD-3
    ## tightness of a fixed-anchor placement).
    @test std(T) > 30.0
end

@testitem "infection_model: two-phase anchor seeding" tags=[:slow] begin
    using Statistics: mean, std
    using Turing: @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: infection_model, nuts_sample

    n = 101
    anchor = 31          # genetic-TMRCA grid day
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
    r0 = vec(Array(chn[:r0]))

    @test all(isfinite, T) && all(isfinite, C_T) && all(C_T .> 0)
    ## T is WIDE (prior-dominated by the m prior), not pinned near τ_obs.
    @test std(T) > 20.0
    @test mean(T) > τ_obs    # origin sits before the anchor (cryptic phase)
    ## The anchor seed grows forward over τ_obs to the prior size scale.
    @test all(sa .> 0)
    @test all(isapprox.(sa .* exp.(r0 .* τ_obs), 2.0 .^ vec(Array(chn[:m]));
        rtol = 1e-6))
end
