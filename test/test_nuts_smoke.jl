## Fast NUTS smoke fits for each latent submodel and every composer.
## All use a small number of samples and a single chain with progress
## disabled to stay quick in CI; each asserts the chain is non-empty and
## the key tracked quantities are finite.
##
## Submodels are wrapped in a tiny parent that re-exposes the checked
## quantity with `:=` so it surfaces as a bare chain key (a nested
## submodel return otherwise surfaces under a dotted prefix).

@testitem "NUTS smoke: rt_walk_model (n=40)" tags = [:slow] begin
    using Turing: @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: rt_walk_model, nuts_sample

    ## `log_R0` is now a DERIVED argument (forward Euler–Lotka from the
    ## sampled growth rate), not sampled in the walk; inject a fixed value.
    @model function _wrap_rt()
        st ~ to_submodel(
            rt_walk_model(40, log(1.6); breakpoint = 30), false)
        R_T := st.Rt[end]
        return st
    end

    chn = nuts_sample(_wrap_rt(); samples = 12, chains = 1, progress = false)
    R_T = vec(Array(chn[:R_T]))
    @test length(R_T) == 12
    @test all(isfinite, R_T) && all(R_T .> 0)
end

@testitem "NUTS smoke: generation_interval_model (nmax=40)" tags = [:slow] begin
    using Turing: @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: generation_interval_model, nuts_sample

    @model function _wrap_gi()
        st ~ to_submodel(generation_interval_model(40), false)
        gi_mean := st.gi_mean
        return st
    end

    chn = nuts_sample(_wrap_gi(); samples = 12, chains = 1, progress = false)
    gi_mean = vec(Array(chn[:gi_mean]))
    @test length(gi_mean) == 12
    @test all(isfinite, gi_mean)
    @test all(gi_mean .> 0)
end

@testitem "NUTS smoke: infection_model (n=40)" tags = [:slow] begin
    using Turing: @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: infection_model, nuts_sample

    @model function _wrap_inf()
        st ~ to_submodel(infection_model(40), false)
        r := st.r
        C_T := st.C_T
        return st
    end

    chn = nuts_sample(_wrap_inf(); samples = 12, chains = 1, progress = false)
    r = vec(Array(chn[:r]))
    C_T = vec(Array(chn[:C_T]))
    @test length(r) == 12
    @test all(isfinite, r)
    @test all(C_T .> 0)
end

@testitem "NUTS smoke: exports_only_model (n=40, obs=2)" tags = [:slow] begin
    import FlexiChains
    using BVDOutbreakSize: exports_only_model, nuts_sample

    chn = nuts_sample(exports_only_model(40, 2);
        samples = 12, chains = 1, progress = false)
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 12
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "NUTS smoke: deaths_only_model (n=40, obs=18)" tags = [:slow] begin
    import FlexiChains
    using BVDOutbreakSize: deaths_only_model, nuts_sample

    chn = nuts_sample(
        deaths_only_model(40, 18;
            deaths_history = (; days = [13, 18, 23], counts = [10, 14, 18]));
        samples = 12, chains = 1, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 12
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "NUTS smoke: cases_only_model (n=40, obs=905)" tags = [:slow] begin
    import FlexiChains
    using BVDOutbreakSize: cases_only_model, nuts_sample

    chn = nuts_sample(
        cases_only_model(40, 905;
            reported_history = (; days = [13, 18, 23],
                counts = [340, 516, 905]));
        samples = 12, chains = 1, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 12
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "NUTS smoke: confirmed_only_model (n=40, obs=27)" tags = [:slow] begin
    import FlexiChains
    using BVDOutbreakSize: confirmed_only_model, nuts_sample

    chn = nuts_sample(
        confirmed_only_model(40, 27;
            confirmed_history = (; days = [13, 18, 23],
                counts = [9, 17, 27]));
        samples = 12, chains = 1, progress = false
    )
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 12
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "NUTS smoke: exports_deaths_only_model (n=40, obs=0)" tags = [:slow] begin
    import FlexiChains
    using BVDOutbreakSize: exports_deaths_only_model, nuts_sample

    ## This composer keeps the deaths and exports submodels only for their
    ## CFR / onset-to-death PMF / export onsets, leaving their own counts
    ## missing; that leaves two redundant sampled discrete draws, so the
    ## model check is disabled (see nuts_sample).
    chn = nuts_sample(exports_deaths_only_model(40, 0);
        samples = 12, chains = 1, progress = false, check_model = false)
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 12
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "NUTS smoke: bvd_joint plain (n=40)" tags = [:slow] begin
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, nuts_sample

    ## Streams are positional: n, exported_cases, total_deaths,
    ## reported_cases, exports_deaths, confirmed_cases. Confirmed deaths
    ## are a keyword count; pass one so the joint has no sampled discrete.
    chn = nuts_sample(
        bvd_joint(40, 2, 18, 905, 0, 27; confirmed_deaths = 5);
        samples = 12, chains = 1, progress = false
    )
    for name in [:C_T, :R_T, :CFR, :k, :p_drc, :p_uganda]
        vals = vec(Array(chn[name]))
        @test length(vals) == 12
        @test all(isfinite, vals)
        @test all(vals .> 0)
    end
end

@testitem "NUTS smoke: bvd_joint composition death link (n=40)" tags = [:slow] begin
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, nuts_sample

    ## The main-like confirmed-death lab process (issue #225): the
    ## composition link turns the suspect-death background on so the
    ## death-pool composition is informative, and shares the case-lab
    ## sensitivity/specificity into the death positivity.
    chn = nuts_sample(
        bvd_joint(40, 2, 18, 905, 0, 27; confirmed_deaths = 5,
            confirmed_positivity_link = :composition,
            confirmed_death_link = :composition);
        samples = 12, chains = 1, progress = false
    )
    for name in [:C_T, :lambda_bg_death, :death_confirmation, :m_death]
        vals = vec(Array(chn[name]))
        @test length(vals) == 12
        @test all(isfinite, vals)
    end
    ## The death background must be strictly positive (composition link
    ## forces it on) so the death-pool composition does not collapse to 1.
    @test all(vec(Array(chn[:lambda_bg_death])) .> 0)
    @test all(0 .< vec(Array(chn[:death_confirmation])) .< 1)
end

@testitem "NUTS smoke: bvd_joint with histories + breakpoint + genetic" tags = [:slow] begin
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, genetic_seeding_model, nuts_sample

    n = 40
    dh = (; days = [13, 18, 23], counts = [10, 14, 18])
    rh = (; days = [13, 18, 23], counts = [340, 516, 905])
    ch = (; days = [13, 18, 23], counts = [9, 17, 27])

    chn = nuts_sample(
        bvd_joint(n, 2, 18, 905, 0, 27;
            confirmed_deaths = 5,
            deaths_history = dh,
            reported_history = rh,
            confirmed_history = ch,
            breakpoint = 30,
            genetic = genetic_seeding_model,
            tmrca_days = 30.0,
            tmrca_days_sd = 15.0);
        samples = 12, chains = 1, progress = false
    )
    ## C_T, r, T and CFR must be finite. doubling_time = log(2) / r is
    ## allowed to be non-finite as r crosses zero (documented), so it is
    ## only required to be present and non-NaN (Inf is a valid limit).
    for name in [:C_T, :r, :T, :CFR]
        vals = vec(Array(chn[name]))
        @test length(vals) == 12
        @test all(isfinite, vals)
    end
    dt = vec(Array(chn[:doubling_time]))
    @test length(dt) == 12
    @test !any(isnan, dt)
    C_T = vec(Array(chn[:C_T]))
    @test all(C_T .> 0)
end
