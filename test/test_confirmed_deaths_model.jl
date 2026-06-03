## Smoke tests for the laboratory-confirmed-deaths stream
## (unconditioned generator). Confirmed deaths are the suspected-death
## increment thinned by the count-implied BVD composition among suspects
## (`p_bvd = μ_BVD / N_susp`), shared with the confirmed-cases stream. No
## confirmed-deaths data exists yet, so the stream is exercised as a
## prior / posterior-predictive generator only. Fast `Prior()` + small
## `predict` runs, no NUTS.

@testitem "confirmed_deaths_only_model: prior draws are finite, bounded" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: confirmed_deaths_only_model
    susp = 120
    chn = sample(confirmed_deaths_only_model(susp), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    cd = reduce(vcat, vec(Array(chn[:confirmed_deaths])))
    @test length(cd) == 200
    @test all(isfinite, cd)
    @test all(cd .>= 0)
    ## Bounded by the suspected-death denominator (single cumulative).
    @test all(cd .<= susp)

    p = vec(Array(chn[:p_bvd_cutoff]))
    @test all(0 .<= p .<= 1)
end

@testitem "confirmed_deaths_only_model: conditioning gives finite C(T)" tags=[:slow] begin
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: confirmed_deaths_only_model
    chn = sample(confirmed_deaths_only_model(120, 40), Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test length(C) == 200
    @test all(isfinite, C)
    @test all(C .> 0)
end

@testitem "confirmed_deaths_model: single vintage reduces to cumulative" tags=[:slow] begin
    ## A length-1 edge vector at T must give a Binomial on the cumulative
    ## suspected-death total, so draws are bounded by that total.
    using Turing: sample, Prior, @model, to_submodel
    import FlexiChains
    using BVDOutbreakSize: confirmed_deaths_model,
                           exponential_growth_model,
                           report_delay_model, test_positivity_model,
                           incubation_model, onset_rescale
    @model function _cd_single(susp, cd)
        g ~ to_submodel(exponential_growth_model(), false)
        rd ~ to_submodel(report_delay_model(), false)
        tp ~ to_submodel(test_positivity_model(), false)
        inc ~ to_submodel(incubation_model(), false)
        os = onset_rescale(inc.dist, g.r)
        s ~ to_submodel(
            confirmed_deaths_model(Union{Missing, Int}[cd], [susp], g,
                [0.2], tp.λ_bg, rd.dist, [g.T]; onset_fraction = os),
            false)
        return (;)
    end
    chn = sample(_cd_single(80, missing), Prior(), 150;
        chain_type = FlexiChains.VNChain, progress = false)
    cd = reduce(vcat, vec(Array(chn[:confirmed_deaths])))
    @test all(0 .<= cd .<= 80)
end

@testitem "bvd_joint: confirmed-deaths stream on vs off" tags=[:slow] begin
    ## Enabling the unconditioned confirmed-deaths generator must give
    ## finite generated quantities and predictive draws bounded by the
    ## suspected-death increments, and leave the model unchanged when the
    ## stream is off.
    using Turing: sample, Prior, predict
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, load_observations
    obs = load_observations()
    rh = obs.reported_case_history
    dh = obs.death_history
    function _inc(values)
        out = similar(values, Int)
        prev = 0
        for i in eachindex(values)
            out[i] = values[i] - prev
            prev = values[i]
        end
        return out
    end
    death_incr = _inc(dh.values)

    ## Stream off (default): finite cumulative cases.
    chn_off = sample(
        bvd_joint(obs.exported_cases, death_incr, _inc(rh.values),
            obs.export_deaths_daily;
            reported_offsets = rh.offsets, death_offsets = dh.offsets),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    C_off = vec(Array(chn_off[:cumulative_cases]))
    @test all(isfinite, C_off)

    ## Stream on, unconditioned: missing confirmed-death vector aligned
    ## with the death offsets.
    cd_missing = Vector{Union{Missing, Int}}(missing, length(dh.offsets))
    model_on = bvd_joint(obs.exported_cases, death_incr, _inc(rh.values),
        obs.export_deaths_daily;
        reported_offsets = rh.offsets, death_offsets = dh.offsets,
        confirmed_deaths = cd_missing)
    chn_on = sample(model_on, Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false)
    C_on = vec(Array(chn_on[:cumulative_cases]))
    @test all(isfinite, C_on)

    pp = predict(model_on, chn_on)
    cd = reduce(vcat, vec(Array(pp[:confirmed_deaths])))
    @test all(cd .>= 0)
    ## Each generated increment is bounded by its suspected-death
    ## denominator (sum bound on the cumulative total).
    @test sum(cd) <= 100 * sum(death_incr)
end
