## Smoke tests for the per-vintage (per-sitrep) confirmed-cases
## likelihood. The base `confirmed_cases_model` is now per-vintage:
## a cumulative-count vector + edge times; a length-1 vector reduces
## to the cumulative single-total likelihood. Tiny fits exercise the
## real `bvd_joint` composer with the full history.

@testitem "confirmed_cases: prior draws are finite, non-negative" tags=[:slow] begin
    using BVDOutbreakSize: bvd_joint, load_observations
    using Turing: sample, Prior
    import FlexiChains

    obs = load_observations()
    rep = obs.reported_case_history
    conf = obs.confirmed_case_history
    dh = obs.death_history
    n_rep = length(rep.values)
    n_conf = length(conf.values)
    n_dh = length(dh.values)
    ## Prior-predictive: pass missing vectors with matching offsets.
    m = bvd_joint(missing,
        fill(missing, n_dh), fill(missing, n_rep);
        reported_offsets = rep.offsets,
        death_offsets = dh.offsets,
        confirmed_cases = fill(missing, n_conf),
        confirmed_offsets = conf.offsets)
    chn = sample(m, Prior(), 200;
        chain_type = FlexiChains.VNChain, progress = false)
    raw = vec(Array(chn[:confirmed_cases]))
    flat = reduce(vcat, raw)
    @test all(isfinite, flat)
    @test all(flat .>= 0)
    ## Per-vintage increments reconstruct a monotone cumulative series.
    @test all(issorted(cumsum(v)) for v in raw)
end

@testitem "confirmed_cases: tiny fit stays positive" tags=[:slow] begin
    using BVDOutbreakSize: bvd_joint, load_observations
    using Turing: sample, Prior
    import FlexiChains

    obs = load_observations()
    rep = obs.reported_case_history
    conf = obs.confirmed_case_history
    dh = obs.death_history
    ## Models observe between-vintage increments, not cumulative totals.
    function _increments(values)
        out = similar(values, Int)
        prev = 0
        for i in eachindex(values)
            out[i] = values[i] - prev
            prev = values[i]
        end
        return out
    end
    m = bvd_joint(obs.exported_cases,
        _increments(dh.values), _increments(rep.values),
        obs.export_deaths_daily;
        reported_offsets = rep.offsets,
        death_offsets = dh.offsets,
        confirmed_cases = _increments(conf.values),
        confirmed_offsets = conf.offsets,
        tests_analysed = obs.cumulative_tests_analysed,
        tests_offset = 0,
        first_export_detection_delta = obs.first_export_detection_delta)
    chn = sample(m, Prior(), 50;
        chain_type = FlexiChains.VNChain, progress = false)
    C = vec(Array(chn[:cumulative_cases]))
    @test length(C) == 50
    @test all(isfinite, C)
    @test all(C .> 0)
end
