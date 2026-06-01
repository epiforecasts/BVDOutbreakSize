## The DRC sitrep streams are proper per-vintage observations scored with
## `~`, so `predict` replicates them and the joint fits under NUTS without
## the `check_model = false` escape. These items exercise both.

@testitem "bvd_joint fits under NUTS without check_model=false" tags = [
    :slow] begin
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, nuts_sample

    n = 40
    dh = (; days = [13, 18, 40], counts = [10, 14, 18])
    rh = (; days = [13, 18, 40], counts = [340, 516, 905])
    ch = (; days = [13, 18, 40], counts = [9, 17, 27])
    ## Default check_model = true: a passing fit proves no stream leaves a
    ## sampled discrete latent.
    chn = nuts_sample(
        bvd_joint(n, 2, 18, 905, 0, 27, 50;
            deaths_history = dh,
            reported_history = rh,
            confirmed_history = ch,
            lab_history = (; days = [18, 40], counts = [30, 50]),
            breakpoint = 30);
        samples = 12, chains = 1, progress = false)
    C_T = vec(Array(chn[:C_T]))
    @test length(C_T) == 12
    @test all(isfinite, C_T)
    @test all(C_T .> 0)
end

@testitem "predict replicates the per-vintage DRC streams" tags = [
    :slow] begin
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, nuts_sample
    using Turing: predict

    n = 40
    dh = (; days = [13, 18, 40], counts = [10, 14, 18])
    rh = (; days = [13, 18, 40], counts = [340, 516, 905])
    ch = (; days = [13, 18, 40], counts = [9, 17, 27])
    fitted = bvd_joint(n, 2, 18, 905, 0, 27;
        deaths_history = dh, reported_history = rh, confirmed_history = ch,
        breakpoint = 30)
    chn = nuts_sample(fitted; samples = 12, chains = 1, progress = false)

    ## Keep each stream's vintage day grid but drop the counts, so the
    ## increments are resampled by `predict` rather than held fixed.
    _days_only(h) = (; days = h.days, counts = Int[])
    gen = bvd_joint(n, missing, missing, missing, missing, missing;
        deaths_history = _days_only(dh),
        reported_history = _days_only(rh),
        confirmed_history = _days_only(ch),
        breakpoint = 30)
    pp = predict(gen, chn)

    ## The indexed per-vintage increment variables must be replicated for
    ## each DRC stream, under the stream's prefixed submodel name.
    for (prefix, m) in (
        ("cases_state.reported_increments", length(rh.days)),
        ("confirmed_state.confirmed_increments", length(ch.days)),
        ("deaths_state.death_increments", length(dh.days)))
        for i in 1:m
            key = Symbol("$prefix.increments[$i]")
            v = vec(Array(pp[key]))
            @test all(isfinite, v)
            @test all(v .>= 0)
        end
    end
end
