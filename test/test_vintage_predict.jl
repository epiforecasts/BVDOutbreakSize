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
            confirmed_deaths = 5,
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
    using BVDOutbreakSize: bvd_joint, nuts_sample, confirmed_positivity_windows
    using Turing: predict, @varname

    n = 40
    dh = (; days = [13, 18, 40], counts = [10, 14, 18])
    rh = (; days = [13, 18, 40], counts = [340, 516, 905])
    ch = (; days = [13, 18, 40], counts = [9, 17, 27])
    lh = (; days = [18, 40], counts = [30, 50])
    ## Daily new-suspect inflow on days after the cumulative series: per-day
    ## counts, scored single-day, not cumulated.
    sdh = (; days = [35, 36, 37], counts = [12, 9, 7])
    fitted = bvd_joint(n, 2, 18, 905, 0, 27;
        confirmed_deaths = 5,
        deaths_history = dh, reported_history = rh, confirmed_history = ch,
        lab_history = lh, suspected_daily_history = sdh, breakpoint = 30)
    chn = nuts_sample(fitted; samples = 12, chains = 1, progress = false)

    ## Keep the reported and death vintage day grids but drop their counts,
    ## so their increments are resampled by `predict`. The confirmed early
    ## counts and observed positives resample when `confirmed_cases` is
    ## `missing`; their denominators (both histories' counts) stay fixed.
    _days_only(h) = (; days = h.days, counts = Int[])
    gen = bvd_joint(n, missing, missing, missing, missing, missing;
        deaths_history = _days_only(dh),
        reported_history = _days_only(rh),
        confirmed_history = ch,
        lab_history = lh,
        suspected_daily_history = _days_only(sdh),
        breakpoint = 30)
    pp = predict(gen, chn)

    ## The daily inflow replicates as one per-day count vector under
    ## `cases_state.suspected_daily.increments`, one entry per report day, all
    ## finite and non-negative — the vector the docs PPC daily panel plots.
    sd = vec(pp[FlexiChains.Parameter(
        @varname(cases_state.suspected_daily.increments))])
    @test !isempty(sd)
    @test all(d -> length(d) == length(sdh.days), sd)
    @test all(d -> all(isfinite, d) && all(>=(0), d), sd)

    ## `predict` replicates each stream's per-vintage variable under its
    ## prefixed submodel name. The confirmed stream splits into early
    ## counts and observed-window positives. Each draw is the full vector;
    ## every replicated count must be non-negative and finite.
    w = confirmed_positivity_windows(ch, lh)
    keys = (
        FlexiChains.Parameter(@varname(cases_state.reported_increments.increments)),
        FlexiChains.Parameter(@varname(confirmed_state.early_increments.increments)),
        FlexiChains.Parameter(@varname(confirmed_state.confirmed_positives.positives)),
        FlexiChains.Parameter(@varname(deaths_state.death_increments.increments)))
    lens = (length(rh.days), length(w.early_days), length(w.obs_analysed),
        length(dh.days))
    for (key, m) in zip(keys, lens)
        draws = vec(pp[key])
        @test !isempty(draws)
        @test all(d -> length(d) == m, draws)
        @test all(d -> all(isfinite, d) && all(>=(0), d), draws)
    end
end

@testitem "predict replicates late windows incl 24h-anchored days" tags = [
    :slow] begin
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, nuts_sample, confirmed_positivity_windows
    using Turing: predict, @varname

    ## A scenario with confirmed vintages after the last cumulative
    ## laboratory date (28) — so there are late windows — one of which (35)
    ## publishes a 24h analysed count and is anchored. The late group is a
    ## single submodel over all late days, so its predict-key vector has one
    ## entry per late day and the docs PPC reconstruction (early ++ obs ++
    ## late) lines up with the window-day grid (regression guard for #235).
    n = 40
    dh = (; days = [13, 18, 40], counts = [10, 14, 18])
    rh = (; days = [13, 18, 40], counts = [340, 516, 905])
    ch = (; days = [13, 18, 28, 35, 40], counts = [9, 17, 27, 33, 40])
    lh = (; days = [18, 28], counts = [30, 50])
    ldh = (; days = [35], counts = [12])
    fitted = bvd_joint(n, 2, 18, 905, 0, 40;
        confirmed_deaths = 5,
        deaths_history = dh, reported_history = rh, confirmed_history = ch,
        lab_history = lh, lab_daily_history = ldh, breakpoint = 30)
    chn = nuts_sample(fitted; samples = 12, chains = 1, progress = false)

    _days_only(h) = (; days = h.days, counts = Int[])
    gen = bvd_joint(n, missing, missing, missing, missing, missing;
        deaths_history = _days_only(dh), reported_history = _days_only(rh),
        confirmed_history = ch, lab_history = lh, lab_daily_history = ldh,
        breakpoint = 30)
    pp = predict(gen, chn)

    w = confirmed_positivity_windows(ch, lh, ldh)
    @test w.late_days == [35, 40]
    @test w.late_analysed == [12, 0]

    late = vec(pp[FlexiChains.Parameter(
        @varname(confirmed_state.late_increments.increments))])
    early = vec(pp[FlexiChains.Parameter(
        @varname(confirmed_state.early_increments.increments))])
    obs = vec(pp[FlexiChains.Parameter(
        @varname(confirmed_state.confirmed_positives.positives))])
    ## One predict entry per late day, and the three groups partition the
    ## confirmed window grid (early ++ obs ++ late) exactly.
    @test all(d -> length(d) == length(w.late_days), late)
    n_windows = length(w.early_days) + length(w.obs_analysed) +
                length(w.late_days)
    @test length(first(early)) + length(first(obs)) + length(first(late)) ==
          n_windows
    @test all(d -> all(isfinite, d) && all(>=(0), d), late)
end
