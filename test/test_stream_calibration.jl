## Tests for the per-stream posterior-predictive calibration metrics:
## `bias_sample` and the `stream_calibration` summary table.

@testitem "bias_sample matches the scoringutils sign convention" begin
    using BVDOutbreakSize: bias_sample

    pred = collect(1.0:100.0)
    ## Observation at the predictive median → zero bias.
    @test abs(bias_sample(50.5, pred)) < 1.0e-9
    ## Observation above the whole sample → the stream is under-predicted.
    @test bias_sample(200.0, pred) == -1.0
    ## Observation below the whole sample → over-predicted.
    @test bias_sample(-5.0, pred) == 1.0
    ## Ties are split half-and-half: a sample all equal to the observation
    ## gives zero bias.
    @test bias_sample(7.0, fill(7.0, 10)) == 0.0
    ## Empty predictive sample is undefined.
    @test isnan(bias_sample(1.0, Float64[]))
end

@testitem "stream_calibration scores a well-calibrated cumulative panel" begin
    using BVDOutbreakSize: stream_calibration
    using DataFrames: DataFrame, nrow
    using Random: MersenneTwister

    rng = MersenneTwister(1)
    ## Observed cumulative counts with a constant increment of 10. Each draw
    ## carries a per-vintage increment vector centred (symmetrically) on 10,
    ## so the conditional cumulative sits with its median on the observed
    ## count at every vintage: bias ≈ 0 and coverage ≈ nominal.
    observed = [10, 20, 30, 40]
    ndraw = 4_000
    replicates = [10 .+ rand(rng, -5:5, length(observed)) for _ in 1:ndraw]
    panel = (;
        title = "Calibrated", observed = observed,
        replicates = replicates,
    )

    df = stream_calibration([panel])
    @test df isa DataFrame
    @test names(df) ==
        ["Stream", "Vintages", "Bias", "50% coverage", "90% coverage"]
    @test nrow(df) == 1
    @test df[1, "Stream"] == "Calibrated"
    @test df[1, "Vintages"] == length(observed)
    @test abs(df[1, "Bias"]) < 0.1
    @test df[1, "90% coverage"] == 1.0
end

@testitem "stream_calibration detects over- and under-prediction" begin
    using BVDOutbreakSize: stream_calibration
    using DataFrames: DataFrame

    observed = [10, 20, 30]
    ## Increments far above the observed increment of 10 → the conditional
    ## cumulative sits above the observed count everywhere → positive bias
    ## (over-prediction) and the observed never falls inside the interval.
    over = (;
        title = "Over", observed = observed,
        replicates = [fill(40.0, length(observed)) for _ in 1:200],
    )
    ## Increments far below → negative bias (under-prediction).
    under = (;
        title = "Under", observed = observed,
        replicates = [fill(1.0, length(observed)) for _ in 1:200],
    )

    df = stream_calibration([over, under])
    @test df[1, "Bias"] == 1.0
    @test df[1, "90% coverage"] == 0.0
    @test df[2, "Bias"] == -1.0
    @test df[2, "90% coverage"] == 0.0
end

@testitem "stream_calibration handles a non-cumulative daily panel" begin
    using BVDOutbreakSize: stream_calibration
    using Random: MersenneTwister

    rng = MersenneTwister(2)
    ## A `cumulative = false` panel has no previous-vintage baseline: each
    ## replicate is its own daily count, scored directly against the observed
    ## daily count. Centre the draws on the observed value for calibration.
    observed = [5, 8, 3, 6]
    replicates = [observed .+ rand(rng, -2:2, length(observed)) for _ in 1:3_000]
    panel = (;
        title = "Daily", observed = observed,
        replicates = replicates, cumulative = false,
    )

    df = stream_calibration([panel])
    @test df[1, "Vintages"] == length(observed)
    @test abs(df[1, "Bias"]) < 0.1
    @test df[1, "90% coverage"] == 1.0
end

@testitem "province_composition_panels: one daily panel per province" begin
    using BVDOutbreakSize: province_composition_panels, stream_calibration

    nd = 400
    shares = [0.6 0.5 0.5; 0.3 0.35 0.35; 0.1 0.15 0.15]
    obs = [60 0 50; 30 0 35; 10 0 15]
    chn = (;
        province_shares = [shares for _ in 1:nd],
        province_composition_rho = fill(0.05, nd),
    )
    panels = province_composition_panels(
        chn; share_key = :province_shares, obs_increments = obs,
        stream = "Confirmed cases", n_patches = 3,
        patch_labels = ["A", "B", "C"]
    )
    @test length(panels) == 3
    @test panels[2].title == "Confirmed cases, B"
    ## The vintage with no observed cases has no split to predict, so it is
    ## dropped rather than scored as a certain zero.
    @test panels[1].observed == [60, 50]
    @test all(length(r) == 2 for r in panels[1].replicates)
    @test length(panels[1].replicates) == nd
    @test !panels[1].cumulative
    ## The replicates are counts that partition each vintage's observed
    ## total, as the fitted composition does.
    for d in 1:nd
        @test sum(p.replicates[d][1] for p in panels) == 100
        @test sum(p.replicates[d][2] for p in panels) == 100
    end
    tbl = stream_calibration(panels)
    @test tbl[!, "Stream"] == [
        "Confirmed cases, A", "Confirmed cases, B",
        "Confirmed cases, C",
    ]
    ## Without the overdispersion there is no predictive to score.
    @test_throws ErrorException province_composition_panels(
        (; province_shares = [shares for _ in 1:nd]);
        share_key = :province_shares, obs_increments = obs,
        stream = "Confirmed cases", n_patches = 3
    )
end

@testitem "_composition_predictive: a fixed seed redraws the same band" begin
    using BVDOutbreakSize: _composition_predictive

    ## Rebuilt reports must redraw the same composition band. The random
    ## stream differs between Julia versions, so the check is that a seed
    ## repeats itself and that another seed does not, not a pinned draw.
    shares = [0.6 0.5; 0.3 0.35; 0.1 0.15]
    ms = [shares for _ in 1:5]
    totals = [200, 300]
    draw(seed) = _composition_predictive(
        ms, fill(0.05, 5), totals, 2; seed
    )
    @test draw(1) == draw(1)
    @test draw(1) != draw(2)
end

@testitem "province_count_panels: splits each draw's national total" begin
    using BVDOutbreakSize: province_count_panels, stream_calibration

    nd = 300
    shares = [0.6 0.5 0.5; 0.3 0.35 0.35; 0.1 0.15 0.15]
    chn = (;
        province_shares = [shares for _ in 1:nd],
        province_composition_rho = fill(0.05, nd),
    )
    ## The national grid carries days the province tables skip (3 and 6),
    ## whose increments merge into the next province vintage. The first
    ## province vintage is the cumulative count to date, so it takes the
    ## baseline as well.
    national_days = [1, 2, 3, 4, 5, 6, 7]
    national = [[d, 1, 2, 3, 0, 4, 5] for d in 1:nd]
    province_days = [2, 4, 7]
    obs = [60 3 5; 30 1 3; 10 1 1]
    panels = province_count_panels(
        chn; share_key = :province_shares, obs_increments = obs,
        province_days, national_days, national_replicates = national,
        baseline = 10, stream = "Confirmed cases", n_patches = 3,
        patch_labels = ["A", "B", "C"]
    )
    @test length(panels) == 3
    @test panels[3].title == "Confirmed cases, C"
    @test panels[1].observed == [60, 3, 5]
    @test length(panels[1].replicates) == nd
    @test !panels[1].cumulative
    ## Each draw's province counts partition that draw's national total over
    ## the province vintage, not the observed one.
    for d in 1:nd
        @test sum(p.replicates[d][1] for p in panels) == 10 + d + 1
        @test sum(p.replicates[d][2] for p in panels) == 5
        @test sum(p.replicates[d][3] for p in panels) == 9
    end
    @test size(stream_calibration(panels), 1) == 3
end

@testitem "province_count_panels: a zero national total gives zeros" begin
    using BVDOutbreakSize: province_count_panels

    nd = 50
    shares = [0.5 0.5; 0.5 0.5]
    chn = (;
        province_death_shares = [shares for _ in 1:nd],
        province_death_composition_rho = fill(0.05, nd),
    )
    panels = province_count_panels(
        chn; share_key = :province_death_shares,
        obs_increments = [1 0; 1 0], province_days = [1, 2],
        national_days = [1, 2],
        national_replicates = [[4, 0] for _ in 1:nd],
        stream = "Confirmed deaths", n_patches = 2,
        patch_labels = ["A", "B"]
    )
    @test all(r[2] == 0 for p in panels for r in p.replicates)
    @test all(sum(p.replicates[d][1] for p in panels) == 4 for d in 1:nd)
end

@testitem "province_count_panels: rejects grids that do not line up" begin
    using BVDOutbreakSize: province_count_panels

    nd = 10
    shares = [0.5; 0.5;;]
    chn = (;
        province_shares = [shares for _ in 1:nd],
        province_composition_rho = fill(0.05, nd),
    )
    kw = (;
        share_key = :province_shares, obs_increments = [1; 1;;],
        stream = "Confirmed cases", n_patches = 2,
        patch_labels = ["A", "B"],
    )
    ## A province vintage off the national grid cannot be binned.
    @test_throws ErrorException province_count_panels(
        chn; kw..., province_days = [5], national_days = [1, 2],
        national_replicates = [[1, 1] for _ in 1:nd]
    )
    ## Every replicate must run over the whole national grid.
    @test_throws ErrorException province_count_panels(
        chn; kw..., province_days = [2], national_days = [1, 2],
        national_replicates = [[1, 1, 1] for _ in 1:nd]
    )
    ## The national replicates must pair one to one with the chain draws.
    @test_throws ErrorException province_count_panels(
        chn; kw..., province_days = [2], national_days = [1, 2],
        national_replicates = [[1, 1] for _ in 1:(nd - 1)]
    )
end
