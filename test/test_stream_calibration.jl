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
