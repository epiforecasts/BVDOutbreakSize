## Smoke tests for predict_no_onward_deaths. A prior-predictive draw from
## the joint model provides the bare `:CFR`, `:C_T`, and
## `:expected_deaths_T` deterministics the function reads (the joint
## composer re-exposes these; the single-stream composers expose only
## `C_T`).

@testitem "predict_no_onward_deaths returns the documented columns" tags=[
    :slow
] begin
    using DataFrames: DataFrame, nrow
    using Turing: sample, Prior
    import FlexiChains
    using BVDOutbreakSize: bvd_joint, predict_no_onward_deaths

    chn = sample(
        bvd_joint(40, missing, missing),
        Prior(), 100;
        chain_type = FlexiChains.VNChain, progress = false
    )
    obs_deaths = 18
    df = predict_no_onward_deaths(chn; obs_deaths = obs_deaths)

    @test df isa DataFrame
    @test sort(names(df)) == sort(["delta_deaths", "total_projected"])
    @test nrow(df) == 100

    @test all(isfinite, df.delta_deaths)
    @test all(isfinite, df.total_projected)
    @test all(df.delta_deaths .>= 0)
    @test all(df.total_projected .>= obs_deaths)
    @test maximum(
        abs.(df.total_projected .- (obs_deaths .+ df.delta_deaths))
    ) < 1e-8
end

## Both plotted quantities have a hard floor: `delta_deaths` is clamped at
## zero, so `total_projected` cannot fall below the deaths already observed.
## A free kernel density spreads mass past the smallest draw, which put a
## visible negative tail on the still-expected-deaths panel.

@testitem "plot_no_onward_deaths keeps each density inside its bound" setup=[
    HeadlessMakie
] begin
    using DataFrames: DataFrame
    using CairoMakie: Axis
    using BVDOutbreakSize: plot_no_onward_deaths

    ## Draws piled against the zero clamp, which is what produced the tail.
    delta = max.(collect(range(-40, 160; length = 500)), 0.0)
    obs_deaths = 2325
    df = DataFrame(delta_deaths = delta,
        total_projected = obs_deaths .+ delta)

    fig = plot_no_onward_deaths(df; obs_deaths = obs_deaths)
    axes = [c for c in fig.content if c isa Axis]
    @test length(axes) == 2

    ## Each axis starts at its bound, so no impossible value is on show: no
    ## negative still-expected deaths, and no projected total below the
    ## deaths already observed.
    @test first(axes[1].limits[])[1] == 0
    @test first(axes[2].limits[])[1] == obs_deaths

    ## The upper limit still clears the draws, so clipping the impossible
    ## side does not also clip the distribution.
    @test last(first(axes[1].limits[])) >= maximum(delta)
    @test last(first(axes[2].limits[])) >= maximum(df.total_projected)
end
