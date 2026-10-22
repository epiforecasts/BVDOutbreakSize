## The onset-to-admission Gamma scale θ (onset-to-report, onset-to-detection
## and the onset-to-admission leg of onset-to-death) carries a LogNormal
## prior, so the sampler meets no hard lower bound near zero (#1087).

@testitem "onset-to-admission scale prior is LogNormal at the reanalysis median" begin
    using BVDOutbreakSize: reported_cases_model, deaths_model, exports_model
    using Turing.DynamicPPL: extract_priors, @varname
    using Distributions

    onsets = fill(5.0, 40)
    empty = (; days = Int[], counts = Int[])
    rep = extract_priors(reported_cases_model(empty, missing, onsets, 5.0, 0.3))
    dth = extract_priors(deaths_model(empty, missing, onsets, 5.0))
    xpt = extract_priors(
        exports_model(missing, onsets, 0.1; incubation_pmf = [0.5, 0.5])
    )
    priors = [
        rep[@varname(report_state.θ)],
        dth[@varname(od_state.oa.θ)],
        xpt[@varname(detect_state.θ)],
    ]
    for d in priors
        @test d isa LogNormal
        @test median(d) ≈ 3.694
        @test d.σ ≈ 0.32
        ## Open at zero: no truncation bound for the sampler to sit on.
        @test minimum(d) == 0
        ## On the log scale the sampler moves in, a scale of 0.13 (an
        ## onset-to-report mean near 0.15 d) is more than 40 nats below the
        ## median. The truncated Normal it replaces put it 4.4 nats below.
        lp(θ) = logpdf(d, θ) + log(θ)
        @test lp(0.13) - lp(3.694) < -40
    end
end

@testitem "AD gradient: reported_cases_model is finite near a small scale (Mooncake)" tags = [
    :ad,
] begin
    using BVDOutbreakSize: reported_cases_model, default_adtype
    using Turing.DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: Xoshiro

    onsets = collect(range(1.0, 30.0; length = 40))
    hist = (; days = [13, 18, 40], counts = [40, 70, 150])
    model = reported_cases_model(hist, 150, onsets, 5.0, 0.3)
    vi0 = DynamicPPL.link(DynamicPPL.VarInfo(Xoshiro(1), model), model)
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi0; adtype = default_adtype()
    )
    for θ in (0.13, 3.694)
        init = DynamicPPL.InitFromParams(
            Dict(@varname(report_state.θ) => θ), DynamicPPL.InitFromPrior()
        )
        vi = DynamicPPL.VarInfo(Xoshiro(1), model, init)
        x = collect(DynamicPPL.link(vi, model)[:])
        @test log(θ) in x
        logp, grad = logdensity_and_gradient(ldf, x)
        @test isfinite(logp)
        @test length(grad) == length(x)
        @test all(isfinite, grad)
        @test any(!iszero, grad)
    end
end
