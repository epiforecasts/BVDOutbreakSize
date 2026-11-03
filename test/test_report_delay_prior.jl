## The onset-to-admission Gamma scale θ (onset-to-report, onset-to-detection
## and the onset-to-admission leg of onset-to-death) carries a LogNormal
## prior that is open at zero.

@testitem "onset-to-admission scale prior is LogNormal at the reanalysis median" begin
    using BVDOutbreakSize: reported_cases_model, deaths_model, exports_model,
        ONSET_ADMISSION_THETA_PRIOR
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
        @test d == ONSET_ADMISSION_THETA_PRIOR
        ## On the log scale the sampler moves in, a scale of 0.13 is more
        ## than 50 nats below the median.
        lp(θ) = logpdf(d, θ) + log(θ)
        @test lp(0.13) - lp(3.694) < -50
    end
end
