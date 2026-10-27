## AD-gradient smoke check: the package default backend (Mooncake) must be
## able to differentiate the models — the one property the NUTS sampler
## actually relies on. A single unconstrained-space log-density gradient per
## component is the minimal, fast way to assert that, avoiding a full NUTS
## fit as a differentiability check: that takes ~25 min and can flake during
## sampler adaptation rather than in the gradient itself.
##
## The components come from `test/ad_fixtures.jl`, which
## `benchmark/benchmarks.jl` includes too. One source means a component
## cannot be benchmarked without also being asserted differentiable, or
## asserted without being timed. Adding a model there adds it to both.
##
## The full `bvd_joint` is not a component here: one gradient sits behind a
## cold Mooncake compile of roughly 18 minutes. Its gradient in the suite is
## the production joint with and without the rules in
## `test/test_mooncake_rules.jl`.
## The components below are the submodels, the latent process and every
## single-stream composer, so the joint is a composition of surfaces each
## asserted differentiable here. Its own gradient is timed
## by the benchmark suite and exercised end to end by the NUTS fits the docs
## build runs.
##
## Tagged `:ad`, the tag the downgrade-compat run skips (AD gradients drift
## below the package's pinned dependency versions), matching the other
## AD-sensitive items. The gradient pattern mirrors the Enzyme-extension
## check in `test/enzyme/runtests.jl`, which validates the same models.

@testitem "AD gradient: every component differentiates (Mooncake)" tags = [
    :ad,
] begin
    using LogDensityProblems: logdensity_and_gradient
    using BVDOutbreakSize: default_adtype
    include(joinpath(@__DIR__, "ad_fixtures.jl"))

    scenarios = ADFixtures.scenarios()
    @test length(scenarios) >= ADFixtures.MIN_SCENARIOS

    for scen in scenarios
        @testset "$(scen.group): $(scen.name)" begin
            ldf, x = ADFixtures.log_density_function(scen, default_adtype())
            logp, grad = logdensity_and_gradient(ldf, x)
            @test isfinite(logp)
            @test length(grad) == length(x)
            @test all(isfinite, grad)
            ## A non-trivial gradient (not identically zero) confirms AD
            ## actually ran.
            @test any(!iszero, grad)
        end
    end
end

@testitem "AD gradient: province_composition_model differentiates (Mooncake)" tags = [:ad] begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: province_composition_model, default_adtype
    using Distributions: truncated, Normal

    ## The stick-breaking BetaBinomial composition is the likelihood the
    ## spatial data enter through, so its gradient has to exist.
    ##
    ## The reshaping that turns the per-province histories into this matrix
    ## (`province_increment_matrix`) looks provinces up by name in a
    ## `Dict{String}`. It must stay OUTSIDE the model body: a string compare
    ## on the tape is a `memcmp` foreigncall that Mooncake has no rule for,
    ## and it aborts the gradient of the whole joint.
    seed!(20260518)
    obs = [853 21 42; 77 2 5; 3 0 0]
    modelled = [800.0 20.0 40.0; 70.0 2.5 4.0; 2.0 0.1 0.2]
    ## The death composition adds the logit severity contrast about the
    ## national case-fatality ratio.
    models = (
        province_composition_model(obs, modelled),
        province_composition_model(
            obs, modelled;
            severity_sd_prior = truncated(Normal(0, 0.15); lower = 0),
            severity_baseline = 0.4
        ),
    )
    for model in models
        vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
        x0 = collect(vi[:])
        ldf = DynamicPPL.LogDensityFunction(
            model, DynamicPPL.getlogjoint, vi; adtype = default_adtype()
        )
        logp, grad = logdensity_and_gradient(ldf, x0)
        @test isfinite(logp)
        @test length(grad) == length(x0)
        @test all(isfinite, grad)
        @test any(!iszero, grad)
    end
end

@testitem "AD gradient: province suspected-case detection differentiates (Mooncake)" tags = [:ad] begin
    using Turing: DynamicPPL, @model, to_submodel
    using LogDensityProblems: logdensity_and_gradient
    using Random: seed!
    using BVDOutbreakSize: province_suspected_ascertainment_model,
        reported_cases_model, province_composition_model,
        _patch_confirmed_increments, _patch_carried, default_adtype
    using Distributions: Beta

    ## The joint's detection path: the patch probabilities about the
    ## national centre weight the patch onsets, whose sum feeds the national
    ## suspected cases and whose rows feed the case composition.
    seed!(20261022)
    n = 40
    onsets_matrix = vcat(
        reshape(collect(1.0:n) .^ 1.5, 1, :),
        reshape(collect(1.0:n), 1, :),
        reshape(fill(2.0, n), 1, :),
    )
    history = (; days = [10, 20, 30, 40], counts = [30, 120, 400, 900])
    prov = [60 200 300; 10 30 60; 2 5 9]
    days = [20, 30, 40]
    @model function detection_path(onsets_matrix, history, prov, days)
        p_drc ~ Beta(15, 5)
        detection_state ~ to_submodel(
            province_suspected_ascertainment_model(size(onsets_matrix, 1), p_drc)
        )
        detected = detection_state.relative .* onsets_matrix
        cases_state ~ to_submodel(
            reported_cases_model(
                history, history.counts[end],
                vec(sum(detected; dims = 1)), 20.0, p_drc
            )
        )
        carried = _patch_carried(detected, cases_state.report_pmf)
        composition_state ~ to_submodel(
            province_composition_model(
                prov, _patch_confirmed_increments(carried, 0.9, days)
            )
        )
    end
    model = detection_path(onsets_matrix, history, prov, days)
    vi = DynamicPPL.link(DynamicPPL.VarInfo(model), model)
    x0 = collect(vi[:])
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype()
    )
    logp, grad = logdensity_and_gradient(ldf, x0)
    @test isfinite(logp)
    @test length(grad) == length(x0)
    @test all(isfinite, grad)
    ## The detection coordinates reach the density.
    names = string.(keys(vi))
    @test any(contains("τ_susp"), names)
    @test any(!iszero, grad)
end

## The health-zone scenario the benchmark suite times with
## by default. It is not in the default scenario list, so it is asserted
## here, at the size it is timed at and with mixing and the zone
## correlation on.
@testitem "AD gradient: the zone benchmark scenario differentiates (Mooncake)" tags = [
    :ad,
] begin
    using LogDensityProblems: logdensity_and_gradient
    using BVDOutbreakSize: default_adtype
    include(joinpath(@__DIR__, "ad_fixtures.jl"))

    scen = only(
        filter(s -> s.group == "Zone", ADFixtures.scenarios(; zone = true))
    )
    @test !any(s -> s.group == "Zone", ADFixtures.scenarios())
    zd = scen.model.args.zd
    @test size(zd.counts) == (60, 95)
    @test length(zd.patch_ranges) == 4
    @test zd.n == 150
    @test zd.mixing !== nothing
    @test !isempty(zd.zone_distances)
    ldf, x = ADFixtures.log_density_function(scen, default_adtype())
    logp, grad = logdensity_and_gradient(ldf, x)
    @test isfinite(logp)
    @test length(grad) == length(x)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
end

## Gradient check far from the prior median: at an onset-to-report scale of
## 0.13 the Mooncake gradient along log θ matches a central finite
## difference of the log density.
@testitem "AD gradient: reported_cases_model matches finite differences at a small scale (Mooncake)" tags = [
    :ad,
] begin
    using LogDensityProblems: logdensity, logdensity_and_gradient
    using BVDOutbreakSize: reported_cases_model, default_adtype
    using Turing: DynamicPPL
    using Random: Xoshiro
    include(joinpath(@__DIR__, "ad_fixtures.jl"))

    onsets = ADFixtures.latent_draw().onsets
    model = reported_cases_model(
        ADFixtures.REPORTED_HISTORY, ADFixtures.REPORTED_TOTAL, onsets, 5.0,
        0.3; suspected_daily_history = ADFixtures.SUSPECTED_DAILY_HISTORY
    )
    θ = 0.13
    init = DynamicPPL.InitFromParams(
        Dict(DynamicPPL.@varname(report_state.θ) => θ),
        DynamicPPL.InitFromPrior()
    )
    vi = DynamicPPL.link(
        DynamicPPL.VarInfo(Xoshiro(ADFixtures.SEED), model, init), model
    )
    x = collect(vi[:])
    idx = findall(≈(log(θ)), x)
    @test length(idx) == 1
    i = only(idx)
    ldf = DynamicPPL.LogDensityFunction(
        model, DynamicPPL.getlogjoint, vi; adtype = default_adtype()
    )
    logp, grad = logdensity_and_gradient(ldf, x)
    @test isfinite(logp)
    h = 1.0e-5
    e = zeros(length(x))
    e[i] = h
    fd = (logdensity(ldf, x .+ e) - logdensity(ldf, x .- e)) / (2h)
    @test grad[i] ≈ fd rtol = 1.0e-4
end

@testitem "AD gradient: cases_only_model with the background walk matches finite differences (Mooncake)" tags = [
    :ad,
] begin
    using LogDensityProblems: logdensity, logdensity_and_gradient
    using BVDOutbreakSize: cases_only_model, reported_cases_model,
        background_walk_model, default_adtype
    using Random: Xoshiro
    include(joinpath(@__DIR__, "ad_fixtures.jl"))

    ## The single-stream composers run the constant background. Injecting
    ## the centred walk, as the production joint does, puts its sum-to-zero
    ## projection on the gradient path.
    walk = (nn; kw...) -> background_walk_model(nn, 0.2; onset = 5, kw...)
    cases = (args...; kw...) -> reported_cases_model(
        args...; background_re = walk, kw...
    )
    model = cases_only_model(
        ADFixtures.N, ADFixtures.REPORTED_TOTAL;
        reported_history = ADFixtures.REPORTED_HISTORY,
        suspected_daily_history = ADFixtures.SUSPECTED_DAILY_HISTORY,
        cases
    )
    scen = ADFixtures.Scenario(
        "cases_only_model (background walk)", "Composer", model,
        ADFixtures.SEED
    )
    ldf, x = ADFixtures.log_density_function(scen, default_adtype())
    logp, grad = logdensity_and_gradient(ldf, x)
    @test isfinite(logp)
    @test all(isfinite, grad)
    @test any(!iszero, grad)
    ## A directional finite difference checks the whole gradient at once.
    v = randn(Xoshiro(1), length(x))
    v ./= sqrt(sum(abs2, v))
    h = 1.0e-5
    fd = (logdensity(ldf, x .+ h .* v) - logdensity(ldf, x .- h .* v)) / (2h)
    @test sum(grad .* v) ≈ fd rtol = 1.0e-4
end
