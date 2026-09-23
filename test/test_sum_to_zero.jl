## Tests for the sum-to-zero helpers in src/sum_to_zero.jl: the basis, the
## loading matrix and its implied moments, each against a naive reference.

@testsnippet SumToZeroReference begin
    using Random: Xoshiro, randn
    using Statistics: mean
    using Distributions: LKJCholesky, Normal, truncated
    using BVDOutbreakSize: sum_to_zero_basis, sum_to_zero_factor,
        sum_to_zero, sum_to_zero_moments

    ## Naive dense references, written out entry by entry.
    matmul(A, B) = [
        sum(A[i, m] * B[m, j] for m in axes(A, 2); init = 0.0)
            for i in axes(A, 1), j in axes(B, 2)
    ]
    diagm_(v) = [i == j ? v[i] : 0.0 for i in eachindex(v), j in eachindex(v)]
    centring(n) = [(i == j) - 1 / n for i in 1:n, j in 1:n]
    covariance(F) = matmul(F, permutedims(F))

    ## The construction the patch model used before: `n` per-patch scales,
    ## an `n × n` correlation factor, then the mean subtracted.
    old_factor(σ, L) = matmul(centring(length(σ)), matmul(diagm_(σ), L))
end

@testitem "sum_to_zero_basis: orthonormal and orthogonal to the ones" setup = [
    SumToZeroReference,
] begin
    for n in 1:7
        Q = sum_to_zero_basis(n)
        @test size(Q) == (n, n - 1)
        QtQ = matmul(permutedims(Q), Q)
        @test all(
            isapprox(QtQ[i, j], i == j; atol = 1.0e-12)
                for i in 1:(n - 1), j in 1:(n - 1)
        )
        @test all(abs(sum(Q[:, j])) < 1.0e-12 for j in 1:(n - 1))
        ## `Q Qᵀ` is the centring projector, so `σ Q z` is a centred
        ## vector of independent draws.
        @test matmul(Q, permutedims(Q)) ≈ centring(n) atol = 1.0e-12
    end
    @test_throws ArgumentError sum_to_zero_basis(0)
end

@testitem "sum_to_zero: every draw sums to zero" setup = [SumToZeroReference] begin
    rng = Xoshiro(1)
    for n in 2:6, _ in 1:200
        Q = sum_to_zero_basis(n)
        s = abs.(randn(rng, n - 1))
        L = n > 2 ? rand(rng, LKJCholesky(n - 1, 2.0)).L : nothing
        z = randn(rng, n - 1)
        δ = sum_to_zero(sum_to_zero_factor(Q, s, L), z)
        @test abs(sum(δ)) < 1.0e-12
        ## The loading matrix is `Q diag(s) L`.
        ref = matmul(Q, isnothing(L) ? diagm_(s) : matmul(diagm_(s), L))
        @test sum_to_zero_factor(Q, s, L) ≈ ref atol = 1.0e-14
        @test δ ≈ vec(matmul(ref, reshape(z, :, 1))) atol = 1.0e-14
    end
    @test_throws DimensionMismatch sum_to_zero(
        sum_to_zero_factor(sum_to_zero_basis(3), 1.0), zeros(3)
    )
end

@testitem "sum_to_zero_factor: a scalar scale is the centred iid vector" setup = [
    SumToZeroReference,
] begin
    ## The exchangeable form has exactly the covariance of `n` independent
    ## `N(0, σ²)` draws with their mean subtracted, which is what the
    ## importation and ascertainment deviations were before, with one draw
    ## fewer.
    for n in 2:6
        σ = 0.37
        F = sum_to_zero_factor(sum_to_zero_basis(n), σ)
        @test covariance(F) ≈ σ^2 .* centring(n) atol = 1.0e-14
        old = old_factor(fill(σ, n), diagm_(ones(n)))
        @test covariance(F) ≈ covariance(old) atol = 1.0e-14
        mom = sum_to_zero_moments(F)
        @test all(≈(σ * sqrt((n - 1) / n)), mom.sd)
        @test all(
            isapprox(mom.cor[i, j], i == j ? 1.0 : -1 / (n - 1); atol = 1.0e-12)
                for i in 1:n, j in 1:n
        )
    end
end

@testitem "sum_to_zero_moments: match the covariance and the draws" setup = [
    SumToZeroReference,
] begin
    rng = Xoshiro(2)
    n = 4
    Q = sum_to_zero_basis(n)
    s = [0.05, 0.03, 0.08]
    L = rand(rng, LKJCholesky(n - 1, 2.0)).L
    F = sum_to_zero_factor(Q, s, L)
    mom = sum_to_zero_moments(F)
    Σ = covariance(F)
    @test mom.sd ≈ sqrt.([Σ[i, i] for i in 1:n])
    @test all(
        mom.cor[i, j] ≈ Σ[i, j] / (mom.sd[i] * mom.sd[j])
            for i in 1:n, j in 1:n
    )
    ## Every row of the covariance sums to zero, so no entry can be
    ## positively correlated with all the others.
    @test all(abs(sum(Σ[i, :])) < 1.0e-14 for i in 1:n)
    @test all(any(mom.cor[i, j] < 0 for j in 1:n if j != i) for i in 1:n)

    ## Empirical moments of the draws agree with the implied ones.
    N = 200_000
    draws = [sum_to_zero(F, randn(rng, n - 1)) for _ in 1:N]
    emp_sd = [sqrt(mean(d[i]^2 for d in draws)) for i in 1:n]
    @test emp_sd ≈ mom.sd rtol = 0.01
    for i in 1:n, j in (i + 1):n
        emp = mean(d[i] * d[j] for d in draws) / (emp_sd[i] * emp_sd[j])
        @test emp ≈ mom.cor[i, j] atol = 0.01
    end

    ## With three entries the standard deviations fix the correlations.
    L3 = [1.0 0.0; 0.3 sqrt(1 - 0.09)]
    F3 = sum_to_zero_factor(sum_to_zero_basis(3), [0.2, 0.5], L3)
    m3 = sum_to_zero_moments(F3)
    sd = m3.sd
    @test m3.cor[1, 2] ≈ (sd[3]^2 - sd[1]^2 - sd[2]^2) / (2 * sd[1] * sd[2])

    ## One entry: the vector is identically zero.
    F1 = sum_to_zero_factor(sum_to_zero_basis(1), 0.5)
    @test size(F1) == (1, 0)
    @test sum_to_zero(F1, Float64[]) == [0.0]
    @test sum_to_zero_moments(F1).cor == ones(1, 1)
end

@testitem "sum_to_zero: prior matches the centred construction it replaces" setup = [
    SumToZeroReference,
] begin
    ## Reference check against the construction the patch model used
    ## before. With `n` per-patch scales `σ_p ~ half-N(0, c)` and an
    ## `n × n` LKJ correlation, then centred, each patch's deviation has
    ## expected variance `c² (n - 1) / n`. The basis form, `n - 1` scales
    ## from the same prior and an `(n - 1) × (n - 1)` LKJ, has the same
    ## expected variance for every patch.
    rng = Xoshiro(3)
    c = 0.05
    sd_prior = truncated(Normal(0, c); lower = 0)
    N = 40_000
    for n in (3, 4)
        Q = sum_to_zero_basis(n)
        old_var = zeros(n)
        new_var = zeros(n)
        for _ in 1:N
            Lo = rand(rng, LKJCholesky(n, 2.0)).L
            Fo = old_factor(rand(rng, sd_prior, n), Lo)
            Ln = rand(rng, LKJCholesky(n - 1, 2.0)).L
            Fn = sum_to_zero_factor(Q, rand(rng, sd_prior, n - 1), Ln)
            old_var .+= sum_to_zero_moments(Fo).sd .^ 2 ./ N
            new_var .+= sum_to_zero_moments(Fn).sd .^ 2 ./ N
        end
        @test all(isapprox.(old_var, c^2 * (n - 1) / n; rtol = 0.03))
        @test all(isapprox.(new_var, c^2 * (n - 1) / n; rtol = 0.03))
    end
end

@testitem "sum_to_zero: Mooncake matches finite differences" tags = [:ad] setup = [
    SumToZeroReference,
] begin
    using Mooncake: Mooncake
    ## No hand-written rule: this checks the rule Mooncake derives for each
    ## helper against finite differences.
    rng = Xoshiro(4)
    Q = sum_to_zero_basis(4)
    s = [0.05, 0.03, 0.08]
    L = rand(rng, LKJCholesky(3, 2.0)).L
    F = sum_to_zero_factor(Q, s, L)
    Mooncake.TestUtils.test_rule(
        rng, sum_to_zero_factor, Q, s, L; is_primitive = false
    )
    Mooncake.TestUtils.test_rule(
        rng, sum_to_zero_factor, Q, 0.4; is_primitive = false
    )
    Mooncake.TestUtils.test_rule(
        rng, sum_to_zero, F, randn(rng, 3); is_primitive = false
    )
    Mooncake.TestUtils.test_rule(
        rng, F -> sum_to_zero_moments(F).cor, F; is_primitive = false
    )
end

@testitem "sum-to-zero submodels: finite log density and Mooncake gradient" tags = [
    :ad,
] begin
    using Turing: DynamicPPL
    using LogDensityProblems: logdensity_and_gradient
    using ADTypes: AutoForwardDiff
    import ForwardDiff
    using Distributions: Normal, truncated
    using Random: Xoshiro
    using BVDOutbreakSize: patch_rt_model, patch_infection_model,
        province_composition_model, default_adtype

    ## Each submodel whose deviations sit on the sum-to-zero basis, at a
    ## prior draw in unconstrained space, against ForwardDiff on the same
    ## log density.
    function check(model)
        vi = DynamicPPL.link(DynamicPPL.VarInfo(Xoshiro(11), model), model)
        x = collect(vi[:])
        grad_with(adtype) = logdensity_and_gradient(
            DynamicPPL.LogDensityFunction(
                model, DynamicPPL.getlogjoint, vi; adtype
            ), x
        )
        logp, grad = grad_with(default_adtype())
        @test isfinite(logp)
        @test all(isfinite, grad)
        @test any(!iszero, grad)
        _, grad_fd = grad_with(AutoForwardDiff())
        @test grad ≈ grad_fd rtol = 1.0e-6
        return nothing
    end

    kernel4 = [p == q ? 0.0 : 1.0e-4 for p in 1:4, q in 1:4]
    obs = [853 21 42; 77 2 5; 3 0 0; 10 1 2]
    modelled = [800.0 20.0 40.0; 70.0 2.5 4.0; 2.0 0.1 0.2; 9.0 1.0 1.5]
    @testset "$name" for (name, model) in (
            (
                "patch_rt_model, two patches",
                patch_rt_model(60, 2, log(1.5); rt_start = 10),
            ),
            (
                "patch_rt_model, four patches",
                patch_rt_model(60, 4, log(1.5); rt_start = 10),
            ),
            (
                "patch_infection_model, coupled",
                patch_infection_model(60, 4; importation_kernel = kernel4),
            ),
            (
                "province_composition_model, covariate and severity",
                province_composition_model(
                    obs, modelled;
                    testing_covariate = [0.6, -0.2, -0.3, -0.1],
                    severity_sd_prior = truncated(Normal(0, 0.3); lower = 0)
                ),
            ),
        )
        check(model)
    end
end
