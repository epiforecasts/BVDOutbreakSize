## Tests for the sum-to-zero helpers in src/sum_to_zero.jl: the basis, the
## loading matrix and its implied moments, each against a naive reference.

@testsnippet SumToZeroReference begin
    using Random: Xoshiro, randn
    using Statistics: mean
    using Distributions: LKJCholesky, Normal, Chi, truncated
    using BVDOutbreakSize: sum_to_zero_basis, sum_to_zero_factor,
        sum_to_zero, sum_to_zero_moments, bartlett_factor

    ## Naive dense references, written out entry by entry.
    matmul(A, B) = [
        sum(A[i, m] * B[m, j] for m in axes(A, 2); init = 0.0)
            for i in axes(A, 1), j in axes(B, 2)
    ]
    diagm_(v) = [i == j ? v[i] : 0.0 for i in eachindex(v), j in eachindex(v)]
    centring(n) = [(i == j) - 1 / n for i in 1:n, j in 1:n]
    covariance(F) = matmul(F, permutedims(F))

    ## A centred per-patch construction: `n` per-patch scales, an `n × n`
    ## correlation factor, then the mean subtracted.
    centred_factor(σ, L) = matmul(centring(length(σ)), matmul(diagm_(σ), L))

    ## A Bartlett factor of a `Wishart(ν, I_k)` draw, written out directly.
    function bartlett_draw(rng, k, ν)
        A = zeros(k, k)
        for i in 1:k
            A[i, i] = sqrt(sum(abs2, randn(rng, ν - i + 1)))
            for j in 1:(i - 1)
                A[i, j] = randn(rng)
            end
        end
        return A
    end
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
    ## `N(0, σ²)` draws with their mean subtracted, from one draw fewer.
    for n in 2:6
        σ = 0.37
        F = sum_to_zero_factor(sum_to_zero_basis(n), σ)
        @test covariance(F) ≈ σ^2 .* centring(n) atol = 1.0e-14
        old = centred_factor(fill(σ, n), diagm_(ones(n)))
        @test covariance(F) ≈ covariance(old) atol = 1.0e-14
        mom = sum_to_zero_moments(F)
        @test all(≈(σ * sqrt((n - 1) / n)), mom.sd)
        @test all(
            isapprox(mom.cor[i, j], i == j ? 1.0 : -1 / (n - 1); atol = 1.0e-12)
                for i in 1:n, j in 1:n
        )
    end
end

@testitem "bartlett_factor: fills a lower-triangular factor row by row" setup = [
    SumToZeroReference,
] begin
    A = bartlett_factor([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    @test A == [1.0 0.0 0.0; 4.0 2.0 0.0; 5.0 6.0 3.0]
    @test bartlett_factor([0.7], Float64[]) == fill(0.7, 1, 1)
    @test_throws DimensionMismatch bartlett_factor([1.0, 2.0], [1.0, 2.0])
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

@testitem "sum_to_zero: prior matches a centred per-patch prior" setup = [
    SumToZeroReference,
] begin
    ## With `n` per-patch scales `σ_p ~ half-N(0, c)` and an
    ## `n × n` LKJ correlation, then centred, each patch's deviation has
    ## expected variance `c² (n - 1) / n`. The basis form,
    ## `σ √((n - 1) / tr(A Aᵀ)) Q A` with `σ ~ half-N(0, c)` and
    ## `A Aᵀ ~ Wishart(ν, I_{n-1})`, has the same expected variance for every
    ## patch, and the same distribution of the sd for every patch, whatever
    ## order the patches come in.
    rng = Xoshiro(3)
    c = 0.05
    sd_prior = truncated(Normal(0, c); lower = 0)
    N = 40_000
    for n in (3, 4)
        Q = sum_to_zero_basis(n)
        ν = n - 1
        old_var = zeros(n)
        new_var = zeros(n)
        new_sd = zeros(N, n)
        for d in 1:N
            Lo = rand(rng, LKJCholesky(n, 2.0)).L
            Fo = centred_factor(rand(rng, sd_prior, n), Lo)
            A = bartlett_draw(rng, n - 1, ν)
            Fn = sum_to_zero_factor(
                Q, rand(rng, sd_prior) * sqrt((n - 1) / sum(abs2, A)), A
            )
            old_var .+= sum_to_zero_moments(Fo).sd .^ 2 ./ N
            sd = sum_to_zero_moments(Fn).sd
            new_var .+= sd .^ 2 ./ N
            new_sd[d, :] = sd
        end
        @test all(isapprox.(old_var, c^2 * (n - 1) / n; rtol = 0.03))
        @test all(isapprox.(new_var, c^2 * (n - 1) / n; rtol = 0.03))
        ## Every patch has the same mass near zero, and at least a fifth of
        ## the prior mass.
        near_zero = [mean(new_sd[:, p] .< 0.4 * c) for p in 1:n]
        @test maximum(near_zero) - minimum(near_zero) < 0.015
        @test minimum(near_zero) > 0.2
    end
end

@testitem "patch_rt_model: the prior treats every patch alike" tags = [
    :slow,
] setup = [SumToZeroReference] begin
    using BVDOutbreakSize: patch_rt_model
    using Turing: sample, Prior
    using Turing.DynamicPPL: returned
    using Statistics: median, quantile

    ## Reordering the patches rotates the sum-to-zero basis, and the Wishart
    ## prior on the basis covariance does not change under a rotation. So
    ## every patch's sd and every pair's correlation has the same prior. A
    ## prior built from independent scales on the basis directions fails
    ## this: the last patch loads on one direction only.
    np = 4
    m = patch_rt_model(60, np, log(1.5); rt_start = 10)
    rets = vec(
        returned(m, sample(Xoshiro(5), m, Prior(), 6000; progress = false))
    )
    sd_near_zero = [mean(r.σ_δ[p] < 0.02 for r in rets) for p in 1:np]
    sd_q90 = [quantile([r.σ_δ[p] for r in rets], 0.9) for p in 1:np]
    @test maximum(sd_near_zero) - minimum(sd_near_zero) < 0.04
    @test maximum(sd_q90) / minimum(sd_q90) < 1.06
    pairs = [(p, q) for p in 1:np for q in (p + 1):np]
    cor_med = [median([r.Ω[p, q] for r in rets]) for (p, q) in pairs]
    cor_q90 = [quantile([r.Ω[p, q] for r in rets], 0.9) for (p, q) in pairs]
    @test maximum(cor_med) - minimum(cor_med) < 0.06
    @test maximum(cor_q90) - minimum(cor_q90) < 0.08
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
        rng, sum_to_zero_factor, Q, 0.4, L; is_primitive = false
    )
    Mooncake.TestUtils.test_rule(
        rng, sum_to_zero, F, randn(rng, 3); is_primitive = false
    )
    Mooncake.TestUtils.test_rule(
        rng, F -> sum_to_zero_moments(F).cor, F; is_primitive = false
    )
    Mooncake.TestUtils.test_rule(
        rng, bartlett_factor, [1.2, 0.8, 0.5], randn(rng, 3);
        is_primitive = false
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
                "patch_rt_model, three patches",
                patch_rt_model(60, 3, log(1.5); rt_start = 10),
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
                "province_composition_model, ascertainment and severity",
                province_composition_model(
                    obs, modelled;
                    severity_sd_prior = truncated(Normal(0, 0.3); lower = 0)
                ),
            ),
        )
        check(model)
    end
end

@testitem "sum-to-zero sites: deviations follow the basis passed in" setup = [
    SumToZeroReference,
] begin
    using BVDOutbreakSize: patch_rt_model, province_composition_model
    using Turing.DynamicPPL: OnlyAccsVarInfo, RawValueAccumulator,
        InitFromPrior, UnlinkAll, init!!, get_raw_values, @varname

    ## A prior draw's return value and its parameter values.
    function prior_draw(model, seed)
        accs = OnlyAccsVarInfo(RawValueAccumulator(false))
        r, vi = init!!(Xoshiro(seed), model, accs, InitFromPrior(), UnlinkAll())
        return r, get_raw_values(vi)
    end

    ## Deviation knots from the formulas: level `s Q B z_level`, innovation
    ## `η_k = c Q B z_k` with `z_k` the `k - 1`-th block of `nd` draws, and
    ## the AR(1) in closed form `δ_k = φ^(k-1) δ_1 + Σ_{j ≤ k} φ^(k-j) η_j`.
    function reference_knots(Q, B, s, c, z_level, z_drift, φ, nb)
        n, nd = size(Q)
        QB = matmul(Q, B)
        apply(v) = vec(matmul(QB, reshape(v, :, 1)))
        level = s .* apply(z_level)
        η(k) = c .* apply(z_drift[((k - 2) * nd + 1):((k - 1) * nd)])
        knots = zeros(n, nb)
        for k in 1:nb
            knots[:, k] = φ^(k - 1) .* level
            for j in 2:k
                knots[:, k] .+= φ^(k - j) .* η(j)
            end
        end
        return knots
    end
    ## A basis of the same subspace in a different orientation: the default
    ## one with its first two columns rotated.
    function rotated(Q)
        size(Q, 2) < 2 && return -Q
        R = diagm_(ones(size(Q, 2)))
        R[1, 1], R[1, 2], R[2, 1], R[2, 2] = 0.6, -0.8, 0.8, 0.6
        return matmul(Q, R)
    end

    n = 60
    for np in (2, 3, 4), seed in 1:3
        nd = np - 1
        Q = sum_to_zero_basis(np)
        for (basis, model) in (
                (
                    Q,
                    patch_rt_model(n, np, log(1.5); rt_start = 10),
                ),
                (
                    rotated(Q),
                    patch_rt_model(
                        n, np, log(1.5); rt_start = 10,
                        basis = rotated(Q)
                    ),
                ),
            )
            r, vi = prior_draw(model, seed)
            φ = exp2(-7 / vi[@varname(δ_halflife)])
            ## The first knot is drawn at the stationary scale.
            σ_level = vi[@varname(σ_drift)] / sqrt(1 - φ^2)
            ## Two patches draw no lower entry.
            d = vi[@varname(bartlett_diag)]
            o = np > 2 ? vi[@varname(bartlett_lower)] : Float64[]
            B = zeros(nd, nd)
            m = 0
            for i in 1:nd
                B[i, i] = d[i]
                for j in 1:(i - 1)
                    m += 1
                    B[i, j] = o[m]
                end
            end
            s = σ_level * sqrt(nd / sum(abs2, B))
            c = vi[@varname(σ_drift)] * sqrt(nd / sum(abs2, B))
            nb = size(r.δ_knots, 2)
            ref = reference_knots(
                basis, B, s, c, vi[@varname(z_level)], vi[@varname(z_drift)],
                φ, nb
            )
            @test r.δ_knots ≈ ref rtol = 1.0e-12 atol = 1.0e-14
            @test r.drift_factor ≈ c .* matmul(basis, B) rtol = 1.0e-12
        end
    end

    ## The composition multipliers on a passed basis: the log ascertainment
    ## and the log severity are each `τ Q z`.
    obs = [853 21 42; 77 2 5; 3 0 0; 10 1 2]
    modelled = [800.0 20.0 40.0; 70.0 2.5 4.0; 2.0 0.1 0.2; 9.0 1.0 1.5]
    Qr = rotated(sum_to_zero_basis(4))
    model = province_composition_model(
        obs, modelled;
        severity_sd_prior = truncated(Normal(0, 0.3); lower = 0), basis = Qr
    )
    for seed in 1:3
        r, vi = prior_draw(model, seed)
        apply(v) = vec(matmul(Qr, reshape(v, :, 1)))
        log_asc = vi[@varname(τ_asc)] .* apply(vi[@varname(z_asc)])
        @test log.(r.province_ascertainment) ≈ log_asc rtol = 1.0e-12
        @test log.(r.province_severity) ≈
            vi[@varname(τ_sev)] .* apply(vi[@varname(z_sev)]) rtol = 1.0e-12
    end
end
