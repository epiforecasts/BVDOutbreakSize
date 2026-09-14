## Absconding competes with the clinical exits rather than adding to them.
## Deaths and recoveries discharge the whole of `A_bvd` and rule-outs the whole
## of `A_bg`, so an unthinned schedule plus an abscond outflow removes more mass
## than was admitted.

@testitem "abscond_thinned leaves a PMF alone at zero hazard" begin
    using BVDOutbreakSize: abscond_thinned

    pmf = [0.5, 0.3, 0.15, 0.05]
    @test abscond_thinned(pmf, 0.0) == pmf
    @test sum(abscond_thinned(pmf, 0.0)) ≈ 1.0
end

@testitem "abscond_thinned discounts by cohort age, not calendar day" begin
    using BVDOutbreakSize: abscond_thinned

    pmf = [0.25, 0.25, 0.25, 0.25]
    κ = 0.1
    out = abscond_thinned(pmf, κ)
    ## Day `d` of a stay carries `(1 - κ)^d`, so the discount compounds with
    ## time in care and not with position on the grid.
    @test out ≈ [0.25 * (1 - κ)^d for d in 0:3]
    @test out[1] == pmf[1]
    @test issorted(out; rev = true)
    ## The mass removed is the abscond share, so what remains is what leaves
    ## clinically.
    @test sum(out) < 1.0
    @test sum(out) ≈ sum(0.25 * (1 - κ)^d for d in 0:3)
end

@testitem "thinning keeps the occupied stock off zero as absconding rises" begin
    ## The pathology the thinning removes: with an unthinned schedule the
    ## running balance is drained below zero on a declining tail and clipped,
    ## which flattens the likelihood. With it, the stock stays positive at
    ## every abscond rate.
    using BVDOutbreakSize: accumulate_occupancy, convolve_delay,
                           discretise_censored, abscond_thinned
    using Distributions: Gamma

    n = 210
    A_bvd = [60.0 * exp(-0.5 * ((t - 90) / 35)^2) for t in 1:n]
    A_bg = [90.0 * exp(-0.5 * ((t - 90) / 35)^2) for t in 1:n]
    mk(m, s) = discretise_censored(Gamma((m / s)^2, s^2 / m), 60)
    pd, pr, pro = mk(9.0, 5.0), mk(14.0, 6.0), mk(4.0, 2.0)
    conf = fill(0.18, n)

    for κ in (0.0, 0.0063, 0.02, 0.05, 0.1)
        acc = accumulate_occupancy(A_bvd, A_bg,
            convolve_delay(0.44 .* A_bvd, abscond_thinned(pd, κ)),
            convolve_delay(0.56 .* A_bvd, abscond_thinned(pr, κ)),
            convolve_delay(A_bg, abscond_thinned(pro, κ)), κ, conf)
        window = findall(t -> A_bvd[t] + A_bg[t] > 0.5, 1:n)
        @test all(acc.demand[t] > 1e-9 for t in window)
        @test all(acc.O_conf .<= acc.O_bvd .+ 1e-8)
        @test all(acc.O_susp .>= -1e-8)
    end
end

@testitem "an unthinned schedule does drain the stock to zero" begin
    ## The control for the item above: without thinning, the same inputs floor
    ## the stock on a declining tail, and more days floor as the abscond rate
    ## rises. This is the behaviour being fixed.
    using BVDOutbreakSize: accumulate_occupancy, convolve_delay,
                           discretise_censored
    using Distributions: Gamma

    n = 210
    A_bvd = [60.0 * exp(-0.5 * ((t - 90) / 35)^2) for t in 1:n]
    A_bg = [90.0 * exp(-0.5 * ((t - 90) / 35)^2) for t in 1:n]
    mk(m, s) = discretise_censored(Gamma((m / s)^2, s^2 / m), 60)
    deaths = convolve_delay(0.44 .* A_bvd, mk(9.0, 5.0))
    recover = convolve_delay(0.56 .* A_bvd, mk(14.0, 6.0))
    ruleout = convolve_delay(A_bg, mk(4.0, 2.0))
    conf = fill(0.18, n)

    floored(κ) = count(<=(1e-12),
        accumulate_occupancy(A_bvd, A_bg, deaths, recover, ruleout, κ,
            conf).demand)

    @test floored(0.0) < floored(0.02)
    @test floored(0.02) <= floored(0.1)
end
