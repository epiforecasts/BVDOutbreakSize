## The partitioned patch renewal: one national renewal at the national
## trend, split across the patches by their force of infection.

@testitem "partitioned_patch_renewal: the patches partition the national renewal" begin
    using BVDOutbreakSize: partitioned_patch_renewal, renewal_infections,
        province_importation_kernel, PROVINCE_POPULATIONS

    g = [0.1, 0.3, 0.3, 0.2, 0.1]
    n = 60
    Rt_national = [1.6 - 0.8 * t / n for t in 1:n]
    ## Deviations far wider than the prior reaches, one patch at e^3 and one
    ## at e^-3, so a share that could go negative or overflow would.
    δ = [3.0, -3.0, 0.4]
    Rt = [Rt_national[t] * exp(δ[p] * sin(t / 9)) for p in 1:3, t in 1:n]
    seeds = [5.0 8.0; 0.0 0.0; 0.0 0.0]
    K = province_importation_kernel(PROVINCE_POPULATIONS[1:3])
    pools = [4.0e6, 7.6e6, 2.0e6]
    for ε in (0.0, 0.01, [0.02 * (1 + t / n) for _ in 1:3, t in 1:n])
        st = partitioned_patch_renewal(Rt, Rt_national, g, seeds, K, ε, pools)
        national = renewal_infections(
            Rt_national, g, vec(sum(seeds; dims = 1)), sum(pools)
        )
        @test size(st.infections) == (3, n)
        @test vec(sum(st.infections; dims = 1)) ≈ national
        @test all(>=(0), st.infections)
        @test all(isfinite, st.infections)
        @test all(>=(0), st.importation)
        @test st.infections[:, 1:2] == seeds
    end
    ## Uncoupled, the unseeded patches have no route to infections.
    st = partitioned_patch_renewal(Rt, Rt_national, g, seeds, 0 * K, 0.0, pools)
    @test all(iszero, st.infections[2:3, :])
end

@testitem "partitioned_patch_renewal: a common factor sets the size, not the split" begin
    using BVDOutbreakSize: partitioned_patch_renewal, free_patch_renewal,
        province_importation_kernel, PROVINCE_POPULATIONS

    g = [0.2, 0.3, 0.3, 0.2]
    n = 50
    Rt_national = fill(1.2, n)
    Rt = [Rt_national[t] * exp([0.3, -0.1, -0.2][p]) for p in 1:3, t in 1:n]
    seeds = [3.0 4.0; 1.0 1.0; 0.5 0.5]
    K = province_importation_kernel(PROVINCE_POPULATIONS[1:3])
    pools = fill(1.0e12, 3)
    st = partitioned_patch_renewal(Rt, Rt_national, g, seeds, K, 0.01, pools)
    ## Scaling every patch's reproduction number by one factor leaves the
    ## shares where they were.
    st2 = partitioned_patch_renewal(
        2.5 .* Rt, Rt_national, g, seeds, K, 0.01, pools
    )
    @test st2.infections ≈ st.infections
    ## The patches ran at a common rescaling of their own reproduction
    ## numbers, one on the seeded days.
    ratio = st.Rt_matrix ./ Rt
    @test ratio[:, 1:2] == ones(3, 2)
    @test all(t -> all(≈(ratio[1, t]), ratio[:, t]), 1:n)
    ## With no deviation and pools too large to deplete, the partition is the
    ## free patch renewal.
    flat = repeat(Rt_national', 3)
    part = partitioned_patch_renewal(
        flat, Rt_national, g, seeds, K, 0.01, pools
    )
    free = free_patch_renewal(flat, Rt_national, g, seeds, K, 0.01, pools)
    @test part.infections ≈ free.infections rtol = 1.0e-8
    @test part.importation ≈ free.importation rtol = 1.0e-8
end

@testitem "partitioned_patch_infection_model: prior draws sum to the single-patch renewal" begin
    using BVDOutbreakSize: partitioned_patch_infection_model, infection_model,
        PROVINCE_POPULATIONS
    using Turing: returned
    using Random: Xoshiro

    n, np, rt_start = 80, 3, 10
    m = partitioned_patch_infection_model(n, np; rt_start)
    pops = float.(PROVINCE_POPULATIONS[1:np])
    single = infection_model(n; rt_start, population = sum(pops))
    for i in 1:5
        draw = rand(Xoshiro(i), m)
        s = returned(m, draw)
        @test vec(sum(s.infections_matrix; dims = 1)) ≈ s.infections_total
        @test all(>=(0), s.infections_matrix)
        ## The single-patch model reads the same generation interval, growth
        ## and national walk from the draw and ignores the rest, so it must
        ## give the same national infections.
        ref = returned(single, draw)
        @test s.infections_total ≈ ref.infections
        @test s.C_T ≈ ref.C_T
    end
end
