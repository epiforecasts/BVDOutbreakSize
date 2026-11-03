using JET: test_opt
using BVDOutbreakSize

# Type-stability spot checks on the discrete-time renewal primitives.
# They run in an isolated environment so JET's version constraints stay out
# of the main test environment.
test_opt(
    BVDOutbreakSize.euler_lotka_r, (Float64, Vector{Float64});
    target_modules = (BVDOutbreakSize,)
)
test_opt(
    BVDOutbreakSize.convolve_delay,
    (Vector{Float64}, Vector{Float64});
    target_modules = (BVDOutbreakSize,)
)
test_opt(
    BVDOutbreakSize.renewal_infections,
    (Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64);
    target_modules = (BVDOutbreakSize,)
)
