using JET: test_opt
using BVDOutbreakSize

# Type-stability spot checks on the discrete-time renewal primitives.
# Kept here (in an isolated env) so JET version constraints do not bleed
# into the main test environment. These pure, allocation-light helpers
# back the generating infection process and the delay convolutions, so
# they must stay type-stable for the model to differentiate cleanly.
test_opt(BVDOutbreakSize.euler_lotka_r, (Float64, Vector{Float64});
    target_modules = (BVDOutbreakSize,))
test_opt(BVDOutbreakSize.convolve_delay,
    (Vector{Float64}, Vector{Float64});
    target_modules = (BVDOutbreakSize,))
test_opt(BVDOutbreakSize.renewal_infections,
    (Vector{Float64}, Vector{Float64}, Vector{Float64});
    target_modules = (BVDOutbreakSize,))
