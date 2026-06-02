using JET: test_opt
using BVDOutbreakSize

# Type-stability spot check on the gamma CDF rrule plumbing. Kept here
# (in an isolated env) so JET version constraints do not bleed into
# the main test environment.
test_opt(BVDOutbreakSize._gamma_cdf, (Float64, Float64, Float64);
    target_modules = (BVDOutbreakSize,))

# Shared Gauss-Legendre `integrate` must specialise on the integrand and
# infer a concrete result; the previous SciML `solve` path returned `Any`.
_jet_integrand(x) = exp(0.013 * x) + 0.5
test_opt(BVDOutbreakSize.integrate,
    (typeof(_jet_integrand), Float64, Float64);
    target_modules = (BVDOutbreakSize,))
test_opt(BVDOutbreakSize.integrate,
    (typeof(_jet_integrand), Float64, Float64, Float64);
    target_modules = (BVDOutbreakSize,))
