# Fixed package constants: published scenarios, Ituri population and
# travel priors, observed counts, and the shared quadrature settings.

"""
    REPORT_SCENARIOS

Published point estimates of cumulative cases `C_T` from McCabe et
al. (Imperial College London, 20 May 2026 update), as `(label, value)`
tuples in the order they appear in Tables 1 and 2.
"""
const REPORT_SCENARIOS = [
    ("Method 1 Ituri, w=10 d", 470),
    ("Method 1 Ituri, w=15 d", 313),
    ("Method 1 Ituri, w=20 d", 235),
    ("Method 1 +N. Kivu, w=10", 617),
    ("Method 1 +N. Kivu, w=15", 412),
    ("Method 1 +N. Kivu, w=20", 309),
    ("Method 2 τ=14 d, CFR 26%", 860),
    ("Method 2 τ=14 d, CFR 33%", 678),
    ("Method 2 τ=14 d, CFR 40%", 559),
    ("Method 2 τ= 7 d, CFR 26%", 1386),
    ("Method 2 τ= 7 d, CFR 33%", 1092),
    ("Method 2 τ= 7 d, CFR 40%", 901),
    ("Method 2 τ=21 d, CFR 26%", 730),
    ("Method 2 τ=21 d, CFR 33%", 575),
    ("Method 2 τ=21 d, CFR 40%", 474)
]

"""
    ITURI_POPULATION

Source population for the Ituri Province (McCabe et al., Table 1).
"""
const ITURI_POPULATION = 4_392_200

"""
    ITURI_DAILY_TRAVEL

Default prior mean for the daily outbound traveller volume from
Ituri Province across seven points of entry.
"""
const ITURI_DAILY_TRAVEL = 1_871

"""
    ITURI_DAILY_TRAVEL_SD

Default prior SD for the daily outbound traveller volume, covering
point-of-entry-to-point-of-entry variation and reporting uncertainty
in the underlying mobility survey.
"""
const ITURI_DAILY_TRAVEL_SD = 200

"""
    GaussLegendre

Fixed Gauss-Legendre quadrature rule on the reference domain `[-1, 1]`,
carrying its `nodes` and `weights` (from
`FastGaussQuadrature.gausslegendre(n)`). Construct with the node count,
`GaussLegendre(; n = 64)`. The shared [`integrate`](@ref) reduces an
integrand against these fields directly; the nodes and weights are the
same rule the previous `Integrals.GaussLegendre` carried, so results are
unchanged.
"""
struct GaussLegendre{N, W}
    nodes::N
    weights::W
end

function GaussLegendre(; n = 250)
    nodes, weights = FastGaussQuadrature.gausslegendre(n)
    return GaussLegendre(nodes, weights)
end

"""
    DEATH_INTEGRAL_ALG

Gauss-Legendre quadrature scheme (`n = 64`) used for the deaths
onset-to-death convolution, the no-onward-transmission counterfactual,
and the forecast deaths integral.
"""
const DEATH_INTEGRAL_ALG = GaussLegendre(; n = 64)

"""
    CUMULATIVE_INTEGRAL_ALG

Gauss-Legendre quadrature scheme (`n = 32`) used for the at-risk
person-time export integral and the deaths-among-exports convolution
(outer and inner integrals).
"""
const CUMULATIVE_INTEGRAL_ALG = GaussLegendre(; n = 32)

"""
    DELAY_SUPPORT_K

Number of standard deviations beyond the mean used as the clustering
scale for a delay distribution in the onset-to-death convolution
integrals. `mean + DELAY_SUPPORT_K · std` is the width near the cut-off
over which the clustered [`integrate`](@ref) packs roughly half its
nodes, so the quadrature tracks the delay's scale as it is sampled.
"""
const DELAY_SUPPORT_K = 10

"""
    EXPORT_DELAY_GRID_POINTS

Number of evenly spaced grid points used to precompute the onset-to-death
CDF in [`ExportDeathDelay`](@ref).
"""
const EXPORT_DELAY_GRID_POINTS = 256
