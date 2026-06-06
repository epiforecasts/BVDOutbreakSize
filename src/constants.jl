# Fixed package constants: published scenarios, Ituri population and
# travel priors.

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
    M_PRIOR_BASE_DATE

Base date for the doubling-count prior centre: McCabe et al.'s first
report (18 May 2026), whose Method 2 central scenario of 501 cases implies
`m ≈ log2(501) ≈ 9`. Used by [`m_prior_centre`](@ref).
"""
const M_PRIOR_BASE_DATE = "2026-05-18"

"""
    M_PRIOR_DOUBLING_DAYS

Central doubling time (days) for the size and growth priors, from
Cuomo-Dannenburg & Ghafari's molecular-clock reanalysis (mean 15.2-24.5 d
across six clock assumptions). The doubling-count prior centre advances by
one doubling per `M_PRIOR_DOUBLING_DAYS` of elapsed time to the cut-off.
"""
const M_PRIOR_DOUBLING_DAYS = 20.0

"""
    M_PRIOR_BASE

Centre of the wide doubling-count prior in
[`exponential_growth_model`](@ref).
The doubling count sets the cryptic-phase length through `T = m·τ`, so the
centre is chosen to place the prior outbreak age near the
genetically-plausible window (about three months) rather than the much
longer span implied by a doubling count tied to the cut-off case total.
A wider doubling count would push the seeding time implausibly early.
"""
const M_PRIOR_BASE = 4.0

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
