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
The doubling count `m` now counts ONLY the CRYPTIC-phase doublings (origin →
anchor): the cryptic duration is `m·τ` and the total outbreak age is
`T = m·τ + τ_obs`, with `τ_obs = n − anchor` the observed window. With the
molecular-clock doubling `τ ≈ 20 d` (from `R0 ≈ 1.6` through Euler–Lotka)
and the anchor a small lead after the genetic TMRCA (`τ_obs ≈ 63 d`), a
centre of 2 places the cryptic duration at `m·τ ≈ 40 d`, so the prior
outbreak age `T ≈ 103 d` lands mid-window in the genetically-plausible span
(origin roughly Feb–Mar, `T ≈ 70–130 d`). The SD-3 spread keeps `T` wide and
prior-dominated; the genetic seeding bound pulls the lower tail to sit at or
before the MRCA. The earlier centre (4) was tied to a `T = m·τ` that
counted the whole age, so it over-counts now that `τ_obs` is added
separately.
"""
const M_PRIOR_BASE = 2.0

"""
    SEEDING_ANCHOR_LEAD

Days the renewal anchor (the day the reproduction-number walk starts, where
the analytic cryptic phase hands off to the recursion) sits AFTER the
genetic TMRCA day, past the TMRCA's molecular-clock uncertainty. Placing the
anchor a small lead after the TMRCA — rather than exactly on it — leaves the
observed span `τ_obs = n − anchor` strictly shorter than `tmrca_days`, so the
genetic censored bound on the total age `T = m·τ + τ_obs` stays informative
(it pulls the origin to sit at or before the MRCA, bounding the cryptic
duration `m·τ` from below). A few days past the TMRCA is where sustained
transmission can be treated as confidently established.
"""
const SEEDING_ANCHOR_LEAD = 7

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
