# Fixed package constants: published scenarios, Ituri population and
# travel priors.

"""
    REPORT_SCENARIOS

Published point estimates of cumulative cases `C_T` from McCabe et
al. (Imperial College London, 20 May 2026 update), as `(label, value)`
tuples in the order they appear in Tables 1 and 2. These are the scenario
means; the matching 95% confidence intervals are carried by
[`REPORT_SCENARIOS_CI`](@ref).
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
    REPORT_SCENARIOS_CI

Published McCabe et al. scenario estimates WITH their reported 95%
confidence intervals, for the two situation-report vintages: the 18 May
2026 report [mccabe2026](@cite) and the 20 May 2026 update
[mccabe2026update](@cite). Each entry is a
`(date, label, mean, lower, upper)` tuple, where `date` is the report's
own cut-off date. Method 1 (geographic spread from exported cases and
travel volume) is unchanged between the two reports, so it is recorded
once under the 20 May vintage. Method 2 (back-calculation from deaths)
differs between the vintages: the 18 May report used 88 deaths and CFR
24/30/40%, the 20 May update used 131 deaths and the corrected CFR
26/33/40%. Confidence intervals are exact negative-binomial (Method 1)
and Poisson likelihood-profile (Method 2), as reported in Tables 1 and 2
of each report.
"""
const REPORT_SCENARIOS_CI = [
    ## Method 1 (geographic spread): identical across both reports.
    ("2026-05-20", "M1 Ituri, w=10 d", 470, 58, 1306),
    ("2026-05-20", "M1 Ituri, w=15 d", 313, 39, 870),
    ("2026-05-20", "M1 Ituri, w=20 d", 235, 29, 652),
    ("2026-05-20", "M1 +N. Kivu, w=10 d", 617, 76, 1718),
    ("2026-05-20", "M1 +N. Kivu, w=15 d", 412, 51, 1145),
    ("2026-05-20", "M1 +N. Kivu, w=20 d", 309, 38, 858),
    ## Method 2 (back-calculation from deaths), 18 May report: 88 deaths,
    ## CFR 24/30/40%.
    ("2026-05-18", "M2 τ=14 d, CFR 24%", 626, 503, 765),
    ("2026-05-18", "M2 τ=14 d, CFR 30%", 501, 402, 612),
    ("2026-05-18", "M2 τ=14 d, CFR 40%", 376, 302, 459),
    ("2026-05-18", "M2 τ= 7 d, CFR 24%", 1008, 808, 1230),
    ("2026-05-18", "M2 τ= 7 d, CFR 30%", 807, 649, 982),
    ("2026-05-18", "M2 τ= 7 d, CFR 40%", 605, 485, 789),
    ("2026-05-18", "M2 τ=21 d, CFR 24%", 531, 425, 645),
    ("2026-05-18", "M2 τ=21 d, CFR 30%", 425, 341, 517),
    ("2026-05-18", "M2 τ=21 d, CFR 40%", 319, 255, 387),
    ## Method 2, 20 May update: 131 deaths, corrected CFR 26/33/40%.
    ("2026-05-20", "M2 τ=14 d, CFR 26%", 860, 721, 1015),
    ("2026-05-20", "M2 τ=14 d, CFR 33%", 678, 568, 800),
    ("2026-05-20", "M2 τ=14 d, CFR 40%", 559, 469, 660),
    ("2026-05-20", "M2 τ= 7 d, CFR 26%", 1386, 1160, 1636),
    ("2026-05-20", "M2 τ= 7 d, CFR 33%", 1092, 914, 1289),
    ("2026-05-20", "M2 τ= 7 d, CFR 40%", 901, 754, 1062),
    ("2026-05-20", "M2 τ=21 d, CFR 26%", 730, 612, 862),
    ("2026-05-20", "M2 τ=21 d, CFR 33%", 575, 482, 679),
    ("2026-05-20", "M2 τ=21 d, CFR 40%", 474, 398, 560)
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
The doubling count `m` counts ONLY the cryptic-phase doublings (the origin to
the renewal-process start): the cryptic duration is `m·τ` and the total
outbreak age is `T = m·τ + τ_obs`, with `τ_obs` the observed window. A centre
of 3 places the prior cryptic phase at `m·τ ≈ 60 d` at the central 20-day
doubling and a prior seed of `2^m = 8` infections at the renewal start. The
genetic seeding bound pulls the lower tail of the outbreak age to sit at or
before the most recent common ancestor.
"""
const M_PRIOR_BASE = 3.0

"""
    RENEWAL_START_LEAD

Days the renewal start (the day the reproduction-number walk starts, where
the analytic cryptic phase hands off to the recursion) sits AFTER the
genetic TMRCA day, past the TMRCA's molecular-clock uncertainty. Placing the
renewal start a 14-day lead after the TMRCA — rather than exactly on it —
leaves the observed span `τ_obs = n − renewal_start` strictly shorter than
`tmrca_days`, so the genetic censored bound on the total age
`T = m·τ + τ_obs` stays informative (it pulls the origin to sit at or before
the MRCA, bounding the cryptic duration `m·τ` from below). Two weeks past the
TMRCA leaves room for the TMRCA's own molecular-clock uncertainty before
sustained transmission is treated as confidently established.
"""
const RENEWAL_START_LEAD = 14

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
