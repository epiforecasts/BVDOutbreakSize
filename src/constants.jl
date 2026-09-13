# Fixed package constants: fit-id labels, published scenarios, Ituri
# population and travel priors.

"""
    JOINT_FIT

Value of a `fit` column identifying a row as the joint model's forecast or
estimate. A [`FROZEN_FIT`](@ref) row is drawn in the same role, being the
same model re-fit at a past cut-off rather than a separate fit.
"""
const JOINT_FIT = "joint"

"""
    BASELINE_FIT

Value of a `fit` column identifying a row as the persistence baseline,
which is no model's output.
"""
const BASELINE_FIT = "baseline"

"""
    FROZEN_FIT

Value of a `fit` column identifying a row as a frozen-fit forecast, the
joint model re-fit and evaluated at a past cut-off
(`forecast_frozen.csv`). That asset carries no `fit` column of its own, so
the label is applied at scoring time.
"""
const FROZEN_FIT = "frozen"

"""
    REPORT_SCENARIOS

Published point estimates of cumulative cases `C_T` from McCabe et
al. (Imperial College London, 20 May 2026 update), as `(label, value)`
tuples in the order they appear in Tables 1 and 2. These are the scenario
means. The matching 95% confidence intervals are carried by
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

Published McCabe et al. scenario estimates with their reported 95%
confidence intervals, for three vintages: the 18 May 2026 report
[mccabe2026](@cite), the 20 May 2026 update [mccabe2026update](@cite) and
the peer-reviewed Lancet Infectious Diseases publication
[mccabe2026lancet](@cite), whose inputs are as of 27 May 2026. Each entry
is a `(date, label, mean, lower, upper)` tuple, where `date` is that
vintage's own cut-off date. Method 1 (geographic spread from exported
cases and travel volume) is unchanged between the two Imperial reports,
so it is recorded once under the 20 May vintage. Method 2
(back-calculation from deaths) differs between the vintages: the 18 May
report used 88 deaths and CFR 24/30/40%, the 20 May update used 131
deaths and the corrected CFR 26/33/40%. Confidence intervals are exact
negative-binomial (Method 1) and Poisson likelihood-profile (Method 2),
as reported in Tables 1 and 2 of each report.

The 27 May Lancet vintage uses 240 deaths and three Uganda imports, and
varies the epidemic doubling time `T_d` (7/10/14 d) for both methods
rather than the earlier geographic window `w` / onset-to-death `τ`. Its
back-calculation fixes the onset-to-death gamma (mean 11.37 d, SD 5.41)
and assumes 30% of deaths are attributable to Ebola. The published paper
swaps the method numbers (its "method 1" is the back-calculation and its
"method 2" the geographic spread). This package keeps its own convention
(M1 = geographic spread, negative-binomial CIs; M2 = back-calculation,
Poisson CIs), confirmed by the paper's reported CI types.
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
    ("2026-05-20", "M2 τ=21 d, CFR 40%", 474, 398, 560),
    ## 27 May 2026 Lancet publication: geographic spread (Method 1 here =
    ## the paper's "method 2"), Table 2; doubling time T_d, source pop.
    ("2026-05-27", "M1 Ituri, T_d=10 d", 945, 196, 2274),
    ("2026-05-27", "M1 Ituri, T_d= 7 d", 1100, 228, 2647),
    ("2026-05-27", "M1 Ituri, T_d=14 d", 847, 176, 2037),
    ("2026-05-27", "M1 +N. Kivu, T_d=10 d", 1164, 241, 2800),
    ("2026-05-27", "M1 +N. Kivu, T_d= 7 d", 1354, 281, 3260),
    ("2026-05-27", "M1 +N. Kivu, T_d=14 d", 1042, 216, 2508),
    ## 27 May 2026 Lancet publication: back-calculation from 240 deaths
    ## (Method 2 here = the paper's "method 1"), Table 1; 30% of deaths
    ## attributable to Ebola, doubling time T_d, CFR 26/33/40%.
    ("2026-05-27", "M2 T_d=10 d, CFR 26%", 573, 502, 648),
    ("2026-05-27", "M2 T_d=10 d, CFR 33%", 451, 396, 511),
    ("2026-05-27", "M2 T_d=10 d, CFR 40%", 372, 327, 421),
    ("2026-05-27", "M2 T_d= 7 d, CFR 26%", 756, 664, 855),
    ("2026-05-27", "M2 T_d= 7 d, CFR 33%", 596, 523, 674),
    ("2026-05-27", "M2 T_d= 7 d, CFR 40%", 491, 432, 556),
    ("2026-05-27", "M2 T_d=14 d, CFR 26%", 471, 413, 533),
    ("2026-05-27", "M2 T_d=14 d, CFR 33%", 371, 326, 420),
    ("2026-05-27", "M2 T_d=14 d, CFR 40%", 306, 269, 346)
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

Median doubling time (days) for the size and growth priors, from the
BEAST X analysis (mbalaplacide2026, Exponential growth model,
11.7 d, 95% HPD 6.8-17.5). The doubling-count prior centre
advances by one doubling per `M_PRIOR_DOUBLING_DAYS` of elapsed time to
the cut-off.
"""
const M_PRIOR_DOUBLING_DAYS = 11.7

"""
    M_PRIOR_BASE

Base centre (at [`M_PRIOR_BASE_DATE`](@ref)) for the advancing doubling-count
prior centre used by the backfill fits via [`m_prior_centre`](@ref):
`m_0 = M_PRIOR_BASE + (as_of − M_PRIOR_BASE_DATE) / M_PRIOR_DOUBLING_DAYS`.
The doubling count `m` counts only the cryptic-phase doublings (the origin to
the renewal-process start): the cryptic duration is `m·τ` and the total
outbreak age is `T = m·τ + τ_obs`, with `τ_obs` the observed window. The
genetic seeding bound pulls the lower tail of the outbreak age to sit at or
before the most recent common ancestor. The main fit's `m` prior centre is
set directly in [`exponential_growth_model`](@ref), from field intelligence
and genetic evidence.
"""
const M_PRIOR_BASE = 3.0

"""
    RENEWAL_START_LEAD

Days the renewal start (the day the reproduction-number walk starts, where
the analytic cryptic phase hands off to the recursion) sits after the
genetic TMRCA day. Placing the renewal start a 14-day lead after the TMRCA,
rather than exactly on it, leaves the observed span
`τ_obs = n − renewal_start` strictly shorter than `tmrca_days`, so the
genetic censored bound on the total age `T = m·τ + τ_obs` stays informative:
it pulls the origin to sit at or before the MRCA, bounding the cryptic
duration `m·τ` from below. The lead accounts for the TMRCA's own
molecular-clock uncertainty before sustained transmission is treated as
confidently established.
"""
const RENEWAL_START_LEAD = 14

"""
    RT_WALK_LEAD

Days before the first situation report (`breakpoint`) at which the
reproduction-number random walk is allowed to start moving, rather than
holding `R_t` flat at `R0` until the report. A month lets the walk capture
transmission dynamics before the outbreak is first reported, while staying
floored at the renewal start so the walk never precedes the seeded
trajectory.
"""
const RT_WALK_LEAD = 28

"""
    RT_INTERVENTION_RAMP

Time scale in days of the logistic ramp over which the outbreak-response
intervention takes effect on `R_t`, centred on the `breakpoint`
([`sigmoid_ramp`](@ref)). A ramped rather than instantaneous step: the
response damps transmission over weeks, not on a single day. This constant
is the single source of truth for the ramp, referenced by both the model
(the `rt_walk_model` prior and `sigmoid_ramp`) and every `reconstruct_rt`
caller that rebuilds the daily `R_t` from the chain, so the reconstruction
cannot drift from the value the model fitted.
"""
const RT_INTERVENTION_RAMP = 21.0

"""
    ITURI_POPULATION

Source population for the Ituri Province (McCabe et al., Table 1).
"""
const ITURI_POPULATION = 4_392_200

"""
    PROVINCE_SOURCE_NAMES

Every province the situation reports carry a per-province row for, in the
order they first appear. These key the per-province blocks of
`data/observations.toml`. They are not the model's patches: see
[`PROVINCE_NAMES`](@ref) for those and [`PROVINCE_MEMBERS`](@ref) for how
the two relate.
"""
const PROVINCE_SOURCE_NAMES = ["ituri", "nord_kivu", "sud_kivu",
    "haut_uele", "tshopo", "bas_uele"]

"""
    PROVINCE_SOURCE_POPULATIONS

Resident population of each province in [`PROVINCE_SOURCE_NAMES`](@ref)
order: 2019 figures from the Democratic Republic of the Congo's Institut
National de la Statistique, *Annuaire statistique RDC 2020* (March 2021),
as tabulated at
<https://en.wikipedia.org/wiki/Provinces_of_the_Democratic_Republic_of_the_Congo>
(accessed 12 September 2026).

One source for all six rather than the best figure for each. Only the
relative sizes enter the model, through the importation kernel and the
per-capita testing covariate, so consistency between provinces matters more
than the accuracy of any one of them.
"""
const PROVINCE_SOURCE_POPULATIONS = [4_008_000, 7_574_000, 6_565_000,
    2_046_000, 2_582_000, 1_250_000]

"""
    PROVINCE_SOURCE_CAPITALS

Latitude and longitude of each province's capital in
[`PROVINCE_SOURCE_NAMES`](@ref) order, as `(latitude, longitude)` in decimal
degrees north and east: Bunia, Goma, Bukavu, Isiro, Kisangani and Buta.
Coordinates from GeoNames (<https://www.geonames.org>), the source for the
distance term in [`province_importation_kernel`](@ref).

The capital stands in for the province. That is coarse, but it is the level
the data are reported at, and the provinces are far enough apart that the
ordering of the distances between them does not depend on the choice of
point within each one.
"""
const PROVINCE_SOURCE_CAPITALS = [
    (1.56667, 30.25000),    # Bunia, Ituri
    (-1.67918, 29.22195),   # Goma, Nord-Kivu
    (-2.50000, 28.86667),   # Bukavu, Sud-Kivu
    (2.77374, 27.61674),    # Isiro, Haut-Uele
    (0.51528, 25.19099),    # Kisangani, Tshopo
    (2.78594, 24.73876)     # Buta, Bas-Uele
]

"""
    PROVINCE_NAMES

The patches of the meta-population model, in patch order. The first entry is
the primary patch: the origin of the outbreak, the reference for the
per-patch reproduction-number deviations in [`patch_rt_model`](@ref), and
the reference for the Uganda export propensities in
[`province_export_pressure_model`](@ref).

Three provinces are patches in their own right and the rest are pooled into
`other`. Ituri, Nord-Kivu and Haut-Uele carry signal: at the cut-off they
hold 5508, 1139 and 264 confirmed cases. Sud-Kivu, Tshopo and Bas-Uele hold
3, 24 and 4 between them, and Sud-Kivu has reported no new confirmed case
since 26 May. Giving each of those its own reproduction number and its own
ascertainment would sample dimensions nothing informs, and their estimates
would be the deviation prior read back. Pooled they are one weak patch,
which is what they are.
"""
const PROVINCE_NAMES = ["ituri", "nord_kivu", "haut_uele", "other"]

"""
    PROVINCE_LABELS

Display names for the patches, in [`PROVINCE_NAMES`](@ref) order, for table
rows and figure panels.
"""
const PROVINCE_LABELS = ["Ituri", "Nord-Kivu", "Haut-Uele",
    "Other provinces"]

"""
    PROVINCE_MEMBERS

Which source provinces each patch pools, keyed by [`PROVINCE_NAMES`](@ref)
and valued in [`PROVINCE_SOURCE_NAMES`](@ref). Every source province belongs
to exactly one patch, so the patches partition the national totals and the
composition likelihoods stay conditional on them.
"""
const PROVINCE_MEMBERS = Dict(
    "ituri" => ["ituri"],
    "nord_kivu" => ["nord_kivu"],
    "haut_uele" => ["haut_uele"],
    "other" => ["sud_kivu", "tshopo", "bas_uele"])

## Index of each patch's member provinces in `PROVINCE_SOURCE_NAMES`, built
## once so the pooled population and capital below, and the pooled increments
## in `province_increment_matrix`, all read the same membership.
const _PROVINCE_MEMBER_IDX = [[findfirst(==(m), PROVINCE_SOURCE_NAMES)
                               for m in PROVINCE_MEMBERS[name]] for name in PROVINCE_NAMES]

"""
    PROVINCE_POPULATIONS

Resident population of each patch, in [`PROVINCE_NAMES`](@ref) order, summed
over the provinces it pools. Used to put the per-province testing effort on
a per-capita scale (the covariate for the provincial ascertainment) and to
weight the between-province importation kernel.
"""
const PROVINCE_POPULATIONS = [sum(PROVINCE_SOURCE_POPULATIONS[idx])
                              for idx in _PROVINCE_MEMBER_IDX]

"""
    PROVINCE_CAPITALS

Representative point of each patch, in [`PROVINCE_NAMES`](@ref) order, as
`(latitude, longitude)` in decimal degrees north and east. A patch holding
one province takes its capital from [`PROVINCE_SOURCE_CAPITALS`](@ref); a
pooled patch takes the population-weighted mean of its members' capitals.

A weighted mean rather than a member's capital, because the pooled patch is
a stand-in for several places at once and the kernel asks where its
population is. The pooled patch here spans Bukavu, Kisangani and Buta, so
its point sits between them and its distance to the epicentre is a weighted
compromise rather than any one province's.
"""
const PROVINCE_CAPITALS = [(
                               sum(PROVINCE_SOURCE_POPULATIONS[i] *
                                   PROVINCE_SOURCE_CAPITALS[i][1]
                               for i in idx) / sum(PROVINCE_SOURCE_POPULATIONS[idx]),
                               sum(PROVINCE_SOURCE_POPULATIONS[i] *
                                   PROVINCE_SOURCE_CAPITALS[i][2]
                               for i in idx) / sum(PROVINCE_SOURCE_POPULATIONS[idx]))
                           for idx in _PROVINCE_MEMBER_IDX]

"""
    PROVINCE_DISTANCE_DECAY

Exponent on the distance term of [`province_importation_kernel`](@ref). One
is the conventional gravity value, and it is fixed rather than sampled: the
importation intensity it scales is already weakly identified against the
secondary provinces' seeds, so a second free parameter on the same term
would not be determined by anything.
"""
const PROVINCE_DISTANCE_DECAY = 1.0

"""
    haversine_km(a, b)

Great-circle distance in kilometres between two `(latitude, longitude)`
points in decimal degrees, on a spherical Earth of radius 6371 km.
"""
function haversine_km(a::Tuple{<:Real, <:Real}, b::Tuple{<:Real, <:Real})
    φ1, λ1 = deg2rad(a[1]), deg2rad(a[2])
    φ2, λ2 = deg2rad(b[1]), deg2rad(b[2])
    h = sin((φ2 - φ1) / 2)^2 +
        cos(φ1) * cos(φ2) * sin((λ2 - λ1) / 2)^2
    return 2 * 6371.0 * asin(sqrt(clamp(h, 0.0, 1.0)))
end

"""
    province_distance_matrix(capitals = PROVINCE_CAPITALS)

Great-circle distances in kilometres between every pair of province
capitals, as a symmetric matrix with a zero diagonal. Built from
[`PROVINCE_CAPITALS`](@ref) by [`haversine_km`](@ref).
"""
function province_distance_matrix(
        capitals::AbstractVector = PROVINCE_CAPITALS)
    np = length(capitals)
    D = zeros(Float64, np, np)
    @inbounds for p in 1:np, q in 1:np

        p == q && continue
        D[p, q] = haversine_km(capitals[p], capitals[q])
    end
    return D
end

"""
    province_importation_kernel(pops = PROVINCE_POPULATIONS; distances, decay)

Between-province importation kernel `K`, where `K[p, q]` is the relative
rate of infectious travel from province `q` into province `p`. Diagonal is
zero (no self-importation). The overall intensity is carried by the sampled
`ε` in [`patch_infection_model`](@ref), so only the relative structure
matters here.

This is a gravity kernel: travel from `q` to `p` scales with the destination
population and falls with the distance between the two capitals,

```math
K_{p,q} \\propto \\frac{N_p}{d_{p,q}^{\\gamma}},
```

with `γ` fixed at [`PROVINCE_DISTANCE_DECAY`](@ref) and `d` from
[`province_distance_matrix`](@ref). Pass `distances` as a zero matrix to
recover the population-only kernel.

Each column is scaled so that its off-diagonal entries sum to `1 - N_q/N`,
the share of the country that is not `q` itself. That is what the kernel
summed to before the distance term was added, so the distance changes where
a province's exported transmission lands without changing how much of it
leaves. `ε` keeps its meaning, and every column's off-diagonal sum stays
below one at any `ε` in `[0, 1]`, which is what stops a province exporting
more transmission than it generates.

There is no origin-destination or mobility data for this outbreak, so the
kernel is still a structural assumption rather than a measurement, and `ε`
is weakly identified against the secondary-patch seeds (both can raise a
secondary province's early incidence). Treat the split between imported and
locally-seeded infections as poorly determined even though their sum is not.
"""
function province_importation_kernel(
        pops::AbstractVector = PROVINCE_POPULATIONS;
        distances::AbstractMatrix = province_distance_matrix(
            PROVINCE_CAPITALS[1:min(length(pops), end)]),
        decay::Real = PROVINCE_DISTANCE_DECAY)
    np = length(pops)
    tot = sum(pops)
    K = zeros(Float64, np, np)
    size(distances) == (np, np) || error(
        "province_importation_kernel: `distances` is $(size(distances)) " *
        "but there are $np provinces.")
    @inbounds for q in 1:np
        ## Relative pull of each destination from origin `q`, by destination
        ## size and by how far it is. A zero distance matrix leaves the
        ## population-only kernel.
        pull = zeros(Float64, np)
        for p in 1:np
            p == q && continue
            d = distances[p, q]
            pull[p] = d > 0 ? pops[p] / d^decay : Float64(pops[p])
        end
        s = sum(pull)
        s > 0 || continue
        ## Hold the column total at the pre-distance value, so the distance
        ## redistributes a province's exports without changing their volume.
        outflow = 1 - pops[q] / tot
        for p in 1:np
            K[p, q] = outflow * pull[p] / s
        end
    end
    return K
end

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
    CHAMLA_CONFIRMED_CENTRAL

Central-scenario cumulative laboratory-confirmed case projection of Chamla et
al. [chamla2026](@cite) (WHO Regional Office for Africa, Lancet Infectious
Diseases, 2026), as `(date, median, lower_90, upper_90)` tuples from their
Table 1 (mean accepted `R₀ = 1.71`). Their stochastic SEIRD ensemble is
calibrated by simulation filtering to 598 cumulative confirmed cases on 8 June
2026, with the reporting fraction fixed at 1.0, so these are projected
confirmed cases (a floor on the true size) rather than the
ascertainment-corrected cumulative cases that this package's `C_T` and McCabe
et al. estimate. The 8/10 June row is their calibration anchor. The dates
from 24 June on are forward projections. The bounds are 90% prediction
intervals.
"""
const CHAMLA_CONFIRMED_CENTRAL = [
    ("2026-05-18", 294, 212, 375),
    ("2026-05-27", 399, 289, 502),
    ("2026-06-10", 648, 470, 812),
    ("2026-06-24", 990, 709, 1293),
    ("2026-07-22", 2114, 1450, 2980),
    ("2026-08-19", 4242, 2748, 6528),
    ("2026-09-16", 8210, 5063, 13498)
]

"""
    CHAMLA_CONFIRMED_W12

Chamla et al. [chamla2026](@cite) week-12 (24 June 2026) cumulative
confirmed-case projection under their three transmissibility scenarios, as
`(scenario, median, lower_90, upper_90)` tuples (low `R₀ = 1.42`, central
`R₀ = 1.71`, high `R₀ = 2.08`). Week 12 is the forward horizon closest to
this analysis's current cut-off, so the scenario spread sits beside the
observed confirmed count and our matched-date projection. The bounds are
90% prediction intervals. The same estimand caveat as
[`CHAMLA_CONFIRMED_CENTRAL`](@ref) applies: these are confirmed cases, not
ascertainment-corrected total cases.
"""
const CHAMLA_CONFIRMED_W12 = [
    ("Chamla low (R₀=1.42)", 870, 641, 1133),
    ("Chamla central (R₀=1.71)", 990, 709, 1293),
    ("Chamla high (R₀=2.08)", 1364, 975, 1807)
]
