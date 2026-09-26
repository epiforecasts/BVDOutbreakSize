# Health-zone composition model, the second stage of a two-stage Markov
# melding. Stage one is the four-patch joint fit. Stage two takes its
# posterior as the prior on the quantity the two stages share and fits the
# per-zone confirmed-case and death tables as a within-patch composition,
# with no feedback to stage one.
# The shared quantity is the province model's weekly infections in each
# patch. The zone stage samples it from a multivariate normal fitted to the
# province model's draws (`zone_meld_block`) and conditions on the daily
# patch trajectory it implies. The province model's uncertainty, and the
# correlation it learns across patches, reach every zone quantity through
# the fit. The functions here build the inputs from a
# parent chain (`zone_fit_inputs`), run the share renewal and its
# observation model as plain functions the model and the render both call,
# define the Turing model (`bvd_zone`) and fit it (`fit_zone`).

## Truncation lags of the stage-1 delay PMFs, the defaults of
## `patch_infection_model` (generation interval, incubation) and
## `lab_delay_model` (receipt), so the PMFs built here are discretised as
## the parent fit discretised them.
const ZONE_GI_NMAX = cdf_nmax(Gamma(2.71, 5.65))
const ZONE_INCUBATION_NMAX = cdf_nmax(lognormal_meansd(6.3, 3.5))
const ZONE_RECEIPT_NMAX = cdf_nmax(lognormal_meansd(4.5, 4.0))

## Chain keys of the parent quantities the cut reads. The joint attaches its
## latent submodel unprefixed, so the generation interval and incubation
## sit under their submodel names, and the receipt delay under the
## confirmed stream's.
const _ZONE_PARENT_KEYS = (
    gi_alpha = Symbol("gi_state.α"),
    gi_theta = Symbol("gi_state.θ"),
    inc_mean = Symbol("inc_state.delay_mean"),
    inc_sd = Symbol("inc_state.delay_sd"),
    receipt_mean = Symbol("confirmed_state.receipt_state.d.delay_mean"),
    receipt_sd = Symbol("confirmed_state.receipt_state.d.delay_sd"),
    infections = :infections_patch,
    death_confirmation = :onset_to_death_confirmation_pmf,
    C_T = :C_T,
    importation_epsilon = :importation_epsilon_patch,
    importation = :importation_patch,
    ascertainment_sd = :province_ascertainment_sd,
    province_ascertainment = :province_ascertainment,
    province_severity = :province_cfr_relative,
    rho = :province_composition_rho,
    rho_death = :province_death_composition_rho,
    drift_sd = :region_drift_sd,
    correlation = :region_corr_primary_secondary,
)

## Every parent key the zone stage and its diagnostics read, so an extract
## of these draws stands in for the full chain.
const _ZONE_EXTRACT_KEYS = (values(_ZONE_PARENT_KEYS)..., :province_shares)

"""
$(TYPEDSIGNATURES)

The draws the zone stage reads from a parent `bvd_joint` chain, as plain
arrays: one iterations-by-chains matrix per key of `_ZONE_EXTRACT_KEYS` the
chain carries (a vector-valued quantity holds one vector per cell), under
`draws`, with the number of grid days `n` the patch infections span (their
length over `n_patches`), the keys present and the `source` label. The extract is a few megabytes where
the chain is tens, serialises without a chains library, and is accepted
wherever the zone stage takes a parent chain ([`zone_fit_inputs`](@ref),
[`zone_forecast_draws`](@ref)). The fit script writes one next to a cached
fit (`docs/fits/summary.jl`).
"""
function zone_parent_extract(
        chn; source::AbstractString = "",
        n_patches::Integer = length(PROVINCE_NAMES)
    )
    draws = Dict{Symbol, Any}()
    for k in _ZONE_EXTRACT_KEYS
        _has_key(chn, k) || continue
        draws[k] = collect(chn[k])
    end
    haskey(draws, _ZONE_PARENT_KEYS.infections) || error(
        "zone_parent_extract: the chain carries no `infections_patch`. It " *
            "must be a `bvd_joint` chain sampled with the patch structure on."
    )
    first_draw = first(vec(draws[_ZONE_PARENT_KEYS.infections]))
    n = length(first_draw) ÷ n_patches
    return (;
        draws, n, keys = sort!(collect(keys(draws)); by = string),
        source = String(source), created = string(now()),
    )
end

## A parent given as an extract reads through its `draws`. A chain passes
## through unchanged.
_zone_parent_chain(x) = x isa NamedTuple && haskey(x, :draws) ? x.draws : x

## --- Pure building blocks ------------------------------------------------

"""
$(TYPEDSIGNATURES)

Initial zone shares at the grid start: a within-patch softmax of the
centred standard-normal draws `z_w` at a fixed `scale`,
`w_z = exp(scale (z_z − mean_p z)) / Σ_p`. The centring removes the flat
direction of the softmax. Returns a vector over the zones, summing to one
within every patch range.
"""
function zone_initial_shares(
        z_w::AbstractVector,
        patch_ranges::AbstractVector{<:UnitRange}, scale::Real
    )
    Tp = promote_type(eltype(z_w), typeof(float(scale)))
    w = zeros(Tp, length(z_w))
    @inbounds for zs in patch_ranges
        isempty(zs) && continue
        m = zero(Tp)
        for z in zs
            m += z_w[z]
        end
        m /= length(zs)
        tot = zero(Tp)
        for z in zs
            w[z] = exp(scale * (z_w[z] - m))
            tot += w[z]
        end
        for z in zs
            w[z] /= tot
        end
    end
    return w
end

"""
$(TYPEDSIGNATURES)

The linear map from knot values to the days `t0 … n`, as an `(n_days ×
n_knots)` weight matrix `W` with `δ_daily = W δ_knotsᵀ`, so the
interpolation of every zone is one matrix product. Each column is
[`interpolate_knots`](@ref) applied to a unit vector. The reverse pass
over a per-day interpolation loop costs far more than the product.
"""
function zone_interpolation_weights(
        knots::AbstractVector{<:Integer},
        t0::Integer, n::Integer
    )
    nb = length(knots)
    nd = n - t0 + 1
    W = zeros(Float64, nd, nb)
    unit = zeros(Float64, nb)
    for k in 1:nb
        fill!(unit, 0.0)
        unit[k] = 1.0
        daily = interpolate_knots(unit, knots, n)
        for j in 1:nd
            W[j, k] = daily[t0 + j - 1]
        end
    end
    return W
end

"""
$(TYPEDSIGNATURES)

The delay convolution over the days `t0 … n` as an `(n_days × n_days)`
lower-triangular Toeplitz matrix `F` with `F[j, j − s] = f_s` (lag 0 on
the diagonal), so the convolution of every zone's infections is one
matrix product `F I`. The days before `t0` enter separately through
[`zone_report_pre_rows`](@ref).
"""
function zone_delay_operator(f::AbstractVector, nd::Integer)
    F = zeros(Float64, nd, nd)
    for j in 1:nd, s in 0:min(length(f) - 1, j - 1)

        F[j, j - s] = f[s + 1]
    end
    return F
end

"""
$(TYPEDSIGNATURES)

The pre-`t0` reporting term of each zone's patch on the days `t0 … n`, as
an `(n_days × n_zones)` matrix so that a zone's daily reports are
`F I[:, z] + w_z(t_0) pre[:, z]`. `report_pre` is the `(n_patches × n)`
term of [`zone_fixed_terms`](@ref).
"""
function zone_report_pre_rows(
        report_pre::AbstractMatrix,
        patch_of_zone::AbstractVector{<:Integer}, t0::Integer, n::Integer
    )
    nd = n - t0 + 1
    out = zeros(Float64, nd, length(patch_of_zone))
    for (z, p) in enumerate(patch_of_zone), j in 1:nd

        out[j, z] = report_pre[p, t0 + j - 1]
    end
    return out
end

"""
$(TYPEDSIGNATURES)

The parts of the zone renewal and its delays that involve only patch
infections before the grid start `t0`, which are constant under the cut.
Every zone holds its initial share over those days, so each term enters a
zone's series multiplied by `w_z(t0)` only.

Returns, with `I_bar` `(n_patches × n)` and PMFs `g` (lag 1) and `f`
(lag 0):

- `force_pre[p, t] = Σ_{s ≥ 1, t − s < t0} g_s Ī_p(t − s)`, the pre-`t0`
  force of infection on every day;
- `report_pre[p, t] = Σ_{s ≥ 0, t − s < t0} f_s Ī_p(t − s)`, the
  pre-`t0` contribution to the daily expected reports;
- `report_pre_cum[p, t] = Σ_{t' ≤ min(t, t0 − 1)} report_pre[p, t']`, the
  expected reports accrued over the days before `t0`, so a vintage window
  opening before the grid start bins them by difference;
- `infections_pre[p] = Σ_{t < t0} Ī_p(t)`, the patch infections before
  the grid start.
"""
function zone_fixed_terms(
        I_bar::AbstractMatrix, g::AbstractVector,
        f::AbstractVector, t0::Integer
    )
    np, n = size(I_bar)
    Tp = promote_type(eltype(I_bar), eltype(g), eltype(f))
    force_pre = zeros(Tp, np, n)
    report_pre = zeros(Tp, np, n)
    report_pre_cum = zeros(Tp, np, n)
    infections_pre = zeros(Tp, np)
    @inbounds for p in 1:np
        for t in 1:n
            acc = zero(Tp)
            for s in 1:min(length(g), t - 1)
                t - s < t0 || continue
                acc += g[s] * I_bar[p, t - s]
            end
            force_pre[p, t] = acc
            acc = zero(Tp)
            for s in 0:min(length(f) - 1, t - 1)
                t - s < t0 || continue
                acc += f[s + 1] * I_bar[p, t - s]
            end
            report_pre[p, t] = acc
        end
        run = zero(Tp)
        for t in 1:n
            t < t0 && (run += report_pre[p, t])
            report_pre_cum[p, t] = run
        end
        for t in 1:min(t0 - 1, n)
            infections_pre[p] += I_bar[p, t]
        end
    end
    return (; force_pre, report_pre, report_pre_cum, infections_pre)
end

## --- The shared quantity -------------------------------------------------

"""
$(TYPEDSIGNATURES)

Midpoint day of each weekly window of `knots`. Window `k` spans
`(knots[k], knots[k + 1]]`, so its midpoint is the mean of the two, and
there are `length(knots) - 1` of them. A window's infections are
attributed to its midpoint when the parent's weekly deformation is
interpolated to a daily curve.
"""
function zone_week_midpoints(knots::AbstractVector{<:Integer})
    return [(knots[k] + knots[k + 1]) ÷ 2 for k in 1:(length(knots) - 1)]
end

"""
$(TYPEDSIGNATURES)

The quantity the two stages share: the province model's weekly infections
in each patch, as one multivariate normal on the log scale fitted to its
draws.

`infections` holds one flattened `(n_patches × n)` infection matrix per
parent draw, `knots` the weekly grid and `I_bar` the parent's mean patch
infections. Window `k` of patch `p` spans `(knots[k], knots[k + 1]]`, and
the pair is kept when the parent's mean infections over that window reach
`min_infections`. That drops the windows in which a patch is not yet
seeded, whose log sums are not normal. Over the kept pairs the log sums
have mean `m` and sample covariance `Σ`. The Cholesky factor `L` is taken
after adding `ridge` of each cell's own variance to the diagonal, so that
`log S = m + L η` with `η ∼ N(0, I_d)` reproduces the parent's joint
posterior on the shared quantity to first order. `Σ` is rank deficient
whenever the parent carries fewer draws than there are cells. The
factorisation then blends toward `diag(Σ)` at the smallest weight that
succeeds ([`_zone_meld_factor`](@ref)). Shrinking toward the diagonal of
`Σ` rather than toward the identity keeps every week's own variance and
gives up only the correlations, which is the direction the parent's
posterior supports least.

Only the deviation `a = L η` is used, never `m`, because the zone stage
multiplies the parent's own mean curve by `exp(a)`. `weights` is the
linear map from `a` to that daily log deviation on every patch and day, an
`(n_patches · n) × d` matrix whose rows are flattened patch-major as the
parent stores `infections_patch`. It interpolates each patch's kept
midpoints ([`zone_week_midpoints`](@ref)) and holds flat outside them.
Every patch's first kept midpoint falls on or after the grid start `t0`,
so the deviation is constant over the days before `t0`, which is what lets
the pre-`t0` terms of [`zone_fixed_terms`](@ref) be rescaled by one factor
per patch rather than rebuilt on every forward pass. A patch whose first
kept midpoint fell before `t0` would break that, and is an error.

Returns `(; weights, L, cells_patch, cells_week, midpoints, d, mean_log,
log_sums)`, `log_sums` holding the per-draw log sums `(n_draws × d)`.
"""
## Cholesky factor of a covariance, conditioned by a ridge of `ridge` of
## each cell's own variance and then blended toward that diagonal by the
## least weight that factorises. A covariance that never factorises falls
## back to its diagonal, which drops the correlations and keeps the
## marginal variances.
function _zone_meld_factor(Σ::AbstractMatrix, d::Integer; ridge::Real = 1.0e-6)
    dg = Diagonal(max.(diag(Σ), eps()))
    λ = 0.0
    while λ <= 1.0
        M = Symmetric((1 - λ) * Σ + λ * Matrix(dg) + ridge * Matrix(dg))
        F = cholesky(M; check = false)
        issuccess(F) && return Matrix(F.L)
        λ = λ == 0.0 ? 1.0e-4 : 10λ
    end
    return Matrix(cholesky(Symmetric(Matrix(dg))).L)
end

function zone_meld_block(
        infections::AbstractVector, np::Integer,
        n::Integer, knots::AbstractVector{<:Integer},
        I_bar::AbstractMatrix, t0::Integer;
        min_infections::Real = 1.0, ridge::Real = 1.0e-6
    )
    mids = zone_week_midpoints(knots)
    cells_patch = Int[]
    cells_week = Int[]
    for p in 1:np, k in eachindex(mids)

        lo, hi = knots[k] + 1, min(knots[k + 1], n)
        lo <= hi || continue
        sum(@view I_bar[p, lo:hi]) >= min_infections || continue
        push!(cells_patch, p)
        push!(cells_week, k)
    end
    d = length(cells_patch)
    d > 0 || return (;
        weights = zeros(Float64, np * n, 0),
        L = zeros(Float64, 0, 0), cells_patch, cells_week,
        midpoints = mids, d = 0, mean_log = Float64[],
        log_sums = zeros(Float64, length(infections), 0),
    )
    S = Matrix{Float64}(undef, length(infections), d)
    for (i, v) in enumerate(infections)
        M = reshape(Float64.(v), np, n)
        for c in 1:d
            k = cells_week[c]
            lo, hi = knots[k] + 1, min(knots[k + 1], n)
            S[i, c] = log(safe_rate(sum(@view M[cells_patch[c], lo:hi])))
        end
    end
    mean_log = vec(mean(S; dims = 1))
    ## The sample covariance is rank deficient whenever the parent carries
    ## fewer draws than there are cells, so the factorisation is pulled
    ## toward the diagonal by the least amount that succeeds. Shrinking
    ## toward the diagonal of the sample covariance rather than toward the
    ## identity keeps every week's own marginal variance and gives up only
    ## the correlations, which is the direction the parent's own posterior
    ## is least able to support.
    L = _zone_meld_factor(cov(S), d; ridge)
    weights = zeros(Float64, np * n, d)
    for p in 1:np
        cs = findall(==(p), cells_patch)
        isempty(cs) && continue
        W = zone_interpolation_weights(mids[cells_week[cs]], 1, n)
        for (j, c) in enumerate(cs), t in 1:n

            weights[(t - 1) * np + p, c] = W[t, j]
        end
        for t in 1:min(t0, n), c in cs

            weights[(t - 1) * np + p, c] == weights[p, c] || error(
                "zone_meld_block: patch $p's weekly deformation is not " *
                    "constant before the grid start (day $t0); its first kept " *
                    "week midpoint must fall on or after it."
            )
        end
    end
    return (;
        weights, L, cells_patch, cells_week, midpoints = mids, d,
        mean_log, log_sums = S,
    )
end

"""
$(TYPEDSIGNATURES)

The fixed inputs [`bvd_zone`](@ref) reads past the cut-off `n` when run with
a [`ForecastHorizon`](@ref) of `horizon` days, built from the fitted
shared quantity `meld` ([`zone_meld_block`](@ref)) and the parent's
posterior-predictive draws `forecast` ([`forecast_draws`](@ref) on the
parent), which carry one draw per parent draw in the same order.

The shared quantity is extended over the forecast weeks. Each parent draw's
log weekly patch infections over the fitted cells of `meld` and over the
future windows `(n, n + 7]`, … (kept where the parent's mean reaches
`min_infections`) are stacked, and their covariance factorised with the
fitted block held at `meld.L`:

```math
L = \\begin{pmatrix} L_{11} & 0 \\\\ L_{21} & L_{22} \\end{pmatrix},
\\qquad L_{21} = Σ_{21} L_{11}^{-\\top},
\\qquad L_{22} L_{22}^\\top = Σ_{22} − L_{21} L_{21}^\\top.
```

A fresh `η_f ∼ N(0, I)` then draws the future weeks from the parent's
posterior conditional on the fitted draw `η`, and the fitted model is
unchanged. The daily deformation keeps the fitted rows up to `n` and
interpolates from day `n` to the future midpoints after it. The mean curve
past `n` is the parent's mean log predicted infections, and every delay
operator and pre-`t0` term is rebuilt on the longer grid. The import
fractions hold their value at `n`.

`totals` holds each parent draw's predicted confirmed cases per patch over
`(n, n + horizon]` (`forecast_province_confirmed`), the counts the zone
forecast splits. `horizon` must be one of the parent's weekly forecast
vintages.
"""
function zone_forecast_block(
        forecast, meld, I_bar::AbstractMatrix, g::AbstractVector,
        f::AbstractVector, death_pmf::AbstractVector, t0::Integer,
        knots::AbstractVector{<:Integer}, patch_of_zone, mixing;
        horizon::Integer = 7, week::Integer = 7,
        min_infections::Real = 1.0, ridge::Real = 1.0e-6
    )
    np, n = size(I_bar)
    H = Int(horizon)
    inf = _draw_vectors(forecast, :forecast_infections_patch)
    Hp = length(first(inf)) ÷ np
    Hp >= H || error(
        "zone_forecast_block: the parent forecast covers $Hp days, fewer " *
            "than the $H-day horizon."
    )
    ndraw = length(inf)
    size(meld.log_sums, 1) == ndraw || error(
        "zone_forecast_block: $ndraw parent forecast draws for " *
            "$(size(meld.log_sums, 1)) fitted draws; draw the forecast from " *
            "the chain the zone inputs were built from."
    )
    future = [reshape(Float64.(v), np, Hp)[:, 1:H] for v in inf]
    I_f = exp.(reduce(+, (log.(safe_rate.(M)) for M in future)) ./ ndraw)
    I_ext = hcat(I_bar, I_f)
    ## Future windows and their midpoints, on the days past the cut-off.
    edges = vcat(n, future_knot_days(n, H; week))
    cells_patch = Int[]
    cells_mid = Int[]
    cells_lohi = UnitRange{Int}[]
    for p in 1:np, k in 1:(length(edges) - 1)

        lo, hi = edges[k] + 1 - n, edges[k + 1] - n
        sum(@view I_f[p, lo:hi]) >= min_infections || continue
        push!(cells_patch, p)
        push!(cells_mid, (edges[k] + edges[k + 1]) ÷ 2)
        push!(cells_lohi, lo:hi)
    end
    d, df = meld.d, length(cells_patch)
    S_f = [
        log(safe_rate(sum(@view future[i][cells_patch[c], cells_lohi[c]])))
            for i in 1:ndraw, c in 1:df
    ]
    Σ = cov(hcat(meld.log_sums, S_f))
    L21 = d == 0 ? zeros(df, 0) :
        Σ[(d + 1):end, 1:d] / transpose(LowerTriangular(meld.L))
    L22 = df == 0 ? zeros(0, 0) :
        _zone_meld_factor(
            Symmetric(Σ[(d + 1):end, (d + 1):end] - L21 * transpose(L21)), df;
            ridge
        )
    L = [meld.L zeros(d, df); L21 L22]
    ## The deformation's daily rows: the fitted rows to `n`, then each patch
    ## interpolated from its value on day `n` to its future midpoints.
    weights = zeros(Float64, np * (n + H), d + df)
    weights[1:(np * n), 1:d] .= meld.weights
    for p in 1:np
        cs = findall(==(p), cells_patch)
        anchors = vcat(n, cells_mid[cs])
        W = zone_interpolation_weights(anchors, n + 1, n + H)
        base = meld.weights[(n - 1) * np + p, :]
        for t in 1:H
            row = (n + t - 1) * np + p
            weights[row, 1:d] .= W[t, 1] .* base
            for (j, c) in enumerate(cs)
                weights[row, d + c] = W[t, j + 1]
            end
        end
    end
    fixed = zone_fixed_terms(I_ext, g, f, t0)
    death_fixed = zone_fixed_terms(I_ext, g, death_pmf, t0)
    nd = n + H - t0 + 1
    knots_ext = vcat(knots, future_knot_days(n, H; week))
    mixing_ext = mixing === nothing ? nothing :
        merge(
            mixing, (;
                import_fraction = hcat(
                    mixing.import_fraction,
                    repeat(mixing.import_fraction[:, n], 1, H)
                ),
            )
        )
    conf = _draw_vectors(forecast, :forecast_province_confirmed)
    j = findfirst(==(H), future_knot_days(0, Hp; week))
    j === nothing && error(
        "zone_forecast_block: the $H-day horizon is not one of the parent's " *
            "weekly forecast vintages."
    )
    totals = [
        sum(reshape(Int.(conf[i]), np, :)[p, 1:j]) for i in 1:ndraw, p in 1:np
    ]
    return (;
        horizon = H, knots = knots_ext,
        n_future_knots = length(knots_ext) - length(knots),
        meld_weights = weights, meld_L = L, meld_d_future = df,
        I_bar = I_ext, fixed.force_pre, fixed.report_pre_cum,
        fixed.infections_pre,
        report_pre_rows = zone_report_pre_rows(
            fixed.report_pre, patch_of_zone, t0, n + H
        ),
        death_pre_cum = death_fixed.report_pre_cum,
        death_pre_rows = zone_report_pre_rows(
            death_fixed.report_pre, patch_of_zone, t0, n + H
        ),
        interp = zone_interpolation_weights(knots_ext, t0, n + H),
        report_matrix = zone_delay_operator(f, nd),
        death_matrix = zone_delay_operator(death_pmf, nd),
        mixing = mixing_ext, totals,
    )
end

"""
$(TYPEDSIGNATURES)

The multiplier the sampled parent draw applies to the province model's
mean patch infections, `exp(a_p(t))` as an `(n_patches × n)` matrix, from
the whitened draw `η`. `weights` and `L` are the corresponding fields of
[`zone_meld_block`](@ref). One triangular matrix-vector product and one
interpolation product per forward pass.
"""
function zone_parent_scale(
        weights::AbstractMatrix, L::AbstractMatrix,
        η::AbstractVector, np::Integer, n::Integer
    )
    return reshape(exp.(weights * (L * η)), np, n)
end

"""
$(TYPEDSIGNATURES)

The patch quantities the zone renewal reads, deformed by the sampled
parent draw. `scale` is the multiplier of [`zone_parent_scale`](@ref), or
`nothing` for the cut, which reads the province model's posterior mean and
carries none of its uncertainty. Because the deformation is constant over
the days before the grid start, every pre-`t0` term of
[`zone_fixed_terms`](@ref) is the fixed one times one factor per patch.
"""
function zone_deformation(zd, ::Nothing)
    base = (;
        zd.I_bar, zd.force_pre, zd.report_pre_cum, zd.infections_pre,
        zd.report_pre_rows,
    )
    hasproperty(zd, :death_pre_cum) || return base
    return merge(base, (; zd.death_pre_cum, zd.death_pre_rows))
end

function zone_deformation(zd, scale::AbstractMatrix)
    pre = view(scale, :, 1)
    poz = hasproperty(zd, :patch_of_zone) ? zd.patch_of_zone :
        _zone_patch_of_zone(zd.patch_ranges)
    zpre = transpose(view(pre, poz))
    base = (;
        I_bar = zd.I_bar .* scale,
        force_pre = zd.force_pre .* pre,
        report_pre_cum = zd.report_pre_cum .* pre,
        infections_pre = zd.infections_pre .* pre,
        report_pre_rows = zd.report_pre_rows .* zpre,
    )
    hasproperty(zd, :death_pre_cum) || return base
    return merge(
        base,
        (;
            death_pre_cum = zd.death_pre_cum .* pre,
            death_pre_rows = zd.death_pre_rows .* zpre,
        )
    )
end

"""
$(TYPEDSIGNATURES)

Lower-triangular correlation factors of the zone deviations, one per patch,
for [`deviation_knots`](@ref). The correlation between two zones decays
with the distance between their centroids on one shared length scale,

```math
C_{zq} = \\exp(-d_{zq} / \\ell),
\\qquad \\ell = -\\bar d / \\log \\rho_{\\text{ref}},
```

so `ρ_ref` is the correlation of two zones a reference distance `d̄` apart.
A `ridge` on the diagonal conditions each factorisation.
`distances` holds one matrix per patch, over its zones for the level and
over its walking zones for the innovations. The meld carries correlation
between patches, through the parent draw every zone of a patch shares.
This carries it within a patch, where the zones multiply one trajectory and
only their deviations separate them, so a cluster of neighbours can share a
local excursion.
"""
function zone_correlation_factors(
        distances::AbstractVector, ℓ::Real;
        ridge::Real = 1.0e-6
    )
    return [_zone_correlation_factor(D, ℓ, ridge) for D in distances]
end

function _zone_correlation_factor(D::AbstractMatrix, ℓ::Real, ridge::Real)
    m = size(D, 1)
    C = exp.(.-D ./ ℓ) + ridge * Matrix{Float64}(I, m, m)
    return Matrix(cholesky(Symmetric(C)).L)
end

## Great-circle distances between the centroids of the zones in `rows`,
## grouped by `ranges`, or an empty list when any zone lacks metadata.
function _zone_distance_blocks(
        zones, zone_province, zone_names, ranges;
        keep = nothing
    )
    zones === nothing && return Matrix{Float64}[]
    lookup = Dict((r.province, r.zone) => r for r in zones)
    rows = [
        get(lookup, (zone_province[z], zone_names[z]), nothing)
            for z in eachindex(zone_names)
    ]
    any(isnothing, rows) && return Matrix{Float64}[]
    coords = [(r.lat, r.lon) for r in rows]
    sel(zs) = keep === nothing ? collect(zs) : [z for z in zs if keep[z]]
    return [province_distance_matrix(coords[sel(zs)]) for zs in ranges]
end

## Draws of `key` flattened, whether the chain stores a scalar or a vector.
function _zone_draw_values(chn, key::Symbol)
    raw = vec(collect(chn[key]))
    return eltype(raw) <: AbstractVector ? Float64.(reduce(vcat, raw)) :
        Float64.(raw)
end

## Log-normal matched to the parent's draws of a positive scale, so the zone
## prior carries the province model's centre and its uncertainty.
function _zone_parent_lognormal(chn, key::Symbol, fallback::NTuple{2, Float64})
    _has_key(chn, key) || return fallback
    l = log.(max.(_zone_draw_values(chn, key), floatmin(Float64)))
    isempty(l) && return fallback
    return (mean(l), max(length(l) > 1 ? std(l) : 0.5, 0.05))
end

## Mean over the parent's draws of the ratio of two vector-valued
## quantities, empty when the chain does not carry both. Taken per draw and
## then averaged, rather than as a ratio of two averages, so it is a share
## the parent actually held in some draw. Floored at zero and capped just
## below one, since the numerator is a component of the denominator and only
## rounding can put the ratio outside that.
function _mean_parent_ratio(
        chn, num::Symbol, den::Symbol; cap::Real = 0.95
    )
    (_has_key(chn, num) && _has_key(chn, den)) || return Float64[]
    ns = _draw_vectors(chn, num)
    ds = _draw_vectors(chn, den)
    (isempty(ns) || length(ns) != length(ds)) && return Float64[]
    length(first(ns)) == length(first(ds)) || return Float64[]
    m = zeros(Float64, length(first(ns)))
    for (a, b) in zip(ns, ds)
        m .+= clamp.(Float64.(a) ./ max.(Float64.(b), eps()), 0.0, cap)
    end
    return m ./ length(ns)
end

## Mean over the parent's draws of a vector-valued quantity, empty when the
## chain does not carry it.
function _mean_parent_vector(chn, key::Symbol)
    _has_key(chn, key) || return Float64[]
    vs = _draw_vectors(chn, key)
    isempty(vs) && return Float64[]
    m = zeros(Float64, length(first(vs)))
    for v in vs
        m .+= Float64.(v)
    end
    return m ./ length(vs)
end

## Beta matched by moments to the parent's draws of a probability. A
## correlation reaches this clamped to `[0, 1]`: the province model's
## deviation correlation has support on `[-1, 1]`, and a draw in which two
## provinces move apart says nothing about how far two neighbouring zones
## move together, so it enters as no correlation rather than a negative one.
function _zone_parent_beta(
        chn, key::Symbol, fallback::NTuple{2, Float64};
        clamp_unit::Bool = false
    )
    _has_key(chn, key) || return fallback
    v = _zone_draw_values(chn, key)
    clamp_unit && (v = clamp.(v, 0.0, 1.0))
    isempty(v) && return fallback
    m = clamp(mean(v), 1.0e-6, 1 - 1.0e-6)
    ν = max(m * (1 - m) / max(var(v), 1.0e-12) - 1, 1.0e-3)
    ## Both shapes are floored above one so the density is bounded at the
    ## ends of the support. Draws piled against zero or one give a moment
    ## match with a shape below one, whose density diverges there and
    ## invites the sampler into a spike at the boundary.
    return (max(m * ν, BETA_SHAPE_FLOOR), max((1 - m) * ν, BETA_SHAPE_FLOOR))
end

## Smallest Beta shape a parent-fitted prior may carry.
const BETA_SHAPE_FLOOR = 1.05

## Mean of a symmetric matrix's off-diagonal entries, or zero when there is
## no pair.
function _mean_offdiagonal(D::AbstractMatrix)
    m = size(D, 1)
    m > 1 || return 0.0
    return (sum(D) - sum(diag(D))) / (m * (m - 1))
end

## Patch index of every zone, from the patch ranges.
function _zone_patch_of_zone(patch_ranges::AbstractVector{<:UnitRange})
    nz = isempty(patch_ranges) ? 0 : last(last(patch_ranges))
    out = zeros(Int, nz)
    for (p, zs) in enumerate(patch_ranges), z in zs

        out[z] = p
    end
    return out
end

"""
$(TYPEDSIGNATURES)

The share renewal over the days `t0 … n`. Each zone's force of infection is
its own past infections through the generation interval `g` (lag 1), with
the days before `t0` entering as the initial share times the patch's
pre-`t0` force (`force_pre`, from [`zone_fixed_terms`](@ref)):

```math
Λ_z(t) = w_z(t_0)\\, Λ^{pre}_p(t)
    + \\sum_{s ≥ 1,\\; t − s ≥ t_0} g_s I_z(t − s),
\\qquad u_z(t) = e^{δ_z(t)} Λ_z(t).
```

`δ_daily` is `(n_days × n_zones)`, the interpolation weights
([`zone_interpolation_weights`](@ref)) times the knots, `w0` the initial
shares and `I_bar` the patch infections the draw's shared quantity implies.

Without mixing a zone takes the share of its patch its own force earns,
`w_z = u_z / \\sum_{z' ∈ p} u_{z'}` and `I_z = Ī_p w_z`.

With `mix`, the blocks of [`zone_importation_blocks`](@ref) and the
province model's own per-origin intensity and import fraction, each day
splits the patch total into what the patch grew and what it received:

```math
v_z = (1 − ε_z) u_z + \\sum_{q ∈ p,\\, q ≠ z} ε_q K^w_{zq} u_q,
\\qquad
M_p(t) = f_p(t)\\, Ī_p(t),
\\qquad
c_p(t) = \\frac{Ī_p(t) − M_p(t)}{\\sum_{z ∈ p} v_z},
```

```math
h_z(t) = \\sum_{q:\\, p(q) ≠ p(z)} \\barε_{p(q)} K^b_{zq} u_q,
\\qquad
\\mathrm{imp}_z(t) = M_p(t)\\, \\frac{h_z(t)}{\\sum_{z' ∈ p} h_{z'}(t)},
\\qquad
I_z(t) = c_p(t)\\, v_z(t) + \\mathrm{imp}_z(t).
```

The within-patch spill is a transfer, so it conserves `\\sum_{z ∈ p} v_z`,
and the imports are allocated to exactly the parent's own arrivals, so
`\\sum_{z ∈ p} I_z(t) = Ī_p(t)` exactly and no second parent term is added.
Before any other patch has infections the import pattern is empty and the
arrivals fall to the zones in proportion to their own force. `ε` is the
sampled per-origin within-patch intensity. Pass `mix = nothing` or
`ε = nothing` to leave the zones unmixed.

Returns `(; shares, forces, infections, imports)`, each `(n_days ×
n_zones)`.
"""
function zone_share_renewal(
        I_bar::AbstractMatrix, g::AbstractVector,
        δ_daily::AbstractMatrix, w0::AbstractVector,
        patch_ranges::AbstractVector{<:UnitRange}, t0::Integer,
        force_pre::AbstractMatrix;
        mix = nothing,
        ε::Union{Nothing, AbstractVector} = nothing
    )
    return zone_share_renewal_kernel(
        I_bar, g, δ_daily, w0, patch_ranges, t0, force_pre, mix, ε
    )
end

## `zone_share_renewal` with `mix` and `ε` positional, the form the reverse
## rule in `mooncake_rules.jl` is written against.
function zone_share_renewal_kernel(
        I_bar::AbstractMatrix, g::AbstractVector,
        δ_daily::AbstractMatrix, w0::AbstractVector,
        patch_ranges::AbstractVector{<:UnitRange}, t0::Integer,
        force_pre::AbstractMatrix, mix, ε
    )
    st = zone_share_renewal_with_state(
        I_bar, g, δ_daily, w0, patch_ranges, t0, force_pre, mix, ε
    )
    return (; st.shares, st.forces, st.infections, st.imports)
end

## `zone_share_renewal` keeping the per-day intermediates its adjoint reads:
## each zone's own force `u`, the mixed force `v`, the import pattern `h`,
## and per patch and day the unclamped denominators `patch_total` (the sum
## of `u` unmixed, of `v` mixed) and `share_total` (the sum of `h`). The
## rule calls this so the recursion is defined once.
function zone_share_renewal_with_state(
        I_bar::AbstractMatrix, g::AbstractVector,
        δ_daily::AbstractMatrix, w0::AbstractVector,
        patch_ranges::AbstractVector{<:UnitRange}, t0::Integer,
        force_pre::AbstractMatrix, mix, ε
    )
    nd, nz = size(δ_daily)
    np = length(patch_ranges)
    on = mix !== nothing && ε !== nothing
    ## The loop reads the per-day terms by absolute day under `@inbounds`,
    ## so one that stops short of the last day is read past its end rather
    ## than raising. Checked once per call, not per day.
    last_day = t0 + nd - 1
    size(I_bar, 2) >= last_day && size(force_pre, 2) >= last_day || throw(
        DimensionMismatch(
            "zone_share_renewal: the patch trajectories and forces must " *
                "reach day $(last_day); got $(size(I_bar, 2)) and " *
                "$(size(force_pre, 2))."
        )
    )
    if on && size(mix.import_fraction, 2) < last_day
        throw(
            DimensionMismatch(
                "zone_share_renewal: the import fractions must reach day " *
                    "$(last_day); got $(size(mix.import_fraction, 2)). " *
                    "A forecast extension has to carry the mixing forward."
            )
        )
    end
    Tp = promote_type(
        eltype(I_bar), eltype(g), eltype(δ_daily),
        eltype(w0), eltype(force_pre),
        ε === nothing ? Float64 : eltype(ε),
        on ? eltype(mix.within) : Float64
    )
    I = zeros(Tp, nd, nz)
    Λ = zeros(Tp, nd, nz)
    W = zeros(Tp, nd, nz)
    imports = zeros(Tp, nd, nz)
    U = zeros(Tp, nd, nz)
    V = zeros(Tp, nd, nz)
    H = zeros(Tp, nd, nz)
    patch_total = zeros(Tp, np, nd)
    share_total = zeros(Tp, np, nd)
    u = zeros(Tp, nz)
    eu = zeros(Tp, nz)
    wu = zeros(Tp, nz)
    L = length(g)
    floor_ = eps(Tp)
    @inbounds for j in 1:nd
        t = t0 + j - 1
        smax = min(L, j - 1)
        ## Every zone's own force first, since the between-patch pattern
        ## reads the zones of the other patches on the same day.
        for p in eachindex(patch_ranges)
            zs = patch_ranges[p]
            fp = force_pre[p, t]
            for z in zs
                acc = w0[z] * fp
                for s in 1:smax
                    acc += g[s] * I[j - s, z]
                end
                Λ[j, z] = acc
                u[z] = exp(δ_daily[j, z]) * acc
                U[j, z] = u[z]
            end
        end
        if !on
            for p in eachindex(patch_ranges)
                zs = patch_ranges[p]
                isempty(zs) && continue
                tot = zero(Tp)
                for z in zs
                    tot += u[z]
                end
                patch_total[p, j] = tot
                tot = max(tot, floor_)
                Ip = I_bar[p, t]
                for z in zs
                    W[j, z] = u[z] / tot
                    I[j, z] = Ip * W[j, z]
                end
            end
            continue
        end
        ## Within-patch spill, per origin as the province model exports:
        ## zone `q` sends `ε_q` of its force through the within block and
        ## keeps the rest, so the patch total is conserved exactly. The
        ## pattern the parent's arrivals land in comes from the zones of the
        ## other patches, weighted by the province model's own per-origin
        ## intensity.
        ##
        ## Both are one matrix-vector product rather than a nested loop. The
        ## within block is zero on the diagonal and across patches and the
        ## between block is zero within a patch, so the whole matrix times
        ## the weighted force is exactly the restricted sum, and the reverse
        ## pass sees one product per day instead of `n_zones^2` scalar
        ## multiplies.
        for z in 1:nz
            eu[z] = ε[z] * u[z]
            wu[z] = mix.origin_weight[z] * u[z]
        end
        spill = mix.within * eu
        h = mix.between * wu
        for z in 1:nz
            V[j, z] = (one(Tp) - ε[z]) * u[z] + spill[z]
            H[j, z] = h[z]
        end
        for p in eachindex(patch_ranges)
            zs = patch_ranges[p]
            isempty(zs) && continue
            Ip = I_bar[p, t]
            Mp = mix.import_fraction[p, t] * Ip
            sv = zero(Tp)
            sh = zero(Tp)
            for z in zs
                sv += V[j, z]
                sh += H[j, z]
            end
            patch_total[p, j] = sv
            share_total[p, j] = sh
            sv = max(sv, floor_)
            cp = (Ip - Mp) / sv
            for z in zs
                share_h = sh > floor_ ? H[j, z] / sh : V[j, z] / sv
                imports[j, z] = Mp * share_h
                I[j, z] = cp * V[j, z] + imports[j, z]
                W[j, z] = I[j, z] / max(Ip, floor_)
            end
        end
    end
    return (;
        shares = W, forces = Λ, infections = I, imports,
        u = U, v = V, h = H, patch_total, share_total,
    )
end

"""
$(TYPEDSIGNATURES)

Bin the daily expected reports `(n_days × n_zones)` of
[`zone_forward`](@ref) into the vintage windows `(d_{v−1}, d_v]`
given by the grid days `days`, the first window opening on day one. Days
before the grid start `t0` contribute the initial share times the patch's
pre-`t0` accrual (`report_pre_cum`, from [`zone_fixed_terms`](@ref)).
Returns the `(n_zones × n_vintages)` expected increments.
"""
function zone_report_increments(
        reports::AbstractMatrix, w0::AbstractVector,
        patch_ranges::AbstractVector{<:UnitRange},
        days::AbstractVector{<:Integer}, t0::Integer,
        report_pre_cum::AbstractMatrix
    )
    nd, nz = size(reports)
    nv = length(days)
    Tp = promote_type(eltype(reports), eltype(w0), eltype(report_pre_cum))
    C = zeros(Tp, nz, nv)
    @inbounds for v in 1:nv
        lo = v == 1 ? 1 : Int(days[v - 1]) + 1
        hi = Int(days[v])
        for p in eachindex(patch_ranges)
            zs = patch_ranges[p]
            ## Accrued before the grid start, by difference of the cumulative.
            pre = zero(Tp)
            if lo < t0
                top = min(hi, t0 - 1)
                pre = report_pre_cum[p, top] -
                    (lo > 1 ? report_pre_cum[p, lo - 1] : zero(Tp))
            end
            jlo = max(lo, t0) - t0 + 1
            jhi = hi - t0 + 1
            for z in zs
                acc = w0[z] * pre
                for j in jlo:min(jhi, nd)
                    acc += reports[j, z]
                end
                C[z, v] = acc
            end
        end
    end
    return C
end

"""
$(TYPEDSIGNATURES)

Log probability mass of counts `y` under a Dirichlet-multinomial with
concentration vector `α`, `Σ y` trials:

```math
\\log \\frac{N!\\,Γ(A)}{Γ(N + A)}
+ \\sum_i \\log \\frac{Γ(y_i + α_i)}{Γ(α_i)\\, y_i!},
\\qquad A = \\sum_i α_i.
```

Plain arithmetic and `loggamma`, so it differentiates under Mooncake.
"""
function dirichlet_multinomial_logpdf(
        y::AbstractVector{<:Integer},
        α::AbstractVector
    )
    Tp = float(eltype(α))
    N = 0
    A = zero(Tp)
    @inbounds for i in eachindex(y, α)
        N += y[i]
        A += α[i]
    end
    lp = loggamma(Tp(N + 1)) + loggamma(A) - loggamma(N + A)
    @inbounds for i in eachindex(y, α)
        lp += loggamma(y[i] + α[i]) - loggamma(α[i]) - loggamma(Tp(y[i] + 1))
    end
    return lp
end

"""
$(TYPEDSIGNATURES)

The composition log likelihood summed over the scored cells: for cell `c`,
patch `cell_patch[c]` at vintage `cell_vintage[c]` with allocated total
`cell_total[c]`, the observed zone counts in `counts` `(n_zones ×
n_vintages)` follow a Dirichlet-multinomial with concentration `κ π`,
where `π` is the modelled expected reports `C` of that patch's zones
normalised within the patch. `cell_const[c]` is the parameter-free part
of the mass, `log N! − Σ log y!`.
"""
function zone_composition_logpdf(
        counts::AbstractMatrix{<:Integer},
        C::AbstractMatrix, cell_patch::AbstractVector{<:Integer},
        cell_vintage::AbstractVector{<:Integer},
        cell_total::AbstractVector{<:Integer},
        cell_const::AbstractVector{<:Real},
        patch_ranges::AbstractVector{<:UnitRange}, κ::Real
    )
    Tp = promote_type(eltype(C), typeof(float(κ)))
    lp = zero(Tp)
    lgκ = loggamma(κ)
    @inbounds for c in eachindex(cell_patch)
        v = cell_vintage[c]
        zs = patch_ranges[cell_patch[c]]
        tot = zero(Tp)
        for z in zs
            tot += safe_rate(C[z, v])
        end
        acc = lgκ - loggamma(cell_total[c] + κ) + cell_const[c]
        for z in zs
            y = counts[z, v]
            y == 0 && continue
            α = κ * safe_rate(C[z, v]) / tot
            acc += loggamma(y + α) - loggamma(α)
        end
        lp += acc
    end
    return lp
end

## `log N! − Σ_z log y_z!` for the counts in one column of `counts` over the
## zones `zs`: the parameter-free part of the Dirichlet-multinomial mass.
function _zone_cell_const(counts::AbstractMatrix{<:Integer}, zs, v::Integer)
    N = 0
    acc = 0.0
    @inbounds for z in zs
        y = counts[z, v]
        N += y
        acc -= loggamma(y + 1.0)
    end
    return acc + loggamma(N + 1.0)
end

"""
$(TYPEDSIGNATURES)

One forward pass of the zone model from its fixed data `zd` (the
`model_data` of [`zone_fit_inputs`](@ref)) and a draw's deviation knots
`(n_zones × n_knots)`, initial shares `w0` and per-patch mixing fractions
`ε` (`nothing` when mixing is off): the daily deviations (the
interpolation weights times the knots), the share renewal and the binned
expected reports (the delay operator times the infections plus the
pre-`t0` rows). `def` is the patch trajectory the draw's shared quantity
implies ([`zone_deformation`](@ref)); its default is the cut, the
province model's posterior mean. Called by the model and, per draw, by the
render. Returns
`(; shares, forces, infections, imports, reports, increments)`.
"""
function zone_forward(
        zd, δ_knots::AbstractMatrix, w0::AbstractVector,
        ε::Union{Nothing, AbstractVector},
        def = zone_deformation(zd, nothing)
    )
    δ_daily = zd.interp * transpose(δ_knots)
    mix = ε === nothing ? nothing : zd.mixing
    st = zone_share_renewal(
        def.I_bar, zd.g, δ_daily, w0, zd.patch_ranges,
        zd.t0, def.force_pre; mix, ε
    )
    reports = zd.report_matrix * st.infections .+
        def.report_pre_rows .* transpose(w0)
    increments = zone_report_increments(
        reports, w0, zd.patch_ranges,
        zd.days, zd.t0, def.report_pre_cum
    )
    return (;
        st.shares, st.forces, st.infections, st.imports, reports,
        increments,
    )
end

## Cumulative zone infections to the last grid day: the pre-`t0` patch
## infections at the initial share plus the renewal's own days.
function _zone_cumulative_infections(
        infections::AbstractMatrix,
        w0::AbstractVector, patch_ranges::AbstractVector{<:UnitRange},
        infections_pre::AbstractVector
    )
    nd, nz = size(infections)
    Tp = promote_type(eltype(infections), eltype(w0), eltype(infections_pre))
    cum = zeros(Tp, nz)
    @inbounds for p in eachindex(patch_ranges)
        zs = patch_ranges[p]
        for z in zs
            acc = w0[z] * infections_pre[p]
            for j in 1:nd
                acc += infections[j, z]
            end
            cum[z] = acc
        end
    end
    return cum
end

## Implied zone reproduction number on grid row `j`, `I_z / Λ_z`, `NaN` where
## the zone's cumulative infections are below `rt_floor`.
function _zone_rt_at(
        infections::AbstractMatrix, forces::AbstractMatrix,
        cum::AbstractVector, j::Integer, rt_floor::Real
    )
    Tp = promote_type(eltype(infections), eltype(forces))
    nz = size(infections, 2)
    out = zeros(Tp, nz)
    nan = convert(Tp, NaN)
    @inbounds for z in 1:nz
        r = infections[j, z] / max(forces[j, z], eps(Tp))
        out[z] = cum[z] >= rt_floor ? r : nan
    end
    return out
end

## Dirichlet-multinomial concentration `(1 − ρ) / ρ` from the composition
## dispersion `ρ`, with `ρ` and `1 − ρ` floored at machine epsilon so a
## proposal at either end of `[0, 1]` gives a finite, positive `κ`.
_zone_kappa(ρ::Real) = safe_rate(one(ρ) - ρ) / safe_rate(ρ)

## Shares at the knot days as an `(n_zones × n_knots)` matrix.
function _zone_shares_at_knots(
        shares::AbstractMatrix,
        knots::AbstractVector{<:Integer}, t0::Integer
    )
    nz = size(shares, 2)
    out = zeros(eltype(shares), nz, length(knots))
    @inbounds for (k, d) in enumerate(knots), z in 1:nz

        out[z, k] = shares[d - t0 + 1, z]
    end
    return out
end

## --- The model -----------------------------------------------------------

## The zone forecast over the horizon. Each zone's expected confirmed
## reports over `(n, n + H]`, times the case composition's multiplier, set
## its share of the patch, and a paired parent draw's predicted patch total
## is split over the zones by the fitted Dirichlet-multinomial, drawn one
## zone at a time as its sequence of Beta-binomials.
@model function _zone_forecast_counts(zd, zf, fw, asc, ρ, nd::Integer)
    nz = length(asc)
    H = zf.horizon
    future = (nd + 1):(nd + H)
    expected = vec(sum(view(fw.reports, future, :); dims = 1))
    Tp = promote_type(eltype(expected), eltype(asc), typeof(float(ρ)))
    π = zeros(Tp, nz)
    for zs in zd.patch_ranges
        isempty(zs) && continue
        tot = sum(asc[z] * safe_rate(expected[z]) for z in zs)
        for z in zs
            π[z] = asc[z] * safe_rate(expected[z]) / tot
        end
    end
    κ = _zone_kappa(ρ)
    forecast_pair ~ DiscreteUniform(1, size(zf.totals, 1))
    N = zf.totals[forecast_pair, :]
    forecast_zone_draw = Vector{Union{Missing, Int}}(missing, nz)
    counts = zeros(Int, nz)
    for (p, zs) in enumerate(zd.patch_ranges)
        isempty(zs) && continue
        left = N[p]
        tail = one(Tp)
        for z in zs
            if z == last(zs)
                counts[z] = left
            else
                forecast_zone_draw[z] ~ BetaBinomial(
                    left, κ * safe_rate(π[z]), κ * safe_rate(tail - π[z])
                )
                counts[z] = forecast_zone_draw[z]
                left -= counts[z]
                tail -= π[z]
            end
        end
    end
    forecast_zone_confirmed := counts
    forecast_patch_confirmed := N
    forecast_zone_share := π
    forecast_zone_concentration := κ
    forecast_zone_infections := vec(sum(view(fw.infections, future, :); dims = 1))
    return (; counts, π)
end


"""
Health-zone composition model, stage two of the melding. `zd` is the
`model_data` of [`zone_fit_inputs`](@ref): the fixed stage-1 inputs (patch
infections `Ī_p`, the generation-interval and infection-to-report PMFs),
the zone count matrix and the scored cells, the knot days and the walking
mask, all plain arrays, plus the softmax scale `share_scale`, the knot
spacing `week` and the reporting floor `rt_floor`. Nothing in it is
sampled here.

### Parameters

```math
\\begin{aligned}
z^w_z &\\sim N(0, 1), &
w_z(t_0) &= \\mathrm{softmax}_p\\bigl(s\\,(z^w_z − \\bar z^w_p)\\bigr) \\\\
σ_L &\\sim N^+(0, 0.3),\\; z^L_z \\sim N(0, 1), &
h &\\sim \\mathrm{LogNormal}(\\log 42, 0.6),\\; φ = 2^{−w/h} \\\\
σ_{δ,p} &\\sim \\mathrm{LogNormal}(\\text{parent}), &
ρ_{\\text{corr}} &\\sim \\mathrm{Beta}(\\text{parent}),\\;
    ℓ = −\\bar d / \\log ρ_{\\text{corr}} \\\\
ρ,\\; ρ_D &\\sim \\mathrm{Beta}(\\text{parent}), &
    κ &= (1 − ρ)/ρ,\\; κ_D = (1 − ρ_D)/ρ_D \\\\
σ_a &\\sim \\mathrm{LogNormal}(\\text{parent}),\\; σ_θ \\sim N^+(0, 0.1), &
    z^a_p,\\; z^θ_p &\\sim N(0, I_{n_p − 1}) \\\\
ε_w &\\sim \\mathrm{Beta}(1, 20),\\; τ \\sim N^+(0, 0.5), &
    z^m_z &\\sim N(0, 1) \\\\
η &\\sim N(0, I_d)
\\end{aligned}
```

`s` is `zd.share_scale` and `w` is `zd.week`, the days between knots. The
mixing block (`ε_w`, `τ`, `z^m`) is sampled only where `zd.mixing` carries
the kernel, the correlation `ρ_corr` only where `zd.zone_distances` does,
and the shared draw `η` only where `zd.meld_d` is positive.
The deviation knots `δ_z(k)` are [`deviation_knots`](@ref), the province
model's own process, called with one group per patch and the zones of a
patch as its units: the level and the innovations are correlated within a
patch by [`zone_correlation_factors`](@ref) and centred within it, so every
patch sums to zero at every knot.

The innovations exist for the walking zones `W_p` only (cumulative
confirmed cases at the cut-off at or above the threshold, and at least two
such zones in the patch). A level-only zone decays along the AR mean path.
Each block is one `~` over a `product_distribution`, and the innovation
block is absent when no zone walks. `κ` floors both `ρ` and `1 − ρ` at
machine epsilon (`safe_rate`), so a proposal at either end of the prior's
support keeps the composition mass finite.

Every scale the two levels share takes its prior from the province
posterior fitted to its draws (`zd.parent_priors`): the drift scale, the
two composition overdispersions, the relative ascertainment and fatality
scales, and the deviation correlation. A zone therefore starts at its
province's estimate and departs only as far as its own counts require.

### The shared quantity

`η ~ N(0, I_d)` is the whitened draw of the province model's weekly
infections in each patch, the one quantity the two stages share
([`zone_meld_block`](@ref)). Its Cholesky factor carries the province
model's joint posterior over all patches and weeks together, so a draw
moves whole patch trajectories, and moves the patches together where the
province model says they move together. The daily patch trajectory the
draw implies is the province model's mean curve times `exp(a_p(t))` with
`a = L η` ([`zone_parent_scale`](@ref)), and the zone model is conditional
on it.

This is the model's only parent term. The zone infections sum to the
sampled patch total by construction. Scoring those sums against the
province model's posterior again would count that posterior twice. With no
kept cell (`zd.meld_d == 0`) the stage reads the province model's mean
curve.

### Likelihood

The share renewal [`zone_share_renewal`](@ref) gives each zone's
infections as its share of the sampled patch infections, the delay
operator ([`zone_delay_operator`](@ref)) and binning
([`zone_report_increments`](@ref)) give
the expected confirmed reports per vintage window, and the observed zone
increments of each patch and vintage follow a Dirichlet-multinomial on the
allocated total with concentration `κ π`
([`zone_composition_logpdf`](@ref)), as one `@addlogprob!` over every
cell.

The allocated zone deaths of every vintage follow a second
Dirichlet-multinomial, on the infection-to-confirmed-death PMF rather than
the case delay, with its own concentration from `ρ_death`. A case table and
a death table do not disperse alike. The two do not share one.

Cases observe incidence times case-finding and deaths incidence times
lethality, so each composition carries its own multiplier over the zones of
a patch, [`relative_multiplier`](@ref): relative ascertainment on the cases
and relative fatality on the deaths. Both are the province model's own
per-province construction, log contrasts summing to zero within the group,
partially pooled at a scale whose prior is the province posterior for the
same scale between provinces. A composition identifies only the product of
a multiplier and the incidence split, so the pooling is what separates
them: as the scale shrinks the shares weight zones by incidence alone.

Both are relative to the zone's own province. A factor common to a patch
cancels in a within-patch composition, so the province level of each
multiplier is the one the province model estimates from the per-province
compositions, and the zone level is estimated here. The two are multiplied
into `zone_ascertainment_national` and `zone_severity_national` for
reporting on one scale.

Two compositions carry three unknowns per zone, so one has to be pinned.
`severity_sd_prior` is tight where the ascertainment scale is the province
posterior's, which asserts that deaths per infection vary little between
the zones of a patch. The death composition then pins the incidence split
and the case composition identifies ascertainment as the residual, as the
province model's own asymmetry does between provinces. Read
`zone_severity_sd` against its prior: a posterior that has not moved says
the assumption is carrying the identification.

Infections cross zone boundaries through the two
blocks of [`zone_importation_blocks`](@ref) as
[`zone_share_renewal`](@ref) applies them. Between patches the flows are
the province model's own, fixed: its per-origin intensity weights the
origin zones and its arrivals into a patch set how many infections land
there, so the zone stage only distributes them over the destination
patch's zones and the same movement is not counted at both levels. Within
a patch the spill is the zone stage's own mechanism and carries its own
intensity, `ε_w ~ Beta(1, 20)` with a pooled per-origin deviation
`τ (z − z̄)` on the logit scale. Mixing is off, and the correlation below with it, when
the health-zone metadata does not cover every zone or the parent chain
carries no between-patch movement (`zd.mixing === nothing`).

### Deterministics

Flattened column-major where a matrix: `delta_knots_zone` `(n_zones ×
n_knots)`, `share_knots_zone` `(n_zones × n_knots)`, `delta_T_zone`,
`share_T_zone`, `share_start_zone` (the initial shares `w_z(t_0)`),
`R_T_zone` (the implied `I_z / Λ_z` at the cut-off, `NaN` below
`zd.rt_floor` cumulative zone infections), `region_sd_zone` (`σ_L`),
`region_drift_sd_zone` (`σ_δ`, one per patch), `region_halflife_zone`
(`h`), `composition_rho_zone` (`ρ`), `composition_rho_death_zone`,
`zone_ascertainment_national` and `zone_severity_national` (the same
multipliers against the national average, the province model's contrast
times the zone's own),
`correlation_reference_zone` (`ρ_corr`) and `correlation_length_zone`
(`ℓ`), `parent_eta_zone` (the whitened shared draw `η`),
`parent_patch_T_zone` (the sampled patch infections on the last grid day)
and, with mixing, `mixing_epsilon_zone`. Daily trajectories are rebuilt
from these by [`zone_forward`](@ref).

### Forecast

With `forecast` a [`ForecastHorizon`](@ref) (see [`with_horizon`](@ref))
the model runs past the cut-off on the inputs of
[`zone_forecast_block`](@ref), which `zd.forecast` must carry for the same
horizon. The future weeks of the shared quantity are `η_future ∼ N(0, I)`
through the extended factor, the future knots take fresh innovations
`z_drift_future` through the same AR(1), and the renewal and delays run on
over the horizon. Each zone's expected confirmed reports over the horizon,
times its relative ascertainment, set its share `forecast_zone_share` of
the patch. `forecast_pair` picks a parent forecast draw, whose confirmed
cases per patch (`forecast_patch_confirmed`) are split over the zones by
the fitted case composition into `forecast_zone_confirmed`. Every fitted
quantity is computed on the fitted days as without a forecast.
"""
@model function bvd_zone(
        zd;
        region_sd_prior = truncated(Normal(0, 0.3); lower = 0),
        region_halflife_prior = LogNormal(log(42), 0.6),
        severity_sd_prior = truncated(Normal(0, 0.1); lower = 0),
        mixing_within_prior = Beta(1, 20),
        mixing_departure_prior = truncated(Normal(0, 0.5); lower = 0),
        offset_prior = Normal(0, 1),
        forecast::Union{Nothing, ForecastHorizon} = nothing
    )
    nz = size(zd.counts, 1)
    np = length(zd.patch_ranges)
    K = length(zd.knots)
    H = horizon_days(forecast)
    zf = get(zd, :forecast, nothing)
    H == 0 || (zf !== nothing && zf.horizon == H) || error(
        "bvd_zone: a $H-day forecast needs zone inputs built with the " *
            "parent's forecast over the same horizon."
    )
    ## Mixing needs the kernel blocks and the province model's own flows;
    ## without them the zones stay inside their own boundaries.
    mix_on = zd.mixing !== nothing
    z_w ~ product_distribution(fill(offset_prior, nz))
    σ_level ~ region_sd_prior
    z_level ~ product_distribution(fill(offset_prior, nz))
    δ_halflife ~ region_halflife_prior
    ## Every scale below takes its prior from the province posterior.
    pp = zd.parent_priors
    ## One drift scale per patch, as the province model gives each patch its
    ## own, on the scale the province model estimated between provinces.
    σ_δ ~ product_distribution(
        fill(LogNormal(pp.drift_sd[1], pp.drift_sd[2]), np)
    )
    ## The death composition carries its own concentration.
    ρ ~ Beta(pp.rho[1], pp.rho[2])
    ρ_death ~ Beta(pp.rho_death[1], pp.rho_death[2])
    ## Relative ascertainment on the cases and relative fatality on the
    ## deaths, each partially pooled within its patch on the scale the
    ## province model estimated between provinces.
    σ_ascertainment ~ LogNormal(pp.ascertainment_sd[1], pp.ascertainment_sd[2])
    n_contrast = relative_multiplier_dims(zd.patch_ranges)
    z_ascertainment ~ product_distribution(fill(offset_prior, n_contrast))
    ## Deaths per infection are taken as near-uniform across the zones of a
    ## patch, which is what lets the death composition pin the incidence
    ## split and the case composition identify ascertainment as the
    ## residual. The prior is tight and fixed rather than inherited, the
    ## asymmetry the province model makes between death ascertainment and
    ## provincial lethality.
    σ_severity ~ severity_sd_prior
    z_severity ~ product_distribution(fill(offset_prior, n_contrast))
    ## Sampled only when used, or they would be prior-only dimensions.
    n_drift = zd.n_walking * (K - 1)
    if n_drift > 0
        z_drift ~ product_distribution(fill(offset_prior, n_drift))
    else
        z_drift = Float64[]
    end
    if mix_on
        ## Within-patch spill only. The province model's own intensity is a
        ## between-province one that its data pulled four orders of
        ## magnitude below one, and a boundary between two neighbouring
        ## zones is not that quantity, so this carries its own prior. The
        ## between-patch flows take the province model's intensity instead,
        ## fixed, inside the kernel's between block. `τ_mix` is how far one
        ## origin may depart from the shared level, on the logit scale so a
        ## departure never leaves the unit interval.
        ε_within ~ mixing_within_prior
        τ_mix ~ mixing_departure_prior
        z_mix ~ product_distribution(fill(offset_prior, nz))
        ε_mix = logistic.(
            logit(ε_within) .+ τ_mix .* (z_mix .- sum(z_mix) / nz)
        )
    else
        ε_mix = nothing
    end
    ## The shared quantity, whitened, as the docstring's "shared quantity"
    ## section describes it.
    n_meld = zd.meld_d
    if n_meld > 0
        η ~ product_distribution(fill(offset_prior, n_meld))
    else
        η = Float64[]
    end
    ## Past the cut-off the grid, the delays and the shared quantity run on
    ## over the horizon; the fitted days are unchanged.
    zx = H == 0 ? zd :
        merge(
            zd, (;
                zf.I_bar, zf.force_pre, zf.report_pre_cum, zf.infections_pre,
                zf.report_pre_rows, zf.death_pre_cum, zf.death_pre_rows,
                zf.interp, zf.report_matrix, zf.death_matrix, zf.mixing,
            )
        )
    if H > 0 && zf.meld_d_future > 0
        η_future ~ product_distribution(fill(offset_prior, zf.meld_d_future))
        scale = zone_parent_scale(
            zf.meld_weights, zf.meld_L, vcat(η, η_future), np, zd.n + H
        )
    elseif H > 0
        scale = n_meld > 0 ?
            zone_parent_scale(zf.meld_weights, zf.meld_L, η, np, zd.n + H) :
            nothing
    else
        scale = n_meld > 0 ?
            zone_parent_scale(zd.meld_weights, zd.meld_L, η, np, zd.n) :
            nothing
    end
    def = zone_deformation(zx, scale)
    ## The correlation of two zones a reference distance apart, the prior
    ## taken from the province model's own learned correlation between its
    ## patches at the distance between their capitals.
    on = !isempty(zd.zone_distances)
    level_factors = Matrix{Float64}[]
    drift_factors = Matrix{Float64}[]
    if on
        ρ_corr ~ Beta(pp.correlation[1], pp.correlation[2])
        ℓ_corr = -zd.correlation_distance / log(safe_rate(ρ_corr))
        correlation_reference_zone := ρ_corr
        correlation_length_zone := ℓ_corr
        level_factors = zone_correlation_factors(zd.zone_distances, ℓ_corr)
        drift_factors = zone_correlation_factors(
            zd.zone_walk_distances, ℓ_corr
        )
    end
    w0 = zone_initial_shares(z_w, zd.patch_ranges, zd.share_scale)
    φ = exp2(-zd.week / δ_halflife)
    ## The province model's deviation process ([`deviation_knots`](@ref)),
    ## called with one group per patch and the zones of a patch as its units.
    Kf = H == 0 ? 0 : zf.n_future_knots
    if Kf > 0 && zd.n_walking > 0
        z_drift_future ~ product_distribution(
            fill(offset_prior, zd.n_walking * Kf)
        )
        z_drift_all = vcat(z_drift, z_drift_future)
    else
        z_drift_all = z_drift
    end
    δ_knots_all = deviation_knots(
        z_level, z_drift_all, σ_level, σ_δ[zd.patch_of_zone], φ,
        zd.patch_ranges, level_factors, drift_factors,
        zd.walking, zd.walk_index, zd.n_walking, K + Kf
    )
    δ_knots = H == 0 ? δ_knots_all : δ_knots_all[:, 1:K]
    fw = zone_forward(zx, δ_knots_all, w0, ε_mix, def)
    asc = relative_multiplier(
        z_ascertainment, σ_ascertainment, zd.patch_ranges
    )
    zone_ascertainment_sd := σ_ascertainment
    zone_ascertainment_relative := asc
    ## The same multiplier against the national average rather than the
    ## zone's own province: the province model's contrast times this one.
    ## It cancels in the composition above and is reported, not fitted.
    zone_ascertainment_national := zd.province_ascertainment .* asc
    @addlogprob! zone_composition_logpdf(
        zd.counts, fw.increments .* asc,
        zd.cell_patch, zd.cell_vintage, zd.cell_total, zd.cell_const,
        zd.patch_ranges, _zone_kappa(ρ)
    )
    ## The allocated deaths of every vintage, through the
    ## infection-to-confirmed-death delay rather than the case delay.
    death_daily = zx.death_matrix * fw.infections .+
        def.death_pre_rows .* transpose(w0)
    D = zone_report_increments(
        death_daily, w0, zd.patch_ranges,
        zd.death_days, zd.t0, def.death_pre_cum
    )
    sev = relative_multiplier(
        z_severity, σ_severity, zd.patch_ranges
    )
    zone_severity_sd := σ_severity
    zone_severity_relative := sev
    zone_severity_national := zd.province_severity .* sev
    @addlogprob! zone_composition_logpdf(
        zd.death_counts, D .* sev,
        zd.death_cell_patch, zd.death_cell_vintage,
        zd.death_cell_total, zd.death_cell_const, zd.patch_ranges,
        _zone_kappa(ρ_death)
    )
    nd = zd.n - zd.t0 + 1
    cum = _zone_cumulative_infections(
        view(fw.infections, 1:nd, :), w0, zd.patch_ranges,
        def.infections_pre
    )
    R_T_zone := _zone_rt_at(fw.infections, fw.forces, cum, nd, zd.rt_floor)
    parent_eta_zone := η
    parent_patch_T_zone := def.I_bar[:, zd.n]
    if H > 0
        forecast_zone ~ to_submodel(
            _zone_forecast_counts(zd, zf, fw, asc, ρ, nd), false
        )
    end
    delta_knots_zone := vec(δ_knots)
    delta_T_zone := δ_knots[:, K]
    share_T_zone := fw.shares[nd, :]
    share_start_zone := w0
    share_knots_zone := vec(_zone_shares_at_knots(fw.shares, zd.knots, zd.t0))
    region_sd_zone := σ_level
    region_drift_sd_zone := σ_δ
    region_halflife_zone := δ_halflife
    composition_rho_zone := ρ
    composition_rho_death_zone := ρ_death
    if mix_on
        mixing_epsilon_zone := ε_mix
        mixing_within_zone := ε_within
        mixing_departure_zone := τ_mix
        import_share_T_zone := [
            fw.imports[nd, z] / max(fw.infections[nd, z], eps(Float64))
                for z in 1:nz
        ]
    end
    return (;
        shares = fw.shares, forces = fw.forces,
        infections = fw.infections, increments = fw.increments,
        δ_knots, w0, cum,
    )
end

## --- Fixed inputs from the parent chain ---------------------------------

## Draws of `key` from the parent chain (per-draw vectors with
## `vectors = true`), with a clear error naming the key when the chain does
## not carry it.
function _zone_parent_draws(chn, key::Symbol; vectors::Bool = false)
    _has_key(chn, key) || error(
        "zone_fit_inputs: the parent chain carries no `$(key)`; it must be " *
            "a `bvd_joint` chain sampled with the patch structure on."
    )
    return vectors ? _draw_vectors(chn, key) : _draws(chn, key)
end

## Generation-interval PMF (lag 1) of one parent draw, discretised as the
## parent did.
function _zone_gi_pmf(α::Real, θ::Real)
    pmf = discretise_censored(Gamma(α, θ), ZONE_GI_NMAX)
    return pmf[2:end] ./ sum(pmf[2:end])
end

## Mean over the parent's draws of a per-draw PMF.
function _zone_mean_pmf(build, ndraws::Integer)
    acc = build(1)
    for i in 2:ndraws
        acc = acc .+ build(i)
    end
    return acc ./ ndraws
end

"""
$(TYPEDSIGNATURES)

The stage-1 quantities the zone model conditions on, read from a
`bvd_joint` parent chain: the patch infections `Ī_p(t)` `(n_patches × n)`,
the generation-interval PMF `g` (lag 1), the infection-to-confirmed-report
PMF `f = incubation ⊛ receipt` (lag 0, the delay the parent's per-province
composition applies to infections) and the infection-to-confirmed-death PMF
`death_pmf = incubation ⊛ onset-to-death ⊛ receipt`.

Each is a posterior mean: `log_infections` is the mean over draws of the
log infections per day, flattened as the parent stores `infections_patch`,
so `Ī_p` is its exponential reshaped to `(n_patches × n)`, and the PMFs are
the mean of the per-draw PMFs. The parent's uncertainty about the patch
trajectory does not enter here: it enters the fit through the shared
quantity of [`zone_meld_block`](@ref), which deforms this mean curve.
Returns `(; log_infections, g, f, death_pmf, log_importation, priors)`.
"""
function zone_parent_inputs(chn)
    chn = _zone_parent_chain(chn)
    keys_ = _ZONE_PARENT_KEYS
    infs = _zone_parent_draws(chn, keys_.infections; vectors = true)
    ndraws = length(infs)
    logI = let acc = log.(safe_rate.(Float64.(infs[1])))
        for i in 2:ndraws
            acc .+= log.(safe_rate.(Float64.(infs[i])))
        end
        acc ./ ndraws
    end
    α = _zone_parent_draws(chn, keys_.gi_alpha)
    θ = _zone_parent_draws(chn, keys_.gi_theta)
    inc_m = _zone_parent_draws(chn, keys_.inc_mean)
    inc_s = _zone_parent_draws(chn, keys_.inc_sd)
    rec_m = _zone_parent_draws(chn, keys_.receipt_mean)
    rec_s = _zone_parent_draws(chn, keys_.receipt_sd)
    g = _zone_mean_pmf(i -> _zone_gi_pmf(α[i], θ[i]), ndraws)
    inc_pmf(i) = discretise_censored(
        lognormal_meansd(inc_m[i], inc_s[i]),
        ZONE_INCUBATION_NMAX
    )
    rec_pmf(i) = discretise_censored(
        lognormal_meansd(rec_m[i], rec_s[i]),
        ZONE_RECEIPT_NMAX
    )
    f = _zone_mean_pmf(i -> convolve_pmf(inc_pmf(i), rec_pmf(i)), ndraws)
    death_pmf = if _has_key(chn, keys_.death_confirmation)
        dc = _draw_vectors(chn, keys_.death_confirmation)
        _zone_mean_pmf(
            i -> convolve_pmf(inc_pmf(i), Float64.(dc[i])), ndraws
        )
    else
        Float64[]
    end
    ## The province model's own between-patch movement, at its posterior
    ## mean: the per-origin export intensity, and the arrivals into each
    ## patch as a fraction of that patch's infections. The zone stage reads
    ## the intensity as a relative weight across origin patches and the
    ## fraction as the share of a patch's infections that came from
    ## elsewhere, so the between-patch flows it applies are the province
    ## model's rather than a second estimate of them.
    origin_epsilon = _mean_parent_vector(chn, keys_.importation_epsilon)
    ## The province model's own relative ascertainment and fatality, at its
    ## posterior mean. A factor common to a patch cancels in a within-patch
    ## composition, so these never reach the likelihood; they carry the
    ## province level of the multipliers into the reported quantities.
    province_ascertainment = _mean_parent_vector(
        chn, keys_.province_ascertainment
    )
    province_severity = _mean_parent_vector(chn, keys_.province_severity)
    import_fraction = _mean_parent_ratio(
        chn, keys_.importation, keys_.infections
    )
    ## Priors the zone stage inherits from the province posterior, fitted
    ## to its draws: a log-normal for each positive scale and a Beta for
    ## each composition overdispersion.
    priors = (;
        ascertainment_sd = _zone_parent_lognormal(
            chn, keys_.ascertainment_sd, (log(0.2), 0.5)
        ),
        drift_sd = _zone_parent_lognormal(chn, keys_.drift_sd, (log(0.05), 0.5)),
        rho = _zone_parent_beta(chn, keys_.rho, (1.0, 24.0)),
        rho_death = _zone_parent_beta(chn, keys_.rho_death, (1.0, 24.0)),
        correlation = _zone_parent_beta(
            chn, keys_.correlation, (1.0, 3.0); clamp_unit = true
        ),
    )
    return (;
        log_infections = logI, g, f, death_pmf,
        origin_epsilon, import_fraction, priors,
        province_ascertainment, province_severity,
    )
end

## Cumulative count of `zone` in `prov` at the last vintage of `history`, or
## zero when the zone is absent from the block.
## One zone's cumulative series, or `nothing` when the table omits it.
function _zone_history_series(
        history, prov::AbstractString,
        zone::AbstractString
    )
    haskey(history, prov) || return nothing
    haskey(history[prov], zone) || return nothing
    return history[prov][zone]
end

function _zone_last_cumulative(
        history, prov::AbstractString,
        zone::AbstractString
    )
    haskey(history, prov) || return 0
    haskey(history[prov], zone) || return 0
    c = history[prov][zone].counts
    return isempty(c) ? 0 : Int(c[end])
end

"""
$(TYPEDSIGNATURES)

The two blocks of the zone importation kernel, both normalisations of one
gravity pull over every zone ([`gravity_pull`](@ref)).

`within[z, q]`, for `z` and `q` in the same patch, is column-stochastic
within that patch, so a zone's within-patch spill is a transfer that
conserves the patch total.

`between[z, q]`, for zones in different patches, carries the province
model's own patch-to-patch flow:

```math
K^b_{zq} = K_{p(z)p(q)}\\,
    \\frac{\\mathrm{pull}_{zq}}{\\sum_{z' \\in p(z)} \\mathrm{pull}_{z'q}},
```

so its column over a destination patch's zones sums to that patch's entry
of `parent_kernel`, the province model's own
[`province_importation_kernel`](@ref). Summed over the zones of a patch,
the zone stage's between-patch flow is the province model's for the same
origin intensity, which is what keeps one movement from being counted at
both levels.
"""
function zone_importation_blocks(
        pops::AbstractVector, coords::AbstractVector,
        patch_of_zone::AbstractVector{<:Integer},
        parent_kernel::AbstractMatrix;
        decay::Real = PROVINCE_DISTANCE_DECAY
    )
    nz = length(pops)
    pull = gravity_pull(
        pops; distances = province_distance_matrix(coords), decay
    )
    np = size(parent_kernel, 1)
    within = zeros(Float64, nz, nz)
    between = zeros(Float64, nz, nz)
    col = zeros(Float64, np)
    @inbounds for q in 1:nz
        pq = patch_of_zone[q]
        ## Column totals of the pull into each patch, so both blocks
        ## normalise over one sweep of the origin's column.
        fill!(col, 0.0)
        for z in 1:nz
            col[patch_of_zone[z]] += pull[z, q]
        end
        for z in 1:nz
            p = patch_of_zone[z]
            col[p] > 0 || continue
            if p == pq
                within[z, q] = pull[z, q] / col[p]
            else
                between[z, q] = parent_kernel[p, pq] * pull[z, q] / col[p]
            end
        end
    end
    return (; within, between)
end

"""
$(TYPEDSIGNATURES)

The observations and health-zone metadata reduced to the `top_zones` zones
of each province with most confirmed cases at the cut-off, the rest of the
province's zones pooled into one `pool_name` zone so every patch total is
still allocated. The pooled zone's cumulative case and death series are
the sums of the pooled zones', its population their sum and its centroid
their population-weighted mean, so mixing and the distance correlation
stay on. A province with `top_zones` zones or fewer is unchanged, and the
`unallocated` rows pass through. For fast iteration on real data rather
than for reporting: [`fit_zone`](@ref) takes `top_zones` and applies this.
Returns `(; obs, zones, kept)`, with `kept` the retained zone keys per
province.
"""
function zone_subset(
        obs; top_zones::Integer, zones = _default_health_zones(),
        pool_name::AbstractString = "rest"
    )
    top_zones >= 1 || throw(ArgumentError("top_zones must be at least one"))
    conf = obs.zone_confirmed_history
    kept = Dict{String, Vector{String}}()
    for (prov, block) in conf
        named = [z for z in keys(block) if z != "unallocated"]
        last_count(z) = isempty(block[z].counts) ? 0 : Int(block[z].counts[end])
        sort!(named; by = z -> (-last_count(z), z))
        kept[prov] = named[1:min(top_zones, end)]
    end
    pool(history) = begin
        out = Dict{String, Dict{String, NamedTuple}}()
        for (prov, block) in history
            keep = get(kept, prov, String[])
            newblock = Dict{String, NamedTuple}()
            rest = nothing
            pooled = 0
            for (z, h) in block
                if z == "unallocated" || z in keep
                    newblock[z] = h
                    continue
                end
                pooled += 1
                if rest === nothing
                    rest = (; days = copy(h.days), counts = collect(Int, h.counts))
                else
                    rest.days == h.days || error(
                        "zone_subset: zone `$(prov).$(z)` is reported on " *
                            "different vintage days to the zones pooled with it."
                    )
                    rest = (; rest.days, counts = rest.counts .+ Int.(h.counts))
                end
            end
            ## A single pooled zone keeps its own name rather than becoming
            ## the rest of nothing.
            if rest !== nothing && pooled == 1
                for (z, h) in block
                    (z == "unallocated" || z in keep) && continue
                    newblock[z] = h
                    push!(kept[prov], z)
                end
            elseif rest !== nothing
                newblock[pool_name] = rest
            end
            out[prov] = newblock
        end
        out
    end
    obs2 = merge(obs, (; zone_confirmed_history = pool(conf)))
    if hasproperty(obs, :zone_death_history)
        obs2 = merge(obs2, (; zone_death_history = pool(obs.zone_death_history)))
    end
    zones2 = zones === nothing ? nothing : _zone_subset_metadata(
            zones, kept, obs2.zone_confirmed_history, pool_name
        )
    return (; obs = obs2, zones = zones2, kept)
end

## Metadata rows for the kept zones plus one pooled row per province that
## has a pooled zone: the population summed, the centroid population
## weighted.
function _zone_subset_metadata(zones, kept, conf, pool_name)
    out = eltype(zones)[]
    for r in zones
        keep = get(kept, r.province, String[])
        r.zone in keep && push!(out, r)
    end
    for (prov, block) in conf
        haskey(block, pool_name) || continue
        keep = get(kept, prov, String[])
        rows = [
            r for r in zones
                if r.province == prov && !(r.zone in keep) && r.zone != "unallocated"
        ]
        isempty(rows) && continue
        pop = sum(r.population for r in rows)
        w = pop > 0 ? [r.population / pop for r in rows] :
            fill(1 / length(rows), length(rows))
        push!(
            out, (;
                zone = String(pool_name),
                label = "Rest of " * titlecase(replace(prov, "_" => " ")),
                province = prov, population = pop,
                lat = sum(w[i] * rows[i].lat for i in eachindex(rows)),
                lon = sum(w[i] * rows[i].lon for i in eachindex(rows)),
                zscode = "",
            )
        )
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Everything the zone model and its render need, built once from the
parent chain and the observations: the fixed stage-1 inputs
([`zone_parent_inputs`](@ref)), the zone units and their data, the grid
and the sampler's starting point.

Zones are the health-zone rows of `obs.zone_confirmed_history`
([`zone_increment_matrix`](@ref)), nested in the stage-1 patches of
[`PROVINCE_NAMES`](@ref) and ordered patch by patch; `unallocated`
pseudo-rows are never units. Per vintage the zone increments are the
consecutive-vintage differences clamped at zero, the first vintage's
increment its cumulative, and the allocated patch total their sum. A
vintage on which a patch's unallocated count falls in either zone block
is a reattribution into the zones and is zeroed for that patch
([`zone_increment_matrix`](@ref)); `excluded` lists those `(patch, date)`
pairs. Cells with a zero allocated total are dropped here, so the model
scores only positive totals. The walking set is the zones whose
cumulative confirmed count at the cut-off is at least `walk_threshold`,
in patches with at least two such zones. The grid starts `lead_days`
before the first vintage and
carries knots every `week` days ([`knot_days`](@ref)) to the cut-off.
`share_scale` is the softmax scale of the initial shares and `rt_floor`
the cumulative zone infections below which `R_T_zone` is undefined; both
are stored in `model_data` for the model and in the returned inputs for
the render, so the two read one value.

The shared quantity's multivariate normal is fitted here too
([`zone_meld_block`](@ref)), over the weekly windows of the same knot grid
and every patch together, and `meld_min_infections` is the parent mean
below which a window is dropped as not yet seeded.

`parent_forecast`, the parent's posterior-predictive draws
([`forecast_draws`](@ref)), adds the inputs the forecast reads
([`zone_forecast_block`](@ref)) over `horizon` days as `model_data.forecast`.
The fitted model does not read them, so a chain fitted without them serves
the forecast model.

`zones` is the metadata table of [`load_health_zones`](@ref), read from the
package data by default, for labels and, with mixing, populations and
centroids. Returns a named tuple whose `model_data` field is the plain-array
input of [`bvd_zone`](@ref), alongside the zone keys and labels, the patch
index and labels, the vintage days and dates, the cumulative counts, the
walking mask, the knots, the grid start, the initial-share start values
and the shared quantity's block.
"""
function zone_fit_inputs(
        parent_chain, obs;
        walk_threshold::Integer = 30,
        lead_days::Integer = 42,
        week::Integer = 7,
        share_scale::Real = 2.0,
        rt_floor::Real = 10.0,
        meld_min_infections::Real = 1.0,
        zones = _default_health_zones(),
        patch_names::AbstractVector = PROVINCE_NAMES,
        patch_labels::AbstractVector = PROVINCE_LABELS,
        parent_forecast = nothing,
        horizon::Integer = 7
    )
    hasproperty(obs, :zone_confirmed_history) || error(
        "zone_fit_inputs: `obs` carries no `zone_confirmed_history`; load " *
            "the observations with a manifest that has the zone tables."
    )
    parent_chain = _zone_parent_chain(parent_chain)
    n = obs.n
    parent = zone_parent_inputs(parent_chain)
    np = length(patch_names)
    length(parent.log_infections) == np * n || error(
        "zone_fit_inputs: the parent's `infections_patch` holds " *
            "$(length(parent.log_infections)) entries but $np patches by $n " *
            "days is $(np * n); the parent must be fitted to the same data " *
            "cut-off."
    )
    I_bar = exp.(reshape(parent.log_infections, np, n))
    ## Zone units and counts, patch by patch. A vintage on which a patch's
    ## unallocated count falls in either block is a reattribution into the
    ## zones and is left out of that patch's composition.
    death_history = hasproperty(obs, :zone_death_history) ?
        obs.zone_death_history : Dict{String, Any}()
    reattribution = mergewith(
        vcat,
        zone_reattribution_days(obs.zone_confirmed_history),
        zone_reattribution_days(death_history)
    )
    per_patch = zone_increment_matrix(
        obs.zone_confirmed_history, patch_names;
        reattribution
    )
    excluded = [
        (;
            patch = String(entry.patch),
            date = obs.seeding + Day(d - 1),
        )
            for entry in per_patch for d in entry.excluded
    ]
    isempty(per_patch) && error(
        "zone_fit_inputs: the observations carry no zone histories."
    )
    zone_province = String[]
    zone_names = String[]
    patch_of_zone = Int[]
    patch_ranges = UnitRange{Int}[]
    days = Int[]
    blocks = Matrix{Int}[]
    for (p, entry) in enumerate(per_patch)
        start = length(zone_names) + 1
        for (prov, zone) in entry.zones
            push!(zone_province, prov)
            push!(zone_names, zone)
            push!(patch_of_zone, p)
        end
        push!(patch_ranges, start:length(zone_names))
        isempty(entry.zones) && continue
        if isempty(days)
            days = copy(entry.days)
        elseif entry.days != days
            error(
                "zone_fit_inputs: patch `$(entry.patch)` is reported on " *
                    "different vintage days to the first patch with zones; " *
                    "the zone tables must share one `dates` array."
            )
        end
        push!(blocks, entry.increments)
    end
    nz = length(zone_names)
    nz > 0 || error("zone_fit_inputs: no health zones in the observations.")
    nv = length(days)
    counts = reduce(vcat, blocks)
    size(counts) == (nz, nv) || error(
        "zone_fit_inputs: the zone count blocks do not stack to " *
            "$(nz) zones by $(nv) vintages."
    )
    zone_keys = [zone_province[z] * "." * zone_names[z] for z in 1:nz]
    labels = _zone_labels(zones, zone_province, zone_names)
    ## Scored cells: every (patch, vintage) with a positive allocated total.
    cell_patch = Int[]
    cell_vintage = Int[]
    cell_total = Int[]
    cell_const = Float64[]
    for v in 1:nv, (p, zs) in enumerate(patch_ranges)

        isempty(zs) && continue
        N = sum(@view counts[zs, v])
        N > 0 || continue
        push!(cell_patch, p)
        push!(cell_vintage, v)
        push!(cell_total, N)
        push!(cell_const, _zone_cell_const(counts, zs, v))
    end
    ## Cumulative confirmed at the cut-off, from the manifest rather than the
    ## clamped increments, and the walking set.
    cumulative = [
        _zone_last_cumulative(
            obs.zone_confirmed_history,
            zone_province[z], zone_names[z]
        ) for z in 1:nz
    ]
    eligible = cumulative .>= walk_threshold
    walking = falses(nz)
    for zs in patch_ranges
        count(eligible[zs]) >= 2 || continue
        walking[zs] .= eligible[zs]
    end
    walk_index = zeros(Int, nz)
    n_walking = 0
    for z in 1:nz
        walking[z] || continue
        n_walking += 1
        walk_index[z] = n_walking
    end
    ## Grid and knots.
    t0 = clamp(days[1] - lead_days, 1, n)
    knots = knot_days(n; week, start = t0)
    fixed = zone_fixed_terms(I_bar, parent.g, parent.f, t0)
    nd = n - t0 + 1
    interp = zone_interpolation_weights(knots, t0, n)
    report_matrix = zone_delay_operator(parent.f, nd)
    report_pre_rows = zone_report_pre_rows(
        fixed.report_pre, patch_of_zone,
        t0, n
    )
    ## Death composition, scored per vintage as the cases are. The rows are
    ## built by looking each zone up by name rather than by restacking the
    ## death table, so they line up with the case rows even where a zone is
    ## missing from one table. A zone absent from the death table counts
    ## zero throughout.
    ##
    ## The exclusions widen the unallocated rule. A revision that moves
    ## deaths out of named zones leaves the unallocated row flat or rising,
    ## so the unallocated rule alone does not see it, and the increment
    ## clamp would absorb the fall silently.
    death_reattribution = zone_reattribution_days(
        death_history;
        include_zone_falls = true
    )
    death_counts = zeros(Int, nz, nv)
    for z in 1:nz
        h = _zone_history_series(
            death_history, zone_province[z],
            zone_names[z]
        )
        h === nothing && continue
        h.days == days || error(
            "zone_fit_inputs: the death table for `$(zone_province[z])." *
                "$(zone_names[z])` is reported on different vintage days to " *
                "the confirmed-case tables; both compositions read one grid."
        )
        prev = 0
        for v in 1:nv
            c = Int(h.counts[v])
            death_counts[z, v] = max(c - prev, 0)
            prev = c
        end
    end
    for (p, zs) in enumerate(patch_ranges)
        isempty(zs) && continue
        provs = unique(zone_province[z] for z in zs)
        for pr in provs, d in get(death_reattribution, pr, Int[])

            v = findfirst(==(d), days)
            v === nothing && continue
            death_counts[zs, v] .= 0
        end
    end
    death_cell_patch = Int[]
    death_cell_vintage = Int[]
    death_cell_total = Int[]
    death_cell_const = Float64[]
    for v in 1:nv, (p, zs) in enumerate(patch_ranges)

        isempty(zs) && continue
        N = sum(@view death_counts[zs, v])
        N > 0 || continue
        push!(death_cell_patch, p)
        push!(death_cell_vintage, v)
        push!(death_cell_total, N)
        push!(death_cell_const, _zone_cell_const(death_counts, zs, v))
    end
    death_pmf = isempty(parent.death_pmf) ? [1.0] : parent.death_pmf
    death_fixed = zone_fixed_terms(I_bar, parent.g, death_pmf, t0)
    death_matrix = zone_delay_operator(death_pmf, nd)
    death_pre_rows = zone_report_pre_rows(
        death_fixed.report_pre,
        patch_of_zone, t0, n
    )
    ## Great-circle distances between the centroids of each patch's zones,
    ## over all of them for the deviation level and over its walking zones
    ## for the innovations.
    zone_distances = _zone_distance_blocks(
        zones, zone_province, zone_names,
        patch_ranges
    )
    zone_walk_distances = _zone_distance_blocks(
        zones, zone_province, zone_names,
        patch_ranges; keep = walking
    )
    ## The distance the inherited deviation correlation refers to: the mean
    ## great-circle distance between the province capitals, which is the
    ## separation the province model's own correlation was learned at.
    correlation_distance = _mean_offdiagonal(
        province_distance_matrix(
            PROVINCE_CAPITALS[1:min(np, length(PROVINCE_CAPITALS))]
        )
    )
    ## The zone mixing structure: the two kernel blocks, the province
    ## model's per-origin intensity as a relative weight across origin
    ## patches, and its arrivals into each patch as a fraction of that
    ## patch's infections. All fixed from the parent posterior. Only the
    ## within-patch intensity is sampled here.
    mixing = _zone_mixing_or_nothing(
        zones, zone_province, zone_names, patch_ranges, patch_of_zone,
        parent, np, n
    )
    ## The shared quantity: the parent's weekly patch infections, over the
    ## same weekly grid the deviations use, as one multivariate normal.
    meld = zone_meld_block(
        _zone_parent_draws(
            parent_chain, _ZONE_PARENT_KEYS.infections;
            vectors = true
        ),
        np, n, knots, I_bar, t0; min_infections = meld_min_infections
    )
    forecast = parent_forecast === nothing ? nothing :
        zone_forecast_block(
            parent_forecast, meld, I_bar, parent.g, parent.f, death_pmf, t0,
            knots, patch_of_zone, mixing; horizon, week,
            min_infections = meld_min_infections
        )
    model_data = (;
        counts, cell_patch, cell_vintage, cell_total, cell_const,
        days, I_bar, g = parent.g, f = parent.f, patch_ranges, patch_of_zone,
        knots, t0, n,
        week, share_scale, rt_floor,
        walking = collect(walking), walk_index, n_walking,
        fixed.force_pre, fixed.report_pre_cum,
        fixed.infections_pre, mixing, interp, report_matrix,
        report_pre_rows,
        death_counts, death_cell_patch, death_cell_vintage,
        death_cell_total, death_cell_const, death_days = days,
        death_pre_cum = death_fixed.report_pre_cum, death_matrix,
        death_pre_rows,
        meld_weights = meld.weights, meld_L = meld.L, meld_d = meld.d,
        zone_distances, zone_walk_distances,
        correlation_distance, parent_priors = parent.priors,
        province_ascertainment = _zone_province_factor(
            parent.province_ascertainment, patch_of_zone
        ),
        province_severity = _zone_province_factor(
            parent.province_severity, patch_of_zone
        ),
        forecast,
    )
    dates = [obs.seeding + Day(d - 1) for d in days]
    return (;
        model_data, zone_keys, zone_labels = labels, zone_province,
        zone_names, patch_of_zone, patch_ranges,
        patch_names = collect(String, patch_names),
        patch_labels = collect(String, patch_labels),
        days, dates, counts, cumulative, excluded, walking = collect(walking),
        walk_threshold, knots, t0, n, week, share_scale, rt_floor,
        seeding = obs.seeding, cutoff = obs.cutoff,
        I_bar, g = parent.g, f = parent.f, death_pmf, meld,
    )
end

## One province factor per zone from a per-patch vector, all ones when the
## parent chain does not carry it.
function _zone_province_factor(
        values::AbstractVector, patch_of_zone::AbstractVector{<:Integer}
    )
    isempty(values) && return ones(Float64, length(patch_of_zone))
    return Float64[values[p] for p in patch_of_zone]
end

## The package's health-zone metadata, or `nothing` when the file is absent.
function _default_health_zones()
    path = joinpath(@__DIR__, "..", "..", "data", "health_zones.csv")
    isfile(path) || return nothing
    return load_health_zones(path)
end

## Display label per zone from the metadata, falling back to the key.
function _zone_labels(
        zones, zone_province::AbstractVector,
        zone_names::AbstractVector
    )
    fallback(z) = titlecase(replace(zone_names[z], "_" => " "))
    zones === nothing && return [fallback(z) for z in eachindex(zone_names)]
    lookup = Dict((r.province, r.zone) => r.label for r in zones)
    return [
        get(lookup, (zone_province[z], zone_names[z]), fallback(z))
            for z in eachindex(zone_names)
    ]
end

## The zone mixing structure, or `nothing` when the metadata misses a zone
## or the parent chain carries no between-patch movement (mixing is then
## refused by `fit_zone`).
function _zone_mixing_or_nothing(
        zones, zone_province, zone_names, patch_ranges, patch_of_zone,
        parent, np::Integer, n::Integer
    )
    nz = length(zone_names)
    zones === nothing && return nothing
    lookup = Dict((r.province, r.zone) => r for r in zones)
    rows = [
        get(lookup, (zone_province[z], zone_names[z]), nothing)
            for z in 1:nz
    ]
    any(isnothing, rows) && return nothing
    length(parent.origin_epsilon) == np || return nothing
    length(parent.import_fraction) == np * n || return nothing
    pops = Float64[r.population for r in rows]
    coords = [(r.lat, r.lon) for r in rows]
    blocks = zone_importation_blocks(
        pops, coords, patch_of_zone,
        province_importation_kernel(
            PROVINCE_POPULATIONS[1:min(np, length(PROVINCE_POPULATIONS))]
        )
    )
    ## The intensity is a relative weight across origin patches inside a
    ## pattern that is normalised over the destination patch, so its overall
    ## scale cancels; it is centred at one so the weights stay near unity,
    ## and carried per zone so the renewal reads it without a patch lookup.
    ε_bar = parent.origin_epsilon ./ max(mean(parent.origin_epsilon), eps())
    origin_weight = ε_bar[patch_of_zone]
    ## Arrivals as a fraction of each patch's own infections, below one so
    ## the locally grown part of a patch's total stays positive under any
    ## deformation of the shared quantity.
    import_fraction = reshape(parent.import_fraction, np, n)
    return (;
        blocks.within, blocks.between, origin_weight,
        import_fraction, patch_of_zone,
    )
end

## --- Fitting -------------------------------------------------------------

## Iterations-by-chains matrix of a sampler statistic, or `nothing`.
function _zone_stat(chn, name::Symbol)
    for e in FlexiChains.extras(chn)
        e.name === name || continue
        return Float64.(coalesce.(chn[e], NaN))
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Fit the health-zone model [`bvd_zone`](@ref) melded from `parent_chain`, a
`bvd_joint` chain with the patch structure on, to the zone tables in `obs`.
Builds the fixed inputs with [`zone_fit_inputs`](@ref) and runs
[`nuts_sample`](@ref) at the settings the caller gives, which the fit
registry sets to the ones every other fit in the report uses. Logs the
adapted step size, the divergences and the fraction of iterations at the
tree-depth cap per chain. The zones mix wherever the health-zone metadata
and the parent chain allow it, and the log says when they do not.
`top_zones` fits the [`zone_subset`](@ref) of the observations, for fast
iteration. `walk_threshold`, `lead_days`, `week`, `share_scale`, `rt_floor`, `zones`,
`patch_names` and `patch_labels` pass through to
[`zone_fit_inputs`](@ref), `model_kwargs` to [`bvd_zone`](@ref) and any
other keyword to [`nuts_sample`](@ref). Priors belong in `model_kwargs`;
passed loose they reach the sampler, which drops what it does not know.
Returns the chain.
"""
function fit_zone(
        parent_chain, obs;
        samples::Integer = 500,
        chains::Integer = 2,
        n_adapts::Integer = min(200, samples ÷ 2),
        target_accept::Real = 0.85,
        max_depth::Integer = 10,
        seed::Integer = 20260518,
        callback = nothing,
        top_zones::Union{Nothing, Integer} = nothing,
        walk_threshold::Integer = 30,
        lead_days::Integer = 42,
        week::Integer = 7,
        share_scale::Real = 2.0,
        rt_floor::Real = 10.0,
        zones = _default_health_zones(),
        patch_names::AbstractVector = PROVINCE_NAMES,
        patch_labels::AbstractVector = PROVINCE_LABELS,
        model_kwargs = (;),
        kwargs...
    )
    if top_zones !== nothing
        sub = zone_subset(obs; top_zones, zones)
        obs, zones = sub.obs, sub.zones
    end
    inputs = zone_fit_inputs(
        parent_chain, obs;
        walk_threshold, lead_days, week, share_scale, rt_floor, zones,
        patch_names, patch_labels
    )
    zd = inputs.model_data
    zd.mixing === nothing && @warn(
        "fit_zone: the zones are not mixed; mixing needs a population and " *
            "centroid for every zone in the health-zone metadata, and a " *
            "parent chain carrying the province model's between-patch " *
            "movement."
    )
    isempty(zd.death_cell_patch) &&
        error(
        "fit_zone: the observations carry no allocated zone deaths, which " *
            "the death composition is fitted to."
    )
    model = bvd_zone(zd; model_kwargs...)
    n_zones = length(inputs.zone_keys)
    n_walking = count(inputs.walking)
    n_knots = length(inputs.knots)
    n_cells = length(zd.cell_patch)
    n_shared = zd.meld_d
    mixing = zd.mixing !== nothing
    correlation = !isempty(zd.zone_distances)
    dimension = length(
        link(VarInfo(model), model)[:]
    )
    @info(
        "fit_zone: starting", n_zones, n_walking, n_knots, n_cells,
        n_shared, mixing, correlation, dimension
    )
    t = time()
    chn = nuts_sample(
        model; samples, chains, n_adapts, target_accept,
        max_depth, seed, callback, kwargs...
    )
    elapsed = time() - t
    steps = _zone_stat(chn, :step_size)
    depth = _zone_stat(chn, :tree_depth)
    div = _zone_stat(chn, :numerical_error)
    ## Per-chain values as one string each: the logger elides vectors in a
    ## CI log, and these are the numbers the log exists to show.
    per_chain(m, f) = m === nothing ? "n/a" :
        join((string(f(view(m, :, c))) for c in 1:size(m, 2)), " / ")
    minutes = round(elapsed / 60; digits = 1)
    step_size = per_chain(steps, x -> round(x[end]; sigdigits = 3))
    depth_cap_fraction = per_chain(
        depth, x -> round(mean(x .>= max_depth); digits = 3)
    )
    divergences = per_chain(div, x -> Int(sum(x)))
    @info(
        "fit_zone: finished", minutes, step_size, depth_cap_fraction,
        divergences
    )
    return chn
end
