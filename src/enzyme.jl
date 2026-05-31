# Enzyme reverse-mode AD glue: the package-default `AutoEnzyme` config
# and the EnzymeRules needed for the gamma-CDF and gamma normalising
# constants. Enzyme is a hard dependency, so these load with the package.

using Enzyme: Enzyme
using Enzyme.EnzymeRules: EnzymeRules
using SpecialFunctions: gamma, digamma

"""
Enzyme reverse-mode AD with a `Duplicated` function annotation (so the
closure over the observed data is differentiated, not treated as
read-only). Runtime activity is not enabled: it was previously needed
because the old `Integrals.solve`-based `integrate` returned `Any`,
defeating Enzyme's static activity analysis through the quadrature. Now
that `integrate` is type-stable, activity is resolved statically and
plain `Enzyme.Reverse` produces correct gradients (validated against
Mooncake and finite differences across the joint and single-stream
models). Dropping runtime activity is faster and allocates less than the
runtime-activity-on config. This is the config the `_gamma_cdf` / `gamma`
rules below are validated against, and the package default returned by
`default_adtype`. Returns an `ADTypes.AutoEnzyme`; pass to
`nuts_sample(model; adtype = ...)`.
"""
function enzyme_adtype()
    return AutoEnzyme(;
        mode = Enzyme.Reverse,
        function_annotation = Enzyme.Duplicated)
end

# `EnzymeRules.@easy_rule` expands into both the reverse-mode
# (`augmented_primal` / `reverse`) and forward-mode (`forward`) rules for
# `_gamma_cdf`. The analytical (dα, dθ, dx) come from `_gamma_cdf_partials`
# in `src/gamma_cdf.jl`, the same helper used by the ChainRules rrule that
# Mooncake/ReverseDiff pick up. Routing `_gamma_cdf` through this rule
# avoids Enzyme differentiating `SpecialFunctions.gamma_inc` directly,
# which it cannot lower (recursive series + DomainError branches), and
# which the `@import_rrule` lift gets wrong on the shape partial.
EnzymeRules.@easy_rule(_gamma_cdf(α::Real, θ::Real, x::Real),
    @setup(_p=_gamma_cdf_partials(α, θ, x),
        dα=_p[1],
        dθ=_p[2],
        dx=_p[3],),
    (dα, dθ, dx))

# Rule for `SpecialFunctions.gamma`, derivative `d/dx Γ(x) = Γ(x) ψ(x)`
# (`Ω` binds to the primal `Γ(x)`). Enzyme's own `EnzymeSpecialFunctionsExt`
# ships no `gamma` rule and instead mis-lowers `gamma(x)` to the `loggamma`
# known-op, returning `ψ(x)` — wrong by a factor of `Γ(x)`. Distribution
# normalising constants (e.g. `Gamma`, `Beta`, `NegativeBinomial`) reach
# `gamma` outside the `_gamma_cdf` rule, so without this their gradients are
# corrupted.
EnzymeRules.@easy_rule(gamma(x::Real), (Ω * digamma(x),))
