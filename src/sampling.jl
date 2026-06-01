# Sampler glue: the package-default AD type and the NUTS driver used
# to fit every Turing model.

"""
Mooncake reverse-mode AD with default `Mooncake.Config()`. Used as
the NUTS `adtype` keyword.
"""
default_adtype() = AutoMooncake(; config = Mooncake.Config())

"""
    enzyme_adtype()

Enzyme reverse-mode AD type, an opt-in alternative to the default
[`default_adtype`](@ref) (Mooncake). Defined by the package's Enzyme
weak-dependency extension (`ext/BVDOutbreakSizeEnzymeExt.jl`); calling it
without `Enzyme` loaded raises a `MethodError`. Load `Enzyme` to activate
the extension, which also installs the `EnzymeRules` for
`SpecialFunctions.gamma` that the Beta and NegativeBinomial normalising
constants reach. Differentiating the full renewal joint under Enzyme is
work in progress; Mooncake remains the package default.
"""
function enzyme_adtype end

"""
NUTS on `model`, parallel chains via `MCMCThreads`. Chains
initialise from the prior (`InitFromPrior()`) to keep the sampler
in regions with reasonable physical interpretation. Pass `init =
Turing.DynamicPPL.InitFromUniform()` to fall back to unconstrained
uniform initialisation.

`check_model = false` disables Turing's pre-sampling model check. A
composer that drops a stream (passes `missing`) leaves that stream's
likelihood as a sampled discrete draw (`Poisson` / `NegativeBinomial`)
whose value feeds nothing downstream. The check rejects any model with
a sampled discrete variable, even a redundant one, so a composer that
conditions on one stream while leaving another's count `missing` (e.g.
[`exports_deaths_only_model`](@ref), which keeps the deaths and exports
submodels only for their `CFR`, onset-to-death PMF and export onsets)
cannot otherwise be fitted. The continuous parameters are unaffected.
"""
function nuts_sample(model;
        samples::Integer = 1_000,
        chains::Integer = 4,
        target_accept::Real = 0.95,
        seed::Integer = 20260518,
        progress::Bool = false,
        adtype = default_adtype(),
        init = InitFromPrior(),
        check_model::Bool = true)
    rng = MersenneTwister(seed)
    return sample(
        rng,
        model,
        NUTS(target_accept; adtype),
        MCMCThreads(),
        samples, chains;
        initial_params = fill(init, chains),
        progress = progress,
        check_model = check_model
    )
end
