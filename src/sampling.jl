# Sampler glue: the package-default AD type and the NUTS driver used
# to fit every Turing model.

"""
Mooncake reverse-mode AD with default `Mooncake.Config()`. Used as
the NUTS `adtype` keyword.
"""
default_adtype() = AutoMooncake(; config = Mooncake.Config())

"""
    enzyme_adtype()

Enzyme reverse-mode AD type, a validated opt-in alternative to the default
[`default_adtype`](@ref) (Mooncake). Defined by the package's Enzyme
weak-dependency extension (`ext/BVDOutbreakSizeEnzymeExt.jl`); calling it
without `Enzyme` loaded raises a `MethodError`. Load `Enzyme` to activate
the extension. The `SpecialFunctions.gamma` `EnzymeRule` that the Beta and
NegativeBinomial normalising constants reach is supplied by
CensoredDistributions' own Enzyme extension. The full renewal joint
differentiates under Enzyme and the gradient matches Mooncake (see
`test/test_enzyme.jl`); Mooncake remains the package default.
"""
function enzyme_adtype end

"""
NUTS on `model`, parallel chains via `MCMCThreads`. Chains
initialise from the prior (`InitFromPrior()`) to keep the sampler
in regions with reasonable physical interpretation. Pass `init =
Turing.DynamicPPL.InitFromUniform()` to fall back to unconstrained
uniform initialisation.

`check_model = false` disables Turing's pre-sampling model check, which
rejects any model with a sampled discrete variable even when its value
feeds nothing downstream. The per-vintage DRC streams are now scored as
observed `~` data, so a composer conditioning on them passes the check
with the default `check_model = true`. The escape is needed only by
[`exports_deaths_only_model`](@ref), which runs the exports submodel in
predictive mode (`exported_cases ~ Poisson` with a `missing` count) purely
for the export onsets, leaving a sampled discrete `Poisson` draw. The
continuous parameters are unaffected.
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
