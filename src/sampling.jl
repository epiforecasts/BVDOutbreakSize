# Sampler glue: the package-default AD type and the NUTS driver used
# to fit every Turing model.

"""
The package-default AD type used as the NUTS `adtype` keyword. Returns
the Enzyme reverse-mode (runtime-activity-off) config from
`enzyme_adtype`. Enzyme is the default because it is validated correct
against Mooncake and finite differences, is ~14% faster per gradient,
and the type-stable Gauss-Legendre `integrate` removed the need for
runtime activity. Mooncake remains available via `mooncake_adtype`.
"""
default_adtype() = enzyme_adtype()

"""
Mooncake reverse-mode AD with default `Mooncake.Config()`. The previous
package default, kept as an opt-in backend; pass to
`nuts_sample(model; adtype = mooncake_adtype())`.
"""
mooncake_adtype() = AutoMooncake(; config = Mooncake.Config())

"""
NUTS on `model`, parallel chains via `MCMCThreads`. Chains
initialise from the prior (`InitFromPrior()`) to keep the sampler
in regions with reasonable physical interpretation. Pass `init =
Turing.DynamicPPL.InitFromUniform()` to fall back to unconstrained
uniform initialisation.
"""
function nuts_sample(model;
        samples::Integer = 1_000,
        chains::Integer = 4,
        target_accept::Real = 0.95,
        seed::Integer = 20260518,
        progress::Bool = false,
        adtype = default_adtype(),
        init = InitFromPrior())
    rng = MersenneTwister(seed)
    return sample(
        rng,
        model,
        NUTS(target_accept; adtype),
        MCMCThreads(),
        samples, chains;
        initial_params = fill(init, chains),
        progress = progress
    )
end
