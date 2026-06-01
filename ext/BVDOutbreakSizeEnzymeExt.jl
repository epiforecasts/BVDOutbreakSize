module BVDOutbreakSizeEnzymeExt

import BVDOutbreakSize
using ADTypes: AutoEnzyme
using Enzyme: Enzyme

# Reverse-mode Enzyme with runtime activity (so per-value activity is
# resolved through the distribution constructors) and a `Duplicated`
# function annotation (so the closure over the observed data is
# differentiated, not treated as read-only). An opt-in alternative to the
# default Mooncake backend; differentiating the full renewal joint under
# Enzyme is still work in progress.
#
# The `SpecialFunctions.gamma` EnzymeRule that the Beta and
# NegativeBinomial normalising constants need is supplied by
# `CensoredDistributions`'s own Enzyme extension (a package dependency),
# so it is not redefined here — doing so would overwrite that method and
# break precompilation.
function BVDOutbreakSize.enzyme_adtype()
    return AutoEnzyme(;
        mode = Enzyme.set_runtime_activity(Enzyme.Reverse),
        function_annotation = Enzyme.Duplicated)
end

end
