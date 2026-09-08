## A `Core.Box` in a model body is a silent, expensive defect, and this is the
## only test that can see it.
##
## Julia boxes any local that a closure captures and something later
## reassigns. Two idioms in the observation models did it: the
## `if eltype(x) === Any; x = convert(...); end` guards, and accumulators
## rebound inside a loop whose comprehensions capture them. A boxed local is
## type-unstable at every use, and Mooncake answers type instability with
## `DynamicDerivedRule` — a dictionary lookup per call site, on every gradient
## evaluation.
##
## Nothing about that is visible to an ordinary test. The gradient is correct,
## the log-density is correct, the parameter count is unchanged and the model
## samples: the only symptom is an AD-to-primal ratio nobody can account for.
## A profiler does not localise it either, because a boxed closure's cost is
## charged to the frame that CAPTURES the variable rather than the loop that
## PAYS it, so the cost appears in the composer rather than in the submodel.
##
## Removing the boxes left the log-density bit-identical and took roughly a
## third off every gradient, so the pattern is worth keeping out. If this test
## fails, look for a variable that is assigned inside an `if` (or rebound in a
## loop) and also read by a comprehension, `map`, or `do` block nearby; the fix
## is to compute into a separate binding and assign the final name once.
@testitem "no boxed captures in the observation models on the AD path" begin
    using BVDOutbreakSize: confirmed_cases_model, confirmed_deaths_model,
                           treatment_flow_model

    for f in (confirmed_cases_model, confirmed_deaths_model,
        treatment_flow_model)
        boxes = sum(Base.code_lowered(f); init = 0) do ci
            count(x -> occursin("Core.Box", string(x)), ci.code)
        end
        @test boxes == 0
    end
end
