@testitem "Run docstring tests" tags = [:quality] begin
    using Documenter, BVDOutbreakSize
    doctest(BVDOutbreakSize)
end
