## Tests for the package AD backends: `default_adtype` is now Enzyme
## reverse-mode, with Mooncake still available via `mooncake_adtype`.

@testitem "default_adtype returns an AutoEnzyme" begin
    using ADTypes: AutoEnzyme
    using Enzyme
    using BVDOutbreakSize: default_adtype
    ad = default_adtype()
    @test ad isa AutoEnzyme
end

@testitem "mooncake_adtype returns an AutoMooncake" begin
    using ADTypes: AutoMooncake
    using BVDOutbreakSize: mooncake_adtype
    ad = mooncake_adtype()
    @test ad isa AutoMooncake
end
