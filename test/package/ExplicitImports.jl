@testitem "ExplicitImports: no stale" tags = [:quality] begin
    using ExplicitImports, BVDOutbreakSize
    @test check_no_stale_explicit_imports(BVDOutbreakSize) === nothing
end

@testitem "ExplicitImports: no implicit" tags = [:quality] begin
    using ExplicitImports, BVDOutbreakSize
    @test check_no_implicit_imports(BVDOutbreakSize) === nothing
end

@testitem "ExplicitImports: via owners" tags = [:quality] begin
    using ExplicitImports, BVDOutbreakSize
    @test check_all_explicit_imports_via_owners(BVDOutbreakSize) === nothing
end
