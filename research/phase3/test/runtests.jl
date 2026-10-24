using Test
include(joinpath(@__DIR__, "..", "harness.jl"))

@testset "Phase 3 research code" begin
    for f in sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(@__DIR__)))
        include(joinpath(@__DIR__, f))
    end
end
