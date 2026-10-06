using Test
include(joinpath(@__DIR__, "..", "harness.jl"))
isfile(joinpath(@__DIR__, "..", "expected.jl")) && include(joinpath(@__DIR__, "..", "expected.jl"))
isfile(joinpath(@__DIR__, "..", "fit.jl")) && include(joinpath(@__DIR__, "..", "fit.jl"))
isfile(joinpath(@__DIR__, "..", "toy.jl")) && include(joinpath(@__DIR__, "..", "toy.jl"))
isfile(joinpath(@__DIR__, "..", "timestep.jl")) && include(joinpath(@__DIR__, "..", "timestep.jl"))

@testset "Phase 2 research code" begin
    for f in sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(@__DIR__)))
        include(joinpath(@__DIR__, f))
    end
end
