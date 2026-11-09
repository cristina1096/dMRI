using Test
include(joinpath(@__DIR__, "..", "harness.jl"))
include(joinpath(@__DIR__, "..", "continuum_radial.jl"))
include(joinpath(@__DIR__, "..", "continuum_annulus.jl"))
include(joinpath(@__DIR__, "..", "analysis.jl"))
include(joinpath(@__DIR__, "..", "p3_7_mc.jl"))
include(joinpath(@__DIR__, "..", "p3_7b_mc.jl"))

@testset "Phase 3 research code" begin
    for f in sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(@__DIR__)))
        include(joinpath(@__DIR__, f))
    end
end
