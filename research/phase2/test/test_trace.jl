include(joinpath(@__DIR__, "..", "p2_2_9_layer_trace_lib.jl"))

@testset "trace: independent layer fraction" begin
    # walls every 2 um, layers [0, h] and [2 - h, 2] in each gap
    @test sampled_layer_fraction(0.1, 0.3, 0.5) ≈ 1.0 atol=1e-4
    @test sampled_layer_fraction(0.6, 1.4, 0.5) == 0.0
    @test sampled_layer_fraction(0.4, 0.6, 0.5) ≈ 0.5 atol=1e-4
    @test sampled_layer_fraction(0.6, -0.6, 0.5) ≈ 1.0 / 1.2 atol=1e-4   # crosses wall 0, both layers
    @test sampled_layer_fraction(0.9, 1.1, 1.5) ≈ 2.0 atol=1e-4          # overlap counts twice
end

@testset "trace: brute-force profile average" begin
    @test sampled_profile_average(0.1, 0.3, 0.5) ≈ sampled_layer_fraction(0.1, 0.3, 0.5)      # step: same as before
    @test sampled_profile_average(0.0, 0.5, 0.5; shape=:linear) ≈ 1.0 atol=1e-3              # g averages to 1 across the layer
    @test sampled_profile_average(0.3, 1.7, 2.0; shape=:linear) ≈ 2.0 atol=1e-12             # h = w: g(lower) + g(upper) = 2 everywhere
end
