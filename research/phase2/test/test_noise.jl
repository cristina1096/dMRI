include(joinpath(@__DIR__, "..", "noise.jl"))

@testset "noise: required spins" begin
    @test required_spins(100_000, 20.0) == 100_000          # already ≥ 10 SEM: keep
    @test required_spins(100_000, 10.0) == 100_000
    @test required_spins(100_000, 5.0) == 400_000           # SEM ∝ 1/√N: half the separation → 4× spins
    @test required_spins(10_000, 2.5; target=10) == 160_000
end
