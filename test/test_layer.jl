@testset "test_layer.jl: segment maths" begin
    Layers = mr.Geometries.Internal.Layers

    @testset "segment_overlap: entry and exit points" begin
        # enters the layer [0, 0.5]: starts at 1.0 (outside), ends at 0.2 (inside)
        (s_in, s_out) = Layers.segment_overlap(1.0, 0.2, 0.0, 0.5)
        @test s_in ≈ 0.625
        @test s_out ≈ 1.0
        # exits the layer
        (s_in, s_out) = Layers.segment_overlap(0.2, 1.0, 0.0, 0.5)
        @test s_in ≈ 0.0
        @test s_out ≈ 0.375
    end

    @testset "segment_fraction" begin
        # neither endpoint inside, the segment passes through [0.25, 0.5]
        @test Layers.segment_fraction(-1.0, 2.0, 0.25, 0.5) ≈ 0.25 / 3
        # no overlap
        @test Layers.segment_fraction(0.6, 1.4, 0.0, 0.5) == 0.0
        # stationary inside / outside
        @test Layers.segment_fraction(0.3, 0.3, 0.0, 0.5) == 1.0
        @test Layers.segment_fraction(0.7, 0.7, 0.0, 0.5) == 0.0
        # stationary exactly on the wall (d = 0) and on the layer edge (d = h): inside
        @test Layers.segment_fraction(0.0, 0.0, 0.0, 0.5) == 1.0
        @test Layers.segment_fraction(0.5, 0.5, 0.0, 0.5) == 1.0
    end

    @testset "side_exponent: step profile" begin
        side = Layers.LayerSide(0.05, 0.5)          # ΔR2(0) = 0.05 / 0.5 = 0.1 /ms
        @test Layers.surface_rate(side) ≈ 0.1
        @test Layers.side_exponent(side, 0.1, 0.3, 0.2) ≈ 0.1 * 0.2         # fully inside
        @test Layers.side_exponent(side, 0.4, 0.6, 1.0) ≈ 0.1 * 1.0 * 0.5   # half inside
        @test Layers.side_exponent(Layers.LayerSide(0.0, 0.0), 0.1, 0.3, 0.2) == 0.0
    end
end
