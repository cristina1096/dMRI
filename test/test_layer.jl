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

@testset "test_layer.jl: wall layer parameters" begin
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry

    @testset "no layer gives nothing" begin
        @test geom(mr.Walls(repeats=2.))[1].layer === nothing
        @test geom(mr.Walls(repeats=2., layer_h=0.5))[1].layer === nothing
    end

    @testset "same layer on both sides" begin
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5))[1].layer[1]
        @test l.position == 0.0
        @test (l.positive.rho, l.positive.h) == (0.05, 0.5)
        @test (l.negative.rho, l.negative.h) == (0.05, 0.5)
    end

    @testset "one side only" begin
        l = geom(mr.Walls(repeats=2., layer_rho_positive=0.05, layer_h_positive=0.5))[1].layer[1]
        @test (l.positive.rho, l.positive.h) == (0.05, 0.5)
        @test l.negative.rho == 0.0
    end

    @testset "side-specific values override both-sides values" begin
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5,
            layer_rho_negative=0.2, layer_h_negative=0.4))[1].layer[1]
        @test (l.positive.rho, l.positive.h) == (0.05, 0.5)
        @test (l.negative.rho, l.negative.h) == (0.2, 0.4)
    end

    @testset "per-wall values reach the right wall" begin
        layers = geom(mr.Walls(position=[0., 1.], layer_rho=[0.05, 0.], layer_h=0.4))[1].layer
        @test [l.position for l in layers] == [0.0, 1.0]
        @test layers[1].positive.rho == 0.05
        @test layers[2].positive.rho == 0.0
    end

    @testset "validation" begin
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=-0.1, layer_h=0.5))
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.1, layer_h=-0.5))
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.1))              # h = 0
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.1, layer_h=2.5)) # h > spacing 2
        @test_throws ErrorException geom(mr.Walls(position=[0., 1.], layer_rho=0.1, layer_h=1.5))
    end

    @testset "warning when spins can stick to the wall" begin
        @test_logs (:warn, r"not applied while a spin is stuck") geom(
            mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, density=1., dwell_time=1.))
    end

    @testset "JSON round trip keeps the layer fields" begin
        walls = mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_rho_negative=0.)
        io = IOBuffer()
        mr.write_geometry(io, walls)
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_rho.value == 0.05
        @test back.layer_h.value == 0.5
        @test back.layer_rho_negative.value == 0.0
        @test isnothing(back.layer_h_positive.value)
    end
end
