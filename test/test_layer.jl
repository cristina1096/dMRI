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

@testset "test_layer.jl: layer exponent on walls" begin
    Layers = mr.Geometries.Internal.Layers
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry
    v(x) = SVector{3, Float64}(x, 0.3, -1.2)
    both = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5))   # ΔR2(0) = 0.1 /ms on both sides

    @testset "inside, outside and partly inside" begin
        @test Layers.layer_exponent(both, v(0.1), v(0.3), 0.2) ≈ 0.1 * 0.2        # positive side of wall 0
        @test Layers.layer_exponent(both, v(1.7), v(1.9), 0.2) ≈ 0.1 * 0.2        # negative side of wall 2
        @test Layers.layer_exponent(both, v(0.6), v(1.4), 1.0) == 0.0             # middle of the gap
        @test Layers.layer_exponent(both, v(0.4), v(0.6), 1.0) ≈ 0.1 * 1.0 * 0.5
    end

    @testset "neither endpoint inside a layer" begin
        # 0.6 → -0.6 crosses both layers of wall 0: 1.0 of 1.2 um inside
        @test Layers.layer_exponent(both, v(0.6), v(-0.6), 1.2) ≈ 0.1 * 1.2 * (1.0 / 1.2)
    end

    @testset "repeated wall images" begin
        @test Layers.layer_exponent(both, v(10.1), v(10.3), 0.2) ≈ 0.1 * 0.2
        @test Layers.layer_exponent(both, v(-7.9), v(-7.7), 0.2) ≈ 0.1 * 0.2
        # segment longer than the spacing: 0.1 → 4.1 is inside layers for 2.0 of its 4.0 um
        @test Layers.layer_exponent(both, v(0.1), v(4.1), 1.0) ≈ 0.1 * 1.0 * 0.5
    end

    @testset "one side only" begin
        pos = geom(mr.Walls(repeats=2., layer_rho_positive=0.05, layer_h_positive=0.5))
        @test Layers.layer_exponent(pos, v(0.1), v(0.3), 0.2) ≈ 0.1 * 0.2
        @test Layers.layer_exponent(pos, v(1.7), v(1.9), 0.2) == 0.0
    end

    @testset "different on each side" begin
        asym = geom(mr.Walls(repeats=2., layer_rho_positive=0.05, layer_h_positive=0.5,
            layer_rho_negative=0.2, layer_h_negative=0.4))                       # 0.1 /ms and 0.5 /ms
        @test Layers.layer_exponent(asym, v(0.1), v(0.3), 1.0) ≈ 0.1
        @test Layers.layer_exponent(asym, v(1.7), v(1.9), 1.0) ≈ 0.5
    end

    @testset "overlapping layers add (Eq. 10)" begin
        wide = geom(mr.Walls(repeats=2., layer_rho=0.15, layer_h=1.5))           # 0.1 /ms, h > w/2
        @test Layers.layer_exponent(wide, v(0.9), v(1.1), 1.0) ≈ 2 * 0.1         # inside both layers
    end

    @testset "h = w/2 covers every point exactly once" begin
        half = geom(mr.Walls(repeats=2., layer_rho=0.1, layer_h=1.0))            # 0.1 /ms
        for (a, b) in ((0.05, 0.3), (0.8, 1.2), (1.0, 1.9), (0.01, 1.99))
            @test Layers.layer_exponent(half, v(a), v(b), 1.0) ≈ 0.1 rtol=1e-14
        end
    end

    @testset "rotation and shifted position" begin
        rot = geom(mr.Walls(repeats=2., rotation=:y, layer_rho=0.05, layer_h=0.5))
        @test Layers.layer_exponent(rot, SVector(5.0, 0.1, 7.0), SVector(9.0, 0.3, 7.0), 0.2) ≈ 0.1 * 0.2
        shifted = geom(mr.Walls(repeats=2., position=0.5, layer_rho=0.05, layer_h=0.5))
        @test Layers.layer_exponent(shifted, v(0.6), v(0.8), 0.2) ≈ 0.1 * 0.2
    end

    @testset "per-wall values" begin
        two = geom(mr.Walls(position=[0., 1.], layer_rho=[0.05, 0.], layer_h=0.4))  # only wall at 0
        @test Layers.layer_exponent(two, v(0.1), v(0.3), 1.0) ≈ 0.05 / 0.4
        @test Layers.layer_exponent(two, v(0.7), v(0.9), 1.0) == 0.0
    end

    @testset "no layer" begin
        none = geom(mr.Walls(repeats=2.))
        @test !Layers.has_layer(none)
        @test Layers.layer_exponent(none, v(0.1), v(0.3), 0.2) == 0.0
        @test Layers.has_layer(both)
    end
end
