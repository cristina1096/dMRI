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

@testset "test_layer.jl: layer applied during the simulation" begin
    per_spin_transverse(snap) = [s.orientations[1].transverse for s in snap.spins]

    @testset "crossing with reflection in one step (P2.2.4)" begin
        # wall at 0, layer [0, 0.5] on the positive side, ΔR2(0) = 0.1 /ms
        walls = mr.Walls(position=0., layer_rho_positive=0.05, layer_h_positive=0.5)
        seq = build_sequence([1.2, :readout])
        sim = mr.Simulation(seq; geometry=walls, diffusivity=3., verbose=false)
        # the first part is the zero-length instant at t = 0; take the 1.2 ms step after it
        part = first(p for p in mr.parts([seq], 0., mr.TimeStep(1.2, Inf)) if p.duration > 0)
        @test part.duration ≈ 1.2
        spin = mr.Spin(position=[0.8, 0., 0.], transverse=1., longitudinal=0.)
        # proposed move 0.8 → -0.4 reflects at 0 and ends at 0.4.
        # in layer: last 0.5 of the 0.8 um before the wall (time 0.5) + all 0.4 um after (time 0.4) = 0.9 of 1.2 ms
        mr.Evolve.draw_step!(spin, sim, part, [3.], [-0.4, 0., 0.])
        @test spin.position[1] ≈ 0.4 atol=1e-12
        @test spin.orientations[1].transverse ≈ exp(-0.1 * 0.9) rtol=1e-12
        @test spin.orientations[1].phase == 0.0
    end

    @testset "no layer: bit-identical to the unmodified code (P2.2.8)" begin
        function run(walls)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., R2=1/80, verbose=false)
            Random.seed!(7)
            snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [50.]; return_snapshot=true)[1]
        end
        a = run(mr.Walls(repeats=2.))
        b = run(mr.Walls(repeats=2., layer_h=0.5))          # ρ = 0: no layer
        @test mr.position.(a) == mr.position.(b)
        @test per_spin_transverse(a) == per_spin_transverse(b)
    end

    @testset "h = w/2: every spin decays at exactly ΔR2(0)" begin
        walls = mr.Walls(repeats=2., layer_rho=0.1, layer_h=1.0)                # ΔR2(0) = 0.1 /ms
        sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., verbose=false)
        Random.seed!(1)
        snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
        res = mr.readout(snap, sim, [10., 50.]; return_snapshot=true)
        for (t, s) in zip((10., 50.), res)
            @test all(isapprox.(per_spin_transverse(s), exp(-0.1 * t); rtol=1e-10))
        end
    end

    @testset "bulk R2 and the layer multiply exactly (R2_total = R2_bulk + ΔR2)" begin
        walls = mr.Walls(repeats=2., layer_rho=0.03, layer_h=0.3)
        function run(R2)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., R2=R2, verbose=false)
            Random.seed!(3)
            snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [50.]; return_snapshot=true)[1]
        end
        layer_only = run(0.)
        with_bulk = run(1/80)
        @test mr.position.(layer_only) == mr.position.(with_bulk)
        @test all(isapprox.(per_spin_transverse(with_bulk), per_spin_transverse(layer_only) .* exp(-50 / 80); rtol=1e-10))
        @test minimum(per_spin_transverse(layer_only)) < 1.0                     # the layer did act
    end

    @testset "positive-only and negative-only layers are mirror images" begin
        function mean_signal(walls)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., verbose=false)
            Random.seed!(11)
            snap = mr.Snapshot(20000, sim, 500; transverse=1., longitudinal=0.)
            m = per_spin_transverse(mr.readout(snap, sim, [50.]; return_snapshot=true)[1])
            return (mean(m), std(m) / sqrt(length(m)))
        end
        (p, ep) = mean_signal(mr.Walls(repeats=2., layer_rho_positive=0.05, layer_h_positive=0.5))
        (n, en) = mean_signal(mr.Walls(repeats=2., layer_rho_negative=0.05, layer_h_negative=0.5))
        @test abs(p - n) < 5 * sqrt(ep^2 + en^2)
        @test p < 0.99                                                          # the layer did act
    end

    @testset "finite RF pulse with a layer is rejected" begin
        seq = mr.SequenceParts.SequenceWaveform(
            (([], []), ([], []), ([], [])),
            [(0., 9., [mr.SequenceParts.ConstantPulse(0.25 / 9, 0., 0.)])],
            [], [10.], 10.,
        )
        sim = mr.Simulation(seq; geometry=mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5), verbose=false)
        @test_throws Exception mr.readout(100, sim, [10.])
    end
end

@testset "test_layer.jl: layer timestep constraint" begin
    Layers = mr.Geometries.Internal.Layers
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry
    @test Layers.max_layer_rate(geom(mr.Walls(repeats=2.))) == 0.0
    @test Layers.max_layer_rate(geom(mr.Walls(repeats=2., layer_h=0.5))) == 0.0             # ρ = 0 → no layer
    @test Layers.max_layer_rate(geom(mr.Walls(repeats=2., layer_rho_positive=0.01, layer_h_positive=0.1,
        layer_rho_negative=0.05, layer_h_negative=0.5))) ≈ 0.1                                # max(0.01/0.1, 0.05/0.5)
    # constraint binds for a strong layer (τ ≤ c/ΔR2(0)), is absent without a layer
    s = mr.Simulation([]; geometry=mr.Walls(repeats=2., layer_rho=0.1, layer_h=0.1), diffusivity=3., verbose=false, timestep=(layer=0.005,))
    @test s.timestep.max_timestep ≈ 0.005 / 1.0
    s0 = mr.Simulation([]; geometry=mr.Walls(repeats=2.), diffusivity=3., verbose=false, timestep=(layer=0.005,))
    @test s0.timestep.max_timestep ≈ 0.03 * 2.0^2 / 3                                       # tortuosity, unchanged
    sd = mr.Simulation([]; geometry=mr.Walls(repeats=2., layer_rho=0.1, layer_h=0.1), diffusivity=3., verbose=false)
    @test sd.timestep.max_timestep ≈ 0.005                                                   # default c = 0.005
    # runtests.jl calls Logging.disable_logging(Logging.Info) globally; lift it for this one check and restore it
    Logging.disable_logging(Logging.Debug)
    try
        @test_logs (:info, r"near-surface layer") match_mode=:any mr.Simulation([]; geometry=mr.Walls(repeats=2., layer_rho=0.1, layer_h=0.1),
            diffusivity=3., verbose=true, timestep=(layer=0.005,))
    finally
        Logging.disable_logging(Logging.Info)
    end
end

@testset "test_layer.jl: linear profile maths" begin
    Layers = mr.Geometries.Internal.Layers

    @testset "profile values and normalisation (P2.4.2)" begin
        @test Layers.profile_g(:linear, 0.0) == 2.0
        @test Layers.profile_g(:linear, 0.5) == 1.0
        @test Layers.profile_g(:linear, 1.0) == 0.0
        @test Layers.profile_g(:linear, 1.5) == 0.0
        @test Layers.profile_g(:step, 0.3) == 1.0
        @test Layers.profile_g(:step, 1.2) == 0.0
        n = 100_000
        for shape in (:step, :linear)
            @test sum(Layers.profile_g(shape, (k - 0.5) / n) for k in 1:n) / n ≈ 1.0 atol=1e-10
        end
    end

    @testset "step is unchanged" begin
        side = Layers.LayerSide(0.05, 0.5)
        @test side.shape == :step
        @test Layers.side_exponent(side, 0.1, 0.3, 0.2) ≈ 0.1 * 0.2
        @test Layers.side_exponent(side, 0.4, 0.6, 1.0) ≈ 0.1 * 1.0 * 0.5
    end

    @testset "linear: surface rate and edges" begin
        side = Layers.LayerSide(0.05, 0.5, :linear)            # ρ/h = 0.1, ΔR2(0) = 0.2
        @test Layers.surface_rate(side) ≈ 0.2
        @test Layers.side_exponent(side, 0.0, 0.0, 1.0) ≈ 0.2  # sitting on the wall
        @test Layers.side_exponent(side, 0.5, 0.5, 1.0) ≈ 0.0  # sitting at d = h
        @test Layers.side_exponent(side, 0.7, 0.7, 1.0) == 0.0 # outside
    end

    @testset "linear: pieces" begin
        side = Layers.LayerSide(0.05, 0.5, :linear)            # ρ/h = 0.1
        # whole layer crossed, 0 → h: average g = 1 → exponent = (ρ/h)·dt
        @test Layers.side_exponent(side, 0.0, 0.5, 1.0) ≈ 0.1
        # 0.8 → 0.2: inside for the second half (d 0.5 → 0.2, mean d = 0.35, g = 2(1 − 0.7) = 0.6)
        @test Layers.side_exponent(side, 0.8, 0.2, 1.0) ≈ 0.1 * 0.6 * 0.5
        # matches brute-force sampling on random pieces
        Random.seed!(5)
        for _ in 1:200
            d0, d1, h = 1.2 * rand() - 0.1, 1.2 * rand() - 0.1, 0.1 + 0.9 * rand()
            s = Layers.LayerSide(0.3 * h, h, :linear)
            n = 200_000
            brute = sum(Layers.profile_g(:linear, (d0 + (d1 - d0) * (k - 0.5) / n) / h) for k in 1:n) / n * (s.rho / h)
            @test Layers.side_exponent(s, d0, d1, 1.0) ≈ brute atol=2e-5
        end
    end

    @testset "timestep constraint uses the profile's ΔR2(0)" begin
        @test Layers.max_layer_rate([Layers.WallLayer(0.0, Layers.LayerSide(0.01, 0.1, :linear), Layers.LayerSide(0.0, 0.0))]) ≈ 0.2
    end
end

@testset "test_layer.jl: layer_shape on Walls" begin
    Layers = mr.Geometries.Internal.Layers
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry
    v(x) = SVector{3, Float64}(x, 0.3, -1.2)

    @testset "resolution" begin
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5))[1].layer[1]
        @test (l.positive.shape, l.negative.shape) == (:step, :step)
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape="linear"))[1].layer[1]
        @test (l.positive.shape, l.negative.shape) == (:linear, :linear)
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape_negative="linear"))[1].layer[1]
        @test (l.positive.shape, l.negative.shape) == (:step, :linear)
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape="gaussian"))
    end

    @testset "mixed shapes use their own formula" begin
        # positive side step, negative side linear; ρ/h = 0.1
        g = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape_negative="linear"))
        @test Layers.layer_exponent(g, v(0.1), v(0.1), 1.0) ≈ 0.1            # positive (step) side of wall 0: g = 1
        @test Layers.layer_exponent(g, v(1.9), v(1.9), 1.0) ≈ 0.16           # negative (linear) side of wall 2, d = 0.1: g = 1.6
        @test Layers.layer_exponent(g, v(0.0), v(0.0), 1.0) ≈ 0.1 + 0.2      # on wall 0: inside both of its layers (step 1 + linear 2)
    end

    @testset "linear, h = w: uniform rate" begin
        g = geom(mr.Walls(repeats=2., layer_rho=0.1, layer_h=2.0, layer_shape="linear"))   # 2ρ/w = 0.1
        for (a, b) in ((0.05, 0.3), (0.8, 1.2), (1.0, 1.9), (0.01, 1.99), (0.3, 0.3))
            @test Layers.layer_exponent(g, v(a), v(b), 1.0) ≈ 0.1 rtol=1e-12
        end
    end

    @testset "linear, h = w: every spin decays at exactly 0.1 /ms" begin
        walls = mr.Walls(repeats=2., layer_rho=0.1, layer_h=2.0, layer_shape="linear")
        sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., verbose=false)
        Random.seed!(1)
        snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
        res = mr.readout(snap, sim, [10., 50.]; return_snapshot=true)
        for (t, s) in zip((10., 50.), res)
            @test all(isapprox.([x.orientations[1].transverse for x in s.spins], exp(-0.1 * t); rtol=1e-10))
        end
    end

    @testset "JSON round trip keeps layer_shape" begin
        io = IOBuffer()
        mr.write_geometry(io, mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape="linear"))
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_shape.value == "linear"
        @test isnothing(back.layer_shape_positive.value)
    end
end

@testset "test_layer.jl: exponential profile maths" begin
    Layers = mr.Geometries.Internal.Layers

    @testset "values and normalisation" begin
        c = 4.0
        N = 1 / (1 - exp(-c))
        @test Layers.profile_g(:exponential, 0.0, c) ≈ N
        @test Layers.profile_g(:exponential, 1.0, c) ≈ N * exp(-1)
        @test Layers.profile_g(:exponential, c + 1e-9, c) == 0.0
        @test Layers.profile_g(:exponential, -0.1, c) == 0.0
        @test Layers.profile_g(:exponential, 2.0) ≈ exp(-2.0)              # c = Inf: untruncated, ∫ = 1
        n = 400_000
        @test sum(Layers.profile_g(:exponential, (k - 0.5) / n * c, c) for k in 1:n) * c / n ≈ 1.0 atol=1e-9
        @test Layers.profile_g(:step, 0.5, 3.0) == 1.0                     # c ignored for step and linear
        @test Layers.profile_g(:linear, 0.5, 3.0) == 1.0
        @test Layers.SHAPES == (:step, :linear, :exponential)
    end

    @testset "constructor, cutoff and surface rate" begin
        @test_throws ErrorException Layers.LayerSide(0.1, 0.2, :exponential)   # needs a cutoff
        s = Layers.LayerSide(0.1, 0.2, :exponential, 2.0)
        @test s.cutoff == 2.0
        @test Layers.surface_rate(s) ≈ 0.1 / 0.2 / (1 - exp(-10))
        @test Layers.LayerSide(0.05, 0.5).cutoff == 0.5
        @test Layers.LayerSide(0.05, 0.5, :linear).cutoff == 0.5
        @test Layers.max_layer_rate([Layers.WallLayer(0.0, s, Layers.LayerSide(0.0, 0.0))]) ≈ 0.5 / (1 - exp(-10))
    end

    @testset "pieces vs brute force" begin
        Random.seed!(9)
        for _ in 1:200
            h = 0.05 + 0.5 * rand()
            cut = h + (2.0 - h) * rand()
            d0, d1 = 2.4 * rand() - 0.2, 2.4 * rand() - 0.2
            s = Layers.LayerSide(0.3 * h, h, :exponential, cut)
            n = 200_000
            brute = sum(Layers.profile_g(:exponential, (d0 + (d1 - d0) * (k - 0.5) / n) / h, cut / h) for k in 1:n) / n * (s.rho / h)
            @test Layers.side_exponent(s, d0, d1, 1.0) ≈ brute rtol=1e-4 atol=5e-6
        end
        s = Layers.LayerSide(0.02, 0.2, :exponential, 2.0)                 # ρ/h = 0.1
        r0 = Layers.surface_rate(s)
        @test Layers.side_exponent(s, 0.3, 0.3, 0.5) ≈ 0.5 * r0 * exp(-1.5)             # stationary
        @test Layers.side_exponent(s, 0.3, 0.3 + 1e-12, 0.5) ≈ 0.5 * r0 * exp(-1.5)     # tiny piece: series branch
        @test Layers.side_exponent(s, 0.0, 2.0, 1.0) ≈ 0.02 / 2.0 rtol=1e-12             # whole support: ρ/cutoff
        @test Layers.side_exponent(s, 2.5, 2.6, 1.0) == 0.0                              # beyond the cutoff
        @test Layers.side_exponent(s, -0.5, -0.1, 1.0) == 0.0                            # other side of the wall
    end
end

@testset "test_layer.jl: exponential on Walls" begin
    Layers = mr.Geometries.Internal.Layers
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry
    v(x) = SVector{3, Float64}(x, 0.3, -1.2)

    @testset "resolution and validation" begin
        l = geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))[1].layer[1]
        @test (l.positive.shape, l.positive.cutoff, l.negative.cutoff) == (:exponential, 2.0, 2.0)
        l = geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=1.0))[1].layer[1]
        @test l.positive.cutoff == 1.0
        l = geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape_negative="exponential", layer_cutoff=1.0))[1].layer[1]
        @test (l.positive.shape, l.positive.cutoff, l.negative.shape, l.negative.cutoff) == (:step, 0.2, :exponential, 1.0)
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=2.5))
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=0.0))
        @test_throws ErrorException geom(mr.Walls(layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))   # single wall: no spacing
        @test geom(mr.Walls(layer_shape="exponential"))[1].layer === nothing                                  # ρ = 0: inert, no error
        l = geom(mr.Walls(layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=3.0))[1].layer[1]
        @test l.positive.cutoff == 3.0
    end

    @testset "reach includes the opposite wall" begin
        g = geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))   # ρ/h = 0.1, cutoff 2
        N = 1 / (1 - exp(-10))
        # stationary at x = 0.5: d = 0.5 from wall 0 (positive side), 1.5 from wall 2 (negative side)
        @test Layers.layer_exponent(g, v(0.5), v(0.5), 1.0) ≈ 0.1 * N * (exp(-2.5) + exp(-7.5)) rtol=1e-12
    end

    @testset "JSON round trip keeps layer_cutoff" begin
        io = IOBuffer()
        mr.write_geometry(io, mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=1.5))
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_shape.value == "exponential"
        @test back.layer_cutoff.value == 1.5
    end
end

@testset "test_layer.jl: exponential thin layer, long piece toward the wall (no NaN)" begin
    Layers = mr.Geometries.Internal.Layers
    s = Layers.LayerSide(1e-4, 0.001, :exponential, 2.0)              # h = 1 nm
    a = Layers.side_exponent(s, 0.8, 0.0, 1.0)                         # moving toward the wall
    b = Layers.side_exponent(s, 0.0, 0.8, 1.0)                         # moving away
    @test isfinite(a)
    @test a ≈ b rtol=1e-12                                             # the integral does not depend on direction
    @test a ≈ s.rho / 0.8 rtol=1e-9                                    # whole layer crossed: ρ/(piece length), tail e^(-800) negligible
end
