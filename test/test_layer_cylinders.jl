@testset "test_layer_cylinders.jl: shell maths" begin
    Layers = mr.Geometries.Internal.Layers
    P(x, y) = SVector{2, Float64}(x, y)
    O = P(0., 0.)

    @testset "normalisation κ and surface rate (D3.1)" begin
        side = Layers.LayerSide(0.02, 0.2)
        @test Layers.centroid(side) == 0.1
        @test Layers.kappa(side, 1.0, true) ≈ 0.9
        @test Layers.kappa(side, 1.0, false) ≈ 1.1
        @test Layers.cylinder_surface_rate(side, 1.0, true) ≈ 0.02 / (0.2 * 0.9)
        @test Layers.cylinder_surface_rate(side, 1.0, false) ≈ 0.02 / (0.2 * 1.1)
        @test Layers.cylinder_surface_rate(Layers.LayerSide(0.0, 0.0), 1.0, true) == 0.0
        @test_throws ErrorException Layers.centroid(Layers.LayerSide(0.02, 0.2, :linear))
    end

    @testset "disk_length" begin
        @test Layers.disk_length(1.0, 0.0, 0.0, 0.25) ≈ 0.5          # r² = s², s ≤ 0.5
        @test Layers.disk_length(1.0, 0.0, 2.0, 1.0) == 0.0           # never inside
        @test Layers.disk_length(0.0, 0.0, 0.5, 1.0) == 1.0           # stationary inside
        @test Layers.disk_length(0.0, 0.0, 1.5, 1.0) == 0.0           # stationary outside
    end

    @testset "shell_fraction: exact cases" begin
        # radial through the centre, shell 0.8 ≤ r ≤ 1: |x| ∈ [0.8, 1] → 0.4 of a length-2 piece
        @test Layers.shell_fraction(P(-1., 0.), P(1., 0.), O, 0.8, 1.0) ≈ 0.2
        # chord at y = 0.5 from x = −2 to 2: inside r ≤ 1 for |x| ≤ √0.75, inside r ≤ 0.8 for |x| ≤ √0.39
        @test Layers.shell_fraction(P(-2., 0.5), P(2., 0.5), O, 0.8, 1.0) ≈ (2sqrt(0.75) - 2sqrt(0.39)) / 4
        # tangent to the outer circle r = 1: measure zero
        @test Layers.shell_fraction(P(-2., 1.0), P(2., 1.0), O, 0.8, 1.0) == 0.0
        # tangent to the inner circle r = 0.8: the whole chord |x| ≤ 0.6 is in the shell
        @test Layers.shell_fraction(P(-2., 0.8), P(2., 0.8), O, 0.8, 1.0) ≈ 0.3
        # outer shell 1 ≤ r ≤ 1.3, radial piece 1 → 2
        @test Layers.shell_fraction(P(1., 0.), P(2., 0.), O, 1.0, 1.3) ≈ 0.3
        # r_lo = 0 (shell = whole disk, h = R)
        @test Layers.shell_fraction(P(-0.5, 0.), P(0.5, 0.), O, 0.0, 1.0) == 1.0
        # stationary
        @test Layers.shell_fraction(P(0.9, 0.), P(0.9, 0.), O, 0.8, 1.0) == 1.0
        @test Layers.shell_fraction(P(0.5, 0.), P(0.5, 0.), O, 0.8, 1.0) == 0.0
        # shifted centre
        @test Layers.shell_fraction(P(4., 3.), P(6., 3.), P(5., 3.), 0.8, 1.0) ≈ 0.2
    end

    @testset "Review Focus 1: nearly axial motion equals stationary" begin
        f = Layers.shell_fraction(P(0.9, 0.), P(0.9, 1e-15), O, 0.8, 1.0)
        @test isfinite(f)
        @test f ≈ 1.0
        @test Layers.shell_fraction(P(0.5, 0.), P(0.5 + 1e-15, 0.), O, 0.8, 1.0) == 0.0
    end

    @testset "Review Focus 2: piece starting exactly on the surface" begin
        # after a reflection at r = R = 1, moving inward to r = 0.7: inner shell [0.8, 1] → 0.2/0.3
        @test Layers.shell_fraction(P(1., 0.), P(0.7, 0.), O, 0.8, 1.0) ≈ 2 / 3
        @test Layers.shell_fraction(P(1., 0.), P(0.7, 0.), O, 1.0, 1.3) ≈ 0.0 atol=1e-12
        # moving outward from the surface to r = 1.5: outer shell [1, 1.3] → 0.3/0.5
        θ = 0.7
        u = P(cos(θ), sin(θ))
        @test Layers.shell_fraction(u, 1.5 .* u, O, 1.0, 1.3) ≈ 0.6
        @test Layers.shell_fraction(u, 1.5 .* u, O, 0.8, 1.0) ≈ 0.0 atol=1e-12
    end

    @testset "shell_fraction vs brute-force sampling" begin
        Random.seed!(2026)
        K = 100_000
        maxerr = 0.0
        for _ in 1:1000
            p0 = P((4 .* rand(2) .- 2)...)
            p1 = P((4 .* rand(2) .- 2)...)
            for (lo, hi) in ((0.8, 1.0), (1.0, 1.3), (0.0, 1.0))
                brute = count(k -> lo <= norm(p0 + (k - 0.5) / K * (p1 - p0)) <= hi, 1:K) / K
                maxerr = max(maxerr, abs(Layers.shell_fraction(p0, p1, O, lo, hi) - brute))
            end
        end
        @test maxerr < 3e-5          # ≤ 4 boundary crossings, each ≤ 1/(2K)
    end

    @testset "cylinder_exponent: both sides" begin
        layer = Layers.CylinderLayer(O, 1.0, Layers.LayerSide(0.02, 0.2), Layers.LayerSide(0.02, 0.3))
        rin = 0.02 / (0.2 * 0.9)
        rout = 0.02 / (0.3 * 1.15)
        # radial piece r = 0.5 → 1.5 of duration 2 ms: 0.2 in the inner shell, 0.3 in the outer shell (of 1.0)
        @test Layers.cylinder_exponent(layer, O, P(0.5, 0.), P(1.5, 0.), 2.0) ≈ 2.0 * (0.2 * rin + 0.3 * rout)
        # inner side off
        layer_out = Layers.CylinderLayer(O, 1.0, Layers.LayerSide(0.0, 0.0), Layers.LayerSide(0.02, 0.3))
        @test Layers.cylinder_exponent(layer_out, O, P(0.5, 0.), P(0.9, 0.), 1.0) == 0.0
    end
end

@testset "test_layer_cylinders.jl: cylinder layer parameters" begin
    geom(c) = mr.Simulation([]; geometry=c, verbose=false).geometry

    @testset "no layer gives nothing" begin
        @test geom(mr.Cylinders(radius=1.))[1].layer === nothing
        @test geom(mr.Cylinders(radius=1., layer_h=0.2))[1].layer === nothing       # ρ = 0
    end

    @testset "same layer on both sides" begin
        l = geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2))[1].layer[1]
        @test l.centre == SVector(0., 0.)
        @test l.radius == 1.0
        @test (l.inner.rho, l.inner.h, l.inner.shape) == (0.02, 0.2, :step)
        @test (l.outer.rho, l.outer.h, l.outer.shape) == (0.02, 0.2, :step)
    end

    @testset "side-specific values override both-sides values" begin
        l = geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_rho_outer=0.05, layer_h_inner=0.4))[1].layer[1]
        @test (l.inner.rho, l.inner.h) == (0.02, 0.4)
        @test (l.outer.rho, l.outer.h) == (0.05, 0.2)
        l = geom(mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2))[1].layer[1]
        @test l.outer.rho == 0.0
    end

    @testset "per-cylinder values and centres" begin
        layers = geom(mr.Cylinders(radius=[1., 0.5], position=[[0., 0.], [5., 0.]], layer_rho=[0.02, 0.], layer_h=0.2))[1].layer
        @test [l.centre for l in layers] == [SVector(0., 0.), SVector(5., 0.)]
        @test [l.radius for l in layers] == [1.0, 0.5]
        @test layers[2].inner.rho == 0.0
    end

    @testset "centres are wrapped into the unit cell as the obstructions are" begin
        l = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], position=[2.0, 0.0], layer_rho=0.02, layer_h=0.2))[1].layer[1]
        @test l.centre ≈ SVector(-0.5, 0.0)
    end

    @testset "h_inner = R is allowed (uniform inner rate)" begin
        l = geom(mr.Cylinders(radius=1., layer_rho_inner=0.05, layer_h_inner=1.0))[1].layer[1]
        @test l.inner.h == 1.0
    end

    @testset "validation" begin
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=-0.1, layer_h=0.2))
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=-0.2))
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1))                              # h = 0
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=0.1, layer_h_inner=1.5))  # h_inner > R
        @test_throws ErrorException geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.1, layer_h=0.6))  # h_outer > gap 0.5
        @test_throws ErrorException geom(mr.Cylinders(radius=1., position=[[0., 0.], [2.4, 0.]], layer_rho=0.1, layer_h=0.5))  # h 0.5 > gap 0.4
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=0.2, layer_shape="linear"))   # step only
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=0.2, layer_shape="bogus"))
        # h_outer = gap exactly is allowed; a side with ρ = 0 is not checked against the gap
        @test geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.1, layer_h=0.5))[1].layer !== nothing
        @test geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho_inner=0.1, layer_h_inner=0.2, layer_h_outer=0.9))[1].layer !== nothing
    end

    @testset "warning when spins can stick to the surface" begin
        @test_logs (:warn, r"not applied while a spin is stuck") geom(
            mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, density=1., dwell_time=1.))
    end

    @testset "JSON round trip keeps the layer fields" begin
        c = mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_h_outer=0.3)
        io = IOBuffer()
        mr.write_geometry(io, c)
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_rho.value == 0.02
        @test back.layer_h.value == 0.2
        @test back.layer_h_outer.value == 0.3
        @test isnothing(back.layer_rho_inner.value)
        @test back.layer_shape.value == "step"
    end
end

@testset "test_layer_cylinders.jl: layer exponent on cylinder groups" begin
    Layers = mr.Geometries.Internal.Layers
    geom(c) = mr.Simulation([]; geometry=c, verbose=false).geometry
    V3(x, y, z=0.) = SVector{3, Float64}(x, y, z)
    rin(ρ, h, R) = ρ / (h * (1 - h / (2R)))           # D3.1, written independently of the package
    rout(ρ, h, R) = ρ / (h * (1 + h / (2R)))

    "Brute-force summed rate at (x, y) over cylinders (centres, radii), images within ±2 cells if repeating."
    function brute_rate(x, y, centres, radii, ρ, h_in, h_out; repeats=nothing)
        cells = isnothing(repeats) ? [(0., 0.)] : [(i * repeats[1], j * repeats[2]) for i in -2:2 for j in -2:2]
        total = 0.0
        for (c, R) in zip(centres, radii), (ox, oy) in cells
            r = hypot(x - c[1] - ox, y - c[2] - oy)
            (R - h_in <= r <= R) && (total += rin(ρ, h_in, R))
            (R <= r <= R + h_out) && (total += rout(ρ, h_out, R))
        end
        return total
    end
    function brute_exponent(p0, p1, dt, args...; kw...)
        K = 20_000
        return dt * sum(brute_rate((p0 + (k - 0.5) / K * (p1 - p0))[1:2]..., args...; kw...) for k in 1:K) / K
    end
    function max_error(g, centres, radii, ρ, h_in, h_out, box; repeats=nothing, n=200)
        Random.seed!(77)
        err = 0.0
        for _ in 1:n
            p0 = V3((box .* (2 .* rand(2) .- 1))..., rand())
            p1 = V3((box .* (2 .* rand(2) .- 1))..., rand())
            err = max(err, abs(Layers.layer_exponent(g, p0, p1, 1.0) - brute_exponent(p0, p1, 1.0, centres, radii, ρ, h_in, h_out; repeats=repeats)))
        end
        return err
    end

    @testset "single cylinder: stationary points (rate) and a radial piece" begin
        g = geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h_inner=0.2, layer_h_outer=0.3))
        @test Layers.layer_exponent(g, V3(0.9, 0.), V3(0.9, 0.), 1.0) ≈ rin(0.02, 0.2, 1.0)
        @test Layers.layer_exponent(g, V3(1.2, 0.), V3(1.2, 0.), 1.0) ≈ rout(0.02, 0.3, 1.0)
        @test Layers.layer_exponent(g, V3(0.5, 0.), V3(0.5, 0.), 1.0) == 0.0
        @test Layers.layer_exponent(g, V3(1.5, 0.), V3(1.5, 0.), 1.0) == 0.0
        # axial motion does not change r
        @test Layers.layer_exponent(g, V3(0.9, 0., 0.), V3(0.9, 0., 5.), 2.0) ≈ 2 * rin(0.02, 0.2, 1.0)
    end

    @testset "Review Focus 3: rotated axis (cylinder along x)" begin
        g = geom(mr.Cylinders(radius=1., rotation=:x, layer_rho_inner=0.02, layer_h_inner=0.2))
        @test Layers.layer_exponent(g, V3(5., 0.9, 0.), V3(5., 0.9, 0.), 1.0) ≈ rin(0.02, 0.2, 1.0)
        @test Layers.layer_exponent(g, V3(0., 0., 0.9), V3(0., 0., 0.9), 1.0) ≈ rin(0.02, 0.2, 1.0)
        @test Layers.layer_exponent(g, V3(0.9, 0., 0.), V3(0.9, 0., 0.), 1.0) == 0.0     # on the axis
        @test Layers.layer_exponent(g, V3(-3., 0.9, 0.), V3(3., 0.9, 0.), 1.0) ≈ rin(0.02, 0.2, 1.0)  # along the axis
    end

    @testset "two overlapping outer layers (non-repeating) vs brute force" begin
        centres, radii = [(0., 0.), (2.2, 0.)], [1.0, 0.8]          # gap 0.4, outer h 0.3 → overlap 0.2 wide
        g = geom(mr.Cylinders(radius=radii, position=[[0., 0.], [2.2, 0.]], layer_rho=0.03, layer_h_inner=0.2, layer_h_outer=0.3))
        @test Layers.layer_exponent(g, V3(1.15, 0.), V3(1.15, 0.), 1.0) ≈ rout(0.03, 0.3, 1.0) + rout(0.03, 0.3, 0.8)
        @test max_error(g, centres, radii, 0.03, 0.2, 0.3, (2.5, 1.5)) < 1e-4
    end

    @testset "four overlapping cylinders (non-repeating) vs brute force" begin
        centres = [(1.2, 1.2), (-1.2, 1.2), (1.2, -1.2), (-1.2, -1.2)]
        g = geom(mr.Cylinders(radius=1., position=[[c...] for c in centres], layer_rho=0.03, layer_h=0.3))
        @test max_error(g, centres, fill(1.0, 4), 0.03, 0.3, 0.3, (2.5, 2.5)) < 1e-4
    end

    @testset "square array (images) vs brute force" begin
        g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.03, layer_h_inner=0.2, layer_h_outer=0.3))
        @test max_error(g, [(0., 0.)], [1.0], 0.03, 0.2, 0.3, (4.0, 4.0); repeats=(2.5, 2.5)) < 1e-4
    end

    @testset "Review Focus 4: two cylinders per cell, rectangular array vs brute force" begin
        centres = [(0., 0.), (1.5, 1.25)]
        g = geom(mr.Cylinders(radius=[0.6, 0.5], position=[[c...] for c in centres], repeats=[3.0, 2.5],
            layer_rho=0.03, layer_h_inner=0.2, layer_h_outer=0.25))
        @test max_error(g, centres, [0.6, 0.5], 0.03, 0.2, 0.25, (4.0, 4.0); repeats=(3.0, 2.5)) < 1e-4
    end

    @testset "timestep: max_layer_rate is the summed-rate bound" begin
        # separated (gap 0.5 > 2h = 0.4): the largest single rate, inside
        g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2))
        @test Layers.max_layer_rate(g) ≈ rin(0.02, 0.2, 1.0)
        # overlapping outer layers (2h = 0.6 > gap 0.5): own + 4 nearest neighbours (diagonals at 3.54 > 2.6 do not reach)
        g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.1, layer_h=0.3))
        @test Layers.max_layer_rate(g) ≈ 5 * rout(0.1, 0.3, 1.0)
        s = mr.Simulation([]; geometry=mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.1, layer_h=0.3), diffusivity=3., verbose=false)
        @test s.timestep.max_timestep ≈ 0.005 / (5 * rout(0.1, 0.3, 1.0))
        # single cylinder
        @test Layers.max_layer_rate(geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h_inner=0.2, layer_h_outer=0.1))) ≈ max(rin(0.02, 0.2, 1.0), rout(0.02, 0.1, 1.0))
        @test Layers.max_layer_rate(geom(mr.Cylinders(radius=1.))) == 0.0
    end
end

@testset "test_layer_cylinders.jl: layer applied during the simulation" begin
    per_spin(snap) = [s.orientations[1].transverse for s in snap.spins]
    inside_mask(snap; R=1.) = [hypot(s.position[1], s.position[2]) < R for s in snap.spins]
    rin(ρ, h, R) = ρ / (h * (1 - h / (2R)))
    rout(ρ, h, R) = ρ / (h * (1 + h / (2R)))
    box = mr.BoundingBox([-1.5, -1.5, -1.], [1.5, 1.5, 1.])

    @testset "radial step with a reflection, inside and outside" begin
        cyl = mr.Cylinders(radius=1., layer_rho=0.02, layer_h_inner=0.2, layer_h_outer=0.3)
        seq = build_sequence([1.2, :readout])
        sim = mr.Simulation(seq; geometry=cyl, diffusivity=3., verbose=false)
        part = first(p for p in mr.parts([seq], 0., mr.TimeStep(1.2, Inf)) if p.duration > 0)
        # inside: 0.5 → 1.3 reflects at 1 and ends at 0.7; path 0.8 um in 1.2 ms; in [0.8, 1]: 0.2 + 0.2 um → 0.6 ms
        spin = mr.Spin(position=[0.5, 0., 0.], transverse=1., longitudinal=0.)
        mr.Evolve.draw_step!(spin, sim, part, [3.], [1.3, 0., 0.])
        @test spin.position[1] ≈ 0.7 atol=1e-12
        @test spin.orientations[1].transverse ≈ exp(-rin(0.02, 0.2, 1.) * 0.6) rtol=1e-12
        # outside: 1.5 → 0.7 reflects at 1 and ends at 1.3; in [1, 1.3]: 0.3 + 0.3 um → 0.9 ms
        spin = mr.Spin(position=[1.5, 0., 0.], transverse=1., longitudinal=0.)
        mr.Evolve.draw_step!(spin, sim, part, [3.], [0.7, 0., 0.])
        @test spin.position[1] ≈ 1.3 atol=1e-12
        @test spin.orientations[1].transverse ≈ exp(-rout(0.02, 0.3, 1.) * 0.9) rtol=1e-12
    end

    @testset "ρ = 0: bit-identical to no layer (P3.2.5)" begin
        function run(c)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=c, diffusivity=3., R2=1/80, verbose=false)
            Random.seed!(7)
            snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [50.]; return_snapshot=true)[1]
        end
        a = run(mr.Cylinders(radius=1., repeats=[2.5, 2.5]))
        b = run(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_h=0.3))
        @test mr.position.(a) == mr.position.(b)
        @test per_spin(a) == per_spin(b)
    end

    @testset "h_inner = R: every inside spin decays at exactly 2ρ/R; outside spins untouched" begin
        cyl = mr.Cylinders(radius=1., layer_rho_inner=0.05, layer_h_inner=1.0)       # κ = 1/2 → ΔR2 = 0.1 /ms
        sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=cyl, diffusivity=3., verbose=false)
        Random.seed!(1)
        snap = mr.Snapshot(2000, sim, box; transverse=1., longitudinal=0.)
        mask = inside_mask(snap)
        res = mr.readout(snap, sim, [10., 50.]; return_snapshot=true)
        for (t, s) in zip((10., 50.), res)
            m = per_spin(s)
            @test all(isapprox.(m[mask], exp(-0.1 * t); rtol=1e-10))
            @test all(m[.!mask] .== 1.0)
        end
    end

    @testset "bulk R2 and the layer multiply exactly" begin
        cyl = mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2)
        function run(R2)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=cyl, diffusivity=3., R2=R2, verbose=false)
            Random.seed!(3)
            snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [50.]; return_snapshot=true)[1]
        end
        a = run(0.)
        b = run(1/80)
        @test mr.position.(a) == mr.position.(b)
        @test all(isapprox.(per_spin(b), per_spin(a) .* exp(-50 / 80); rtol=1e-10))
        @test minimum(per_spin(a)) < 1.0
    end

    @testset "side-specific: each side acts only in its own compartment (P3.4.7)" begin
        function final(c)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=c, diffusivity=3., verbose=false)
            Random.seed!(5)
            snap = mr.Snapshot(2000, sim, box; transverse=1., longitudinal=0.)
            mask = inside_mask(snap)
            return (per_spin(mr.readout(snap, sim, [50.]; return_snapshot=true)[1]), mask)
        end
        (m, mask) = final(mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2))
        @test all(m[.!mask] .== 1.0)
        @test minimum(m[mask]) < 1.0
        (m, mask) = final(mr.Cylinders(radius=1., layer_rho_outer=0.02, layer_h_outer=0.2))
        # reflection points land a roundoff distance outside r = R, inside the outer shell's boundary (measure zero)
        @test all(1 .- m[mask] .<= 1e-12)
        @test minimum(m[.!mask]) < 1.0
    end

    @testset "Review Focus 5: permeable cylinder with a layer runs, M⊥ finite and ≤ 1" begin
        cyl = mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, permeability=1.)
        sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=cyl, diffusivity=3., verbose=false)
        Random.seed!(9)
        snap = mr.Snapshot(1000, sim, box; transverse=1., longitudinal=0.)
        m = per_spin(mr.readout(snap, sim, [20.]; return_snapshot=true)[1])
        @test all(isfinite.(m))
        @test all(0 .< m .<= 1.0)
        @test minimum(m) < 1.0
    end

    @testset "finite RF pulse with a cylinder layer is rejected" begin
        seq = mr.SequenceParts.SequenceWaveform(
            (([], []), ([], []), ([], [])),
            [(0., 9., [mr.SequenceParts.ConstantPulse(0.25 / 9, 0., 0.)])],
            [], [10.], 10.,
        )
        sim = mr.Simulation(seq; geometry=mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2), verbose=false)
        @test_throws Exception mr.readout(100, sim, [10.])
    end
end
