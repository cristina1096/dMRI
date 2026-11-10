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
        @test Layers.centroid(Layers.LayerSide(0.02, 0.3, :linear)) ≈ 0.1          # linear now supported (plan 2026-11-06)
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
        @test geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=0.2, layer_shape="linear"))[1].layer[1].inner.shape == :linear   # supported since plan 2026-11-06
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=0.2, layer_shape="bogus"))
        # h_outer = gap exactly is allowed; a side with ρ = 0 is not checked against the gap
        @test geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.1, layer_h=0.5))[1].layer !== nothing
        # h_outer = gap where the gap is not exact in binary (1 − 0.8 = 0.19999999999999996) is still allowed
        @test geom(mr.Cylinders(radius=0.4, repeats=[1., 1.], layer_rho_outer=0.1, layer_h_outer=0.2))[1].layer !== nothing
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

@testset "test_layer_cylinders.jl: non-step profiles along a piece (P3.3)" begin
    Layers = mr.Geometries.Internal.Layers
    P(x, y) = SVector{2, Float64}(x, y)
    O = P(0., 0.)
    expside(h, L) = Layers.LayerSide(0.02, h, :exponential, L)

    @testset "centroid and κ for all profiles (D3.1)" begin
        @test Layers.centroid(Layers.LayerSide(0.02, 0.3, :linear)) ≈ 0.1
        c = 2.0                                           # support 0.4, h = 0.2
        @test Layers.centroid(expside(0.2, 0.4)) ≈ 0.2 * (1 - (c + 1) * exp(-c)) / (1 - exp(-c))
        n = 200_000                                       # κ by midpoint quadrature, exponential outside, R = 1
        κq = sum(exp(-(k - 0.5) / n * c) / (1 - exp(-c)) * (1 + (k - 0.5) / n * 0.4 / 1.0) for k in 1:n) / n * c
        @test Layers.kappa(expside(0.2, 0.4), 1.0, false) ≈ κq rtol=1e-8
    end

    @testset "quad_adaptive on polynomials, a kink and a jump-free exponential" begin
        @test Layers.quad_adaptive(s -> s^7, 0.0, 1.0, 1e-12) ≈ 1 / 8 rtol=1e-12
        @test Layers.quad_adaptive(s -> exp(-3s), 0.2, 0.9, 1e-12) ≈ (exp(-0.6) - exp(-2.7)) / 3 rtol=1e-12
        @test Layers.quad_adaptive(s -> sqrt(abs(s - 0.3)), 0.0, 0.3, 1e-10) ≈ 2 / 3 * 0.3^1.5 rtol=1e-8
    end

    @testset "disk_interval and disk_length (refactor keeps disk_length)" begin
        @test Layers.disk_interval(1.0, 0.0, 0.0, 0.25) == (0.0, 0.5)
        @test Layers.disk_length(1.0, 0.0, 0.0, 0.25) ≈ 0.5
        (lo, hi) = Layers.disk_interval(1.0, 0.0, 2.0, 1.0)
        @test hi <= lo
    end

    "∫₀¹ √(a s² + b s + q) ds in closed form (independent of the package; a > 0)."
    function sqrtquad(a, b, q, s0, s1)
        F(s) = (Q = a * s^2 + b * s + q; (2a * s + b) * sqrt(Q) / (4a) + (4a * q - b^2) / (8a^1.5) * log(abs(2sqrt(a) * sqrt(Q) + 2a * s + b)))
        return F(s1) - F(s0)
    end

    @testset "linear profile vs closed form (P3.3.2)" begin
        # inside, h = 0.3, support [0.7, 1]: g = 2(1 − (R − r)/h) is linear in r, so ∫g ds = ∫(2 − 2(R − r)/h) ds over the
        # shell parts, with ∫ r ds from the √(quadratic) closed form on each part.
        side = Layers.LayerSide(0.02, 0.3, :linear)
        R, h = 1.0, 0.3
        for (p0, p1) in ((P(-2., 0.5), P(2., 0.5)), (P(-0.9, -0.2), P(0.95, 0.4)), (P(0.75, 0.), P(0.95, 0.3)))
            Δ = p1 - p0; a = Δ ⋅ Δ; b = 2(p0 ⋅ Δ); q = p0 ⋅ p0
            (h0, h1) = Layers.disk_interval(a, b, q, R^2)
            (l0, l1) = Layers.disk_interval(a, b, q, (R - h)^2)
            parts = l1 > l0 ? [(h0, l0), (l1, h1)] : [(h0, h1)]
            closed = sum((x1 - x0) * (2 - 2R / h) + 2 / h * sqrtquad(a, b, q, x0, x1) for (x0, x1) in parts if x1 > x0; init=0.0)
            @test Layers.shell_profile_integral(side, R, true, O, p0, p1) ≈ closed atol=1e-10
        end
    end

    @testset "wall regression: linear d(s) reproduces side_exponent (P3.3.3)" begin
        # along a radial piece r(s) is linear in s, so d(s) is exactly linear, as on a wall
        R = 2.0
        for shape in (:linear, :exponential)
            side = shape == :linear ? Layers.LayerSide(0.02, 0.3, :linear) : Layers.LayerSide(0.02, 0.3, :exponential, 2.0)
            d0, d1 = 0.05, 0.8
            p0, p1 = P(R - d0, 0.), P(R - d1, 0.)          # radial, inside: d = R − r, exactly linear in s
            wall = Layers.side_exponent(side, d0, d1, 1.0) / (side.rho / side.h)
            @test Layers.shell_profile_integral(side, R, true, O, p0, p1) ≈ wall rtol=1e-9
        end
    end

    @testset "brute force, inside and outside, linear and exponential (P3.3.1)" begin
        Random.seed!(11)
        K = 200_000
        g(shape, u, c) = Layers.profile_g(shape, u, c)
        maxerr = 0.0
        for _ in 1:300, (shape, inside, L, h) in ((:linear, true, 0.3, 0.3), (:linear, false, 0.4, 0.4),
                                                    (:exponential, true, 1.0, 0.2), (:exponential, false, 0.5, 0.15))
            side = shape == :linear ? Layers.LayerSide(0.02, h, :linear) : Layers.LayerSide(0.02, h, :exponential, L)
            p0 = P((3 .* rand(2) .- 1.5)...); p1 = P((3 .* rand(2) .- 1.5)...)
            dist(r) = inside ? 1.0 - r : r - 1.0
            brute = sum(k -> (r = norm(p0 + (k - 0.5) / K * (p1 - p0)); d = dist(r); 0 <= d <= L ? g(shape, d / h, L / h) : 0.0), 1:K) / K
            maxerr = max(maxerr, abs(Layers.shell_profile_integral(side, 1.0, inside, O, p0, p1) - brute))
        end
        @test maxerr < 2e-4            # exponential truncation jump: ≤ g(c)·(crossings)/(2K) + smooth part O(K⁻²)
    end

    @testset "Review Focus 1: piece through the axis (inner exponential reaches r = 0)" begin
        side = Layers.LayerSide(0.02, 0.2, :exponential, 1.0)
        # radial through the centre, from (−1, 0) to (1, 0): ds = dx/2, d = 1 − |x|, so
        # ∫₀¹ g(d/h) ds = ∫₀¹ g(d/h) dd = h ∫₀^c g(u) du = h (g normalised over its support)
        expected = 0.2
        @test Layers.shell_profile_integral(side, 1.0, true, O, P(-1., 0.), P(1., 0.)) ≈ expected rtol=1e-9
        @test isfinite(Layers.shell_profile_integral(side, 1.0, true, O, P(-1., 1e-14), P(1., 1e-14)))
    end

    @testset "Review Focus 2: tangent and grazing pieces" begin
        lin = Layers.LayerSide(0.02, 0.2, :linear)
        @test Layers.shell_profile_integral(lin, 1.0, true, O, P(-2., 1.0), P(2., 1.0)) == 0.0          # tangent to r = R
        v = Layers.shell_profile_integral(lin, 1.0, true, O, P(-2., 0.8), P(2., 0.8))                    # tangent to r = R − h
        @test isfinite(v) && v > 0
        ex = Layers.LayerSide(0.02, 0.1, :exponential, 0.3)
        @test isfinite(Layers.shell_profile_integral(ex, 1.0, false, O, P(-2., 1.3), P(2., 1.3)))        # tangent to the support edge
    end

    @testset "stationary points" begin
        lin = Layers.LayerSide(0.02, 0.2, :linear)
        @test Layers.shell_profile_integral(lin, 1.0, true, O, P(0.9, 0.), P(0.9, 0.)) ≈ Layers.profile_g(:linear, 0.5)
        @test Layers.shell_profile_integral(lin, 1.0, true, O, P(0.5, 0.), P(0.5, 0.)) == 0.0
    end
end

@testset "test_layer_cylinders.jl: linear and exponential on cylinders (D3.5)" begin
    Layers = mr.Geometries.Internal.Layers
    geom(c) = mr.Simulation([]; geometry=c, verbose=false).geometry
    V3(x, y, z=0.) = SVector{3, Float64}(x, y, z)
    per_spin(snap) = [s.orientations[1].transverse for s in snap.spins]

    @testset "resolution: shapes and supports" begin
        l = geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.3, layer_shape="linear"))[1].layer[1]
        @test (l.inner.shape, l.inner.cutoff, l.outer.shape, l.outer.cutoff) == (:linear, 0.3, :linear, 0.3)
        l = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))[1].layer[1]
        @test l.inner.cutoff == 1.0                      # inside: support to the axis (D3.5)
        @test l.outer.cutoff ≈ 0.5                       # outside: default = smallest surface gap
        l = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=0.3))[1].layer[1]
        @test l.outer.cutoff == 0.3
        @test l.inner.cutoff == 1.0                      # layer_cutoff does not change the inner support
    end

    @testset "validation" begin
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))   # Review Focus 3
        @test geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=2.))[1].layer !== nothing
        @test_throws ErrorException geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2,
            layer_shape="exponential", layer_cutoff=0.6))                                                                       # > gap
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=1.5, layer_shape="linear"))          # h_inner > R
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=-1.))
        # inert outer side (ρ = 0) needs no cutoff
        @test geom(mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2, layer_shape="exponential"))[1].layer !== nothing
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_shape="bogus"))
    end

    @testset "Review Focus 4: step results bit-identical" begin
        function run(c)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=c, diffusivity=3., verbose=false)
            Random.seed!(21)
            snap = mr.Snapshot(1000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [20.]; return_snapshot=true)[1]
        end
        a = run(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2))
        b = run(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape="step"))
        @test per_spin(a) == per_spin(b)
        # reference values recorded before this plan (Task 2 Step 2 writes them; see the brief)
        @test sum(per_spin(a)) == 447.8892123934866
    end

    @testset "radial step with a reflection: linear and exponential, hand-computed" begin
        seq = build_sequence([1.2, :readout])
        part(sim) = first(p for p in mr.parts([seq], 0., mr.TimeStep(1.2, Inf)) if p.duration > 0)
        n = 200_000
        κq(gfun, umax, h, R) = sum(gfun((k - 0.5) / n * umax) * (1 - (k - 0.5) / n * umax * h / R) for k in 1:n) / n * umax
        # linear inside h = 0.2: radial path 0.5 → 1 → 0.7 (1.5 ms/um); in the shell 0.4 um → 0.6 ms; mean g over a radial
        # traverse of the whole shell is 1, so E = (ρ/(hκ))·0.6 with κ = 1 − h/(3R)
        sim = mr.Simulation(seq; geometry=mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2, layer_shape="linear"),
            diffusivity=3., verbose=false)
        spin = mr.Spin(position=[0.5, 0., 0.], transverse=1., longitudinal=0.)
        mr.Evolve.draw_step!(spin, sim, part(sim), [3.], [1.3, 0., 0.])
        κ = κq(u -> 2(1 - u), 1.0, 0.2, 1.0)
        @test spin.orientations[1].transverse ≈ exp(-0.02 / (0.2 * κ) * 0.6) rtol=1e-9
        # exponential inside h = 0.2 (support R = 1, c = 5): d runs 0.5 → 0 then 0 → 0.3, dt = 1.5 dd, so
        # ∫g dt = 1.5 h [(1 − e^(−2.5)) + (1 − e^(−1.5))]/(1 − e^(−5))
        sim = mr.Simulation(seq; geometry=mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2, layer_shape="exponential"),
            diffusivity=3., verbose=false)
        spin = mr.Spin(position=[0.5, 0., 0.], transverse=1., longitudinal=0.)
        mr.Evolve.draw_step!(spin, sim, part(sim), [3.], [1.3, 0., 0.])
        κ = κq(u -> exp(-u) / (1 - exp(-5)), 5.0, 0.2, 1.0)
        I = 1.5 * 0.2 * ((1 - exp(-2.5)) + (1 - exp(-1.5))) / (1 - exp(-5))
        @test spin.orientations[1].transverse ≈ exp(-0.02 / (0.2 * κ) * I) rtol=1e-8
    end

    @testset "Review Focus 5: reach and timestep use the support" begin
        g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho_outer=0.02, layer_h_outer=0.05, layer_shape="exponential"))
        # point on the diagonal at r = 1.4 (0.4 from this surface: beyond h = 0.05, inside the support 0.5); its distances to
        # the neighbouring surfaces are ≥ 0.80, outside their supports, so only this cylinder contributes
        c = 0.5 / 0.05
        κ = Layers.kappa(Layers.LayerSide(0.02, 0.05, :exponential, 0.5), 1.0, false)
        p = V3(1.4 / sqrt(2), 1.4 / sqrt(2))
        @test Layers.layer_exponent(g, p, p, 1.0) ≈ 0.02 / (0.05 * κ) * exp(-0.4 / 0.05) / (1 - exp(-c)) rtol=1e-10
        # overlap bound counts the 4 neighbours whose outer support (R + 0.5) intersects: 2(R + 0.5) = 3 > 2.5
        @test Layers.max_layer_rate(g) ≈ 5 * Layers.cylinder_surface_rate(Layers.LayerSide(0.02, 0.05, :exponential, 0.5), 1.0, false)
    end

    @testset "layer exponent vs brute force, linear and exponential, square array with overlap" begin
        for (shape, kw) in ((:linear, (layer_h_inner=0.3, layer_h_outer=0.4)), (:exponential, (layer_h_inner=0.2, layer_h_outer=0.15)))
            g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.03, layer_shape=String(shape); kw...))
            l = g[1].layer[1]
            pre_in = l.inner.rho / l.inner.h / Layers.kappa(l.inner, 1.0, true)
            pre_out = l.outer.rho / l.outer.h / Layers.kappa(l.outer, 1.0, false)
            function rate(x, y)
                t = 0.0
                for i in -2:2, j in -2:2
                    r = hypot(x - i * 2.5, y - j * 2.5)
                    d = 1.0 - r
                    0 <= d <= l.inner.cutoff && (t += pre_in * Layers.profile_g(shape, d / l.inner.h, l.inner.cutoff / l.inner.h))
                    d = r - 1.0
                    0 <= d <= l.outer.cutoff && (t += pre_out * Layers.profile_g(shape, d / l.outer.h, l.outer.cutoff / l.outer.h))
                end
                return t
            end
            Random.seed!(5)
            K = 20_000
            err = 0.0
            for _ in 1:100
                p0 = V3((2.5 .* (2 .* rand(2) .- 1))..., 0.); p1 = V3((2.5 .* (2 .* rand(2) .- 1))..., 0.)
                brute = sum(k -> rate((p0 + (k - 0.5) / K * (p1 - p0))[1:2]...), 1:K) / K
                err = max(err, abs(Layers.layer_exponent(g, p0, p1, 1.0) - brute))
            end
            @test err < 2e-5                     # truncation jumps ≤ 0.007/(2K) per crossing, smooth parts O(K⁻²)
        end
    end

    @testset "simulation: side-specific and additivity hold for linear and exponential" begin
        for shape in ("linear", "exponential")
            c = mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape=shape)
            function run(R2)
                sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=c, diffusivity=3., R2=R2, verbose=false)
                Random.seed!(3)
                snap = mr.Snapshot(1000, sim, 500; transverse=1., longitudinal=0.)
                return mr.readout(snap, sim, [20.]; return_snapshot=true)[1]
            end
            a = run(0.); b = run(1 / 80)
            @test mr.position.(a) == mr.position.(b)
            @test all(isapprox.(per_spin(b), per_spin(a) .* exp(-20 / 80); rtol=1e-10))
            @test minimum(per_spin(a)) < 1.0
        end
    end

    @testset "JSON round trip keeps layer_cutoff and shape" begin
        c = mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=0.4)
        io = IOBuffer()
        mr.write_geometry(io, c)
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_cutoff.value == 0.4
        @test back.layer_shape.value == "exponential"
    end
end

# measured through functions that receive everything as arguments: a closure over the non-const `Layers` binding inside a
# @testset allocates by itself (dynamic dispatch), which is not what is being tested
measure_spi(spi, side, O, A, B) = (spi(side, 1.0, true, O, A, B); @allocated spi(side, 1.0, true, O, A, B))
measure_quad(quad, g) = (quad(g, 0.0, 1.0, 1e-10); @allocated quad(g, 0.0, 1.0, 1e-10))
quad_test_integrand(s) = exp(-s)

@testset "test_layer_cylinders.jl: the quadrature path does not allocate (Group 1 review I1)" begin
    Layers = mr.Geometries.Internal.Layers
    side = Layers.LayerSide(0.02, 0.2, :exponential, 1.0)
    O, A, B = SVector(0., 0.), SVector(0.95, 0.1), SVector(0.6, -0.3)
    @test measure_spi(Layers.shell_profile_integral, side, O, A, B) == 0
    @test measure_quad(Layers.quad_adaptive, quad_test_integrand) == 0
end

@testset "test_layer_cylinders.jl: overlapping outer layers, exact summation (P3.5.3)" begin
    Layers = mr.Geometries.Internal.Layers
    geom(c) = mr.Simulation([]; geometry=c, verbose=false).geometry
    V3(x, y, z=0.) = SVector{3, Float64}(x, y, z)
    "Brute-force summed outer rate over images (±3 cells), support inclusive, written from the layer's own parameters."
    function brute_rate(g, x, y, sx, sy, shape)
        l = g[1].layer[1]
        pre = l.outer.rho / l.outer.h / Layers.kappa(l.outer, 1.0, false)
        t = 0.0
        for i in -3:3, j in -3:3
            d = hypot(x - i * sx, y - j * sy) - 1.0
            0 <= d <= l.outer.cutoff && (t += pre * Layers.profile_g(shape, d / l.outer.h, l.outer.cutoff / l.outer.h))
        end
        return t
    end
    "Brute-force ∫ rate dt along p0 → p1 (dt = 1): split at every crossing of r = R and r = R + L of every image, 64-point
    Gauss–Legendre on each smooth piece."
    function brute_exponent(g, p0, p1, sx, sy, shape)
        L = g[1].layer[1].outer.cutoff
        a, b = p0[1:2], p1[1:2] .- p0[1:2]
        cuts = [0.0, 1.0]
        for i in -3:3, j in -3:3, rad in (1.0, 1.0 + L)
            c = [i * sx, j * sy]
            q = a .- c
            A, B, C = b' * b, 2 * (q' * b), q' * q - rad^2
            disc = B^2 - 4A * C
            A > 0 && disc > 0 && append!(cuts, filter(s -> 0 < s < 1, [(-B - sqrt(disc)) / (2A), (-B + sqrt(disc)) / (2A)]))
        end
        sort!(cuts)
        (xg, wg) = Layers.gauss_legendre(64)
        total = 0.0
        for k in 1:length(cuts) - 1
            lo, hi = cuts[k], cuts[k + 1]
            hi - lo < 1e-15 && continue
            for (xi, wi) in zip(xg, wg)
                s = lo + (xi + 1) / 2 * (hi - lo)
                total += wi * (hi - lo) / 2 * brute_rate(g, (a .+ s .* b)..., sx, sy, shape)
            end
        end
        return total
    end
    "Random straight pieces with both ends in the fluid (outside every cylinder) inside the box `centre ± half`."
    function pieces(sx, sy, centre, half, maxlen; n=200, seed=11)
        Random.seed!(seed)
        fluid(p) = all(hypot(p[1] - i * sx, p[2] - j * sy) > 1.0 for i in -3:3, j in -3:3)
        out = Tuple{SVector{3, Float64}, SVector{3, Float64}}[]
        while length(out) < n
            p0 = centre .+ half .* (2 .* rand(2) .- 1)
            dir = 2π * rand()
            p1 = p0 .+ maxlen * rand() .* [cos(dir), sin(dir)]
            (fluid(p0) && fluid(p1)) && push!(out, (V3(p0...), V3(p1...)))
        end
        return out
    end
    relerr(e, b) = abs(e - b) / max(abs(b), 1e-300)

    @testset "(a) rectangular array, one-direction overlap: linear h = 0.4 and exponential (cutoff = gap 0.5)" begin
        for (shape, h) in ((:linear, 0.4), (:exponential, 0.2))
            g = geom(mr.Cylinders(radius=1., repeats=[2.5, 4.0], layer_rho_outer=0.03, layer_h_outer=h, layer_shape=String(shape)))
            @test g[1].layer[1].outer.cutoff ≈ (shape == :linear ? 0.4 : 0.5)
            pts = pieces(2.5, 4.0, [1.25, 0.0], [0.6, 1.5], 0.5)
            @test maximum(relerr(Layers.layer_exponent(g, p, p, 1.0), brute_rate(g, p[1], p[2], 2.5, 4.0, shape)) for (p, _) in pts) <= 1e-12
            @test maximum(relerr(Layers.layer_exponent(g, p0, p1, 1.0), brute_exponent(g, p0, p1, 2.5, 4.0, shape)) for (p0, p1) in pts) <= 1e-10
        end
    end

    @testset "(b) 4-fold interstitial overlap: square 4.0, exponential with default cutoff 2.0" begin
        g = geom(mr.Cylinders(radius=1., repeats=[4.0, 4.0], layer_rho_outer=0.03, layer_h_outer=0.3, layer_shape="exponential"))
        @test g[1].layer[1].outer.cutoff ≈ 2.0
        @test brute_rate(g, 2.0, 2.0, 4.0, 4.0, :exponential) > 0          # the corner point is reached by all four
        pts = pieces(4.0, 4.0, [2.0, 2.0], [0.5, 0.5], 0.5)
        @test maximum(relerr(Layers.layer_exponent(g, p, p, 1.0), brute_rate(g, p[1], p[2], 4.0, 4.0, :exponential)) for (p, _) in pts) <= 1e-12
        @test maximum(relerr(Layers.layer_exponent(g, p0, p1, 1.0), brute_exponent(g, p0, p1, 4.0, 4.0, :exponential)) for (p0, p1) in pts) <= 1e-10
    end

    @testset "(c) h = gap: step and linear at square 2.5, gap midline and facing surfaces" begin
        for shape in (:step, :linear)
            g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho_outer=0.03, layer_h_outer=0.5, layer_shape=String(shape)))
            for p in (V3(1.25, 0.), V3(1.25, 0.3), V3(1.0 + 1e-9, 0.), V3(1.5 - 1e-9, 0.), V3(0., 1.25))
                @test Layers.layer_exponent(g, p, p, 1.0) ≈ brute_rate(g, p[1], p[2], 2.5, 2.5, shape) rtol=1e-12
            end
        end
    end

    @testset "(d) max_layer_rate bounds the summed rate (grid maximum)" begin
        for (rep, kw, shape) in (([2.5, 4.0], (layer_h_outer=0.4, layer_shape="linear"), :linear),
                                 ([2.5, 4.0], (layer_h_outer=0.2, layer_shape="exponential"), :exponential),
                                 ([4.0, 4.0], (layer_h_outer=0.3, layer_shape="exponential"), :exponential),
                                 ([2.5, 2.5], (layer_h_outer=0.5, layer_shape="step"), :step))
            g = geom(mr.Cylinders(radius=1., repeats=rep, layer_rho_outer=0.03; kw...))
            gridmax = maximum(brute_rate(g, x, y, rep[1], rep[2], shape) for x in range(-rep[1] / 2, rep[1] / 2, length=301),
                              y in range(-rep[2] / 2, rep[2] / 2, length=301))
            @test Layers.max_layer_rate(g) >= gridmax * (1 - 1e-12)
        end
    end
end
