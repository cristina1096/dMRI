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
