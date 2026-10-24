@testset "harness: geometry formulas" begin
    @test sv_inside(1.0) == 2.0
    @test sv_outside(1.0, 2.5, 2.5) ≈ 2π / (6.25 - π)
    @test shell_area_inside(1.0, 1.0) == 1.0
    @test shell_area_inside(1.0, 0.2) ≈ 0.36
    @test shell_area_outside(1.0, 0.2, 2.5, 2.5) ≈ π * 0.44 / (6.25 - π)
    @test rate_inside(0.02, 0.2, 1.0) ≈ 0.02 / 0.18
    @test rate_outside(0.02, 0.2, 1.0) ≈ 0.02 / 0.22
    @test radial_distance([3.0, 0.2, 7.0]; sx=2.5, sy=2.5) ≈ hypot(0.5, 0.2)
    @test radial_distance([3.0, 4.0, 0.0]) == 5.0
    @test roundoff_bound([0.0, 50.0], 1e-2, 100_000) ≈ (5000 + 100_000) * eps()
    @test expected_neglog_layer(0.02, 2.0, 50.0) ≈ 2.0
    @test expected_neglog_conventional(0.01, 3.0, 2.0, 50.0) ≈ 0.01 * sqrt(3 / π) * 100
end

@testset "harness: brute_rate" begin
    kw = (R=1.0, s=2.5, h_in=0.2, h_out=0.3, rate_in=1.0, rate_out=10.0)
    @test brute_rate(0.9, 0.0; kw...) == 1.0
    @test brute_rate(1.1, 0.0; kw...) == 10.0
    @test brute_rate(1.25, 0.0; kw...) == 20.0        # in the outer shells of two neighbours
    @test brute_rate(0.5, 0.0; kw...) == 0.0
    @test brute_rate(1.25 + 2.5 * 7, 0.0; kw...) == 20.0
end

@testset "harness: run_cylinders" begin
    R = 1.0
    bb = mr.BoundingBox([-R, -R, -R], [R, R, R])
    r = run_cylinders(Cylinders(radius=R); region=:inside, bounding_box=bb, times=[0.0, 1.0], seeds=1:2, nspins=4000)
    @test size(r.S) == (2, 2)
    @test all(r.S .== 1.0)
    @test all(r.neglog .== 0.0)
    @test all(abs.(r.nregion ./ 4000 .- π / 4) .< 5 * sqrt(π / 4 * (1 - π / 4) / 4000))
    # uniform inner layer h = R: −ln M⊥ = 0.1·t exactly for every inside spin
    u = run_cylinders(Cylinders(radius=R, layer_rho_inner=0.05, layer_h_inner=R); region=:inside, bounding_box=bb,
        times=[0.0, 5.0], seeds=1:2, nspins=2000)
    @test all(isapprox.(u.neglog[:, 2], 0.5; rtol=1e-9))
    o = run_cylinders(Cylinders(radius=R, repeats=[S_REF, S_REF]); region=:outside, sx=S_REF, sy=S_REF,
        times=[0.0, 1.0], seeds=1:1, nspins=4000)
    @test abs(o.nregion[1] / 4000 - (1 - π / S_REF^2)) < 5 * sqrt(0.5 * 0.5 / 4000)
end
