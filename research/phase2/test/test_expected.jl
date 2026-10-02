@testset "expected" begin
    t = [0.0, 10.0, 50.0]

    @testset "exact cases" begin
        @test expected_null(t) == [1.0, 1.0, 1.0]
        @test expected_bulk(t, 0.1) ≈ exp.(-0.1 .* t)
        @test expected_half_gap(t, 0.1) ≈ exp.(-0.1 .* t)
        @test expected_half_gap(t, 0.1; R2_bulk=1/80) ≈ exp.(-(0.1 + 1/80) .* t)
        @test expected_additive([1.0, 0.5, 0.2], t, 0.02) ≈ [1.0, 0.5, 0.2] .* exp.(-0.02 .* t)
        @test roundoff_bound([0.0, 50.0], 0.04, 100_000) ≈ (1250 + 100_000) * eps()
    end

    ref = load_b4()

    @testset "load_b4" begin
        @test sort(collect(keys(ref))) == [0.0, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05]
        @test ref[0.01].times == TIMES_REF
        @test ref[0.0].mean == ones(11)
        @test ref[0.01].mean[end] ≈ 0.61379 atol=1e-5          # baseline_record.md B4 table
        @test ref[0.05].mean[end] ≈ 0.08802 atol=1e-5
        @test all(ref[0.01].sem .<= ref[0.01].sigma)
    end

    @testset "theta_for_signal" begin
        @test theta_for_signal(ref, ref[0.01].mean[end]) ≈ 0.01 rtol=1e-12
        θ = theta_for_signal(ref, (ref[0.01].mean[end] + ref[0.015].mean[end]) / 2)
        @test 0.01 < θ < 0.015
        @test isnan(theta_for_signal(ref, 0.99))      # above S(50) of the smallest θ > 0 (0.952)
        @test isnan(theta_for_signal(ref, 0.01))      # below S(50) of the largest θ (0.088)
        @test theta_for_signal(ref, ref[0.005].mean[6]; t=25.0) ≈ 0.005 rtol=1e-12
    end

    @testset "reference_at" begin
        @test reference_at(ref, 0.01) === ref[0.01].mean                 # grid point: stored curve unchanged
        mid = reference_at(ref, sqrt(0.01 * 0.015))                       # geometric midpoint
        @test all(ref[0.015].mean .<= mid .<= ref[0.01].mean)
        @test_throws ErrorException reference_at(ref, 0.2)                # outside the grid
    end

    @testset "harness reproduces one B4 curve exactly (P2.1.3 criterion)" begin
        r = run_walls(Walls(repeats=2., surface_relaxation=0.01); timestep=1e-2)  # seeds 1–10, 1e5 spins
        @test all(isapprox.(r.mean, ref[0.01].mean; rtol=1e-8))                     # CSV stores 10 significant digits
    end
end
