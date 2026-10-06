@testset "harness" begin
    @testset "rho_for" begin
        @test rho_for(0.1, 0.5) == 0.05
        @test rho_for(0.1, 0.0) == 0.0
    end

    @testset "layer_walls sides" begin
        g(w) = Simulation([]; geometry=w, verbose=false).geometry[1].layer
        @test g(layer_walls()) === nothing
        @test g(layer_walls(h=0.5, rate=0.0)) === nothing
        l = g(layer_walls(h=0.5, rate=0.1))[1]
        @test (l.positive.rho, l.positive.h, l.negative.rho, l.negative.h) == (0.05, 0.5, 0.05, 0.5)
        l = g(layer_walls(h=0.5, rate=0.1, side=:positive))[1]
        @test (l.positive.rho, l.negative.rho) == (0.05, 0.0)
        l = g(layer_walls(h=0.5, rate=0.1, side=:negative))[1]
        @test (l.positive.rho, l.negative.rho) == (0.0, 0.05)
        l = g(layer_walls(h=0.2, rate=0.1, side=:asymmetric, h_neg=0.4, rate_neg=0.05))[1]
        @test (l.positive.rho, l.positive.h, l.negative.rho, l.negative.h) == (rho_for(0.1, 0.2), 0.2, rho_for(0.05, 0.4), 0.4)
        @test_throws ErrorException layer_walls(side=:top)
    end

    @testset "run_walls reproduces the exact cases" begin
        t = [0.0, 10.0, 50.0]
        r = run_walls(layer_walls(); times=t, seeds=1:2, nspins=1000)            # B1: nothing relaxes
        @test size(r.S) == (2, 3)
        @test all(r.S .== 1.0)
        r = run_walls(layer_walls(); R2_bulk=1/80, times=t, seeds=1:2, nspins=1000)  # B2: bulk only
        @test all(isapprox.(r.S, exp.(-t' ./ 80); rtol=1e-11))
        r = run_walls(layer_walls(h=1.0, rate=0.1); times=t, seeds=1:2, nspins=1000) # h = w/2
        @test all(isapprox.(r.S, exp.(-0.1 .* t'); rtol=1e-10))
        @test r.tau == TAU_REF
        @test r.sem ≈ r.sigma ./ sqrt(2)
    end

    @testset "snapshot path agrees with the summed path" begin
        t = [0.0, 20.0]
        a = run_walls(layer_walls(h=0.3, rate=0.1); times=t, seeds=1:2, nspins=2000)
        b = run_walls(layer_walls(h=0.3, rate=0.1); times=t, seeds=1:2, nspins=2000, snapshots=true)
        @test all(isapprox.(a.S, b.S; rtol=1e-12))
        @test length(b.snapshots) == 2
        @test length(b.snapshots[1].spins) == 2000
    end

    @testset "occupancy" begin
        snap = Snapshot([[0.1, 0., 0.], [1.0, 0., 0.], [1.95, 0., 0.], [2.3, 0., 0.]])
        @test occupancy(snap, 0.2) == 0.5          # 0.1 and 1.95 are within 0.2 of a face
    end
end

@testset "harness: shapes" begin
    @test rho_for(0.1, 0.5; shape="linear") == 0.025
    @test rho_for(0.1, 0.5) == 0.05
    l = Simulation([]; geometry=layer_walls(h=0.5, rate=0.1, shape="linear"), verbose=false).geometry[1].layer[1]
    @test (l.positive.shape, l.positive.rho) == (:linear, 0.025)
    @test H_GRID_LINEAR[end] == 2.0
    r = run_walls(layer_walls(h=2.0, rate=0.1, shape="linear"); times=[0.0, 50.0], seeds=1:2, nspins=1000)
    @test all(isapprox.(r.S, exp.(-0.1 .* [0.0 50.0]); rtol=1e-10))            # exact h = w case
end
