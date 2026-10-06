@testset "timestep tools" begin
    @testset "log_taus" begin
        t = log_taus(1e-5, 1e-1)
        @test t[1] ≈ 1e-5 && t[end] ≈ 1e-1
        @test length(t) == 13
        @test all(diff(log10.(t)) .≈ 1 / 3)
    end

    @testset "toy_rate uses the actual simulated time" begin
        # h = w/2: every spin decays at exactly `rate`, so R = rate for any τ, even when τ does not divide T
        r = toy_rate(; h=1.0, rate=0.1, D=3.0, tau=0.003, T=1.0, nspins=50, seeds=1:2)
        @test r.T_eff ≈ round(1.0 / 0.003) * 0.003
        @test all(isapprox.(r.per_seed, 0.1; rtol=1e-10))
    end

    @testset "tau_conv on synthetic data" begin
        taus = [1e-4, 1e-3, 1e-2, 1e-1]
        R = [1.000, 1.001, 1.02, 1.2]                 # relative errors 0, 0.1%, 2%, 20%
        g = tau_conv_grid(taus, R, 1.0; tol=0.01)
        @test g.tau == 1e-3 && !g.all_pass
        i = tau_conv_interp(taus, R, 1.0; tol=0.01)
        @test 1e-3 < i.tau < 1e-2 && !i.lower_bound
        # where |error| crosses 1% between 0.1% (at 1e-3) and 2% (at 1e-2), linearly in log τ
        @test log10(i.tau) ≈ -3 + (0.01 - 0.001) / (0.02 - 0.001) atol=1e-12
        n = tau_conv_interp(taus, [1.0, 1.0, 1.0, 1.005], 1.0; tol=0.01)
        @test n.tau == 1e-1 && n.lower_bound      # never fails: lower bound only
    end

    @testset "noise_limited and fit_c_h" begin
        @test noise_limited(1.0, 0.001; tol=0.001)
        @test !noise_limited(1.0, 0.0001; tol=0.001)
        f = fit_c_h([0.01, 0.04], [0.1, 0.2], [3.0, 3.0])
        @test f.c_h ≈ 3.0 && f.spread ≈ 0.0
    end
end
