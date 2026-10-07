@testset "toy: exact vs endpoint layer rule" begin
    # h = w/2: every point is in exactly one layer, so both rules give exp(-rate·T) for every spin
    for m in (:exact, :endpoint)
        s = toy_layer_signal(; h=1.0, rate=0.1, tau=0.01, T=5.0, nspins=200, seed=1, method=m)
        @test all(isapprox.(s, exp(-0.1 * 5.0); rtol=1e-10))
    end
    # no layer
    @test all(toy_layer_signal(; h=0.0, rate=0.1, tau=0.01, T=5.0, nspins=50, seed=1, method=:exact) .== 1.0)
    # reproducible for a given seed
    a = toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=100, seed=3, method=:exact)
    b = toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=100, seed=3, method=:exact)
    @test a == b
    @test_throws ErrorException toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=10, seed=1, method=:midpoint)
end

@testset "toy: linear profile (P2.6.4)" begin
    # linear, h = w on both faces: ΔR2 = 2ρ/w = ΔR2(0) everywhere, so every spin decays at exactly `rate`
    for m in (:exact, :endpoint)
        s = toy_layer_signal(; h=2.0, rate=0.1, tau=0.01, T=5.0, nspins=200, seed=1, method=m, shape="linear")
        @test all(isapprox.(s, exp(-0.1 * 5.0); rtol=1e-10))
    end
    # step is the default and unchanged
    a = toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=100, seed=3, method=:exact)
    b = toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=100, seed=3, method=:exact, shape="step")
    @test a == b
    # linear at h attenuates less than step at h (same ΔR2(0), half the ρ), and the shape is used
    l = toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=100, seed=3, method=:exact, shape="linear")
    @test mean(l) > mean(a)
    @test_throws ErrorException toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=10, seed=1, method=:exact, shape="gaussian")
    # toy_rate passes the shape through
    r = toy_rate(; h=2.0, rate=0.1, D=3.0, tau=0.003, T=1.0, nspins=50, seeds=1:2, shape="linear")
    @test all(isapprox.(r.per_seed, 0.1; rtol=1e-10))
end
