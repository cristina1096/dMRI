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
