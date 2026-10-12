@testset "continuum benchmark (1D, walls)" begin
    t = collect(0.0:10.0:200.0)
    # eigen-decomposition roundoff (≈ ε·‖A‖ in the rates) limits S to ~1e-8 absolute: tolerance 1e-7
    # no relaxation: S = 1
    @test all(isapprox.(continuum_signal(w=2.0, D=3.0, times=t), 1.0; atol=1e-7))
    # bulk only
    @test all(isapprox.(continuum_signal(w=2.0, D=3.0, times=t, R2_bulk=1 / 80), exp.(-t ./ 80); atol=1e-7))
    # step with h = w/2: every point in exactly one layer → uniform rate ρ/h
    @test all(isapprox.(continuum_signal(w=2.0, D=3.0, times=t, layer=("step", 0.01, 1.0)), exp.(-0.01 .* t); atol=1e-7))
    # linear with h = w: overlapping layers sum to the uniform rate 2ρ/w
    @test all(isapprox.(continuum_signal(w=2.0, D=3.0, times=t, layer=("linear", 0.01, 2.0)), exp.(-0.01 .* t); atol=1e-7))
    # every model starts with slope −2ρ/w (mean-exponent identity at t → 0)
    δ = [0.0, 1e-4]
    for m in ((layer=("exponential", 0.01, 0.2),), (layer=("step", 0.01, 0.05),), (surface=0.01,))
        S = continuum_signal(; w=2.0, D=3.0, times=δ, m...)
        @test (1 - S[2]) / 1e-4 ≈ 2 * 0.01 / 2.0 rtol=1e-3
    end
    # conventional model: lowest decay rate = D·k², k the smallest root of k·tan(k·w/2) = ρ/D (surface sink on both faces)
    ρ, w, D = 0.01, 2.0, 3.0
    f(k) = k * tan(k * w / 2) - ρ / D
    lo, hi = 1e-9, pi / w - 1e-9
    for _ in 1:200
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    k = (lo + hi) / 2
    @test continuum_rates(w=w, D=D, surface=ρ)[1] ≈ D * k^2 rtol=1e-3
    # grid convergence at the P2.5 reference point
    a = continuum_signal(w=2.0, D=3.0, times=[90.0], layer=("exponential", 0.01, 0.2), N=3000)[1]
    b = continuum_signal(w=2.0, D=3.0, times=[90.0], layer=("exponential", 0.01, 0.2), N=6000)[1]
    c = continuum_signal(w=2.0, D=3.0, times=[90.0], surface=0.01, N=3000)[1]
    d = continuum_signal(w=2.0, D=3.0, times=[90.0], surface=0.01, N=6000)[1]
    @test abs((a - c) - (b - d)) <= 0.02 * abs(b - d)
    @test_throws ErrorException continuum_signal(w=2.0, D=3.0, times=t, layer=("gaussian", 0.01, 0.2))
end
