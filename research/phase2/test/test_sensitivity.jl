@testset "sensitivity tools" begin
    Sa = [1.0 0.9 0.8; 1.0 0.92 0.82]
    Sb = [1.0 0.95 0.9; 1.0 0.95 0.9]
    d = paired_difference(Sa, Sb)
    @test d.mean ≈ [0.0, -0.04, -0.09]
    @test d.sem ≈ [0.0, std([-0.05, -0.03]) / sqrt(2), std([-0.1, -0.08]) / sqrt(2)]
    p = peak([0.0, 10.0, 20.0], d.mean)
    @test (p.index, p.t, p.value) == (3, 20.0, d.mean[3])
    t = collect(0.0:10.0:100.0)
    f = monoexp_residual(t, exp.(-0.02 .* t))
    @test f.R ≈ 0.02 && f.residual < 1e-12
    g = monoexp_residual(t, 0.5 .* exp.(-0.01 .* t) .+ 0.5 .* exp.(-0.05 .* t))   # biexponential: residual > 0
    @test g.residual > 1e-3
    @test centroid("step", 0.2) ≈ 0.1
    @test centroid("linear", 0.3) ≈ 0.1
    @test centroid("exponential", 0.5; cutoff=2.0) ≈ 0.5 * (1 - 5 * exp(-4)) / (1 - exp(-4))
    @test centroid("exponential", 0.01; cutoff=2.0) ≈ 0.01 rtol=1e-10                # untruncated limit: h
end
