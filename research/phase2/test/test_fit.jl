@testset "fit: curvature and monotonicity" begin
    t = collect(0.0:5.0:50.0)
    @test curvature(t, exp.(-0.1 .* t)) ≈ 0.0 atol=1e-12
    @test curvature(t, exp.(-0.1 .* t .+ 0.001 .* t .^ 2)) ≈ 0.001 rtol=1e-8
    @test increasing_steps([0.1, 0.2, 0.2001], [0.001, 0.001, 0.001]) == [true, false]
    @test sweep_nspins() isa Int
end

@testset "fit: h* fit" begin
    t = collect(0.0:5.0:50.0)
    hgrid = [0.0, 0.5, 1.0]
    curves = [exp.(-k .* t) for k in (0.0, 0.05, 0.1)]
    sem = [vcat(0.0, fill(1e-3, 10)) for _ in hgrid]          # t = 0: zero uncertainty
    ref = interp_curve(hgrid, curves, 0.3)
    ref_sem = vcat(0.0, fill(1e-3, 10))

    @test interp_curve(hgrid, curves, 0.0) == curves[1]
    @test interp_curve(hgrid, curves, 1.5) == curves[end]      # clamped
    f = fit_h(hgrid, curves, sem, ref, ref_sem)
    @test f.h ≈ 0.3 atol=1e-4
    @test f.chi2 ≈ 0.0 atol=1e-6
    @test f.dof == 10 - 1                                      # 10 readouts with σ > 0, one fitted parameter
    @test !f.at_edge
    @test all(isfinite, f.residuals)

    far = exp.(-0.5 .* t)                                      # decays faster than any model curve
    g = fit_h(hgrid, curves, sem, far, ref_sem)
    @test g.h == 1.0
    @test g.at_edge

    S_by_h = [repeat(c', 10) for c in curves]                  # 10 identical seeds → zero jackknife spread
    j = jackknife_h(hgrid, S_by_h, ref, ref_sem)
    @test j.h ≈ 0.3 atol=1e-4
    @test j.sigma ≈ 0.0 atol=1e-6
end

@testset "fit: χ² and parabola refinement" begin
    t = collect(0.0:5.0:50.0)
    ref = exp.(-0.02 .* t); ref_sem = vcat(0.0, fill(1e-3, 10))
    @test chi2_of(ref, ref_sem, ref, ref_sem) == 0.0
    @test chi2_of(ref .+ vcat(0.0, fill(2e-3, 10)), fill(0.0, 11), ref, ref_sem) ≈ 10 * 4.0
    @test parabola_vertex([0.2, 0.3, 0.45], [(h - 0.3)^2 + 1 for h in (0.2, 0.3, 0.45)]) ≈ 0.3
    @test parabola_vertex([1.0, 2.0, 4.0], [3 * (h - 2.5)^2 for h in (1.0, 2.0, 4.0)]) ≈ 2.5
end

@testset "fit: Gauss–Newton step and sweep slope" begin
    t = collect(0.0:5.0:50.0)
    S(h) = exp.(-0.1 .* h .* t)                     # model: S depends on h through exp(-k h t)
    ref = S(0.30); ref_sem = vcat(0.0, fill(1e-5, 10))
    dSdh(h) = -0.1 .* t .* S(h)
    h1 = 0.29 + gauss_newton_step(S(0.29), zeros(11), ref, ref_sem, dSdh(0.29))
    @test abs(h1 - 0.30) < abs(0.29 - 0.30) / 50      # one step reduces the error by > 50×
    hgrid = [0.0, 0.2, 0.4]
    curves = [S(h) for h in hgrid]
    @test sweep_slope(hgrid, curves, 0.3) ≈ dSdh(0.3) rtol=1e-10   # ln S exactly linear in h here
end
