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

@testset "sweep configurations" begin
    for s in ("h", "w", "D", "rho")
        c = sweep_configs(s)
        @test !isempty(c)
        @test all(x -> !(x.shape == "step" && x.h > x.w / 2), c)               # step excluded for h > w/2
        @test all(x -> x.shape in ("step", "linear", "exponential"), c)
    end
    h = sweep_configs("h")
    @test sort(unique([x.h for x in h])) == [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
    @test count(x -> x.shape == "step", h) == 5                                 # h = 2.0 excluded for step
    @test sort(unique([x.w for x in sweep_configs("w")])) == [1.0, 2.0, 5.0, 10.0]
    @test sort(unique([x.D for x in sweep_configs("D")])) == [1.0, 3.0, 6.0]
    @test sort(unique([x.rho for x in sweep_configs("rho")])) == [0.003, 0.01, 0.03, 0.1]
    @test_throws ErrorException sweep_configs("T")
end

@testset "apparent T2 tools" begin
    @test apparent_t2(1.0, exp(-50 / 70), 50.0) ≈ 70.0
    r = exp(-50 / 70)
    @test t2_noise_rel(r, 100) ≈ sqrt(1 + r^2) / (100 * r * abs(log(r)))
    # per-seed T2, paired: layer decays at 1/69, conventional at 1/70 (two seeds, slightly different)
    times = [0.0, 25.0, 50.0]
    Sa = [exp.(-times ./ 69)'; exp.(-times ./ 68.5)']
    Sb = [exp.(-times ./ 70)'; exp.(-times ./ 70.5)']
    p = paired_t2(Sa, Sb, times, 50.0)
    @test p.T2a ≈ (69 + 68.5) / 2 && p.T2b ≈ (70 + 70.5) / 2
    @test p.dT2 ≈ ((69 - 70) + (68.5 - 70.5)) / 2
    @test p.dT2_sem ≈ std([69 - 70, 68.5 - 70.5]) / sqrt(2)
    @test p.rel ≈ mean([(69 - 70) / 70, (68.5 - 70.5) / 70.5])
    @test_throws ErrorException paired_t2(Sa, Sb, times, 40.0)              # TE must be a readout
    # timestep: tortuosity cap at small w
    @test tau_t2(1 / 70, 0.2, "step"; w=0.5, D=3.0) ≈ 0.03 * 0.25 / 3
    @test tau_t2(1 / 70, 0.2, "step"; w=2.0, D=3.0) == 1e-2
    @test tau_t2(0.06, 0.2, "linear"; w=2.0, D=3.0) ≈ 0.005 / (2 * 0.06 / 0.2)
    @test t2_noise_rel(exp(-5 / 70), 100) > 5 * t2_noise_rel(exp(-50 / 70), 100)    # short TE: much noisier T2 (ratio 6.4)
end

@testset "apparent T2 configurations" begin
    for s in ("w", "h", "D", "rho")
        c = t2_configs(s)
        @test all(x -> !(x.shape == "step" && x.h > x.w / 2), c)
        @test all(x -> x.shape in ("step", "linear", "exponential"), c)
    end
    @test sort(unique([x.w for x in t2_configs("w")])) == [0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
    @test sort(unique([x.h for x in t2_configs("h")])) == [0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
    @test count(x -> x.shape == "step", t2_configs("h")) == 6
    @test sort(unique([x.D for x in t2_configs("D")])) == [0.5, 1.0, 2.0, 3.0]
    @test sort(unique([x.rho for x in t2_configs("rho")])) == [0.005, 1 / 70, 0.03, 0.06]
    @test all(x -> x.rho == 1 / 70, t2_configs("w"))
    te = t2_configs("te")                                                        # Fig. 8A: w ∈ {1, 2, 5} at h = 0.2, h ∈ {0.05, 0.5} at w = 2
    @test sort(unique([(x.w, x.h) for x in te])) == [(1.0, 0.2), (2.0, 0.05), (2.0, 0.2), (2.0, 0.5), (5.0, 0.2)]
    @test length(te) == 15 && all(x -> x.D == 3.0 && x.rho == 1 / 70, te)
    @test_throws ErrorException t2_configs("angle")
end

@testset "same-centroid matching" begin
    @test h_for_centroid("step", 0.1; w=2.0) ≈ 0.2
    @test h_for_centroid("linear", 0.1; w=2.0) ≈ 0.3
    for w in (0.5, 2.0, 10.0), cbar in (0.01, 0.1, 0.2)
        cbar < w / 2 || continue
        h = h_for_centroid("exponential", cbar; w=w)
        @test centroid("exponential", h; cutoff=w) ≈ cbar rtol=1e-9
    end
    @test isnothing(h_for_centroid("step", 0.6; w=2.0))            # step would need h = 1.2 > w/2
    @test isnothing(h_for_centroid("linear", 0.7; w=2.0))          # linear would need h = 2.1 > w
    @test isnothing(h_for_centroid("exponential", 1.0; w=2.0))     # truncated exponential: centroid < w/2
    c = t2_configs("c"; match=:centroid)
    @test sort(unique([x.cbar for x in c])) == [0.01, 0.025, 0.05, 0.1, 0.25, 0.5]
    @test all(x -> centroid(x.shape, x.h; cutoff=x.w) ≈ x.cbar, c)
    @test all(x -> x.rho == 1 / 70, c)
    @test length(t2_configs("te"; match=:centroid)) == 15
    @test all(x -> centroid(x.shape, x.h; cutoff=x.w) ≈ 0.1, t2_configs("w"; match=:centroid))
    @test_throws ErrorException t2_configs("h"; match=:centroid)   # the centroid sweep "c" replaces the h sweep
    @test_throws ErrorException t2_configs("c")                     # "c" needs match = :centroid
end
