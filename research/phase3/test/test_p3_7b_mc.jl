@testset "p3_7b: outside centroid matching (Review Focus 2)" begin
    h = h_for_cbar_out("exponential", 0.1, 0.5)
    @test outer_centroid("exponential", h, 0.5) ≈ 0.1 atol = 1e-10
    @test h_for_cbar_out("step", 0.1, 0.5) == 0.2
    @test h_for_cbar_out("linear", 0.1, 0.5) ≈ 0.3
    @test_throws ErrorException h_for_cbar_out("exponential", 0.25, 0.5)   # c̄ ≥ gap/2 unreachable
    @test_throws ErrorException h_for_cbar_out("linear", 0.1, 0.25)        # support 0.3 > gap 0.25
end

@testset "p3_7b: run_condition_outside (tiny)" begin
    s, ρ, t = 2.5, 0.02, [0.0, 10.0]
    models = [("conventional", "conventional", 0.0), ("step_h0.2", "step", 0.2)]
    r = run_condition_outside(s=s, D=3.0, rho=ρ, R2_bulk=0.0, models=models, times=t, tau=1e-2, nspins=500, seeds=1:2)
    @test r.identical_trajectories
    @test r.S_mean["step_h0.2"][1] == 1.0 && r.S_mean["conventional"][1] == 1.0
    @test r.dS_mean["step_h0.2"][1] == 0.0
    @test r.S_ann["step_h0.2"] ≈ annulus_signal(R=1.0, R_out=krogh_radius(s), D=3.0, times=t, layer=("step", ρ, 0.2, 0.2))
    @test r.S_ann["conventional"] ≈ annulus_signal(R=1.0, R_out=krogh_radius(s), D=3.0, times=t, surface=ρ)
    @test r.kept["step_h0.2"] == 1.0
    @test size(r.S["step_h0.2"]) == (2, 2)
    @test r.S_mean["step_h0.2"][2] < 1.0
end

@testset "p3_7b: sum_compartments (Review Focus 3)" begin
    t = [0.0, 1.0]
    inside = (times=t, S_mean=Dict("m" => [1.0, 0.8]), S_sem=Dict("m" => [0.0, 0.03]),
              dS_mean=Dict("m" => [0.0, -0.01]), dS_sem=Dict("m" => [0.0, 0.004]))
    outside = (times=t, S_mean=Dict("m" => [1.0, 0.6]), S_sem=Dict("m" => [0.0, 0.04]),
               dS_mean=Dict("m" => [0.0, -0.02]), dS_sem=Dict("m" => [0.0, 0.003]))
    c = sum_compartments(0.25, inside, outside, "m")
    @test c.S ≈ [1.0, 0.25 * 0.8 + 0.75 * 0.6]
    @test c.S_sem ≈ [0.0, sqrt(0.25^2 * 0.03^2 + 0.75^2 * 0.04^2)]
    @test c.dS ≈ [0.0, 0.25 * -0.01 + 0.75 * -0.02]
    @test c.dS_sem ≈ [0.0, sqrt(0.25^2 * 0.004^2 + 0.75^2 * 0.003^2)]
    bad = merge(outside, (times=[0.0, 2.0],))
    @test_throws ErrorException sum_compartments(0.25, inside, bad, "m")
end

@testset "p3_7b: atomic writes and safe resume (Review Focus 5)" begin
    dir = mktempdir()
    cf, kf = joinpath(dir, "conditions.csv"), joinpath(dir, "curves.csv")
    write_atomic(kf, [(label="a", t=0.0), (label="b", t=0.0)])
    write_atomic(cf, [(label="a", x=1.0)])                      # "b" written to curves only: a crash between the writes
    @test !isfile(cf * ".tmp") && !isfile(kf * ".tmp")
    @test read_done_both(cf, kf; cheader="label,x", kheader="label,t").labels == Set(["a"])
    @test isempty(read_done_both(joinpath(dir, "none.csv"), kf; cheader="label,x", kheader="label,t").labels)
    @test_throws ErrorException read_done_both(cf, kf; cheader="label,y", kheader="label,t")
end

@testset "p3_7b: resume — label only in conditions.csv is not done (plan Task 2 case)" begin
    dir = mktempdir()
    cf, kf = joinpath(dir, "conditions.csv"), joinpath(dir, "curves.csv")
    write_atomic(kf, [(label="a", t=0.0)])
    write_atomic(cf, [(label="a", x=1.0), (label="b", x=2.0)])
    @test read_done_both(cf, kf; cheader="label,x", kheader="label,t").labels == Set(["a"])
end
