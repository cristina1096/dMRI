@testset "p3_7 paired MC runner" begin
    models = [("conventional", "conventional", 0.0), ("step_h0.2", "step", 0.2)]
    r = run_condition(R=1.0, D=3.0, rho=0.02, R2_bulk=0.0, models=models, times=[0.0, 10.0], tau=1e-2, nspins=500, seeds=1:2)
    @test r.times == [0.0, 10.0]
    @test r.identical_trajectories
    @test all(r.S_mean[m][1] == 1.0 for (m, _, _) in models)
    @test size(r.S["step_h0.2"]) == (2, 2)
    @test r.S_cont["step_h0.2"] ≈ radial_signal(R=1.0, D=3.0, times=[0.0, 10.0], layer=("step", 0.02, 0.2)) rtol=1e-12
    @test r.S_cont["conventional"] ≈ radial_signal(R=1.0, D=3.0, times=[0.0, 10.0], surface=0.02) rtol=1e-12
    @test r.dS_mean["step_h0.2"][1] == 0.0
    @test r.dS_cont["step_h0.2"] ≈ r.S_cont["step_h0.2"] .- r.S_cont["conventional"]
    @test r.S_mean["step_h0.2"][2] < 1.0
    @test r.tau == 1e-2
end
