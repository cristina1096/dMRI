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

@testset "p3_7 resume: read_done / write_merged (Group A review M3)" begin
    dir = mktempdir()
    f = joinpath(dir, "conditions.csv")
    @test read_done(f).labels == Set{String}()                 # no file yet: nothing done
    write_csv(f, [(label="a", x=1.0), (label="b", x=2.0)])
    prev = read_done(f)
    @test prev.labels == Set(["a", "b"])
    write_merged(f, prev, [(label="c", x=3.0)])                # previous rows kept, new row appended
    @test readlines(f) == ["label,x", "a,1", "b,2", "c,3"]
    write_merged(f, read_done(f), NamedTuple[])                # nothing new: file unchanged
    @test readlines(f) == ["label,x", "a,1", "b,2", "c,3"]
    @test_throws ErrorException write_merged(f, read_done(f), [(label="d", y=1.0)])   # header mismatch
end
