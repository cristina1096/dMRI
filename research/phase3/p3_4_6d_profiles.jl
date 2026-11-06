# P3.4.6 for linear and exponential: Monte Carlo vs the radial continuum inside a single cylinder R = 1, D = 3, ρ = 0.02.
# Models: conventional (θ = ρ√(π/D)), linear h = 0.3, 0.6 (centroids 0.1, 0.2), exponential h = 0.1, 0.2 (support R;
# centroids ≈ 0.1, 0.2). No thin layers (user request 2026-10-25): ΔR2(0) ≤ 0.23 /ms, so τ_a = 1e-2 and τ_b = 2.5e-3.
# All models share τ (identical trajectories, checked), ΔS = S_layer − S_conv paired per seed.
# Criterion (tracker P3.4.6): at τ_b, max over readouts |z| ≤ 3 for S and for ΔS vs the continuum (3, not 2: 11
# correlated readouts). Seeds 1–10, NSPINS (default 20000) in the box [−R, R]³, readouts 0:5:50.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_6d_profiles.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
const OUT = phase3_outdir("p3_4_6d_profiles")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO = 0.02
const θ = RHO * sqrt(π / D_REF)
const TAUS = [1e-2, 2.5e-3]
bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
models = [("conventional", Cylinders(radius=R_REF, surface_relaxation=θ), nothing),
          ("linear_h0.3", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.3, layer_shape="linear"), ("linear", 0.3)),
          ("linear_h0.6", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.6, layer_shape="linear"), ("linear", 0.6)),
          ("exponential_h0.1", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.1, layer_shape="exponential"), ("exponential", 0.1)),
          ("exponential_h0.2", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.2, layer_shape="exponential"), ("exponential", 0.2))]
cont = Dict(n => (isnothing(p) ? radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, surface=RHO) :
                                 radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, layer=(p[1], RHO, p[2]))) for (n, _, p) in models)
runs = Dict((n, τ) => run_cylinders(g; region=:inside, bounding_box=bb, timestep=τ, nspins=NSPINS) for τ in TAUS for (n, g, _) in models)
identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[("conventional", τ)].finals, runs[(n, τ)].finals))
                for τ in TAUS for (n, _, _) in models)
@printf("identical trajectories across models (per τ): %s\n", identical)
rows = NamedTuple[]
for τ in TAUS, (n, _, p) in models
    r = runs[(n, τ)]
    zS = (r.mean .- cont[n]) ./ r.sem
    d = r.S .- runs[("conventional", τ)].S
    dm = vec(mean(d, dims=1)); ds = vec(std(d, dims=1)) ./ sqrt(size(d, 1))
    dc = cont[n] .- cont["conventional"]
    zD = isnothing(p) ? zeros(length(dm)) : (dm .- dc) ./ ds
    push!(rows, (model=n, tau=τ, S50_mc=r.mean[end], S50_sem=r.sem[end], S50_cont=cont[n][end], maxabs_zS=maximum(abs.(zS[2:end])),
        dS50_mc=dm[end], dS50_sem=ds[end], dS50_cont=dc[end], maxabs_zdS=maximum(abs.(zD[2:end]))))
    @printf("%-18s τ = %.1e: S(50) %.6f ± %.1e vs %.6f (max|z| %.2f); ΔS(50) %+.2e ± %.1e vs %+.2e (max|z| %.2f)\n",
        n, τ, r.mean[end], r.sem[end], cont[n][end], maximum(abs.(zS[2:end])), dm[end], ds[end], dc[end], maximum(abs.(zD[2:end])))
end
final = filter(r -> r.tau == TAUS[end], rows)
pass = identical && all(r -> r.maxabs_zS <= 3 && r.maxabs_zdS <= 3, final)
write_csv(joinpath(OUT, "comparison.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "pass" => pass, "identical_trajectories" => identical,
    "taus" => TAUS, "nspins" => NSPINS, "criterion" => "tau_b, max|z| <= 3 for S and paired dS"))
println(pass ? "PASS" : "FAIL")
