# P3.4.6 follow-up (plan Task 2 rule): the thin layers failed the τ_b = 2.5e-3 check (paired ΔS vs continuum: h = 0.05
# z 6.6, h = 0.1 z 4.5; h ≥ 0.2 pass) — the timestep bias of the fixed-ρ difference found on walls in Review C1. Rerun the
# conventional model and h = 0.05, 0.1 at τ_c = τ_a/16 = 6.25e-4 ms (same seeds, spins, readouts; identical trajectories
# across the three models, asserted) → comparison_tau16.csv. The verdict over both passes is p3_4_6c_verdict.jl.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_6b_tau16.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
const OUT = phase3_outdir("p3_4_6_continuum_check")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO = 0.02
const θ = RHO * sqrt(π / D_REF)
const τc = 1e-2 / 16
bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
models = [("conventional", Cylinders(radius=R_REF, surface_relaxation=θ), nothing),
          ("step_h0.05", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.05), 0.05),
          ("step_h0.1", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.1), 0.1)]
cont = Dict(name => (isnothing(h) ? radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, surface=RHO) :
                                    radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, layer=("step", RHO, h))) for (name, _, h) in models)
runs = Dict(name => run_cylinders(g; region=:inside, bounding_box=bb, timestep=τc, nspins=NSPINS) for (name, g, _) in models)
identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs["conventional"].finals, runs[name].finals)) for (name, _, _) in models)
@printf("identical trajectories: %s\n", identical)
rows = NamedTuple[]
for (name, _, h) in models
    r = runs[name]
    zS = (r.mean .- cont[name]) ./ r.sem
    if isnothing(h)
        zD = zeros(length(TIMES_REF)); dm = ds = dc = 0.0
    else
        d = r.S .- runs["conventional"].S
        dmean = vec(mean(d, dims=1)); dsem = vec(std(d, dims=1)) ./ sqrt(size(d, 1))
        dcont = cont[name] .- cont["conventional"]
        zD = (dmean .- dcont) ./ dsem
        dm, ds, dc = dmean[end], dsem[end], dcont[end]
    end
    push!(rows, (model=name, h=isnothing(h) ? 0.0 : h, tau=τc, S50_mc=r.mean[end], S50_sem=r.sem[end], S50_cont=cont[name][end],
        maxabs_zS=maximum(abs.(zS[2:end])), dS50_mc=dm, dS50_sem=ds, dS50_cont=dc, maxabs_zdS=maximum(abs.(zD[2:end])),
        pass=maximum(abs.(zS[2:end])) <= 3 && maximum(abs.(zD[2:end])) <= 3))
    @printf("%-14s τ = %.2e: S(50) MC %.6f ± %.1e vs %.6f (max|z| %.2f); ΔS(50) MC %+.2e ± %.1e vs %+.2e (max|z| %.2f)\n",
        name, τc, r.mean[end], r.sem[end], cont[name][end], maximum(abs.(zS[2:end])), dm, ds, dc, maximum(abs.(zD[2:end])))
end
write_csv(joinpath(OUT, "comparison_tau16.csv"), rows)
println("saved to $OUT (verdict: p3_4_6c_verdict.jl)")
