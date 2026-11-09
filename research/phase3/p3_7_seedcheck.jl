# Phase 3.7a fresh-seed check (Group B review I1): the fig9a reference R1.0_c0.1 (R = 1, D = 3, ρ = 1/140, bulk 1/2000,
# same centroid c̄ = 0.1, τ = 1.35e-3, readouts 0:TE*/10:1.5 TE*) rerun with seeds 11–20, independent of the seeds 1–10
# shared by every fig9a/params condition. Question: is ΔS_mc/ΔS_cont ≈ 0.93 a Monte Carlo bias or shared seed noise?
# Output: results/phase3/p3_7_seedcheck/{conditions,curves}.csv (p3_7_run.jl format; analysed with at_te_rows).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_7_seedcheck.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_mc.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
const OUT = phase3_outdir("p3_7_seedcheck")
const SEEDS = 11:20
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
R, D, ρ, bulk, cb = 1.0, 3.0, 1 / 140, 1 / 2000, 0.1
models = vcat([("conventional", "conventional", 0.0)], [("$(s)_c$(cb)", s, h_for_cbar(s, cb, R)) for s in ("step", "linear", "exponential")])
τ = 0.135 * cb^2
te = te_star(t -> radial_signal(R=R, D=D, times=[t], R2_bulk=bulk, surface=ρ)[1]; t_hi=50 / (bulk + 2ρ / R))
times = collect(range(0.0, 1.5te, length=16))
t0 = time()
r = run_condition(R=R, D=D, rho=ρ, R2_bulk=bulk, models=models, times=times, tau=τ, seeds=SEEDS, nspins=NSPINS)
label = "R1.0_c0.1_seeds11-20"
conds = [(task="seedcheck", label=label, R=R, D=D, rho=ρ, R2_bulk=bulk, tau=r.tau, cbar_min=cb, te=te, t_end=times[end],
    identical_trajectories=r.identical_trajectories, runtime_s=time() - t0, seeds="11-20")]
curves = NamedTuple[]
for (name, shape, h) in models, (k, t) in enumerate(times)
    haslayer = shape != "conventional"
    push!(curves, (label=label, model=name, shape=shape, h=h, cbar=haslayer ? cbar(shape, h, R) : 0.0, t=t,
        S_mean=r.S_mean[name][k], S_sem=r.S_sem[name][k], S_cont=r.S_cont[name][k],
        dS_mean=haslayer ? r.dS_mean[name][k] : NaN, dS_sem=haslayer ? r.dS_sem[name][k] : NaN, dS_cont=haslayer ? r.dS_cont[name][k] : NaN))
end
write_csv(joinpath(OUT, "conditions.csv"), conds)
write_csv(joinpath(OUT, "curves.csv"), curves)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "seeds" => collect(SEEDS), "all_identical" => r.identical_trajectories))
for row in dS_ratio_rows(at_te_rows(OUT))
    @printf("%-20s ΔS_mc/ΔS_cont %.3f ± %.3f\n", row.model, row.ratio, row.ratio_sigma)
end
println("saved to $OUT")
