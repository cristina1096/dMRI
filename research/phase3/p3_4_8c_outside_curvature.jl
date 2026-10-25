# P3.4.8 outside, completion (final review I2): the tracker criterion includes the curvature of ln S, which the main
# outside pass did not save. Per seed, least squares ln S(t) = a + b t + c t² over readouts 0:5:50 (c = 0 for a single
# exponential), mean ± SEM over seeds, for the conventional model and every h of both arrays (s = 2.5, 3.0), outer layer
# only, ρ = 0.02, τ = 1e-2, seeds 1–10, 2e4 spins per seed. τ = 1e-2 is enough here: the timestep bias of ΔS (≈ 1e-4 in
# S, nearly constant in t) barely changes the shape of ln S; the curvature is reported as found, no pass/fail.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_8c_outside_curvature.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_4_8_h_sweep_outside")
const RHO = 0.02
const θ = RHO * sqrt(π / D_REF)
const ARRAYS = [(2.5, [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]), (3.0, [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0])]
curv(ts, S) = (X = hcat(ones(length(ts)), ts, ts .^ 2); (X \ log.(S))[3])
rows = NamedTuple[]
for (s, HS) in ARRAYS
    models = vcat([("conventional", Cylinders(radius=R_REF, repeats=[s, s], surface_relaxation=θ))],
                  [("step_h$(h)", Cylinders(radius=R_REF, repeats=[s, s], layer_rho_outer=RHO, layer_h_outer=h)) for h in HS])
    for (name, g) in models
        r = run_cylinders(g; region=:outside, sx=s, sy=s, timestep=1e-2, nspins=20_000)
        c = [curv(r.times, r.S[i, :]) for i in 1:size(r.S, 1)]
        push!(rows, (s=s, model=name, curvature=mean(c), curvature_sem=std(c) / sqrt(length(c)), z=mean(c) / (std(c) / sqrt(length(c)))))
        @printf("s = %.1f %-14s curvature %+.2e ± %.1e ms⁻² (z %+.2f)\n", s, name, mean(c), std(c) / sqrt(length(c)), mean(c) / (std(c) / sqrt(length(c))))
    end
end
write_csv(joinpath(OUT, "curvature.csv"), rows)
write_json(joinpath(OUT, "summary_curvature.json"), Dict("provenance" => provenance(), "rows" => rows,
    "max_abs_curvature" => maximum(abs(r.curvature) for r in rows), "max_abs_z" => maximum(abs(r.z) for r in rows)))
println("saved to $OUT")
