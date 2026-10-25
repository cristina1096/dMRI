# P3.4.8 outside, follow-up: third timestep τ_c = τ_a/16 = 6.25e-4 ms for the thin end of the sweep, because P3.4.6
# showed that at fixed ρ the paired difference of thin layers (h ≤ 0.1) is timestep-biased at τ_a and τ_b (Review C1).
# For each array (s = 2.5, 3.0): conventional and outer-only step layers h = 0.05, 0.1, 0.2, ρ = 0.02, seeds 1–10,
# NSPINS (default 10000) spins per seed (half the main pass; cost ∝ steps), identical trajectories asserted. Paired
# neighbour differences at τ_c, and their shift from τ_b (read from steps.csv of p3_4_8_h_sweep_outside.jl).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_8b_outside_tau16.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_4_8_h_sweep_outside")
const NSPINS = parse(Int, get(ENV, "NSPINS", "10000"))
const RHO = 0.02
const θ = RHO * sqrt(π / D_REF)
const τc = 1e-2 / 16
const HS = [0.05, 0.1, 0.2]
paired(a, b) = (d = a.S[:, end] .- b.S[:, end]; (mean(d), std(d) / sqrt(length(d))))
main = [Dict(zip(split(readlines(joinpath(OUT, "steps.csv"))[1], ','), split(l, ','))) for l in readlines(joinpath(OUT, "steps.csv"))[2:end]]
rows = NamedTuple[]
allidentical = true
for s in (2.5, 3.0)
    models = vcat([("conventional", Cylinders(radius=R_REF, repeats=[s, s], surface_relaxation=θ))],
                  [("step_h$(h)", Cylinders(radius=R_REF, repeats=[s, s], layer_rho_outer=RHO, layer_h_outer=h)) for h in HS])
    runs = Dict(name => run_cylinders(g; region=:outside, sx=s, sy=s, timestep=τc, nspins=NSPINS) for (name, g) in models)
    identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs["conventional"].finals, runs[n].finals)) for (n, _) in models)
    global allidentical &= identical
    names = [n for (n, _) in models]
    for k in 1:length(names) - 1
        (mc, sc) = paired(runs[names[k + 1]], runs[names[k]])
        b = only(r for r in main if parse(Float64, r["s"]) == s && r["from"] == names[k] && r["to"] == names[k + 1])
        mb, sb = parse(Float64, b["d_tau_b"]), parse(Float64, b["sem_tau_b"])
        push!(rows, (s=s, from=names[k], to=names[k + 1], d_tau_c=mc, sem_tau_c=sc, z_c=mc / sc, resolved=abs(mc / sc) > 3,
            sign=sign(mc), d_tau_b=mb, sem_tau_b=sb, z_shift_b_c=(mb - mc) / sqrt(sb^2 + sc^2)))
        @printf("s = %.1f %-14s → %-12s: ΔS(50) τ_c %+.2e ± %.1e (z %+.1f); τ_b %+.2e ± %.1e; shift z %+.2f\n",
            s, names[k], names[k + 1], mc, sc, mc / sc, mb, sb, (mb - mc) / sqrt(sb^2 + sc^2))
    end
end
write_csv(joinpath(OUT, "steps_tau16.csv"), rows)
write_json(joinpath(OUT, "summary_tau16.json"), Dict("provenance" => provenance(), "tau_c" => τc, "nspins" => NSPINS,
    "identical_trajectories" => allidentical, "max_abs_shift_z" => maximum(abs(r.z_shift_b_c) for r in rows)))
println("saved to $OUT")
