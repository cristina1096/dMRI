# P3.1.6: does the conventional model's relaxivity ρ = θ√(D/π) hold on curved surfaces? The per-collision factor is
# e^(−θ√τ), so −ln M⊥ = θ√τ·(number of collisions) for every path, and the test is the collision rate:
#     E[−ln M⊥(t)] = θ√(D/π)·(S/V)·t        (exact for walls, P2.5.0)
# On a cylinder the finite step √(2Dτ) can bias the collision rate by O(√(Dτ)/R). Measured at t = 50 ms for single
# cylinders (inside) R = 0.5, 1, 2 and the square array (inside, outside); D = 1, 3, 6; τ = 1e-2 and 2.5e-3 ms;
# θ = 0.01; seeds 1–10, NSPINS (default 20000). z = (mean − expected)/SEM; |z| ≤ 3 passes, otherwise the measured
# ratio is recorded as the mapping (tracker P3.1.6).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_6_mapping.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_1_6_mapping")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const θ = 0.01
const T = 50.0

configs = [("single_R0.5", 0.5, Inf, :inside), ("single_R1", 1.0, Inf, :inside), ("single_R2", 2.0, Inf, :inside),
           ("array_inside", R_REF, S_REF, :inside), ("array_outside", R_REF, S_REF, :outside)]
rows = NamedTuple[]
for (name, R, s, region) in configs, D in (1.0, 3.0, 6.0), τ in (1e-2, 2.5e-3)
    (name != "single_R1" && name != "array_outside" && D != 3.0) && continue     # D scan on two configurations
    g = isfinite(s) ? Cylinders(radius=R, repeats=[s, s], surface_relaxation=θ) : Cylinders(radius=R, surface_relaxation=θ)
    bb = isfinite(s) ? 500 : mr.BoundingBox([-R, -R, -R], [R, R, R])
    r = run_cylinders(g; region=region, R=R, sx=s, sy=s, bounding_box=bb, timestep=τ, D=D, times=[0.0, T], nspins=NSPINS)
    sv = region == :inside ? sv_inside(R) : sv_outside(R, s, s)
    expct = expected_neglog_conventional(θ, D, sv, T)
    m, e = r.neglog_mean[end], r.neglog_sem[end]
    z = (m - expct) / e
    push!(rows, (config=name, R=R, D=D, tau=τ, sqrt_Dtau_over_R=sqrt(D * τ) / R, expected=expct, measured=m, sem=e, z=z, ratio=m / expct))
    @printf("%-14s D = %.0f τ = %.1e: E[−ln M] %.5f vs %.5f (ratio %.4f, z = %+.2f)\n", name, D, τ, m, expct, m / expct, z)
end
write_csv(joinpath(OUT, "mapping.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_pass" => count(r -> abs(r.z) <= 3, rows),
    "n" => length(rows), "max_abs_z" => maximum(abs(r.z) for r in rows), "theta" => θ, "nspins" => NSPINS))
fig = Figure(size=(600, 420))
ax = Axis(fig[1, 1], xlabel="√(Dτ)/R", ylabel="measured / θ√(D/π)(S/V)t", title="conventional mapping on cylinders")
for (name, mk) in zip(unique(r.config for r in rows), (:circle, :rect, :diamond, :utriangle, :dtriangle))
    q = filter(r -> r.config == name, rows)
    errorbars!(ax, [r.sqrt_Dtau_over_R for r in q], [r.ratio for r in q], [r.sem / r.expected for r in q])
    scatter!(ax, [r.sqrt_Dtau_over_R for r in q], [r.ratio for r in q], marker=mk, label=name)
end
hlines!(ax, [1.0], color=:black)
axislegend(ax, position=:lb, labelsize=9)
save(joinpath(OUT, "p3_1_6_mapping.png"), fig)
println("saved to $OUT")
