# P3.1.4 C-B4: conventional surface relaxation (MCMR θ) reference curves on cylinders. A reference, not the truth.
# Configurations: single cylinder inside R = 0.5, 1, 2; square array s = 2.5, R = 1 inside and outside.
# θ = 0.003, 0.01, 0.03; D = 3; readouts 0:5:50; seeds 1–10, NSPINS (default 20000) spins per seed.
# Timesteps τ_a = min(1e-2, 0.03R²/D) and τ_b = τ_a/4. Convergence: |S_a(50) − S_b(50)|/√(SEM_a² + SEM_b²) ≤ 2
# (one test per curve; the 11 readouts are correlated, so they are reported, not each graded).
# Bulk factorisation (θ = 0.01, single R = 1, τ_a): same seeds with R2 = 1/80 → per spin ratio e^(−t/80) within 1e-10.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_4_theta_reference.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_1_4_theta_reference")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const THETAS = [0.003, 0.01, 0.03]

configs = [
    ("single_R0.5", 0.5, Inf, :inside), ("single_R1", 1.0, Inf, :inside), ("single_R2", 2.0, Inf, :inside),
    ("array_inside", R_REF, S_REF, :inside), ("array_outside", R_REF, S_REF, :outside),
]
geom(R, s, θ) = isfinite(s) ? Cylinders(radius=R, repeats=[s, s], surface_relaxation=θ) : Cylinders(radius=R, surface_relaxation=θ)
bbox(R, s) = isfinite(s) ? 500 : mr.BoundingBox([-R, -R, -R], [R, R, R])

rows = NamedTuple[]
conv = NamedTuple[]
for (name, R, s, region) in configs, θ in THETAS
    τa = min(1e-2, 0.03 * R^2 / D_REF)
    for τ in (τa, τa / 4)
        r = run_cylinders(geom(R, s, θ); region=region, R=R, sx=s, sy=s, bounding_box=bbox(R, s), timestep=τ, nspins=NSPINS)
        append!(rows, [(config=name, R=R, s=s, region=String(region), theta=θ, tau=τ, t=t, mean=m, sem=e)
                       for (t, m, e) in zip(r.times, r.mean, r.sem)])
    end
    a = filter(x -> x.config == name && x.theta == θ && x.tau == τa && x.t == 50.0, rows)[1]
    b = filter(x -> x.config == name && x.theta == θ && x.tau == τa / 4 && x.t == 50.0, rows)[1]
    z = (a.mean - b.mean) / sqrt(a.sem^2 + b.sem^2)
    push!(conv, (config=name, theta=θ, tau_a=τa, S50_a=a.mean, S50_b=b.mean, z=z, pass=abs(z) <= 2))
    @printf("%-14s θ = %.3f: S(50) %.5f (τ %.2e) vs %.5f (τ/4), z = %+.2f\n", name, θ, a.mean, τa, b.mean, z)
end
# bulk factorisation
g = geom(1.0, Inf, 0.01)
kw = (region=:inside, bounding_box=bbox(1.0, Inf), timestep=1e-2, nspins=NSPINS, seeds=1:3)
a = run_cylinders(g; kw...)
b = run_cylinders(g; R2_bulk=1 / 80, kw...)
fac = maximum(maximum(abs(sb.orientations[1].transverse / (sa.orientations[1].transverse * exp(-50 / 80)) - 1)
                      for (sa, sb) in zip(fa.spins, fb.spins)) for (fa, fb) in zip(a.finals, b.finals))
@printf("bulk factorisation: max per-spin deviation %.1e\n", fac)
write_csv(joinpath(OUT, "reference_curves.csv"), rows)
write_csv(joinpath(OUT, "convergence.csv"), conv)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_converged" => count(c -> c.pass, conv),
    "n_curves" => length(conv), "bulk_factorisation_max_dev" => fac, "nspins" => NSPINS))
println("saved to $OUT")
