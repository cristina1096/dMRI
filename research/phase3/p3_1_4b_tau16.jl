# P3.1.4 follow-up (plan Task 8 rule): the one C-B4 curve with |z| > 2 between τ_a and τ_a/4 (array outside, θ = 0.03,
# z = +2.06) is rerun at τ_a/16 = 6.25e-4 ms; if τ_a/4 and τ_a/16 agree within 2 SEM, τ_a/4 is the reference for it.
# Same settings as p3_1_4_theta_reference.jl (seeds 1–10, NSPINS = 20000, readouts 0:5:50).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_4b_tau16.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_1_4_theta_reference")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
lines = readlines(joinpath(OUT, "reference_curves.csv"))
hdr = split(lines[1], ',')
rows = [Dict(zip(hdr, split(l, ','))) for l in lines[2:end]]
b = only(r for r in rows if r["config"] == "array_outside" && parse(Float64, r["theta"]) == 0.03 &&
         parse(Float64, r["tau"]) == 2.5e-3 && parse(Float64, r["t"]) == 50.0)
g = Cylinders(radius=R_REF, repeats=[S_REF, S_REF], surface_relaxation=0.03)
r = run_cylinders(g; region=:outside, sx=S_REF, sy=S_REF, timestep=6.25e-4, nspins=NSPINS)
mb, eb = parse(Float64, b["mean"]), parse(Float64, b["sem"])
z = (mb - r.mean[end]) / sqrt(eb^2 + r.sem[end]^2)
@printf("array_outside θ = 0.03: S(50) %.5f (τ 2.5e-3) vs %.5f ± %.1e (τ 6.25e-4), z = %+.2f\n", mb, r.mean[end], r.sem[end], z)
write_json(joinpath(OUT, "tau16_check.json"), Dict("provenance" => provenance(), "S50_tau_a4" => mb, "sem_tau_a4" => eb,
    "S50_tau_a16" => r.mean[end], "sem_tau_a16" => r.sem[end], "z" => z, "pass" => abs(z) <= 2,
    "curve_tau_a16" => Dict("t" => r.times, "mean" => r.mean, "sem" => r.sem)))
