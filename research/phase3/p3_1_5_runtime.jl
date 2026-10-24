# P3.1.5 C-B5: runtime per spin-step on cylinders vs τ, without the layer (old code) and with the step layer
# (ρ = 0.02, h = 0.2 both sides). 1 seed, 10000 spins, 1000 steps per run. Single runs on a shared machine: ratios carry
# a load uncertainty of order 10–20 %.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_5_runtime.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_1_5_runtime")
const N = 10_000
rows = NamedTuple[]
for (name, s, region) in (("single", Inf, :inside), ("array", S_REF, :outside)), layer in (false, true), τ in (0.1, 1e-2, 1e-3, 1e-4)
    kw = layer ? (layer_rho=0.02, layer_h=0.2) : (;)
    g = isfinite(s) ? Cylinders(radius=R_REF, repeats=[s, s]; kw...) : Cylinders(radius=R_REF; kw...)
    bb = isfinite(s) ? 500 : mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
    run_cylinders(g; region=region, sx=s, sy=s, bounding_box=bb, timestep=τ, times=[0.0, 10τ], seeds=1:1, nspins=100)   # warm-up
    r = run_cylinders(g; region=region, sx=s, sy=s, bounding_box=bb, timestep=τ, times=[0.0, 1000τ], seeds=1:1, nspins=N)
    us = r.runtime / (N * 1000) * 1e6
    push!(rows, (config=name, layer=layer, tau=τ, us_per_spin_step=us))
    @printf("%-6s layer %-5s τ = %.0e: %.4f µs per spin-step\n", name, layer, τ, us)
end
write_csv(joinpath(OUT, "runtime.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows))
println("saved to $OUT")
