# P3.1.5 follow-up (final review finding): runtime per spin-step with the step layer vs the number of cylinders per
# unit cell. The layer exponent loops over every cylinder of the group on every straight piece, so its cost grows with
# the total number of cylinders per cell, not only with those within reach. n × n cylinders (R = 0.7, centre spacing 2,
# gap 0.6) in a cell of side 2n, so the packing density is the same for every n; layer ρ = 0.02, h = 0.1 both sides.
# 2000 spins, 1000 steps of τ = 1e-2 ms, 1 seed; single runs on a shared machine (±10–20 %).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_5b_runtime_many.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_1_5_runtime")
const N, τ = 2000, 1e-2
rows = NamedTuple[]
for n in (1, 3, 10, 20)
    L = 2.0 * n
    pos = [[-L / 2 + 1 + 2i, -L / 2 + 1 + 2j] for i in 0:n-1 for j in 0:n-1]
    t = Dict{Bool, Float64}()
    for layer in (false, true)
        kw = layer ? (layer_rho=0.02, layer_h=0.1) : (;)
        g = Cylinders(radius=0.7, position=pos, repeats=[L, L]; kw...)
        run_cylinders(g; region=:outside, R=0.7, sx=L, sy=L, timestep=τ, times=[0.0, 10τ], seeds=1:1, nspins=100)   # warm-up
        r = run_cylinders(g; region=:outside, R=0.7, sx=L, sy=L, timestep=τ, times=[0.0, 1000τ], seeds=1:1, nspins=N)
        t[layer] = r.runtime / (N * 1000) * 1e6
    end
    push!(rows, (cylinders_per_cell=n^2, us_no_layer=t[false], us_layer=t[true], ratio=t[true] / t[false]))
    @printf("%4d cylinders per cell: %.3f µs without, %.3f µs with the layer (×%.1f)\n", n^2, t[false], t[true], t[true] / t[false])
end
write_csv(joinpath(OUT, "runtime_many.csv"), rows)
write_json(joinpath(OUT, "summary_many.json"), Dict("provenance" => provenance(), "rows" => rows))
println("saved to $OUT")
