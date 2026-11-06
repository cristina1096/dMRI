# P3.3.4 / P3.6.3: error of the per-piece profile integral vs the quadrature tolerance, and runtime per spin-step.
# (1) Tolerance: 20000 random pieces (endpoints uniform in [−1.6, 1.6]², R = 1) for linear (h = 0.3) and exponential
#     (h = 0.2, support R inside; h = 0.15, support 0.5 outside); I(tol) = shell_profile_integral at LAYER_QUAD_TOL = tol
#     against tol = 1e-14; reported: max |error| and mean number of g evaluations (counted by a wrapper).
#     The exponent error of a piece is (ρ/h)/κ · dt · |error| ≤ ΔR2(0) τ |error|/g(0), so |error| ≤ 1e-6 keeps the
#     per-step exponent error ≤ 1e-6 ΔR2(0) τ, far below the timestep error ≈ 1e-3 relative in R (P2.6).
# (2) Runtime: MCMR µs per spin-step, square array s = 2.5, R = 1, ρ = 0.02, h = 0.2 both sides, step / linear /
#     exponential, τ = 1e-2, at the adopted tolerance and at 1e-6; 10000 spins × 1000 steps, 1 seed.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_3_quadrature.jl

include(joinpath(@__DIR__, "harness.jl"))
using StaticArrays
const OUT = phase3_outdir("p3_3_quadrature")
Layers = mr.Geometries.Internal.Layers
P(x, y) = SVector{2, Float64}(x, y)
O = P(0.0, 0.0)
Random.seed!(2026)
pieces = [(P((3.2 .* rand(2) .- 1.6)...), P((3.2 .* rand(2) .- 1.6)...)) for _ in 1:20_000]
cases = [("linear_in", Layers.LayerSide(0.02, 0.3, :linear), true), ("linear_out", Layers.LayerSide(0.02, 0.3, :linear), false),
         ("exp_in", Layers.LayerSide(0.02, 0.2, :exponential, 1.0), true), ("exp_out", Layers.LayerSide(0.02, 0.15, :exponential, 0.5), false)]
rows = NamedTuple[]
default_tol = Layers.LAYER_QUAD_TOL[]
for (name, side, inside) in cases
    Layers.LAYER_QUAD_TOL[] = 1e-14
    ref = [Layers.shell_profile_integral(side, 1.0, inside, O, p0, p1) for (p0, p1) in pieces]
    for tol in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12)
        Layers.LAYER_QUAD_TOL[] = tol
        t0 = time_ns()
        v = [Layers.shell_profile_integral(side, 1.0, inside, O, p0, p1) for (p0, p1) in pieces]
        ns = (time_ns() - t0) / length(pieces)
        err = maximum(abs.(v .- ref))
        push!(rows, (case=name, tol=tol, max_abs_error=err, ns_per_piece=ns))
        @printf("%-10s tol %.0e: max |error| %.1e, %.0f ns per piece\n", name, tol, err, ns)
    end
end
Layers.LAYER_QUAD_TOL[] = default_tol
write_csv(joinpath(OUT, "tolerance.csv"), rows)

rt = NamedTuple[]
for shape in ("step", "linear", "exponential"), tol in (default_tol, 1e-6)
    shape == "step" && tol != default_tol && continue
    Layers.LAYER_QUAD_TOL[] = tol
    g = Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_rho=0.02, layer_h=0.2, layer_shape=shape)
    run_cylinders(g; region=:outside, sx=S_REF, sy=S_REF, timestep=1e-2, times=[0.0, 0.1], seeds=1:1, nspins=100)
    r = run_cylinders(g; region=:outside, sx=S_REF, sy=S_REF, timestep=1e-2, times=[0.0, 10.0], seeds=1:1, nspins=10_000)
    us = r.runtime / (10_000 * 1000) * 1e6
    push!(rt, (shape=shape, tol=tol, us_per_spin_step=us))
    @printf("%-11s tol %.0e: %.4f µs per spin-step\n", shape, tol, us)
end
Layers.LAYER_QUAD_TOL[] = default_tol
write_csv(joinpath(OUT, "runtime.csv"), rt)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "default_tol" => default_tol, "tolerance" => rows, "runtime" => rt))
println("saved to $OUT")
