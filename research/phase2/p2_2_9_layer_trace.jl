# P2.2.9 — end-to-end single-spin trace with the layer on (as P1.2.6, now with ΔR2(d)).
#
# Each step: the script draws the proposed displacement, `draw_step!` returns the reflected piecewise path,
# and the step's M⊥ change is predicted independently: time per piece ∝ 3D piece length; layer time per piece
# by brute-force sampling (sampled_layer_fraction); M⊥ *= exp(-ΔR2(0)·Σ dt·fraction)·exp(-R2_bulk·τ).
#
# Variables: walls every 2 µm, step layer on both faces h = 0.5 µm, ΔR2(0) = 0.1 /ms, R2_bulk = 1/80 /ms,
#            D = 3 µm²/ms, τ = 0.01 ms, 100 spins × 500 steps, Random.seed!(20261003).
# Pass: max over all steps of |M⊥_sim / M⊥_pred − 1| ≤ 1e-6, and at least one step with a reflection
#       during which the spin was inside the layer.
#
# Run: julia --project=research/baseline -t 1 research/phase2/p2_2_9_layer_trace.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "p2_2_9_layer_trace_lib.jl"))
using StaticArrays

const SHAPE = get(ENV, "SHAPE", "step")
const OUT = phase2_outdir(SHAPE == "step" ? "p2_2_9_layer_trace" : "p2_2_9_layer_trace_" * SHAPE)
const TAU = 0.01
const NSTEPS = 500
const NTRACE = 100
const H = 0.5
const RATE = 0.1
const R2B = 1 / 80
const TOL = 1e-6

seq = mr.SequenceParts.SequenceWaveform((([], []), ([], []), ([], [])), [], [], [NSTEPS * TAU], NSTEPS * TAU)
sim = Simulation(seq; geometry=layer_walls(h=H, rate=RATE, shape=SHAPE), diffusivity=D_REF, R2=R2B, verbose=false)
const RHO_OVER_H = rho_for(RATE, H; shape=SHAPE) / H
# the first part is the zero-length instant at t = 0; take the τ-long step after it
part = first(p for p in mr.parts([seq], 0., mr.TimeStep(TAU, Inf)) if p.duration > 0)
@assert part.duration ≈ TAU

Random.seed!(20261003)
rows = NamedTuple[]
max_err = 0.0
n_reflect_in_layer = 0
for i in 1:NTRACE
    spin = Spin(position=[rand() * W_REF, rand(), rand()], transverse=1.0, longitudinal=0.0)
    for step in 1:NSTEPS
        m0 = spin.orientations[1].transverse
        proposed = spin.position .+ randn(SVector{3, Float64}) .* sqrt(2 * D_REF * TAU)
        path = mr.Evolve.draw_step!(spin, sim, part, [3.], proposed)
        lens = [norm(path[k + 1] - path[k]) for k in 1:length(path) - 1]
        total_len = sum(lens)
        expo = 0.0
        for k in eachindex(lens)
            iszero(lens[k]) && continue
            dt = TAU * lens[k] / total_len
            expo += RHO_OVER_H * dt * sampled_profile_average(path[k][1], path[k + 1][1], H; shape=Symbol(SHAPE))
        end
        predicted = m0 * exp(-expo) * exp(-R2B * TAU)
        m1 = spin.orientations[1].transverse
        err = abs(m1 / predicted - 1)
        global max_err = max(max_err, err)
        reflected = length(path) > 2
        reflected && expo > 0 && (global n_reflect_in_layer += 1)
        if i == 1
            push!(rows, (step=step, x_start=path[1][1], x_end=path[end][1], n_pieces=length(lens),
                layer_exponent=expo, transverse=m1, predicted=predicted, rel_error=err))
        end
    end
end
pass = max_err <= TOL && n_reflect_in_layer > 0
@printf("max |M_sim/M_pred - 1| = %.2e (tol %.0e); steps with a reflection inside the layer: %d → %s\n",
    max_err, TOL, n_reflect_in_layer, pass ? "PASS" : "FAIL")
write_csv(joinpath(OUT, "trace_spin1.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(),
    "params" => Dict("tau_ms" => TAU, "nsteps" => NSTEPS, "nspins" => NTRACE, "h_um" => H, "rate_per_ms" => RATE,
        "R2_bulk" => R2B, "D" => D_REF, "w_um" => W_REF, "seed" => 20261003, "samples_per_piece" => 10_000, "shape" => SHAPE),
    "max_rel_error" => max_err, "tolerance" => TOL, "n_steps_with_reflection_in_layer" => n_reflect_in_layer, "pass" => pass))
println("saved to $OUT")
