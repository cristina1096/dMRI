# P3.4.1: single-spin trace on a square array with overlapping outer layers. Each step's M⊥ change is predicted from
# the simulator's own reflected path, the layer exponent of every straight piece by brute-force sampling of the summed
# rate (brute_rate, 10⁴ points per piece, independent of the package's Layers code):
#     M⊥,pred = M⊥(previous) · exp(−R2_bulk τ − Σ_j dt_j · mean_k brute_rate(p_j(s_k))),  dt_j = τ ℓ_j/Σℓ.
# Square array s = 2.5, R = 1, step layer h = 0.3 both sides (outer shells overlap: gap 0.5 < 2h), ρ = 0.03,
# R2_bulk = 1/80, D = 3, τ = 0.01, 100 spins × 500 steps from uniform positions in a cell, Random.seed!(20261024).
# Pass: max |M⊥,sim/M⊥,pred − 1| ≤ 1e-6, ≥ 1 step with a reflection inside a shell, ≥ 1 step through an overlap.
# Sampling error per boundary crossing ≤ rate·τ/(2·10⁴) ≈ 6e-8.
# Run: julia --project=research/baseline -t 1 research/phase3/p3_4_1_trace.jl

include(joinpath(@__DIR__, "harness.jl"))
using StaticArrays, LinearAlgebra
const OUT = phase3_outdir("p3_4_1_trace")
const TAU, NSTEPS, NTRACE, H, RHO, R2B, TOL, K = 0.01, 500, 100, 0.3, 0.03, 1 / 80, 1e-6, 10_000
const RIN, ROUT = rate_inside(RHO, H, R_REF), rate_outside(RHO, H, R_REF)
rate(x, y) = brute_rate(x, y; R=R_REF, s=S_REF, h_in=H, h_out=H, rate_in=RIN, rate_out=ROUT)

seq = mr.SequenceParts.SequenceWaveform((([], []), ([], []), ([], [])), [], [], [NSTEPS * TAU], NSTEPS * TAU)
sim = Simulation(seq; geometry=Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_rho=RHO, layer_h=H), diffusivity=D_REF, R2=R2B, verbose=false)
part = first(p for p in mr.parts([seq], 0., mr.TimeStep(TAU, Inf)) if p.duration > 0)
@assert part.duration ≈ TAU

Random.seed!(20261024)
max_err, n_reflect_shell, n_overlap = 0.0, 0, 0
rows = NamedTuple[]
for i in 1:NTRACE
    spin = Spin(position=[(rand() - 0.5) * S_REF, (rand() - 0.5) * S_REF, rand()], transverse=1.0, longitudinal=0.0)
    for step in 1:NSTEPS
        m0 = spin.orientations[1].transverse
        proposed = spin.position .+ randn(SVector{3, Float64}) .* sqrt(2 * D_REF * TAU)
        path = mr.Evolve.draw_step!(spin, sim, part, [3.], proposed)
        lens = [norm(path[k + 1] - path[k]) for k in 1:length(path) - 1]
        total = sum(lens)
        expo, overlap = 0.0, false
        for k in eachindex(lens)
            iszero(lens[k]) && continue
            rates = [rate((path[k] + (j - 0.5) / K * (path[k + 1] - path[k]))[1:2]...) for j in 1:K]
            overlap |= any(r -> r > ROUT * 1.5 && r < RIN, rates) || any(r -> r ≈ 2ROUT, rates)
            expo += TAU * lens[k] / total * sum(rates) / K
        end
        pred = m0 * exp(-expo) * exp(-R2B * TAU)
        err = abs(spin.orientations[1].transverse / pred - 1)
        global max_err = max(max_err, err)
        (length(path) > 2 && expo > 0) && (global n_reflect_shell += 1)
        overlap && (global n_overlap += 1)
        i == 1 && push!(rows, (step=step, n_pieces=length(lens), layer_exponent=expo, transverse=spin.orientations[1].transverse, rel_error=err))
    end
end
pass = max_err <= TOL && n_reflect_shell > 0 && n_overlap > 0
@printf("max error %.2e (tol %.0e); steps with a reflection in a shell %d; through an overlap %d → %s\n",
    max_err, TOL, n_reflect_shell, n_overlap, pass ? "PASS" : "FAIL")
write_csv(joinpath(OUT, "trace_spin1.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "max_rel_error" => max_err, "tolerance" => TOL,
    "n_steps_reflection_in_shell" => n_reflect_shell, "n_steps_overlap" => n_overlap, "pass" => pass,
    "params" => Dict("tau" => TAU, "nsteps" => NSTEPS, "nspins" => NTRACE, "h" => H, "rho" => RHO, "R" => R_REF, "s" => S_REF,
        "rate_in" => RIN, "rate_out" => ROUT, "samples_per_piece" => K)))
println("saved to $OUT")
