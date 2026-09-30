# B1 — Phase 2 reference configuration, null run (target for P2.3.1 V0b and P2.2.8)
#
# The Phase 2 reference configuration with nothing that relaxes. The unmodified simulator must
# return M⊥(t)/M⊥(0) = 1 and phase = 0 at every readout, to floating-point precision. The run
# also records the τ_max the simulator picks and the runtime per spin per ms of sequence, which
# P2.2.8 (λ = 0 short-circuit: no slowdown) compares against.
#
# Variables:
#   geometry : Walls(repeats=w), w = 2 µm, normal along x (both faces, no surface parameters)
#   D = 3 µm²/ms; R1 = R2 = 0; θ_relax = 0; permeability = 0; off-resonance = 0; no MT/sticking
#   sequence : none (empty sequence, no RF, no gradients). Spins start with transverse = 1,
#              longitudinal = 0, phase = 0 (Snapshot keywords), i.e. free decay from t = 0
#   readouts : t = 0, 10, 20, 30, 40, 50 ms
#   timestep : simulator default (printed and recorded); also a fixed τ = 1e-3 ms for runtime scaling
#   N_spins  : NSPINS (default 100_000) per seed, uniform in a ±500 µm box; seeds 1…NSEEDS (default 10)
#
# Pass: |M⊥(t)/M⊥(0) − 1| ≤ 1e-12 and |phase| ≤ 1e-9° at every readout, every seed.
#
# Run: julia --project=research/baseline -t 8 research/baseline/b1_walls_null.jl

include(joinpath(@__DIR__, "common.jl"))

const OUT = outdir("b1_walls_null")
const NSPINS = parse(Int, get(ENV, "NSPINS", "100000"))
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const NSPINS_SMALL_TAU = parse(Int, get(ENV, "NSPINS_SMALL_TAU", "10000"))
const D = 3.0
const W = 2.0
const TIMES = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]
const TOL_SIGNAL = 1e-12
const TOL_PHASE = 1e-9

geometry = Walls(repeats=W)

function make_sim(timestep)
    kw = isnothing(timestep) ? (;) : (; timestep=timestep)
    Simulation([mr.SequenceParts.empty_sequence()]; geometry=geometry, diffusivity=D, verbose=false, kw...)
end

function run_null(label, sim, nspins, nseeds, io)
    ts = sim.timestep.max_timestep
    tee(io, "\n== $label: τ_max = $ts ms, N_spins = $nspins, seeds = $nseeds")
    rows = NamedTuple[]
    max_dev_signal = 0.0
    max_dev_phase = 0.0
    runtimes = Float64[]
    for seed in 1:nseeds
        Random.seed!(seed)
        snap = Snapshot(nspins, sim, 500; transverse=1.0, longitudinal=0.0)
        t0 = time()
        res = readout(snap, sim, TIMES)
        dt = time() - t0
        push!(runtimes, dt)
        for (t, r) in zip(TIMES, res)
            s = per_spin(r)
            ph = mr.phase(r)
            max_dev_signal = max(max_dev_signal, abs(s - 1))
            max_dev_phase = max(max_dev_phase, abs(ph))
            push!(rows, (config=label, seed=seed, t_ms=t, signal=s, phase_deg=ph, nspins=r.nspins))
        end
        @printf(io, "  seed %2d: M⊥(50)/M⊥(0) − 1 = %.2e   phase(50) = %.2e°   runtime %.2f s\n", seed, per_spin(res[end]) - 1, mr.phase(res[end]), dt)
        @printf("  seed %2d: runtime %.2f s\n", seed, dt)
    end
    # first seed includes compilation; exclude it when more than one seed
    rt = nseeds > 1 ? runtimes[2:end] : runtimes
    us_per_spin_ms = mean(rt) / (nspins * TIMES[end]) * 1e6
    pass = max_dev_signal <= TOL_SIGNAL && max_dev_phase <= TOL_PHASE
    tee(io, @sprintf("  max |M⊥/M⊥0 − 1| = %.2e, max |phase| = %.2e°  → %s", max_dev_signal, max_dev_phase, pass ? "PASS" : "FAIL"))
    tee(io, @sprintf("  runtime: %.2f s per seed (excl. first), %.4f µs per spin per ms of sequence, %.4f µs per spin-step",
        mean(rt), us_per_spin_ms, us_per_spin_ms * ts))
    write_csv(joinpath(OUT, "signal_$(label).csv"), rows)
    return Dict(
        "config" => label, "tau_max_ms" => ts, "nspins_per_seed" => nspins, "nseeds" => nseeds,
        "readout_times_ms" => TIMES, "max_abs_signal_deviation" => max_dev_signal,
        "max_abs_phase_deg" => max_dev_phase, "pass" => pass,
        "runtime_per_seed_s" => runtimes, "runtime_mean_excl_first_s" => mean(rt),
        "us_per_spin_per_ms" => us_per_spin_ms, "us_per_spin_step" => us_per_spin_ms * ts,
    )
end

results = Any[]
open(joinpath(OUT, "run.log"), "w") do io
    sim_default = make_sim(nothing)
    tee(io, "size_scale = $(mr.Geometries.Internal.size_scale(sim_default.geometry)) µm; default timestep = $(sim_default.timestep)")
    push!(results, run_null("default_tau", sim_default, NSPINS, NSEEDS, io))
    push!(results, run_null("tau_1e-3", make_sim(1e-3), NSPINS_SMALL_TAU, 3, io))
end

write_json(joinpath(OUT, "summary.json"), Dict(
    "provenance" => provenance(), "results" => results,
    "geometry" => "Walls(repeats=$W): planes x = k·$W, both faces, no surface parameters",
    "physics" => "D=$D, R1=R2=0, θ_relax=0, permeability=0, off-resonance=0",
    "sequence" => "empty_sequence() — no RF, no gradients; spins start transverse=1, longitudinal=0, phase=0",
    "binding_constraint_default" => "tortuosity 0.03·size_scale²/D with size_scale = w",
    "pass_criterion" => "|M⊥/M⊥0 − 1| ≤ $TOL_SIGNAL and |phase| ≤ $(TOL_PHASE)° at every readout",
))
println("saved to $OUT")
