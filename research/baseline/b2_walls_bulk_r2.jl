# B2 — bulk R2 in the wall reference configuration (baseline for the R2_bulk term of R2_total = R2_bulk + ΔR2(d))
#
# The proposed model adds ΔR2(d) to the existing bulk rate. This run fixes how the UNMODIFIED
# simulator applies R2_bulk alone, so that Phase 2 can check the new code keeps it intact:
#   (1) M⊥(t)/M⊥(0) = exp(−R2_bulk·t) and phase = 0 at every readout (floating-point precision)
#   (2) same seed ⇒ spin positions bit-identical to the R2_bulk = 0 run (bulk R2 does not touch
#       the random trajectories), so later bulk-on vs bulk-off comparisons can be paired exactly
#   (3) runtime equal to the R2_bulk = 0 run (B1)
#
# Variables:
#   geometry : Walls(repeats=w), w = 2 µm, normal along x, no surface parameters
#   D = 3 µm²/ms; R1 = 0; R2 (global, = R2_bulk) ∈ R2_VALUES; θ_relax = 0; no permeability/off-resonance/MT
#   sequence : none (empty sequence, no RF); spins start transverse = 1, longitudinal = 0, phase = 0
#   readouts : t = 0, 10, 20, 30, 40, 50 ms
#   timestep : simulator default (0.04 ms, tortuosity; bulk R2 does not enter τ_max)
#   N_spins  : NSPINS (default 100_000) per seed; seeds 1…NSEEDS (default 10)
#
# Pass: per spin |m/exp(−R2·t) − 1| ≤ 1e-13 (seed 1, t = 50 ms); ensemble |S/exp(−R2·t) − 1| ≤ (n_steps + N)·eps
#       (floating-point bound: n_steps repeated multiplications per spin + summation over N spins, ≈ 2.2e-11);
#       |phase| ≤ 1e-9°; positions identical (==) to R2 = 0 at t = 50 ms.
# (A first run used an ensemble tolerance of 1e-12. MCMR's own ensemble summation reaches ~2.5e-12, so that
#  tolerance was below roundoff; the per-spin check shows the relaxation itself is exact to ~2e-14.)
#
# Run: julia --project=research/baseline -t 8 research/baseline/b2_walls_bulk_r2.jl

include(joinpath(@__DIR__, "common.jl"))

const OUT = outdir("b2_walls_bulk_r2")
const NSPINS = parse(Int, get(ENV, "NSPINS", "100000"))
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const D = 3.0
const W = 2.0
const TIMES = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]
const R2_VALUES = [1 / 80, 0.025, 0.05, 0.1]
const TOL_PER_SPIN = 1e-13
const TOL_PHASE = 1e-9

geometry = Walls(repeats=W)
const TOL_SIGNAL = (TIMES[end] / 0.04 + NSPINS) * eps()   # roundoff bound, default τ = 0.04 ms
make_sim(R2) = Simulation([mr.SequenceParts.empty_sequence()]; geometry=geometry, diffusivity=D, R2=R2, verbose=false)
start(seed, sim) = (Random.seed!(seed); Snapshot(NSPINS, sim, 500; transverse=1.0, longitudinal=0.0))

results = Any[]
rows = NamedTuple[]
open(joinpath(OUT, "run.log"), "w") do io
    # reference positions and runtime at R2 = 0 (same as B1)
    sim0 = make_sim(0.0)
    tee(io, "default timestep = $(sim0.timestep)")
    ref_final = readout(start(1, sim0), sim0, [TIMES[end]]; return_snapshot=true)[1]
    ref_pos = mr.position.(ref_final)
    rt0 = Float64[]
    for seed in 1:NSEEDS
        snap = start(seed, sim0)
        t0 = time(); readout(snap, sim0, TIMES); push!(rt0, time() - t0)
    end
    tee(io, @sprintf("R2 = 0: runtime %.2f s per seed (mean over seeds 2…%d)", mean(rt0[2:end]), NSEEDS))

    for R2 in R2_VALUES
        sim = make_sim(R2)
        tee(io, @sprintf("\n== R2_bulk = %.6g ms⁻¹ (T2 = %.1f ms), τ_max = %g ms", R2, 1 / R2, sim.timestep.max_timestep))
        max_dev = 0.0
        max_ph = 0.0
        rts = Float64[]
        for seed in 1:NSEEDS
            snap = start(seed, sim)
            t0 = time()
            res = readout(snap, sim, TIMES)
            push!(rts, time() - t0)
            for (t, r) in zip(TIMES, res)
                s = per_spin(r)
                dev = s / exp(-R2 * t) - 1
                max_dev = max(max_dev, abs(dev))
                max_ph = max(max_ph, abs(mr.phase(r)))
                push!(rows, (R2_bulk=R2, seed=seed, t_ms=t, signal=s, expected=exp(-R2 * t), rel_dev=dev, phase_deg=mr.phase(r)))
            end
        end
        final = readout(start(1, sim), sim, [TIMES[end]]; return_snapshot=true)[1]
        same_pos = mr.position.(final) == ref_pos
        per_spin_dev = maximum(abs(s.orientations[1].transverse / exp(-R2 * TIMES[end]) - 1) for s in final.spins)
        pass = max_dev <= TOL_SIGNAL && per_spin_dev <= TOL_PER_SPIN && max_ph <= TOL_PHASE && same_pos
        rt_ratio = mean(rts[2:end]) / mean(rt0[2:end])
        tee(io, @sprintf("  ensemble max |S/exp(−R2 t) − 1| = %.2e (tol %.1e); per spin max = %.2e (tol %.0e); max |phase| = %.2e°; positions identical to R2 = 0: %s",
            max_dev, TOL_SIGNAL, per_spin_dev, TOL_PER_SPIN, max_ph, same_pos))
        tee(io, @sprintf("  runtime %.2f s per seed, ratio to R2 = 0: %.3f  → %s", mean(rts[2:end]), rt_ratio, pass ? "PASS" : "FAIL"))
        push!(results, Dict("R2_bulk" => R2, "tau_max_ms" => sim.timestep.max_timestep, "max_abs_rel_dev" => max_dev, "max_abs_rel_dev_per_spin" => per_spin_dev,
            "max_abs_phase_deg" => max_ph, "positions_identical_to_R2_0" => same_pos, "pass" => pass,
            "runtime_per_seed_s" => rts, "runtime_ratio_to_R2_0" => rt_ratio))
    end
end
write_csv(joinpath(OUT, "signal.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict(
    "provenance" => provenance(), "results" => results,
    "geometry" => "Walls(repeats=$W): planes x = k·$W, no surface parameters",
    "physics" => "D=$D, R1=0, global R2 = R2_bulk, θ_relax=0, permeability=0, off-resonance=0",
    "sequence" => "empty_sequence() — no RF; spins start transverse=1, longitudinal=0, phase=0",
    "readout_times_ms" => TIMES, "nspins_per_seed" => NSPINS, "nseeds" => NSEEDS,
    "pass_criterion" => "ensemble |S/exp(−R2 t) − 1| ≤ $TOL_SIGNAL ((n_steps+N)·eps); per spin ≤ $TOL_PER_SPIN; |phase| ≤ $(TOL_PHASE)°; positions == R2=0 run (seed 1, t=50 ms)",
    "criterion_note" => "first run used ensemble tol 1e-12, below MCMR summation roundoff (~2.5e-12); replaced by the roundoff bound plus a per-spin check",
))
println("saved to $OUT")
