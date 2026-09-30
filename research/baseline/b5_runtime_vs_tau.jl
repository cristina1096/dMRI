# B5 (slim) — runtime of the unmodified simulator vs timestep τ on walls (cost reference for P2.6)
#
# The R₂(d) timestep constraint τ ≤ c_h·h²/D (proposal Eq. 11) can require very small τ for thin
# layers. This run measures the per-step cost of the OLD code over τ = 1e-1 … 1e-5 ms, with and without
# surface relaxation, so that (a) the feasible τ range of the Phase 2 timestep study can be planned, and
# (b) the cost of the new code can later be quoted relative to the old one at the same τ.
# Runtime only: signal accuracy vs τ is covered by B4.
#
# Variables:
#   geometry : Walls(repeats=w, surface_relaxation=θ), w = 2 µm, θ ∈ {0, 0.01}
#   D = 3 µm²/ms; R1 = R2 = 0; no RF (empty sequence); spins start transverse = 1
#   τ        : 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 ms (fixed, user-set)
#   workload : NSPINS (default 10_000) spins × NSTEPS (default 20_000) steps per spin, i.e. simulated
#              time T = NSTEPS·τ (2000 ms at τ = 1e-1 … 0.2 ms at τ = 1e-5); one readout at T
#   repeats  : NREP (default 3) timed runs per configuration after one untimed warm-up; median reported
#   threads  : julia -t 8
#
# Reported: µs per spin-step, µs per spin per ms of sequence, and the projected wall time of one seed of
# the Phase 2 reference run (1e5 spins, 50 ms) at each τ.
#
# Run: julia --project=research/baseline -t 8 research/baseline/b5_runtime_vs_tau.jl

include(joinpath(@__DIR__, "common.jl"))

const OUT = outdir("b5_runtime_vs_tau")
const NSPINS = parse(Int, get(ENV, "NSPINS", "10000"))
const NSTEPS = parse(Int, get(ENV, "NSTEPS", "20000"))
const NREP = parse(Int, get(ENV, "NREP", "3"))
const D = 3.0
const W = 2.0
const TAUS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
const THETAS = [0.0, 0.01]

function timed(sim, T)
    Random.seed!(1)
    snap = Snapshot(NSPINS, sim, 500; transverse=1.0, longitudinal=0.0)
    t0 = time()
    readout(snap, sim, [T])
    return time() - t0
end

rows = NamedTuple[]
open(joinpath(OUT, "run.log"), "w") do io
    tee(io, "N_spins = $NSPINS, steps per spin = $NSTEPS, repeats = $NREP, threads = $(Threads.nthreads())")
    for θ in THETAS, τ in TAUS
        sim = Simulation([mr.SequenceParts.empty_sequence()]; geometry=Walls(repeats=W, surface_relaxation=θ),
            diffusivity=D, timestep=τ, verbose=false)
        T = NSTEPS * τ
        timed(Simulation([mr.SequenceParts.empty_sequence()]; geometry=Walls(repeats=W, surface_relaxation=θ),
            diffusivity=D, timestep=τ, verbose=false), 20τ)   # warm-up / compilation
        ts = [timed(sim, T) for _ in 1:NREP]
        med = median(ts)
        us_step = med / (NSPINS * NSTEPS) * 1e6
        us_ms = us_step / τ
        proj = us_ms * 1e5 * 50 / 1e6   # seconds per seed, 1e5 spins × 50 ms
        push!(rows, (theta=θ, tau_ms=τ, sim_time_ms=T, nspins=NSPINS, nsteps=NSTEPS, runtime_median_s=med,
            runtime_min_s=minimum(ts), runtime_max_s=maximum(ts), us_per_spin_step=us_step,
            us_per_spin_per_ms=us_ms, projected_s_per_seed_1e5_50ms=proj))
        tee(io, @sprintf("θ=%-5g τ=%-6g T=%-7g ms  median %.2f s (min %.2f, max %.2f)  %.4f µs/spin-step  %.3g µs/spin/ms  → 1e5 spins × 50 ms ≈ %s per seed",
            θ, τ, T, med, minimum(ts), maximum(ts), us_step, us_ms,
            proj < 120 ? @sprintf("%.0f s", proj) : proj < 7200 ? @sprintf("%.1f min", proj / 60) : @sprintf("%.1f h", proj / 3600)))
    end
end
write_csv(joinpath(OUT, "runtime.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows,
    "geometry" => "Walls(repeats=$W, surface_relaxation=θ)", "physics" => "D=$D, R1=R2=0, no RF",
    "workload" => "$NSPINS spins × $NSTEPS steps, median of $NREP timed runs after warm-up"))
println("saved to $OUT")
