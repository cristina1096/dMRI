# Phase 2 test harness (P2.1.3): the reference wall configuration, run over seeds.
#
# Reference configuration (tracker, "Reference configuration"): Walls(repeats=2) → w = 2 µm, D = 3 µm²/ms,
# no RF, spins start transverse in a ±500 µm box, readouts 0:5:50 ms, τ = 1e-2 ms, seeds 1–10,
# 1e5 spins per seed, layer on both faces, R2_bulk = 0 unless stated, ΔR2(0) = 0.1 /ms.
# Uncertainty of a 10-seed mean = SEM = σ/√n_seeds.

include(joinpath(@__DIR__, "..", "baseline", "common.jl"))

const W_REF = 2.0
const D_REF = 3.0
const TAU_REF = 1e-2
const TIMES_REF = collect(0.0:5.0:50.0)
const SEEDS_REF = 1:10
const NSPINS_REF = 100_000
const RATE_REF = 0.1
"Layer thicknesses (µm) for the P2.3.4(b) sweep and the V4 fit: log-spaced at small h, up to h = w/2."
const H_GRID = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0]

"research/results/phase2/<name>/, created on demand."
function phase2_outdir(name::AbstractString)
    d = joinpath(REPO_ROOT, "research", "results", "phase2", name)
    mkpath(d)
    return d
end

"Layer thicknesses (µm) for the linear-profile sweep: as H_GRID, extended to h = w (the linear exact case)."
const H_GRID_LINEAR = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]

"g(0) of a profile for thickness `h` (µm); the exponential is truncated at `w` (the default `layer_cutoff` = spacing)."
profile_g0(shape, h; w=W_REF) = mr.Geometries.Internal.Layers.profile_g(Symbol(shape), 0.0, w / h)

"ρ (µm/ms) giving surface excess rate ΔR2(0) = `rate` (1/ms) for thickness `h` (µm): ρ = rate·h / g(0)."
rho_for(rate, h; shape="step", w=W_REF) = iszero(h) ? 0.0 : Float64(rate * h / profile_g0(shape, h; w=w))

"""
    layer_walls(; h=0.0, rate=0.0, side=:both, h_neg=h, rate_neg=rate, shape="step", w=W_REF)

Walls every `w` µm with a layer of profile `shape` and surface excess rate `rate`.
`side` is `:both`, `:positive`, `:negative`, or `:asymmetric` (positive side `h`, `rate`; negative side `h_neg`, `rate_neg`).
"""
function layer_walls(; h=0.0, rate=0.0, side=:both, h_neg=h, rate_neg=rate, shape="step", w=W_REF)
    ρ(r, hh) = rho_for(r, hh; shape=shape, w=w)
    side == :both && return Walls(repeats=w, layer_rho=ρ(rate, h), layer_h=h, layer_shape=shape)
    side == :positive && return Walls(repeats=w, layer_rho_positive=ρ(rate, h), layer_h_positive=h, layer_shape=shape)
    side == :negative && return Walls(repeats=w, layer_rho_negative=ρ(rate, h), layer_h_negative=h, layer_shape=shape)
    side == :asymmetric && return Walls(repeats=w,
        layer_rho_positive=ρ(rate, h), layer_h_positive=h,
        layer_rho_negative=ρ(rate_neg, h_neg), layer_h_negative=h_neg, layer_shape=shape)
    error("Unknown side $side; use :both, :positive, :negative or :asymmetric.")
end

"Conventional surface relaxation with the same relaxivity ρ (µm/ms): θ = ρ·√(π/D) (MCMR loses θ·√τ per collision)."
theta_for(ρ; D=D_REF) = ρ * sqrt(pi / D)

"Walls every `w` µm with conventional surface relaxation θ on both faces."
baseline_walls(; theta, w=W_REF) = Walls(repeats=w, surface_relaxation=theta)

"Walls every `w` µm with a layer of integrated relaxivity `rho` (µm/ms), thickness `h` and profile `shape`, both faces."
rho_walls(; rho, h, shape, w=W_REF) = Walls(repeats=w, layer_rho=rho, layer_h=h, layer_shape=shape)

"Timestep for a configuration: the adopted layer constraint τ ≤ 0.005/ΔR2(0), capped at TAU_REF."
tau_for(rho, h, shape; w=W_REF) = min(TAU_REF, mr.TimeSteps.DEFAULT_LAYER_SCALING / (rho / h * profile_g0(shape, h; w=w)))

"Mean |Mxy| over the spins of a snapshot (equals per_spin of the summed signal while all phases are 0)."
mean_transverse(snap) = mean(s.orientations[1].transverse for s in snap.spins)

"Fraction of spins within `h` of either wall face (walls at multiples of W_REF along x; h ≤ W_REF/2)."
function occupancy(snap, h)
    u = [mod(p[1], W_REF) for p in mr.position.(snap)]
    return count(x -> x <= h || W_REF - x <= h, u) / length(u)
end

"""
    run_walls(geometry; R2_bulk=0.0, timestep=TAU_REF, times=TIMES_REF, seeds=SEEDS_REF, nspins=NSPINS_REF, snapshots=false)

Runs the reference configuration with `geometry` for every seed (`Random.seed!(seed)` before each snapshot,
so seed k gives the same initial positions in every run). `timestep=nothing` uses the simulator default.
Returns `(S, mean, sigma, sem, times, tau, runtime, snapshots, sim)`; `S[i, j]` is the per-spin signal of
seed i at times[j]; with `snapshots=true`, `snapshots[i]` is seed i's final Snapshot.
"""
function run_walls(geometry; R2_bulk=0.0, timestep=TAU_REF, times=TIMES_REF, seeds=SEEDS_REF, nspins=NSPINS_REF, snapshots=false, D=D_REF)
    kw = isnothing(timestep) ? (;) : (; timestep=timestep)
    sim = Simulation([mr.SequenceParts.empty_sequence()]; geometry=geometry, diffusivity=D, R2=R2_bulk, verbose=false, kw...)
    times = collect(Float64, times)
    S = zeros(length(seeds), length(times))
    snaps = Snapshot[]
    runtime = 0.0
    for (i, seed) in enumerate(seeds)
        Random.seed!(seed)
        snap = Snapshot(nspins, sim, 500; transverse=1.0, longitudinal=0.0)
        t0 = time()
        res = readout(snap, sim, times; return_snapshot=snapshots)
        runtime += time() - t0
        S[i, :] = snapshots ? [mean_transverse(r) for r in res] : [per_spin(r) for r in res]
        snapshots && push!(snaps, res[end])
    end
    m = vec(mean(S, dims=1))
    s = length(seeds) > 1 ? vec(std(S, dims=1)) : zeros(length(times))
    return (S=S, mean=m, sigma=s, sem=s ./ sqrt(length(seeds)), times=times,
        tau=sim.timestep.max_timestep, runtime=runtime, snapshots=snaps, sim=sim)
end
