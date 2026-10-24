# Phase 3 research harness: near-surface layer on cylinders (single and arrays), no RF.
# Settings, formulas and decisions: research/phase3_record.md (D3.1–D3.6). Baselines (MCMR θ) are references, not truth.

include(joinpath(@__DIR__, "..", "baseline", "common.jl"))
using Statistics, Random

const R_REF = 1.0                       # cylinder radius (um)
const S_REF = 2.5                       # square-array centre spacing (um): surface gap 0.5 um
const D_REF = 3.0                       # um²/ms
const TAU_REF = 1e-2                    # ms
const TIMES_REF = collect(0.0:5.0:50.0)
const SEEDS_REF = 1:10
const NSPINS_REF = 100_000

phase3_outdir(name::AbstractString) = mkpath(joinpath(REPO_ROOT, "research", "results", "phase3", name))

"Surface-to-volume ratio inside a cylinder: 2/R."
sv_inside(R) = 2 / R
"Surface-to-volume ratio outside, one cylinder per sx × sy cell: 2πR/(sx·sy − πR²)."
sv_outside(R, sx, sy) = 2π * R / (sx * sy - π * R^2)
"Area fraction of the inside within h of the surface: (R² − (R − h)²)/R²."
shell_area_inside(R, h) = (R^2 - (R - h)^2) / R^2
"Area fraction of the outside within h of the surface (shells not overlapping): π((R + h)² − R²)/(sx·sy − πR²)."
shell_area_outside(R, h, sx, sy) = π * ((R + h)^2 - R^2) / (sx * sy - π * R^2)
"Step layer rate inside (D3.1): ρ/(h(1 − h/2R)); written independently of the package."
rate_inside(ρ, h, R) = ρ / (h * (1 - h / (2R)))
"Step layer rate outside (D3.1): ρ/(h(1 + h/2R))."
rate_outside(ρ, h, R) = ρ / (h * (1 + h / (2R)))

"Distance from the nearest cylinder axis (one cylinder per cell at the origin; sx = sy = Inf for a single cylinder)."
function radial_distance(pos; sx=Inf, sy=Inf)
    x = isfinite(sx) ? mod(pos[1] + sx / 2, sx) - sx / 2 : pos[1]
    y = isfinite(sy) ? mod(pos[2] + sy / 2, sy) - sy / 2 : pos[2]
    return hypot(x, y)
end

"""
    brute_rate(x, y; R, s, h_in, h_out, rate_in, rate_out)

Summed step-layer ΔR2 at (x, y) for a square array of spacing `s` (one cylinder per cell at the origin; `s = Inf` for
a single cylinder), from the 5 × 5 nearest images. Independent of the package's `Layers` code.
"""
function brute_rate(x, y; R, s, h_in, h_out, rate_in, rate_out)
    cells = isfinite(s) ? (-2:2) : (0:0)
    cx = isfinite(s) ? round(x / s) * s : 0.0
    cy = isfinite(s) ? round(y / s) * s : 0.0
    total = 0.0
    for i in cells, j in cells
        r = hypot(x - cx - i * (isfinite(s) ? s : 0.0), y - cy - j * (isfinite(s) ? s : 0.0))
        (R - h_in <= r <= R) && (total += rate_in)
        (R <= r <= R + h_out) && (total += rate_out)
    end
    return total
end

"Floating-point bound for an ensemble mean after maximum(times)/τ multiplications per spin and a sum over nspins."
roundoff_bound(times, tau, nspins) = (maximum(times) / tau + nspins) * eps()

"Conventional model (MCMR θ): E[−ln M⊥] = θ√(D/π)·(S/V)·t, if the wall-collision rate holds on curved surfaces (P3.1.6)."
expected_neglog_conventional(θ, D, sv, t) = θ * sqrt(D / π) * sv * t
"Layer, ρ per unit surface area (D3.1): E[−ln M⊥] = ρ·(S/V)·t (mean-exponent identity, P3.4.3)."
expected_neglog_layer(ρ, sv, t) = ρ * sv * t

"""
    run_cylinders(geometry; region, R=R_REF, sx=Inf, sy=Inf, bounding_box=500, R2_bulk=0.0, timestep=TAU_REF,
                  times=TIMES_REF, seeds=SEEDS_REF, nspins=NSPINS_REF, D=D_REF)

Free decay (no RF) for every seed (`Random.seed!(seed)` before each snapshot, so seed k gives the same initial
positions in every run). `region` (:inside or :outside) selects spins by their distance from the nearest axis at t = 0;
surfaces are impermeable, so the set does not change. `timestep=nothing` uses the simulator default.
Returns `S[i, j]` (mean M⊥ of seed i's region spins at times[j]) and `neglog[i, j]` (mean −ln M⊥ over them), with
mean, sigma and SEM over seeds; `nregion[i]`, `masks[i]` and `finals[i]` (final Snapshot) per seed.
"""
function run_cylinders(geometry; region::Symbol, R=R_REF, sx=Inf, sy=Inf, bounding_box=500, R2_bulk=0.0,
                       timestep=TAU_REF, times=TIMES_REF, seeds=SEEDS_REF, nspins=NSPINS_REF, D=D_REF)
    region in (:inside, :outside) || error("region must be :inside or :outside, got $region")
    kw = isnothing(timestep) ? (;) : (; timestep=timestep)
    sim = Simulation([mr.SequenceParts.empty_sequence()]; geometry=geometry, diffusivity=D, R2=R2_bulk, verbose=false, kw...)
    times = collect(Float64, times)
    S = zeros(length(seeds), length(times))
    L = zeros(length(seeds), length(times))
    nregion = zeros(Int, length(seeds))
    masks = BitVector[]
    finals = Snapshot[]
    runtime = 0.0
    for (i, seed) in enumerate(seeds)
        Random.seed!(seed)
        snap = Snapshot(nspins, sim, bounding_box; transverse=1.0, longitudinal=0.0)
        r0 = [radial_distance(s.position; sx=sx, sy=sy) for s in snap.spins]
        mask = BitVector(region == :inside ? r0 .< R : r0 .> R)
        t0 = time()
        res = readout(snap, sim, times; return_snapshot=true)
        runtime += time() - t0
        for (j, r) in enumerate(res)
            m = [r.spins[k].orientations[1].transverse for k in findall(mask)]
            S[i, j] = mean(m)
            L[i, j] = mean(-log.(m))
        end
        nregion[i] = count(mask)
        push!(masks, mask)
        push!(finals, res[end])
    end
    n = length(seeds)
    sd(X) = n > 1 ? vec(std(X, dims=1)) : zeros(size(X, 2))
    return (S=S, mean=vec(mean(S, dims=1)), sigma=sd(S), sem=sd(S) ./ sqrt(n), neglog=L,
        neglog_mean=vec(mean(L, dims=1)), neglog_sem=sd(L) ./ sqrt(n), times=times,
        tau=sim.timestep.max_timestep, runtime=runtime, nregion=nregion, masks=masks, finals=finals, sim=sim)
end
