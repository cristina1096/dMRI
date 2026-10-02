# Phase 2 Harness, Expectations, Noise Floor and Layer Trace Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the research tools every Phase 2 verification run uses: a test harness (P2.1.3), the exact expectations and B4 reference interpolation (P2.1.5), the noise floor and spin count (P2.1.4), and the end-to-end single-spin trace with the layer on (P2.2.9).

**Architecture:** Research code lives in `research/phase2/`, next to the Phase 1 scripts in `research/baseline/`, and runs in the same `research/baseline` Julia environment. `harness.jl` wraps "build the reference wall configuration, run it over seeds, return mean/σ/SEM" so every later script uses one code path. `expected.jl` holds the exact expectations and the B4 reference loader and interpolation. Library functions are unit-tested in `research/phase2/test/`. Run scripts write to `research/results/phase2/<run>/`, and results are recorded in `research/phase2_record.md`.

**Tech Stack:** Julia 1.12, MCMRSimulator (this repo, with the walls layer from plan `2026-10-01-near-surface-layer-walls.md`), CairoMakie and JSON (already in `research/baseline/Project.toml`), `Test` stdlib.

**Spec:** `research/project-progress-tracker.md` rows P2.1.3, P2.1.4, P2.1.5, P2.2.9 and the "Reference configuration" and "Exact expectations" tables; proposal (Y. Shi, 2026-09-14) §3.2; the user's framing in memory: the baseline θ_relax model is a reference to fit against, not the truth, and no fast-diffusion arguments are used.

## Global Constraints

- Work only on branch `cc/near-surface-r2`.
- Reference configuration: `Walls(repeats=2)` (w = 2 µm), D = 3 µm²/ms, no RF, spins start transverse (`Snapshot(...; transverse=1, longitudinal=0)`) in a ±500 µm box, readouts t = 0, 5, …, 50 ms, timestep τ = 1e-2 ms, seeds 1–10, 10⁵ spins per seed, layer on both faces, R₂_bulk = 0 unless stated.
- Surface excess rate ΔR₂(0) = 0.1 ms⁻¹ unless stated; ρ = ΔR₂(0)·h for the step profile.
- Uncertainty of a 10-seed mean is the SEM = σ/√10. Never use the single-seed σ as the uncertainty of a mean (lesson from Rep1/Rep3).
- No fast-diffusion or Brownstein–Tarr formulas anywhere in Phase 2 code or records.
- Run research tests with: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `research/phase2/harness.jl` | Create | constants of the reference configuration, `layer_walls`, `rho_for`, `run_walls`, `occupancy`, `phase2_outdir` |
| `research/phase2/expected.jl` | Create | exact expectations, `load_b4`, `theta_for_signal`, `reference_at` |
| `research/phase2/test/runtests.jl` | Create | runs every `test_*.jl` in the folder |
| `research/phase2/test/test_harness.jl` | Create | harness unit tests |
| `research/phase2/test/test_expected.jl` | Create | expectations / B4 tests |
| `research/phase2/p2_2_9_layer_trace_lib.jl` | Create | `sampled_layer_fraction`: brute-force layer time, independent of `Layers` |
| `research/phase2/test/test_trace.jl` | Create | tests for `sampled_layer_fraction` |
| `research/phase2/p2_2_9_layer_trace.jl` | Create | P2.2.9 trace script |
| `research/phase2/noise.jl` | Create | `required_spins` |
| `research/phase2/test/test_noise.jl` | Create | tests for `required_spins` |
| `research/phase2/p2_1_4_noise_floor.jl` | Create | P2.1.4 noise-floor script |
| `research/phase2_record.md` | Create | Phase 2 results record |
| `research/project-progress-tracker.md` | Modify | tick P2.1.3, P2.1.4, P2.1.5, P2.2.9 |

## Review Focus

1. **Signal from snapshots vs from sums.** `run_walls(...; snapshots=true)` averages per-spin |Mxy|, while `snapshots=false` uses MCMR's summed vector. They agree only while every phase is 0, which holds here because there's no off-resonance. Test: Task 1 checks that both paths agree to 1e-12.
2. **Interpolation outside the B4 range** must return `NaN`, never extrapolate. Test: Task 2.
3. **A θ exactly on the B4 grid** must return the stored curve unchanged, not exp(log(x)). Test: Task 2.
4. **The trace must include reflections inside the layer**, or P2.2.9 proves nothing about P2.2.4. Test: Task 3 asserts at least one step with a reflection and non-zero layer time.
5. **The noise-floor rule** must use the SEM of the 10-seed mean. Test: Task 4's `required_spins` unit test.

---

### Task 1: Harness (P2.1.3)

**Files:**
- Create: `research/phase2/harness.jl`, `research/phase2/test/runtests.jl`, `research/phase2/test/test_harness.jl`, `research/phase2_record.md`

**Interfaces:**
- Consumes: `research/baseline/common.jl` (`REPO_ROOT`, `per_spin`, `provenance`, `write_json`, `write_csv`); MCMR `Walls` layer keywords.
- Produces:
  - constants `W_REF = 2.0`, `D_REF = 3.0`, `TAU_REF = 1e-2`, `TIMES_REF = 0:5:50` (Vector{Float64}), `SEEDS_REF = 1:10`, `NSPINS_REF = 100_000`, `RATE_REF = 0.1`, `H_GRID`
  - `phase2_outdir(name)::String`
  - `rho_for(rate, h)::Float64`
  - `layer_walls(; h=0.0, rate=0.0, side=:both, h_neg=h, rate_neg=rate)`, which returns a `Walls`
  - `run_walls(geometry; R2_bulk=0.0, timestep=TAU_REF, times=TIMES_REF, seeds=SEEDS_REF, nspins=NSPINS_REF, snapshots=false)`, which returns a NamedTuple `(S, mean, sigma, sem, times, tau, runtime, snapshots, sim)`. `S` is seeds × times; `snapshots` is a Vector of final Snapshots (empty unless `snapshots=true`).
  - `mean_transverse(snap)::Float64`
  - `occupancy(snap, h)::Float64`

- [ ] **Step 1: Write the failing tests**

`research/phase2/test/runtests.jl`:

```julia
using Test
include(joinpath(@__DIR__, "..", "harness.jl"))
isfile(joinpath(@__DIR__, "..", "expected.jl")) && include(joinpath(@__DIR__, "..", "expected.jl"))
isfile(joinpath(@__DIR__, "..", "fit.jl")) && include(joinpath(@__DIR__, "..", "fit.jl"))
isfile(joinpath(@__DIR__, "..", "toy.jl")) && include(joinpath(@__DIR__, "..", "toy.jl"))

@testset "Phase 2 research code" begin
    for f in sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(@__DIR__)))
        include(joinpath(@__DIR__, f))
    end
end
```

`research/phase2/test/test_harness.jl`:

```julia
@testset "harness" begin
    @testset "rho_for" begin
        @test rho_for(0.1, 0.5) == 0.05
        @test rho_for(0.1, 0.0) == 0.0
    end

    @testset "layer_walls sides" begin
        g(w) = Simulation([]; geometry=w, verbose=false).geometry[1].layer
        @test g(layer_walls()) === nothing
        @test g(layer_walls(h=0.5, rate=0.0)) === nothing
        l = g(layer_walls(h=0.5, rate=0.1))[1]
        @test (l.positive.rho, l.positive.h, l.negative.rho, l.negative.h) == (0.05, 0.5, 0.05, 0.5)
        l = g(layer_walls(h=0.5, rate=0.1, side=:positive))[1]
        @test (l.positive.rho, l.negative.rho) == (0.05, 0.0)
        l = g(layer_walls(h=0.5, rate=0.1, side=:negative))[1]
        @test (l.positive.rho, l.negative.rho) == (0.0, 0.05)
        l = g(layer_walls(h=0.2, rate=0.1, side=:asymmetric, h_neg=0.4, rate_neg=0.05))[1]
        @test (l.positive.rho, l.positive.h, l.negative.rho, l.negative.h) == (rho_for(0.1, 0.2), 0.2, rho_for(0.05, 0.4), 0.4)
        @test_throws ErrorException layer_walls(side=:top)
    end

    @testset "run_walls reproduces the exact cases" begin
        t = [0.0, 10.0, 50.0]
        r = run_walls(layer_walls(); times=t, seeds=1:2, nspins=1000)            # B1: nothing relaxes
        @test size(r.S) == (2, 3)
        @test all(r.S .== 1.0)
        r = run_walls(layer_walls(); R2_bulk=1/80, times=t, seeds=1:2, nspins=1000)  # B2: bulk only
        @test all(isapprox.(r.S, exp.(-t' ./ 80); rtol=1e-11))
        r = run_walls(layer_walls(h=1.0, rate=0.1); times=t, seeds=1:2, nspins=1000) # h = w/2
        @test all(isapprox.(r.S, exp.(-0.1 .* t'); rtol=1e-10))
        @test r.tau == TAU_REF
        @test r.sem ≈ r.sigma ./ sqrt(2)
    end

    @testset "snapshot path agrees with the summed path" begin
        t = [0.0, 20.0]
        a = run_walls(layer_walls(h=0.3, rate=0.1); times=t, seeds=1:2, nspins=2000)
        b = run_walls(layer_walls(h=0.3, rate=0.1); times=t, seeds=1:2, nspins=2000, snapshots=true)
        @test all(isapprox.(a.S, b.S; rtol=1e-12))
        @test length(b.snapshots) == 2
        @test length(b.snapshots[1].spins) == 2000
    end

    @testset "occupancy" begin
        snap = Snapshot([[0.1, 0., 0.], [1.0, 0., 0.], [1.95, 0., 0.], [2.3, 0., 0.]])
        @test occupancy(snap, 0.2) == 0.5          # 0.1 and 1.95 are within 0.2 of a face
    end
end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `could not open file .../research/phase2/harness.jl`.

- [ ] **Step 3: Create `research/phase2/harness.jl`**

```julia
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

"ρ (µm/ms) of a step-profile layer with surface excess rate `rate` (1/ms) and thickness `h` (µm): ρ = rate·h."
rho_for(rate, h) = Float64(rate * h)

"""
    layer_walls(; h=0.0, rate=0.0, side=:both, h_neg=h, rate_neg=rate)

Reference walls with a step-profile layer. `side` is `:both`, `:positive`, `:negative`, or `:asymmetric`
(positive side `h`, `rate`; negative side `h_neg`, `rate_neg`). `rate = 0` or `h = 0` gives no layer.
"""
function layer_walls(; h=0.0, rate=0.0, side=:both, h_neg=h, rate_neg=rate)
    side == :both && return Walls(repeats=W_REF, layer_rho=rho_for(rate, h), layer_h=h)
    side == :positive && return Walls(repeats=W_REF, layer_rho_positive=rho_for(rate, h), layer_h_positive=h)
    side == :negative && return Walls(repeats=W_REF, layer_rho_negative=rho_for(rate, h), layer_h_negative=h)
    side == :asymmetric && return Walls(repeats=W_REF,
        layer_rho_positive=rho_for(rate, h), layer_h_positive=h,
        layer_rho_negative=rho_for(rate_neg, h_neg), layer_h_negative=h_neg)
    error("Unknown side $side; use :both, :positive, :negative or :asymmetric.")
end

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
function run_walls(geometry; R2_bulk=0.0, timestep=TAU_REF, times=TIMES_REF, seeds=SEEDS_REF, nspins=NSPINS_REF, snapshots=false)
    kw = isnothing(timestep) ? (;) : (; timestep=timestep)
    sim = Simulation([mr.SequenceParts.empty_sequence()]; geometry=geometry, diffusivity=D_REF, R2=R2_bulk, verbose=false, kw...)
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
```

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS (testset "harness").

- [ ] **Step 5: Create `research/phase2_record.md`**

```markdown
# Phase 2 · Record (R₂(d) on walls, step profile, no RF)

Results of the modified simulator (branch `cc/near-surface-r2`, layer implemented by plan `docs/superpowers/plans/2026-10-01-near-surface-layer-walls.md`).
Baselines are in [baseline_record.md](baseline_record.md). The baseline θ_relax model is a **reference**, not the true answer.

- Code: [research/phase2/](phase2/) (harness, expectations, scripts); tests: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
- Raw outputs: [research/results/phase2/](results/phase2/), each with `summary.json` (provenance: git commit, versions, hardware)

**Reference configuration:** `Walls(repeats=2)` (w = 2 µm), D = 3 µm²/ms, no RF, spins start transverse in a ±500 µm box, readouts 0, 5, …, 50 ms, τ = 1e-2 ms, seeds 1–10, 10⁵ spins per seed, step-profile layer on both faces, ΔR₂(0) = 0.1 ms⁻¹ (ρ = ΔR₂(0)·h), R₂_bulk = 0 unless stated. Uncertainty of a mean = SEM over seeds.

## Summary

| ID | Test | Criterion | Outcome |
|---|---|---|---|
```

- [ ] **Step 6: Commit**

```bash
git add research/phase2/harness.jl research/phase2/test/runtests.jl research/phase2/test/test_harness.jl research/phase2_record.md
git commit -m "Add Phase 2 research harness

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Exact expectations and B4 reference interpolation (P2.1.5)

**Files:**
- Create: `research/phase2/expected.jl`, `research/phase2/test/test_expected.jl`

**Interfaces:**
- Consumes: `REPO_ROOT`, `run_walls`, `TIMES_REF` (Task 1); `research/results/baseline/b4_theta_reference/reference_curves.csv` (columns `theta, config, tau_ms, R2_bulk, t_ms, mean_signal, sigma, sem, nseeds, nspins`).
- Produces:
  - `expected_null(times)`, `expected_bulk(times, R2_bulk)`, `expected_half_gap(times, rate; R2_bulk=0.0)`, `expected_additive(S_layer, times, R2_bulk)`, `roundoff_bound(times, tau, nspins)`
  - `load_b4(; config="tau_1e-2")`, which returns `Dict{Float64, NamedTuple{(:times, :mean, :sigma, :sem)}}` keyed by θ
  - `theta_for_signal(ref, S_target; t=50.0)::Float64` (`NaN` outside the grid)
  - `reference_at(ref, θ)::Vector{Float64}` (S at every readout; interpolates log S linearly in log θ)

- [ ] **Step 1: Write the failing tests**

`research/phase2/test/test_expected.jl`:

```julia
@testset "expected" begin
    t = [0.0, 10.0, 50.0]

    @testset "exact cases" begin
        @test expected_null(t) == [1.0, 1.0, 1.0]
        @test expected_bulk(t, 0.1) ≈ exp.(-0.1 .* t)
        @test expected_half_gap(t, 0.1) ≈ exp.(-0.1 .* t)
        @test expected_half_gap(t, 0.1; R2_bulk=1/80) ≈ exp.(-(0.1 + 1/80) .* t)
        @test expected_additive([1.0, 0.5, 0.2], t, 0.02) ≈ [1.0, 0.5, 0.2] .* exp.(-0.02 .* t)
        @test roundoff_bound([0.0, 50.0], 0.04, 100_000) ≈ (1250 + 100_000) * eps()
    end

    ref = load_b4()

    @testset "load_b4" begin
        @test sort(collect(keys(ref))) == [0.0, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05]
        @test ref[0.01].times == TIMES_REF
        @test ref[0.0].mean == ones(11)
        @test ref[0.01].mean[end] ≈ 0.61379 atol=1e-5          # baseline_record.md B4 table
        @test ref[0.05].mean[end] ≈ 0.08802 atol=1e-5
        @test all(ref[0.01].sem .<= ref[0.01].sigma)
    end

    @testset "theta_for_signal" begin
        @test theta_for_signal(ref, ref[0.01].mean[end]) ≈ 0.01 rtol=1e-12
        θ = theta_for_signal(ref, (ref[0.01].mean[end] + ref[0.015].mean[end]) / 2)
        @test 0.01 < θ < 0.015
        @test isnan(theta_for_signal(ref, 0.99))      # above S(50) of the smallest θ > 0 (0.952)
        @test isnan(theta_for_signal(ref, 0.01))      # below S(50) of the largest θ (0.088)
        @test theta_for_signal(ref, ref[0.005].mean[6]; t=25.0) ≈ 0.005 rtol=1e-12
    end

    @testset "reference_at" begin
        @test reference_at(ref, 0.01) === ref[0.01].mean                 # grid point: stored curve unchanged
        mid = reference_at(ref, sqrt(0.01 * 0.015))                       # geometric midpoint
        @test all(ref[0.015].mean .<= mid .<= ref[0.01].mean)
        @test_throws ErrorException reference_at(ref, 0.2)                # outside the grid
    end

    @testset "harness reproduces one B4 curve exactly (P2.1.3 criterion)" begin
        r = run_walls(Walls(repeats=2., surface_relaxation=0.01); timestep=1e-2)  # seeds 1–10, 1e5 spins
        @test all(isapprox.(r.mean, ref[0.01].mean; rtol=1e-8))                     # CSV stores 10 significant digits
    end
end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `UndefVarError: expected_null not defined`.

- [ ] **Step 3: Create `research/phase2/expected.jl`**

```julia
# Phase 2 exact expectations and the B4 reference curves (P2.1.5).
#
# Exact cases (tracker, "Exact expectations"): they hold in any diffusion regime.
# Everything else is measured and compared with the B4 baseline θ_relax curves, which are a reference, not the truth.

"h = 0 (no layer), R2_bulk = 0: the signal stays 1."
expected_null(times) = ones(length(times))

"Bulk R2 only: exp(-R2_bulk·t)."
expected_bulk(times, R2_bulk) = exp.(-R2_bulk .* times)

"Step profile, h = w/2, both faces: every point of the gap is in exactly one layer, so S = exp(-(ΔR2(0) + R2_bulk)·t)."
expected_half_gap(times, rate; R2_bulk=0.0) = exp.(-(rate + R2_bulk) .* times)

"Bulk on top of a layer (same seeds): S_bulk+layer = S_layer·exp(-R2_bulk·t)."
expected_additive(S_layer, times, R2_bulk) = S_layer .* exp.(-R2_bulk .* times)

"Floating-point bound for an ensemble mean after maximum(times)/τ multiplications per spin and a sum over nspins."
roundoff_bound(times, tau, nspins) = (maximum(times) / tau + nspins) * eps()

const B4_CSV = joinpath(REPO_ROOT, "research", "results", "baseline", "b4_theta_reference", "reference_curves.csv")

"""
    load_b4(; config="tau_1e-2")

B4 reference curves (bulk off) for one timestep configuration ("default", "tau_1e-2", "tau_1e-3").
Returns Dict θ => (times, mean, sigma, sem), readouts sorted by time.
"""
function load_b4(; config="tau_1e-2")
    lines = readlines(B4_CSV)
    header = split(lines[1], ",")
    col(name) = findfirst(==(name), header)
    rows = [split(l, ",") for l in lines[2:end]]
    num(r, name) = parse(Float64, r[col(name)])
    rows = filter(r -> r[col("config")] == config && num(r, "R2_bulk") == 0.0, rows)
    out = Dict{Float64, NamedTuple{(:times, :mean, :sigma, :sem), NTuple{4, Vector{Float64}}}}()
    for θ in sort(unique(num.(rows, "theta")))
        rs = sort(filter(r -> num(r, "theta") == θ, rows), by=r -> num(r, "t_ms"))
        f(name) = [num(r, name) for r in rs]
        out[θ] = (times=f("t_ms"), mean=f("mean_signal"), sigma=f("sigma"), sem=f("sem"))
    end
    return out
end

"""
    theta_for_signal(ref, S_target; t=50.0)

θ_relax of the baseline whose signal at time `t` equals `S_target`: log S interpolated linearly in log θ
between grid values (θ > 0). Returns NaN if `S_target` is outside the range spanned by the grid.
"""
function theta_for_signal(ref, S_target; t=50.0)
    θs = sort([θ for θ in keys(ref) if θ > 0])
    k = findfirst(==(t), ref[θs[1]].times)
    Ss = [ref[θ].mean[k] for θ in θs]                       # decreasing in θ
    (S_target > Ss[1] || S_target < Ss[end]) && return NaN
    j = findfirst(i -> Ss[i] >= S_target >= Ss[i + 1], 1:length(θs) - 1)
    S_target == Ss[j] && return θs[j]
    f = (log(S_target) - log(Ss[j])) / (log(Ss[j + 1]) - log(Ss[j]))
    return exp(log(θs[j]) + f * (log(θs[j + 1]) - log(θs[j])))
end

"""
    reference_at(ref, θ)

Baseline curve S_θ(t) at every readout for any θ inside the grid: the stored curve on a grid point,
otherwise log S interpolated linearly in log θ between the two neighbouring grid values.
"""
function reference_at(ref, θ)
    haskey(ref, θ) && return ref[θ].mean
    θs = sort([x for x in keys(ref) if x > 0])
    (θ < θs[1] || θ > θs[end]) && error("θ = $θ is outside the B4 grid [$(θs[1]), $(θs[end])].")
    j = searchsortedlast(θs, θ)
    f = (log(θ) - log(θs[j])) / (log(θs[j + 1]) - log(θs[j]))
    return exp.((1 - f) .* log.(ref[θs[j]].mean) .+ f .* log.(ref[θs[j + 1]].mean))
end
```

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS (testsets "harness" and "expected"). The last testset takes about 1 minute.

- [ ] **Step 5: Commit**

```bash
git add research/phase2/expected.jl research/phase2/test/test_expected.jl
git commit -m "Add Phase 2 exact expectations and B4 reference interpolation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: End-to-end single-spin trace with the layer on (P2.2.9)

**Files:**
- Create: `research/phase2/p2_2_9_layer_trace.jl`, `research/phase2/test/test_trace.jl`
- Modify: `research/phase2_record.md` (append)

**Interfaces:**
- Consumes: `layer_walls`, `W_REF`, `D_REF`, `phase2_outdir` (Task 1); `mr.Evolve.draw_step!(spin, sim, part, B0s, proposed)`, which returns the list of positions (start, every collision point, end) when given a proposed position.
- Produces: `sampled_layer_fraction(a, b, h; n=10_000)` (in the script; also unit-tested), the CSV `trace_spin1.csv`, and `summary.json` with `max_rel_error`, `n_steps_with_reflection_in_layer`, `pass`.

**Method.** For each step: draw the proposed displacement in the script, call `draw_step!` with it, and get the piecewise path back. The time on each straight piece is proportional to its 3D length (speed is constant along a reflected path). The layer time on each piece is computed **independently of `Layers`** by sampling 10⁴ points along it. The prediction for the step is

  M⊥(after) = M⊥(before) · exp(−Σ ΔR₂(0)·dt_piece·fraction_piece) · exp(−R₂_bulk·τ)

compared with the simulator's M⊥(after).

- [ ] **Step 1: Write the failing test for the independent integrator**

`research/phase2/test/test_trace.jl`:

```julia
include(joinpath(@__DIR__, "..", "p2_2_9_layer_trace_lib.jl"))

@testset "trace: independent layer fraction" begin
    # walls every 2 um, layers [0, h] and [2 - h, 2] in each gap
    @test sampled_layer_fraction(0.1, 0.3, 0.5) ≈ 1.0 atol=1e-4
    @test sampled_layer_fraction(0.6, 1.4, 0.5) == 0.0
    @test sampled_layer_fraction(0.4, 0.6, 0.5) ≈ 0.5 atol=1e-4
    @test sampled_layer_fraction(0.6, -0.6, 0.5) ≈ 1.0 / 1.2 atol=1e-4   # crosses wall 0, both layers
    @test sampled_layer_fraction(0.9, 1.1, 1.5) ≈ 2.0 atol=1e-4          # overlap counts twice
end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `could not open file .../research/phase2/p2_2_9_layer_trace_lib.jl`.

- [ ] **Step 3: Create the library file `research/phase2/p2_2_9_layer_trace_lib.jl`**

```julia
"""
    sampled_layer_fraction(a, b, h; n=10_000)

Time-averaged number of layers a straight piece from x = a to x = b sits in, by brute-force sampling of
`n` midpoints (walls at multiples of W_REF, layer of thickness `h` on both faces). Deliberately independent
of MCMR's `Layers` code, so it can check it. Sampling error ≤ (number of boundary crossings)/n.
"""
function sampled_layer_fraction(a::Float64, b::Float64, h::Float64; n=10_000)
    total = 0
    for k in 1:n
        u = mod(a + (b - a) * (k - 0.5) / n, W_REF)
        total += (u <= h) + (W_REF - u <= h)
    end
    return total / n
end
```

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS (testset "trace: independent layer fraction" included).

- [ ] **Step 5: Create the script `research/phase2/p2_2_9_layer_trace.jl`**

```julia
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

const OUT = phase2_outdir("p2_2_9_layer_trace")
const TAU = 0.01
const NSTEPS = 500
const NTRACE = 100
const H = 0.5
const RATE = 0.1
const R2B = 1 / 80
const TOL = 1e-6

seq = mr.SequenceParts.SequenceWaveform((([], []), ([], []), ([], [])), [], [], [NSTEPS * TAU], NSTEPS * TAU)
sim = Simulation(seq; geometry=layer_walls(h=H, rate=RATE), diffusivity=D_REF, R2=R2B, verbose=false)
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
            expo += RATE * dt * sampled_layer_fraction(path[k][1], path[k + 1][1], H)
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
        "R2_bulk" => R2B, "D" => D_REF, "w_um" => W_REF, "seed" => 20261003, "samples_per_piece" => 10_000),
    "max_rel_error" => max_err, "tolerance" => TOL, "n_steps_with_reflection_in_layer" => n_reflect_in_layer, "pass" => pass))
println("saved to $OUT")
```

- [ ] **Step 6: Run the script and check it passes**

Run: `julia --project=research/baseline -t 1 research/phase2/p2_2_9_layer_trace.jl`
Expected: a line ending `→ PASS`, with max error ≤ 1e-6 and a positive count of steps with a reflection inside the layer. If it prints `FAIL`, stop and use superpowers:systematic-debugging. A failure means the layer code and the independent integration disagree, which is exactly what P2.2.9 is meant to catch.

- [ ] **Step 7: Append the record**

Append to `research/phase2_record.md`. Add a summary-table row, and add the section below the summary, using the numbers from `research/results/phase2/p2_2_9_layer_trace/summary.json`:

```markdown
| P2.2.9 | Single-spin trace with the layer on | every step: \|M⊥_sim/M⊥_pred − 1\| ≤ 1e-6; ≥ 1 reflection inside the layer | <PASS/FAIL>: max error <max_rel_error>; <n_steps_with_reflection_in_layer> steps with a reflection inside the layer |
```

```markdown
## P2.2.9: single-spin trace with the layer on

Script [p2_2_9_layer_trace.jl](phase2/p2_2_9_layer_trace.jl) → [results/phase2/p2_2_9_layer_trace/](results/phase2/p2_2_9_layer_trace/)

**What / why.** The P1.2.6 trace repeated with the layer on. Each step's change in M⊥ is predicted from the simulator's own reflected path, with the layer time on every straight piece computed by brute-force sampling, independently of the `Layers` code. Agreement shows that distance, crossings, reflection handling, segmentation and the M⊥ update (P2.2.1–P2.2.7) work together.

| Field | Value |
|---|---|
| Configuration | walls every 2 µm; step layer on both faces, h = 0.5 µm, ΔR₂(0) = 0.1 ms⁻¹; R₂_bulk = 1/80 ms⁻¹; D = 3; τ = 0.01 ms |
| Spins / steps | 100 spins × 500 steps; `Random.seed!(20261003)` |
| Max relative error | <max_rel_error> (tolerance 1e-6; sampling error ≤ 1e-7 per step) |
| Steps with a reflection inside the layer | <n_steps_with_reflection_in_layer> |
| Outcome | <PASS/FAIL> |
```

Fill each `<…>` with the value from `summary.json` (the "Outcome" and summary entries from its `pass` field).

- [ ] **Step 8: Commit**

```bash
git add research/phase2/p2_2_9_layer_trace.jl research/phase2/p2_2_9_layer_trace_lib.jl research/phase2/test/test_trace.jl research/results/phase2/p2_2_9_layer_trace research/phase2_record.md
git commit -m "Add P2.2.9 single-spin layer trace and record it

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Noise floor and spin count (P2.1.4)

**Files:**
- Create: `research/phase2/noise.jl`, `research/phase2/test/test_noise.jl`, `research/phase2/p2_1_4_noise_floor.jl`
- Modify: `research/phase2_record.md`, `research/project-progress-tracker.md`

**Interfaces:**
- Consumes: `run_walls`, `layer_walls`, `H_GRID`, `RATE_REF`, `NSPINS_REF` (Task 1).
- Produces: `required_spins(nspins, separation_in_sem; target=10)::Int`, which gives the spins per seed needed so that the closest pair of neighbouring h values is separated by `target` SEM (SEM ∝ 1/√N). Also `summary.json` with `sigma_per_readout` for h = 0.1 and 0.5, the per-pair separations, and `recommended_nspins`.

- [ ] **Step 1: Write the failing test**

`research/phase2/test/test_noise.jl`:

```julia
include(joinpath(@__DIR__, "..", "noise.jl"))

@testset "noise: required spins" begin
    @test required_spins(100_000, 20.0) == 100_000          # already ≥ 10 SEM: keep
    @test required_spins(100_000, 10.0) == 100_000
    @test required_spins(100_000, 5.0) == 400_000           # SEM ∝ 1/√N: half the separation → 4× spins
    @test required_spins(10_000, 2.5; target=10) == 160_000
end
```

`runtests.jl` includes every `test_*.jl`, so no change is needed there.

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `could not open file .../research/phase2/noise.jl`.

- [ ] **Step 3: Create `research/phase2/noise.jl`**

```julia
"""
    required_spins(nspins, separation_in_sem; target=10)

Spins per seed needed so that a separation currently equal to `separation_in_sem` SEM (measured with
`nspins` spins per seed) becomes at least `target` SEM. SEM ∝ 1/√N, so N_new = N·(target/separation)².
Never recommends fewer spins than `nspins`.
"""
function required_spins(nspins::Integer, separation_in_sem::Real; target=10)
    separation_in_sem >= target && return Int(nspins)
    return ceil(Int, nspins * (target / separation_in_sem)^2)
end
```

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS.

- [ ] **Step 5: Create the script `research/phase2/p2_1_4_noise_floor.jl`**

```julia
# P2.1.4 — noise floor of the layer model and the spin count for the h sweep.
#
# (1) σ and SEM of S(t) at every readout for h = 0.1 and 0.5 µm (ΔR2(0) = 0.1 /ms), 10 seeds × NSPINS.
#     At h = 0 σ is exactly 0 (B1), so it is not run here.
# (2) Pilot of the P2.3.4(b) sweep over H_GRID with NSEEDS_PILOT seeds: for each neighbouring pair of h, the
#     separation of 1 − S(50) in units of the combined SEM of 10-seed means (pilot σ / √10).
#     The tracker asks for > 10 SEM between neighbours; recommended_nspins = required_spins(...) for the worst pair.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_1_4_noise_floor.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "noise.jl"))

const OUT = phase2_outdir("p2_1_4_noise_floor")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const NSEEDS_PILOT = parse(Int, get(ENV, "NSEEDS_PILOT", "3"))

floor_runs = Dict{String, Any}()
for h in (0.1, 0.5)
    r = run_walls(layer_walls(h=h, rate=RATE_REF); nspins=NSPINS)
    floor_runs[string(h)] = Dict("mean" => r.mean, "sigma" => r.sigma, "sem" => r.sem, "runtime_s" => r.runtime)
    @printf("h = %.2f: S(50) = %.5f, σ(50) = %.2e, SEM(50) = %.2e, max σ over t = %.2e, runtime %.0f s\n",
        h, r.mean[end], r.sigma[end], r.sem[end], maximum(r.sigma), r.runtime)
end

pilot = [run_walls(layer_walls(h=h, rate=RATE_REF); nspins=NSPINS, seeds=1:NSEEDS_PILOT) for h in H_GRID]
att = [1 - p.mean[end] for p in pilot]
sem10 = [p.sigma[end] / sqrt(10) for p in pilot]
pairs = NamedTuple[]
for k in 1:length(H_GRID) - 1
    sep = abs(att[k + 1] - att[k]) / sqrt(sem10[k]^2 + sem10[k + 1]^2 + eps())
    push!(pairs, (h_low=H_GRID[k], h_high=H_GRID[k + 1], attenuation_low=att[k], attenuation_high=att[k + 1], separation_in_sem=sep))
    @printf("  h %.3f → %.3f: 1-S(50) %.5f → %.5f, separation %.1f SEM\n", H_GRID[k], H_GRID[k + 1], att[k], att[k + 1], sep)
end
worst = minimum(p.separation_in_sem for p in pairs)
recommended = required_spins(NSPINS, worst)
@printf("worst neighbouring separation %.1f SEM → recommended spins per seed: %d\n", worst, recommended)
write_csv(joinpath(OUT, "pilot_pairs.csv"), pairs)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS,
    "nseeds_pilot" => NSEEDS_PILOT, "h_grid" => H_GRID, "rate" => RATE_REF, "times" => TIMES_REF,
    "floor" => floor_runs, "pilot_S50" => [p.mean[end] for p in pilot],
    "worst_separation_in_sem" => worst, "recommended_nspins" => recommended))
println("saved to $OUT")
```

- [ ] **Step 6: Run the script**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_1_4_noise_floor.jl`
Expected: two noise-floor lines, 13 pair lines, and a recommended spin count. This takes about 15–20 min. The pilot also shows whether 1 − S(50) over `H_GRID` covers the B4 range 0.048–0.912: the h = 1.0 entry should give S(50) ≈ e⁻⁵ = 0.0067, and the smallest h > 0 should give S(50) above 0.952. If it doesn't cover that range, record it in the ledger as a ruling and extend `H_GRID` before Plan B.

- [ ] **Step 7: Record and tick the tracker**

Append to `research/phase2_record.md` (summary row + section), using the values from `summary.json`:

```markdown
| P2.1.4 | Noise floor and spin count | neighbouring h in the sweep differ by > 10 SEM | σ(50) = <floor.0.1.sigma[end]> (h = 0.1), <floor.0.5.sigma[end]> (h = 0.5); worst neighbour separation <worst_separation_in_sem> SEM → <recommended_nspins> spins per seed |
```

```markdown
## P2.1.4: noise floor and spin count

Script [p2_1_4_noise_floor.jl](phase2/p2_1_4_noise_floor.jl) → [results/phase2/p2_1_4_noise_floor/](results/phase2/p2_1_4_noise_floor/)

| h (µm) | S(50) | σ(50) (one seed) | SEM(50) (10 seeds) | max σ over readouts |
|---|---|---|---|---|
| 0.1 | <…> | <…> | <…> | <…> |
| 0.5 | <…> | <…> | <…> | <…> |

Pilot sweep (`NSEEDS_PILOT` seeds): the worst separation between neighbouring h values was <worst_separation_in_sem> SEM (pair <h_low>–<h_high> from `pilot_pairs.csv`). Recommended spins per seed for P2.3.4(b) and V4: **<recommended_nspins>**. H_GRID S(50) range: <min> – <max> (B4 range 0.088–0.952).
```

Fill each `<…>` from `summary.json` / `pilot_pairs.csv`. Then, in `research/project-progress-tracker.md`, set P2.1.3, P2.1.4, P2.1.5 and P2.2.9 to ☑ and add the outputs (`research/phase2/harness.jl`, `expected.jl`, the two run folders). Increase the Phase 2 "Done" count by 4.

- [ ] **Step 8: Commit**

```bash
git add research/phase2/noise.jl research/phase2/test/test_noise.jl research/phase2/p2_1_4_noise_floor.jl research/results/phase2/p2_1_4_noise_floor research/phase2_record.md research/project-progress-tracker.md
git commit -m "Add P2.1.4 noise floor and record harness tasks

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
