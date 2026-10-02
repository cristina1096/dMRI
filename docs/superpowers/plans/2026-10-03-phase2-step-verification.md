# Phase 2 Step-Profile Verification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Verify and characterise the step-profile layer on walls at full size (tracker section 2.3: P2.3.1–P2.3.7, P2.3.9, P2.3.11, P2.3.12), fit h\*(θ) to the B4 baseline reference, and tag `rd-step-v1`.

**Architecture:** Each tracker test is one script in `research/phase2/`, built on the Plan A harness (`run_walls`, `layer_walls`, `load_b4`, …). Library functions that need logic (curve interpolation, the h\* fit, curvature, the toy random walk) go in `fit.jl` and `toy.jl` with unit tests. Scripts write to `research/results/phase2/<run>/` and print PASS/FAIL against the tracker criterion. Results go into `research/phase2_record.md`.

**Tech Stack:** Julia 1.12, MCMRSimulator (this repo, layer from plan `2026-10-01-near-surface-layer-walls.md`), CairoMakie, JSON, `Test` stdlib.

**Spec:** `research/project-progress-tracker.md` section 2.3 (rows P2.3.1–P2.3.12, the Gate), the "Reference configuration" and "Exact expectations" tables; proposal §3.2 VI (monotonicity); the user's framing: the baseline is a reference to fit, not the truth, and there are no fast-diffusion arguments.

**Depends on:** plan `2026-10-03-phase2-harness-and-trace.md` (all four tasks done). That plan provides `harness.jl`, `expected.jl`, the P2.2.9 trace result, and `research/results/phase2/p2_1_4_noise_floor/summary.json` with `recommended_nspins`.

## Global Constraints

- Work only on branch `cc/near-surface-r2`.
- Reference configuration and ΔR₂(0) = 0.1 ms⁻¹ as in Plan A's Global Constraints; τ = 1e-2 ms unless a task says otherwise.
- Uncertainty of a 10-seed mean is the SEM. Comparisons between two means use √(SEM₁² + SEM₂²).
- The baseline is a reference, not the truth: V4 residuals are reported, never pass/fail.
- No fast-diffusion or Brownstein–Tarr formulas in code or records. Do not assume the sign of the endpoint-method bias (V11); measure it.
- Run research tests with: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Results at τ = 1e-2 ms are provisional

The layer's own timestep convergence (V1, P2.6) is not done yet. Every result here is labelled "τ = 1e-2 ms". Task 4 also measures h\* at τ = 1e-3 ms for three θ, to show how much the fit moves with the timestep.

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `research/phase2/fit.jl` | Create | `interp_curve`, `fit_h`, `jackknife_h`, `curvature`, `increasing_steps`, `sweep_nspins` |
| `research/phase2/toy.jl` | Create | `toy_layer_signal`: 1D walk between walls, exact-overlap or endpoint layer rule |
| `research/phase2/test/test_fit.jl`, `test_toy.jl` | Create | unit tests |
| `research/phase2/p2_3_exact_checks.jl` | Create | P2.3.1 V0b, P2.3.4(a), P2.3.11 |
| `research/phase2/p2_3_2_density_occupancy.jl` | Create | P2.3.2 V9, P2.3.3 occupancy (with the layer on) |
| `research/phase2/p2_3_4_h_sweep.jl` | Create | P2.3.4(b), P2.3.5, P2.3.6 |
| `research/phase2/p2_3_9_fit_h.jl` | Create | P2.3.9 V4 h\*(θ) fit |
| `research/phase2/p2_3_7_endpoint.jl` | Create | P2.3.7 V11 |
| `research/phase2/p2_3_12_sides.jl` | Create | P2.3.12 side-specific |
| `research/phase2_record.md`, `research/project-progress-tracker.md` | Modify | records, ticks |

## Review Focus

1. **h\* at the edge of the grid:** a θ whose best h is 0 or 1.0 µm is not "found"; it means the grid or ΔR₂(0) can't reach it. Test: Task 4's `fit_h` test with a target outside the model range expects `at_edge = true`.
2. **Readouts with zero uncertainty** (t = 0, where S = 1 exactly in both models) must be excluded from χ², not divide by zero. Test: Task 4's `fit_h` test.
3. **The paired-seed caveat:** model and reference share initial positions (seeds 1–10), so their errors are correlated and √(SEM₁² + SEM₂²) overestimates the noise of the difference. Records must say so wherever it matters (V4, V11, mirror). There's no test; it's a record requirement checked in Task 7.
4. **The toy walk must reproduce MCMR's 1D motion:** same D, w, reflection and timestep. Test: Task 5 compares the toy's exact method with MCMR's sweep at h = 0.1 µm, τ = 1e-2 ms, within 2 combined SEM.
5. **Do not tag `rd-step-v1` if a gate item failed:** Task 7 checks the three gate results from their `summary.json` files before tagging.

---

### Task 1: Exact checks at full size (P2.3.1 V0b, P2.3.4(a), P2.3.11)

**Files:**
- Create: `research/phase2/p2_3_exact_checks.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `run_walls(...; snapshots=true)`, `layer_walls`, `expected_half_gap`, `roundoff_bound` (Plan A).
- Produces: `research/results/phase2/p2_3_exact_checks/summary.json` with entries `V0b_bulk_off`, `V0b_bulk_on`, `half_gap_bulk_off`, `half_gap_bulk_on`, `additivity`, each `{pass::Bool, ...numbers}`.

- [ ] **Step 1: Create the script**

```julia
# P2.3.1 V0b, P2.3.4(a) exact end point, P2.3.11 additivity — full size (seeds 1–10, NSPINS spins).
#
# V0b       : rho = 0 (layer_h = 0.5) vs plain Walls(repeats=2), default τ (as B1/B2), R2_bulk = 0 and 1/80:
#             final positions and per-spin M⊥ bit-identical (==), every seed.
# h = w/2   : h = 1.0 µm, ΔR2(0) = 0.1, τ = 1e-2, R2_bulk = 0 and 1/80: per spin |M⊥/exp(-(ΔR2(0)+R2_bulk)·50) − 1| ≤ 1e-10;
#             ensemble at every readout within roundoff_bound.
# additivity: h = 0.2 µm, ΔR2(0) = 0.1, τ = 1e-2, R2_bulk = 0 vs 1/80, same seeds: positions ==, per spin
#             |M⊥_bulk / (M⊥_layer·exp(-50/80)) − 1| ≤ 1e-10; ensemble within roundoff_bound.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_exact_checks.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "expected.jl"))

const OUT = phase2_outdir("p2_3_exact_checks")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const TOL_SPIN = 1e-10
per_spin_final(r) = [[s.orientations[1].transverse for s in snap.spins] for snap in r.snapshots]
positions_final(r) = [mr.position.(snap) for snap in r.snapshots]

results = Dict{String, Any}()
function report(name, pass, extra)
    results[name] = merge(Dict{String, Any}("pass" => pass), extra)
    @printf("%-18s → %s  %s\n", name, pass ? "PASS" : "FAIL", string(extra))
end

for (label, R2) in (("bulk_off", 0.0), ("bulk_on", 1 / 80))
    a = run_walls(Walls(repeats=W_REF); R2_bulk=R2, timestep=nothing, nspins=NSPINS, snapshots=true)
    b = run_walls(layer_walls(h=0.5, rate=0.0); R2_bulk=R2, timestep=nothing, nspins=NSPINS, snapshots=true)
    same = positions_final(a) == positions_final(b) && per_spin_final(a) == per_spin_final(b) && a.S == b.S
    report("V0b_" * label, same, Dict("tau_ms" => a.tau, "R2_bulk" => R2))
end

bound = roundoff_bound(TIMES_REF, TAU_REF, NSPINS)
for (label, R2) in (("bulk_off", 0.0), ("bulk_on", 1 / 80))
    r = run_walls(layer_walls(h=1.0, rate=RATE_REF); R2_bulk=R2, nspins=NSPINS, snapshots=true)
    target = expected_half_gap(TIMES_REF, RATE_REF; R2_bulk=R2)
    spin_dev = maximum(maximum(abs.(m ./ target[end] .- 1)) for m in per_spin_final(r))
    ens_dev = maximum(abs.(r.S ./ target' .- 1))
    report("half_gap_" * label, spin_dev <= TOL_SPIN && ens_dev <= bound,
        Dict("max_rel_dev_per_spin" => spin_dev, "max_rel_dev_ensemble" => ens_dev, "bound" => bound, "R2_bulk" => R2))
end

layer_only = run_walls(layer_walls(h=0.2, rate=RATE_REF); nspins=NSPINS, snapshots=true)
with_bulk = run_walls(layer_walls(h=0.2, rate=RATE_REF); R2_bulk=1 / 80, nspins=NSPINS, snapshots=true)
same_pos = positions_final(layer_only) == positions_final(with_bulk)
spin_dev = maximum(maximum(abs.(mb ./ (ml .* exp(-TIMES_REF[end] / 80)) .- 1))
    for (mb, ml) in zip(per_spin_final(with_bulk), per_spin_final(layer_only)))
ens_dev = maximum(abs.(with_bulk.S ./ (layer_only.S .* exp.(-TIMES_REF' ./ 80)) .- 1))
report("additivity", same_pos && spin_dev <= TOL_SPIN && ens_dev <= bound,
    Dict("positions_identical" => same_pos, "max_rel_dev_per_spin" => spin_dev, "max_rel_dev_ensemble" => ens_dev,
        "bound" => bound, "S50_layer_only" => layer_only.mean[end]))

write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "results" => results))
println("saved to $OUT")
```

- [ ] **Step 2: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_3_exact_checks.jl`
Expected: five lines (`V0b_bulk_off`, `V0b_bulk_on`, `half_gap_bulk_off`, `half_gap_bulk_on`, `additivity`) all `→ PASS`. Runtime is about 15 min. Any FAIL: stop and use superpowers:systematic-debugging. These are exact identities, so a failure is a code bug, not noise.

- [ ] **Step 3: Record**

Append three summary rows and one section to `research/phase2_record.md`, filled from `results/phase2/p2_3_exact_checks/summary.json`:

```markdown
| P2.3.1 | V0b: ρ = 0 vs no layer (bulk off / on), default τ | bit-identical positions and M⊥ | <PASS/FAIL> / <PASS/FAIL> |
| P2.3.4(a) | Exact end point h = w/2 (bulk off / on) | per spin ≤ 1e-10; ensemble ≤ roundoff bound | <PASS/FAIL>: per spin <max_rel_dev_per_spin>, ensemble <max_rel_dev_ensemble> |
| P2.3.11 | Additivity R₂_bulk + ΔR₂ (h = 0.2 µm) | same positions; per spin ≤ 1e-10; ensemble ≤ roundoff bound | <PASS/FAIL>: per spin <…>, ensemble <…> |
```

```markdown
## P2.3.1, P2.3.4(a), P2.3.11: exact checks at full size

Script [p2_3_exact_checks.jl](phase2/p2_3_exact_checks.jl) → [results/phase2/p2_3_exact_checks/](results/phase2/p2_3_exact_checks/)

These hold in any diffusion regime, so they are pass/fail to floating-point precision. Same seeds (1–10) and <nspins> spins per seed throughout.

| Check | Configuration | Max deviation per spin | Max deviation ensemble | Bound | Outcome |
|---|---|---|---|---|---|
| V0b, bulk off | ρ = 0, h = 0.5 vs no layer; default τ | identical (==) | identical | — | <…> |
| V0b, bulk on (1/80) | same | identical | identical | — | <…> |
| h = w/2, bulk off | h = 1.0, ΔR₂(0) = 0.1, τ = 1e-2 | <…> | <…> | <bound> | <…> |
| h = w/2, bulk on | same + R₂_bulk = 1/80 | <…> | <…> | <bound> | <…> |
| additivity | h = 0.2, ΔR₂(0) = 0.1, R₂_bulk 0 vs 1/80 | <…> | <…> | <bound> | <…> |
```

- [ ] **Step 4: Commit**

```bash
git add research/phase2/p2_3_exact_checks.jl research/results/phase2/p2_3_exact_checks research/phase2_record.md
git commit -m "Add Phase 2 exact checks (V0b, h = w/2, additivity) and record them

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Density and occupancy with the layer on (P2.3.2 V9, P2.3.3) — gate

**Files:**
- Create: `research/phase2/p2_3_2_density_occupancy.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `run_walls(...; snapshots=true)`, `layer_walls`, `occupancy` (Plan A).
- Produces: `results/phase2/p2_3_2_density_occupancy/summary.json` with, per τ, `positions_identical_to_no_layer`, `histogram_z_chi2`, `occupancy` rows, and `pass`.

**Method.** The layer only changes |Mxy|, never the path. So with the same seeds, positions with the layer on must equal positions without it, which makes the density and occupancy identical to B3. The script checks that equality first, then recomputes the B3 statistics (histogram flatness and occupancy = 2h/w) on the layer-on positions at t = 50 ms.

- [ ] **Step 1: Create the script**

```julia
# P2.3.2 V9 equilibrium density and P2.3.3 layer occupancy, with the layer on (h = 0.5 µm, ΔR2(0) = 0.1).
#
# For τ = default (0.04 ms) and 1e-2 ms, seeds 1–10, NSPINS spins:
#  (1) final positions with the layer == final positions without it (same seeds)  [layer never moves spins]
#  (2) histogram of x mod w (100 bins of 0.02 µm) pooled over seeds: χ² vs uniform, z = (χ² − dof)/√(2 dof), |z| ≤ 3
#  (3) occupancy within h of either face for h ∈ HS: z = (n − N·2h/w)/√(N·p(1−p)), |z| ≤ 2 for every h
# Pass (gate P2.3.3): (1) and (2) and (3) for both τ.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_2_density_occupancy.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_3_2_density_occupancy")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const NBINS = 100
const HS = [0.05, 0.1, 0.2, 0.4, 0.5]

summary = Dict{String, Any}()
rows = NamedTuple[]
all_pass = true
for (label, τ) in (("default", nothing), ("tau_1e-2", 1e-2))
    with = run_walls(layer_walls(h=0.5, rate=RATE_REF); timestep=τ, times=[0.0, 50.0], nspins=NSPINS, snapshots=true)
    without = run_walls(Walls(repeats=W_REF); timestep=τ, times=[0.0, 50.0], nspins=NSPINS, snapshots=true)
    identical = [mr.position.(a) for a in with.snapshots] == [mr.position.(b) for b in without.snapshots]
    u = vcat([[mod(p[1], W_REF) for p in mr.position.(s)] for s in with.snapshots]...)
    ntot = length(u)
    hist = zeros(Int, NBINS)
    for x in u
        hist[min(NBINS, floor(Int, x / W_REF * NBINS) + 1)] += 1
    end
    e = ntot / NBINS
    χ2 = sum((hist .- e) .^ 2 ./ e)
    zχ = (χ2 - (NBINS - 1)) / sqrt(2 * (NBINS - 1))
    occ_ok = true
    for h in HS
        p = 2h / W_REF
        n = count(x -> x <= h || W_REF - x <= h, u)
        z = (n - ntot * p) / sqrt(ntot * p * (1 - p))
        occ_ok &= abs(z) <= 2
        push!(rows, (config=label, tau_ms=with.tau, h_um=h, expected_fraction=p, measured_fraction=n / ntot, z=z))
        @printf("  %-8s h=%.2f: fraction %.5f vs %.5f (z %+.2f)\n", label, h, n / ntot, p, z)
    end
    pass = identical && abs(zχ) <= 3 && occ_ok
    global all_pass &= pass
    summary[label] = Dict("tau_ms" => with.tau, "positions_identical_to_no_layer" => identical,
        "histogram_chi2_per_dof" => χ2 / (NBINS - 1), "histogram_z_chi2" => zχ, "occupancy_all_within_2sigma" => occ_ok, "pass" => pass)
    @printf("%-8s: positions identical %s, χ²/dof %.3f (z %+.2f), occupancy ok %s → %s\n",
        label, identical, χ2 / (NBINS - 1), zχ, occ_ok, pass ? "PASS" : "FAIL")
end
write_csv(joinpath(OUT, "occupancy.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "nbins" => NBINS,
    "h_values" => HS, "results" => summary, "pass" => all_pass))
println("saved to $OUT")
```

- [ ] **Step 2: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_3_2_density_occupancy.jl`
Expected: both `default` and `tau_1e-2` lines `→ PASS`, with positions identical. About 5 min. A FAIL on "positions identical" is a code bug: the layer must never change trajectories. Stop and debug.

- [ ] **Step 3: Record**

Summary rows and section in `research/phase2_record.md`, from `summary.json` and `occupancy.csv`:

```markdown
| P2.3.2 | V9 density with the layer on (τ = 0.04, 1e-2) | positions == no-layer run; histogram \|z_χ²\| ≤ 3 | <PASS/FAIL> |
| P2.3.3 | Occupancy within h of either face (**gate**) | fraction = 2h/w, \|z\| ≤ 2 for h = 0.05–0.5 | <PASS/FAIL> |
```

```markdown
## P2.3.2, P2.3.3: density and occupancy with the layer on

Script [p2_3_2_density_occupancy.jl](phase2/p2_3_2_density_occupancy.jl) → [results/phase2/p2_3_2_density_occupancy/](results/phase2/p2_3_2_density_occupancy/)

The layer changes |Mxy| only, so final positions with the layer on are compared with the no-layer run (same seeds); the B3 statistics are then recomputed on the layer-on positions at t = 50 ms (10⁶ spins pooled).

| τ (ms) | Positions identical to no layer | χ²/dof (z) | Occupancy z for h = 0.05 / 0.1 / 0.2 / 0.4 / 0.5 | Outcome |
|---|---|---|---|---|
| 0.04 (default) | <…> | <…> (<…>) | <…> | <…> |
| 1e-2 | <…> | <…> (<…>) | <…> | <…> |
```

- [ ] **Step 4: Commit**

```bash
git add research/phase2/p2_3_2_density_occupancy.jl research/results/phase2/p2_3_2_density_occupancy research/phase2_record.md
git commit -m "Add P2.3.2/P2.3.3 density and occupancy with the layer on

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: h sweep, monotonicity and decay shape (P2.3.4(b), P2.3.5, P2.3.6)

**Files:**
- Create: `research/phase2/fit.jl` (first part), `research/phase2/test/test_fit.jl`, `research/phase2/p2_3_4_h_sweep.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `run_walls`, `layer_walls`, `H_GRID`, `RATE_REF`, `phase2_outdir` (Plan A); `results/phase2/p2_1_4_noise_floor/summary.json` field `recommended_nspins`.
- Produces (in `fit.jl`):
  - `sweep_nspins()::Int`: `recommended_nspins` from P2.1.4, or `NSPINS_REF` if that file is missing
  - `curvature(times, S)::Float64`: c from the least-squares fit ln S = a + b·t + c·t²
  - `increasing_steps(att, sem)::Vector{Bool}`: entry k is true if att[k+1] − att[k] > 2·√(sem_k² + sem_{k+1}²)
- Produces (script): `results/phase2/p2_3_4_h_sweep/sweep.csv` with rows `(h_um, seed, t_ms, signal)`, and `summary.json` with per h: mean, sem, attenuation at 50 ms, curvature mean and SEM.

- [ ] **Step 1: Write the failing tests**

`research/phase2/test/test_fit.jl`:

```julia
@testset "fit: curvature and monotonicity" begin
    t = collect(0.0:5.0:50.0)
    @test curvature(t, exp.(-0.1 .* t)) ≈ 0.0 atol=1e-12
    @test curvature(t, exp.(-0.1 .* t .+ 0.001 .* t .^ 2)) ≈ 0.001 rtol=1e-8
    @test increasing_steps([0.1, 0.2, 0.2001], [0.001, 0.001, 0.001]) == [true, false]
    @test sweep_nspins() isa Int
end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `UndefVarError: curvature not defined` (`runtests.jl` includes `fit.jl` only once it exists).

- [ ] **Step 3: Create `research/phase2/fit.jl`**

```julia
# Phase 2 analysis helpers: sweep spin count, decay-shape curvature, monotonicity, h* fit (Task 4 adds the fit).

"Spins per seed for the h sweep: `recommended_nspins` from P2.1.4, else NSPINS_REF."
function sweep_nspins()
    f = joinpath(REPO_ROOT, "research", "results", "phase2", "p2_1_4_noise_floor", "summary.json")
    isfile(f) || return NSPINS_REF
    return Int(JSON.parsefile(f)["recommended_nspins"])
end

"c in the least-squares fit ln S = a + b·t + c·t² (curvature of the log-decay; 0 for a single exponential)."
function curvature(times, S)
    X = hcat(ones(length(times)), times, times .^ 2)
    return (X \ log.(S))[3]
end

"Entry k: attenuation rises from step k to k+1 by more than 2 combined SEM."
increasing_steps(att, sem) = [att[k + 1] - att[k] > 2 * sqrt(sem[k]^2 + sem[k + 1]^2) for k in 1:length(att) - 1]
```

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS.

- [ ] **Step 5: Create the sweep script `research/phase2/p2_3_4_h_sweep.jl`**

```julia
# P2.3.4(b) attenuation vs h, P2.3.5 monotonicity, P2.3.6 decay shape.
#
# Step profile, both faces, ΔR2(0) = 0.1 /ms, h ∈ H_GRID, τ = 1e-2 ms (provisional, see plan), seeds 1–10,
# NSPINS = sweep_nspins() spins per seed. Records S(t) per seed (sweep.csv; also the input of the V4 fit).
# P2.3.5 pass: attenuation 1 − S(50) strictly increasing with h, each step by > 2 combined SEM.
# P2.3.6: curvature c of ln S(t) per seed (mean ± SEM) vs h — characterisation, no pass/fail.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_4_h_sweep.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "fit.jl"))

const OUT = phase2_outdir("p2_3_4_h_sweep")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(sweep_nspins())))

rows = NamedTuple[]
per_h = Dict{String, Any}[]
for h in H_GRID
    r = run_walls(layer_walls(h=h, rate=RATE_REF); nspins=NSPINS)
    for s in 1:size(r.S, 1), (j, t) in enumerate(r.times)
        push!(rows, (h_um=h, seed=s, t_ms=t, signal=r.S[s, j]))
    end
    curv = [curvature(r.times, r.S[s, :]) for s in 1:size(r.S, 1)]
    push!(per_h, Dict("h_um" => h, "mean" => r.mean, "sem" => r.sem, "sigma" => r.sigma,
        "attenuation_50" => 1 - r.mean[end], "attenuation_50_sem" => r.sem[end],
        "curvature_mean" => mean(curv), "curvature_sem" => std(curv) / sqrt(length(curv)), "runtime_s" => r.runtime))
    @printf("h = %.3f: S(50) = %.5f ± %.5f (SEM), curvature %.2e ± %.1e, %.0f s\n",
        h, r.mean[end], r.sem[end], mean(curv), std(curv) / sqrt(length(curv)), r.runtime)
end
att = [p["attenuation_50"] for p in per_h]
sem = [p["attenuation_50_sem"] for p in per_h]
steps = increasing_steps(att, sem)
monotone = all(steps)
@printf("P2.3.5 monotonicity: %d of %d steps increase by > 2 SEM → %s\n", count(steps), length(steps), monotone ? "PASS" : "FAIL")
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "rate" => RATE_REF,
    "tau_ms" => TAU_REF, "times" => TIMES_REF, "h_grid" => H_GRID, "per_h" => per_h,
    "monotone_steps" => steps, "monotone_pass" => monotone))

using CairoMakie
fig = Figure(size=(1100, 340))
ax1 = Axis(fig[1, 1], xlabel="t (ms)", ylabel="S(t)", yscale=log10, title="S(t) per h (ΔR₂(0) = 0.1 /ms)")
for p in per_h
    lines!(ax1, TIMES_REF, max.(p["mean"], 1e-6), label="h = $(p["h_um"])")
end
Legend(fig[1, 2], ax1, labelsize=8, rowgap=0)
ax2 = Axis(fig[1, 3], xlabel="h (µm)", ylabel="1 − S(50)", title="P2.3.5 attenuation vs h")
errorbars!(ax2, H_GRID, att, 2 .* sem); scatterlines!(ax2, H_GRID, att)
ax3 = Axis(fig[1, 4], xlabel="h (µm)", ylabel="curvature c of ln S", title="P2.3.6 decay shape")
cm = [p["curvature_mean"] for p in per_h]; cs = [p["curvature_sem"] for p in per_h]
errorbars!(ax3, H_GRID, cm, 2 .* cs); scatterlines!(ax3, H_GRID, cm); hlines!(ax3, [0.0], color=:gray, linestyle=:dash)
save(joinpath(OUT, "p2_3_4_h_sweep.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 6: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_3_4_h_sweep.jl`
Expected: 14 `h = …` lines with S(50) decreasing from 1.0 (h = 0) to about e⁻⁵ = 0.0067 (h = 1.0), and a `P2.3.5 monotonicity … → PASS` line. About 30–60 min, depending on the spin count. If monotonicity FAILs at the smallest h steps only, record which steps and why in the record. Don't change the criterion.

- [ ] **Step 7: Record**

Summary rows and section in `research/phase2_record.md`, from `summary.json`:

```markdown
| P2.3.4(b) | Attenuation vs h (characterisation), τ = 1e-2 | recorded with SEM | S(50) from 1 (h = 0) to <…> (h = 1.0); <nspins> spins × 10 seeds |
| P2.3.5 | Monotonicity in h (proposal §3.2 VI) | every step > 2 combined SEM | <PASS/FAIL>: <count> of 13 steps |
| P2.3.6 | Decay shape (curvature of ln S) | characterisation | curvature from <…> to <…> (see table) |
```

```markdown
## P2.3.4(b), P2.3.5, P2.3.6: h sweep

Script [p2_3_4_h_sweep.jl](phase2/p2_3_4_h_sweep.jl) → [results/phase2/p2_3_4_h_sweep/](results/phase2/p2_3_4_h_sweep/) (`sweep.csv` per seed, figure `p2_3_4_h_sweep.png`). τ = 1e-2 ms (provisional until P2.6).

| h (µm) | S(25) | S(50) ± SEM | 1 − S(50) | curvature c ± SEM (ms⁻²) |
|---|---|---|---|---|
| <one row per H_GRID entry from per_h> |

**Interpretation.** <one paragraph: monotone or not, where the curve is steep/flat, whether ln S(t) is straight (c ≈ 0) or curved and how that changes with h. Numbers only from the table; no formula-based expectation.>
```

- [ ] **Step 8: Commit**

```bash
git add research/phase2/fit.jl research/phase2/test/test_fit.jl research/phase2/p2_3_4_h_sweep.jl research/results/phase2/p2_3_4_h_sweep research/phase2_record.md
git commit -m "Add P2.3.4(b)-P2.3.6 h sweep, monotonicity and decay shape

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: V4: fit h\*(θ) to the B4 baseline reference (P2.3.9)

**Files:**
- Modify: `research/phase2/fit.jl` (append), `research/phase2/test/test_fit.jl` (append)
- Create: `research/phase2/p2_3_9_fit_h.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `sweep.csv` (Task 3), `load_b4`, `reference_at` (Plan A), `run_walls`, `layer_walls`.
- Produces (in `fit.jl`):
  - `interp_curve(hgrid, curves, h)::Vector{Float64}`: per-readout linear interpolation in h, clamped to the grid
  - `fit_h(hgrid, model_mean, model_sem, ref_mean, ref_sem; npts=10_001)`, which returns `(h, chi2, dof, residuals, model, at_edge)`. Readouts where the combined σ is 0 are excluded.
  - `jackknife_h(hgrid, S_by_h, ref_mean, ref_sem)`, which returns `(h, sigma)`. `S_by_h[i]` is a seeds × times matrix for `hgrid[i]`.
  - `read_sweep(path)`, which returns `(hgrid, S_by_h)`

- [ ] **Step 1: Append the failing tests to `research/phase2/test/test_fit.jl`**

```julia
@testset "fit: h* fit" begin
    t = collect(0.0:5.0:50.0)
    hgrid = [0.0, 0.5, 1.0]
    curves = [exp.(-k .* t) for k in (0.0, 0.05, 0.1)]
    sem = [vcat(0.0, fill(1e-3, 10)) for _ in hgrid]          # t = 0: zero uncertainty
    ref = interp_curve(hgrid, curves, 0.3)
    ref_sem = vcat(0.0, fill(1e-3, 10))

    @test interp_curve(hgrid, curves, 0.0) == curves[1]
    @test interp_curve(hgrid, curves, 1.5) == curves[end]      # clamped
    f = fit_h(hgrid, curves, sem, ref, ref_sem)
    @test f.h ≈ 0.3 atol=1e-4
    @test f.chi2 ≈ 0.0 atol=1e-6
    @test f.dof == 10 - 1                                      # 10 readouts with σ > 0, one fitted parameter
    @test !f.at_edge
    @test all(isfinite, f.residuals)

    far = exp.(-0.5 .* t)                                      # decays faster than any model curve
    g = fit_h(hgrid, curves, sem, far, ref_sem)
    @test g.h == 1.0
    @test g.at_edge

    S_by_h = [repeat(c', 10) for c in curves]                  # 10 identical seeds → zero jackknife spread
    j = jackknife_h(hgrid, S_by_h, ref, ref_sem)
    @test j.h ≈ 0.3 atol=1e-4
    @test j.sigma ≈ 0.0 atol=1e-6
end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `UndefVarError: interp_curve not defined`.

- [ ] **Step 3: Append to `research/phase2/fit.jl`**

```julia
"Per-readout linear interpolation in h between grid curves; clamped to the grid ends."
function interp_curve(hgrid, curves, h)
    h <= hgrid[1] && return curves[1]
    h >= hgrid[end] && return curves[end]
    j = searchsortedlast(hgrid, h)
    f = (h - hgrid[j]) / (hgrid[j + 1] - hgrid[j])
    return (1 - f) .* curves[j] .+ f .* curves[j + 1]
end

"""
    fit_h(hgrid, model_mean, model_sem, ref_mean, ref_sem; npts=10_001)

h minimising χ²(h) = Σ_t [(S_model(t; h) − S_ref(t)) / σ_t]², σ_t = √(sem_model² + sem_ref²), over a dense
h grid with the model curve interpolated in h. Readouts with σ_t = 0 are excluded. `at_edge` is true when
the best h is the first or last grid value (the reference is not reached inside the grid).
"""
function fit_h(hgrid, model_mean, model_sem, ref_mean, ref_sem; npts=10_001)
    best_h, best_chi2 = NaN, Inf
    for h in range(hgrid[1], hgrid[end], length=npts)
        m = interp_curve(hgrid, model_mean, h)
        σ = sqrt.(interp_curve(hgrid, model_sem, h) .^ 2 .+ ref_sem .^ 2)
        ok = σ .> 0
        χ2 = sum(((m .- ref_mean)[ok] ./ σ[ok]) .^ 2)
        χ2 < best_chi2 && ((best_h, best_chi2) = (h, χ2))
    end
    m = interp_curve(hgrid, model_mean, best_h)
    σ = sqrt.(interp_curve(hgrid, model_sem, best_h) .^ 2 .+ ref_sem .^ 2)
    ok = σ .> 0
    residuals = [ok[i] ? (m[i] - ref_mean[i]) / σ[i] : 0.0 for i in eachindex(m)]
    return (h=best_h, chi2=best_chi2, dof=count(ok) - 1, residuals=residuals, model=m,
        at_edge=(best_h == hgrid[1] || best_h == hgrid[end]))
end

"""
    jackknife_h(hgrid, S_by_h, ref_mean, ref_sem)

h* from the mean over all seeds, and its jackknife standard error: refit leaving out one seed at a time,
σ = √((n−1)/n · Σ (h₋ₛ − mean h₋)²). `S_by_h[i]` is the seeds × times matrix for hgrid[i].
"""
function jackknife_h(hgrid, S_by_h, ref_mean, ref_sem)
    n = size(S_by_h[1], 1)
    fit_rows(rows) = fit_h(hgrid,
        [vec(mean(S[rows, :], dims=1)) for S in S_by_h],
        [vec(std(S[rows, :], dims=1)) ./ sqrt(length(rows)) for S in S_by_h],
        ref_mean, ref_sem).h
    full = fit_rows(1:n)
    loo = [fit_rows(setdiff(1:n, s)) for s in 1:n]
    return (h=full, sigma=sqrt((n - 1) / n * sum((loo .- mean(loo)) .^ 2)))
end

"Reads sweep.csv (h_um, seed, t_ms, signal) into (hgrid, S_by_h)."
function read_sweep(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    col(name) = findfirst(==(name), header)
    data = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    hgrid = sort(unique(r[col("h_um")] for r in data))
    S_by_h = map(hgrid) do h
        rs = filter(r -> r[col("h_um")] == h, data)
        seeds = sort(unique(Int(r[col("seed")]) for r in rs))
        times = sort(unique(r[col("t_ms")] for r in rs))
        M = zeros(length(seeds), length(times))
        for r in rs
            M[findfirst(==(Int(r[col("seed")])), seeds), findfirst(==(r[col("t_ms")]), times)] = r[col("signal")]
        end
        M
    end
    return (hgrid, S_by_h)
end
```

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS (testsets "fit: curvature and monotonicity" and "fit: h* fit").

- [ ] **Step 5: Create the script `research/phase2/p2_3_9_fit_h.jl`**

```julia
# P2.3.9 V4 — fit h*(θ) of the step-profile layer (ΔR2(0) = 0.1 /ms, both faces) to the B4 baseline reference.
#
# (1) For each θ > 0 in B4 (τ = 1e-2 curves): h* from the P2.3.4(b) sweep (interpolated in h), jackknife σ over seeds,
#     χ²/dof and per-readout residuals. at_edge = reference not reachable inside the h grid.
# (2) Confirmation: direct run at h*(θ) (no interpolation), residuals vs the reference.
# (3) Timestep: direct run at h*(θ) with τ = 1e-3 for θ ∈ TAU_CHECK; difference of S(50) from the τ = 1e-2 run, in SEM.
# (4) Bulk on: θ = 0.01, direct run at h* with R2_bulk = 1/80 vs reference × exp(−t/80) (B4 showed this factorises).
# The baseline is a reference, not the truth: residuals are reported, never pass/fail.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_9_fit_h.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "expected.jl"))
include(joinpath(@__DIR__, "fit.jl"))

const OUT = phase2_outdir("p2_3_9_fit_h")
const SWEEP = joinpath(REPO_ROOT, "research", "results", "phase2", "p2_3_4_h_sweep", "sweep.csv")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(sweep_nspins())))
const TAU_CHECK = [0.002, 0.01, 0.05]

ref = load_b4()
(hgrid, S_by_h) = read_sweep(SWEEP)
model_mean = [vec(mean(S, dims=1)) for S in S_by_h]
model_sem = [vec(std(S, dims=1)) ./ sqrt(size(S, 1)) for S in S_by_h]

fits = Dict{String, Any}[]
rows = NamedTuple[]
for θ in sort([x for x in keys(ref) if x > 0])
    f = fit_h(hgrid, model_mean, model_sem, ref[θ].mean, ref[θ].sem)
    j = jackknife_h(hgrid, S_by_h, ref[θ].mean, ref[θ].sem)
    direct = run_walls(layer_walls(h=f.h, rate=RATE_REF); nspins=NSPINS)
    σd = sqrt.(direct.sem .^ 2 .+ ref[θ].sem .^ 2)
    res_direct = [σd[i] > 0 ? (direct.mean[i] - ref[θ].mean[i]) / σd[i] : 0.0 for i in eachindex(σd)]
    entry = Dict{String, Any}("theta" => θ, "h_star" => f.h, "h_star_jackknife_sigma" => j.sigma, "at_edge" => f.at_edge,
        "chi2" => f.chi2, "dof" => f.dof, "residuals_interp" => f.residuals, "residuals_direct" => res_direct,
        "S_direct" => direct.mean, "S_direct_sem" => direct.sem, "S_ref" => ref[θ].mean)
    if θ in TAU_CHECK
        fine = run_walls(layer_walls(h=f.h, rate=RATE_REF); timestep=1e-3, nspins=NSPINS)
        entry["S50_tau_1e-3"] = fine.mean[end]
        entry["S50_tau_1e-3_sem"] = fine.sem[end]
        entry["tau_shift_in_sem"] = (fine.mean[end] - direct.mean[end]) / sqrt(fine.sem[end]^2 + direct.sem[end]^2)
    end
    if θ == 0.01
        bulk = run_walls(layer_walls(h=f.h, rate=RATE_REF); R2_bulk=1 / 80, nspins=NSPINS)
        target = ref[θ].mean .* exp.(-TIMES_REF ./ 80)
        σb = sqrt.(bulk.sem .^ 2 .+ (ref[θ].sem .* exp.(-TIMES_REF ./ 80)) .^ 2)
        entry["residuals_bulk_on"] = [σb[i] > 0 ? (bulk.mean[i] - target[i]) / σb[i] : 0.0 for i in eachindex(σb)]
    end
    push!(fits, entry)
    push!(rows, (theta=θ, h_star=f.h, h_star_sigma=j.sigma, at_edge=f.at_edge, chi2_per_dof=f.chi2 / max(f.dof, 1),
        max_abs_residual_direct=maximum(abs.(res_direct)), S50_ref=ref[θ].mean[end], S50_model=direct.mean[end]))
    @printf("θ = %.4f: h* = %.4f ± %.4f µm%s, χ²/dof = %.2f, max |residual| (direct) = %.2f\n",
        θ, f.h, j.sigma, f.at_edge ? " (AT EDGE)" : "", f.chi2 / max(f.dof, 1), maximum(abs.(res_direct)))
end
write_csv(joinpath(OUT, "h_star.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rate" => RATE_REF, "tau_ms" => TAU_REF,
    "nspins" => NSPINS, "h_grid" => hgrid, "fits" => fits,
    "note" => "model and reference share seeds 1–10 (same initial positions): combined SEM overestimates the noise of the difference"))

using CairoMakie
fig = Figure(size=(1000, 360))
ax1 = Axis(fig[1, 1], xlabel="θ_relax", ylabel="h* (µm)", xscale=log10, title="V4: h*(θ), ΔR₂(0) = 0.1 /ms, τ = 1e-2")
θs = [r.theta for r in rows]
errorbars!(ax1, θs, [r.h_star for r in rows], 2 .* [r.h_star_sigma for r in rows]); scatterlines!(ax1, θs, [r.h_star for r in rows])
ax2 = Axis(fig[1, 2], xlabel="t (ms)", ylabel="(model − reference) / σ", title="residuals of direct runs")
for e in fits
    lines!(ax2, TIMES_REF, e["residuals_direct"], label="θ = $(e["theta"])")
end
hlines!(ax2, [-2.0, 2.0], color=:gray, linestyle=:dash)
Legend(fig[1, 3], ax2, labelsize=8, rowgap=0)
save(joinpath(OUT, "p2_3_9_fit_h.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 6: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_3_9_fit_h.jl`
Expected: 12 `θ = …` lines, each with an h\* and a jackknife σ; none `(AT EDGE)` if the P2.1.4 pilot confirmed the grid covers the B4 range. About 1–1.5 h (12 direct runs plus 3 at τ = 1e-3 plus 1 with bulk). An `AT EDGE` line is a result, not a failure: record it (that θ isn't reachable with ΔR₂(0) = 0.1 inside h ≤ 1 µm).

- [ ] **Step 7: Record**

Summary row and section in `research/phase2_record.md`, from `summary.json` and `h_star.csv`:

```markdown
| P2.3.9 | V4: h*(θ) fit to the B4 reference (τ = 1e-2, ΔR₂(0) = 0.1) | h* found for every θ, with uncertainty and residuals (no pass/fail) | <n found> of 12 θ inside the grid; h* from <…> to <…> µm |
```

```markdown
## P2.3.9: V4, h*(θ) fitted to the baseline reference

Script [p2_3_9_fit_h.jl](phase2/p2_3_9_fit_h.jl) → [results/phase2/p2_3_9_fit_h/](results/phase2/p2_3_9_fit_h/) (`h_star.csv`, figure `p2_3_9_fit_h.png`). Step profile, both faces, ΔR₂(0) = 0.1 ms⁻¹, τ = 1e-2 ms (provisional until P2.6). The baseline is a reference, not the truth: residuals show where the two models differ.

| θ_relax | h* ± σ_jack (µm) | χ²/dof | max \|residual\| (direct run) | S(50) reference | S(50) model at h* |
|---|---|---|---|---|---|
| <one row per h_star.csv row> |

| Timestep check (h* fixed) | θ = 0.002 | θ = 0.01 | θ = 0.05 |
|---|---|---|---|
| S(50) shift τ = 1e-2 → 1e-3, in SEM | <…> | <…> | <…> |

Bulk on (θ = 0.01, R₂_bulk = 1/80): max |residual| <…> (compared with reference × e^(−t/80)).

Caveat: model and reference share seeds 1–10 (same initial positions), so √(SEM₁² + SEM₂²) overestimates the noise of their difference; residuals are conservative.

**Interpretation.** <one paragraph: does one h reproduce the whole curve (residuals flat, χ²/dof ≈ 1) or only part of it (residuals with a trend in t); how h* grows with θ; how much h* depends on τ. Numbers only.>
```

- [ ] **Step 8: Commit**

```bash
git add research/phase2/fit.jl research/phase2/test/test_fit.jl research/phase2/p2_3_9_fit_h.jl research/results/phase2/p2_3_9_fit_h research/phase2_record.md
git commit -m "Add P2.3.9 V4 fit of h*(theta) to the B4 reference

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: V11: exact overlap vs endpoint rule (P2.3.7)

**Files:**
- Create: `research/phase2/toy.jl`, `research/phase2/test/test_toy.jl`, `research/phase2/p2_3_7_endpoint.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `W_REF`, `D_REF`, `phase2_outdir` (Plan A); `mr.Geometries.Internal.Layers.segment_fraction`; `sweep.csv` (Task 3) for the cross-check at h = 0.1.
- Produces: `toy_layer_signal(; h, rate, tau, T=50.0, nspins, seed, method)::Vector{Float64}`, the per-spin final signal of a 1D walk between walls at 0 and w (reflecting). `method = :exact` uses the exact overlap on every reflected piece; `method = :endpoint` adds rate·τ·(number of layers containing the step's end point).

- [ ] **Step 1: Write the failing tests**

`research/phase2/test/test_toy.jl`:

```julia
@testset "toy: exact vs endpoint layer rule" begin
    # h = w/2: every point is in exactly one layer, so both rules give exp(-rate·T) for every spin
    for m in (:exact, :endpoint)
        s = toy_layer_signal(; h=1.0, rate=0.1, tau=0.01, T=5.0, nspins=200, seed=1, method=m)
        @test all(isapprox.(s, exp(-0.1 * 5.0); rtol=1e-10))
    end
    # no layer
    @test all(toy_layer_signal(; h=0.0, rate=0.1, tau=0.01, T=5.0, nspins=50, seed=1, method=:exact) .== 1.0)
    # reproducible for a given seed
    a = toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=100, seed=3, method=:exact)
    b = toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=100, seed=3, method=:exact)
    @test a == b
    @test_throws ErrorException toy_layer_signal(; h=0.2, rate=0.1, tau=0.01, T=5.0, nspins=10, seed=1, method=:midpoint)
end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `UndefVarError: toy_layer_signal not defined`.

- [ ] **Step 3: Create `research/phase2/toy.jl`**

```julia
# 1D toy random walk between reflecting walls at 0 and w, with a step-profile layer of thickness h on both
# faces, for the V11 comparison (P2.3.7). The x-motion of MCMR's 3D walk between Walls(repeats=w) is exactly
# this 1D walk (step ~ N(0, 2Dτ), mirror reflection), so the toy's :exact rule must agree with MCMR.

const _Layers = mr.Geometries.Internal.Layers

"""
    toy_layer_signal(; h, rate, tau, T=50.0, nspins, seed, method)

Per-spin final M⊥ after time T. `method = :exact`: on every reflected straight piece, add
rate·dt_piece·(time fraction within h of each face) using `segment_fraction`; `:endpoint`: add
rate·τ·(number of layers containing the end point of the step). Spin i uses its own RNG Xoshiro(seed·10⁷ + i),
so results do not depend on the thread count.
"""
function toy_layer_signal(; h, rate, tau, T=50.0, nspins, seed, method, w=W_REF, D=D_REF)
    method in (:exact, :endpoint) || error("method must be :exact or :endpoint, got $method")
    nsteps = round(Int, T / tau)
    σ = sqrt(2 * D * tau)
    out = ones(nspins)
    iszero(h) && return out
    Threads.@threads for i in 1:nspins
        rng = Random.Xoshiro(seed * 10_000_000 + i)
        x = rand(rng) * w
        expo = 0.0
        pts = Float64[]
        for _ in 1:nsteps
            step = randn(rng) * σ
            empty!(pts); push!(pts, x)
            while true
                if x + step < 0
                    step = -(step + x); x = 0.0; push!(pts, x)
                elseif x + step > w
                    step = -(step - (w - x)); x = w; push!(pts, x)
                else
                    x += step; push!(pts, x); break
                end
            end
            if method == :exact
                total = sum(abs(pts[k + 1] - pts[k]) for k in 1:length(pts) - 1)
                total > 0 || continue
                for k in 1:length(pts) - 1
                    dt = tau * abs(pts[k + 1] - pts[k]) / total
                    expo += rate * dt * (_Layers.segment_fraction(pts[k], pts[k + 1], 0.0, h) +
                                         _Layers.segment_fraction(w - pts[k], w - pts[k + 1], 0.0, h))
                end
            else
                expo += rate * tau * ((x <= h) + (w - x <= h))
            end
        end
        out[i] = exp(-expo)
    end
    return out
end
```

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS (testset "toy: exact vs endpoint layer rule").

- [ ] **Step 5: Create the script `research/phase2/p2_3_7_endpoint.jl`**

```julia
# P2.3.7 V11 — exact segment overlap vs endpoint rule, measured (no assumed sign of the bias).
#
# 1D toy (toy.jl), h = 0.1 µm, ΔR2(0) = 0.1 /ms, T = 50 ms, τ ∈ TAUS, seeds 1–10.
# (1) bias = S_endpoint − S_exact at T, in combined SEM, per τ (characterisation).
# (2) cross-check: toy :exact at τ = 1e-2 vs MCMR sweep (P2.3.4(b)) at h = 0.1, S(50) within 2 combined SEM.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_7_endpoint.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "fit.jl"))
include(joinpath(@__DIR__, "toy.jl"))

const OUT = phase2_outdir("p2_3_7_endpoint")
const H = 0.1
const TAUS = [1e-3, 1e-2, 4e-2]
const NSPINS = Dict(1e-3 => parse(Int, get(ENV, "NSPINS_1E3", "20000")), 1e-2 => 100_000, 4e-2 => 100_000)

seedmeans(τ, m) = [mean(toy_layer_signal(; h=H, rate=RATE_REF, tau=τ, nspins=NSPINS[τ], seed=s, method=m)) for s in 1:10]
rows = NamedTuple[]
for τ in TAUS
    ex = seedmeans(τ, :exact); ep = seedmeans(τ, :endpoint)
    se(v) = std(v) / sqrt(length(v))
    z = (mean(ep) - mean(ex)) / sqrt(se(ex)^2 + se(ep)^2)
    push!(rows, (tau_ms=τ, nspins=NSPINS[τ], S_exact=mean(ex), S_exact_sem=se(ex), S_endpoint=mean(ep), S_endpoint_sem=se(ep), bias=mean(ep) - mean(ex), bias_in_sem=z))
    @printf("τ = %.0e: S_exact %.5f ± %.5f, S_endpoint %.5f ± %.5f, bias %+.2e (%+.1f SEM)\n", τ, mean(ex), se(ex), mean(ep), se(ep), mean(ep) - mean(ex), z)
end
(hgrid, S_by_h) = read_sweep(joinpath(REPO_ROOT, "research", "results", "phase2", "p2_3_4_h_sweep", "sweep.csv"))
k = findfirst(==(H), hgrid)
mc = S_by_h[k][:, end]
toy = rows[findfirst(r -> r.tau_ms == 1e-2, rows)]
zc = (toy.S_exact - mean(mc)) / sqrt(toy.S_exact_sem^2 + (std(mc) / sqrt(length(mc)))^2)
cross_ok = abs(zc) <= 2
@printf("cross-check toy :exact vs MCMR at h = %.1f, τ = 1e-2: %.5f vs %.5f (%+.1f SEM) → %s\n", H, toy.S_exact, mean(mc), zc, cross_ok ? "PASS" : "FAIL")
write_csv(joinpath(OUT, "endpoint_vs_exact.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "h_um" => H, "rate" => RATE_REF, "T_ms" => 50.0,
    "rows" => rows, "cross_check_z" => zc, "cross_check_pass" => cross_ok,
    "note" => "toy and MCMR use different random streams; the cross-check is statistical"))
println("saved to $OUT")
```

- [ ] **Step 6: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_3_7_endpoint.jl`
Expected: three `τ = …` lines with a measured bias, and a cross-check line `→ PASS`. A cross-check FAIL means the toy and MCMR disagree on the same physics: debug before recording. About 10–20 min.

- [ ] **Step 7: Record**

```markdown
| P2.3.7 | V11: endpoint rule vs exact overlap (toy, h = 0.1) | bias quantified (sign not assumed); toy exact = MCMR within 2 SEM | bias at τ = 1e-3 / 1e-2 / 4e-2: <…> / <…> / <…> SEM; cross-check <PASS/FAIL> (<z> SEM) |
```

```markdown
## P2.3.7: V11, endpoint rule vs exact overlap

Script [p2_3_7_endpoint.jl](phase2/p2_3_7_endpoint.jl) → [results/phase2/p2_3_7_endpoint/](results/phase2/p2_3_7_endpoint/). 1D toy walk identical to MCMR's x-motion between walls (toy.jl); h = 0.1 µm, ΔR₂(0) = 0.1 ms⁻¹, T = 50 ms, seeds 1–10.

| τ (ms) | S(50) exact ± SEM | S(50) endpoint ± SEM | bias (endpoint − exact) | bias in SEM |
|---|---|---|---|---|
| <one row per τ from endpoint_vs_exact.csv> |

Cross-check: toy exact vs MCMR (P2.3.4(b), h = 0.1, τ = 1e-2): <…> vs <…> (<z> SEM) → <PASS/FAIL>.

**Interpretation.** <one paragraph: size and sign of the endpoint bias and how it changes with τ; what that means for using endpoint sampling instead of the exact overlap. Numbers only.>
```

- [ ] **Step 8: Commit**

```bash
git add research/phase2/toy.jl research/phase2/test/test_toy.jl research/phase2/p2_3_7_endpoint.jl research/results/phase2/p2_3_7_endpoint research/phase2_record.md
git commit -m "Add P2.3.7 V11 endpoint vs exact comparison

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Side-specific layer (P2.3.12)

**Files:**
- Create: `research/phase2/p2_3_12_sides.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `run_walls`, `layer_walls(...; side=...)`, `RATE_REF` (Plan A); `sweep_nspins` (Task 3).
- Produces: `results/phase2/p2_3_12_sides/summary.json` with the mirror check (`max_abs_z`, `pass`) and S(t) ± SEM for `:both`, `:positive`, `:negative` and `:asymmetric`.

- [ ] **Step 1: Create the script**

```julia
# P2.3.12 — side-specific layer.
# (a) Mirror check: positive-only vs negative-only (h = 0.2 µm, ΔR2(0) = 0.1), same seeds: |z| ≤ 2 at every readout,
#     z = ΔS / √(SEM₊² + SEM₋²) (conservative: same initial positions).
# (b) Characterisation: both sides, one side, and asymmetric (positive h = 0.2, ΔR2(0) = 0.1; negative h = 0.4,
#     ΔR2(0) = 0.05 → same ρ = 0.02 µm/ms on both sides) — S(t) with SEM. No baseline exists for asymmetric surfaces.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_12_sides.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "fit.jl"))

const OUT = phase2_outdir("p2_3_12_sides")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(sweep_nspins())))

cases = Dict(
    "both" => layer_walls(h=0.2, rate=RATE_REF),
    "positive" => layer_walls(h=0.2, rate=RATE_REF, side=:positive),
    "negative" => layer_walls(h=0.2, rate=RATE_REF, side=:negative),
    "asymmetric" => layer_walls(h=0.2, rate=RATE_REF, side=:asymmetric, h_neg=0.4, rate_neg=0.05),
)
runs = Dict(k => run_walls(v; nspins=NSPINS) for (k, v) in cases)
p, n = runs["positive"], runs["negative"]
σ = sqrt.(p.sem .^ 2 .+ n.sem .^ 2)
z = [σ[i] > 0 ? (p.mean[i] - n.mean[i]) / σ[i] : 0.0 for i in eachindex(σ)]
mirror_ok = all(abs.(z) .<= 2)
@printf("mirror: max |z| = %.2f → %s\n", maximum(abs.(z)), mirror_ok ? "PASS" : "FAIL")
for k in ("both", "positive", "negative", "asymmetric")
    @printf("%-10s S(50) = %.5f ± %.5f\n", k, runs[k].mean[end], runs[k].sem[end])
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "times" => TIMES_REF,
    "mirror_z" => z, "mirror_max_abs_z" => maximum(abs.(z)), "mirror_pass" => mirror_ok,
    "cases" => Dict(k => Dict("mean" => r.mean, "sem" => r.sem) for (k, r) in runs)))
println("saved to $OUT")
```

- [ ] **Step 2: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_3_12_sides.jl`
Expected: `mirror: … → PASS` and four S(50) lines. "positive" and "negative" should match each other, and "both" should attenuate more than either. About 10 min.

- [ ] **Step 3: Record**

```markdown
| P2.3.12 | Side-specific layer: mirror check + characterisation | positive-only = negative-only within 2 SEM at every readout | <PASS/FAIL> (max \|z\| <…>); S(50): both <…>, one side <…>, asymmetric <…> |
```

```markdown
## P2.3.12: side-specific layer

Script [p2_3_12_sides.jl](phase2/p2_3_12_sides.jl) → [results/phase2/p2_3_12_sides/](results/phase2/p2_3_12_sides/). h = 0.2 µm, ΔR₂(0) = 0.1 ms⁻¹ unless stated; asymmetric: positive h = 0.2 / ΔR₂(0) = 0.1, negative h = 0.4 / ΔR₂(0) = 0.05 (same ρ on both sides).

| Case | S(25) ± SEM | S(50) ± SEM |
|---|---|---|
| both sides | <…> | <…> |
| positive only | <…> | <…> |
| negative only | <…> | <…> |
| asymmetric | <…> | <…> |

Mirror check (positive vs negative): max |z| = <…> → <PASS/FAIL> (conservative: shared seeds). No baseline exists for asymmetric surfaces, because MCMR's θ_relax is the same on both sides.
```

- [ ] **Step 4: Commit**

```bash
git add research/phase2/p2_3_12_sides.jl research/results/phase2/p2_3_12_sides research/phase2_record.md
git commit -m "Add P2.3.12 side-specific layer check

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Gate, tracker and tag `rd-step-v1` (P2.3.10)

**Files:**
- Modify: `research/project-progress-tracker.md`, `research/phase2_record.md`

**Interfaces:**
- Consumes: `summary.json` of `p2_2_9_layer_trace` (`pass`), `p2_3_2_density_occupancy` (`pass`), `p2_3_exact_checks` (`results.half_gap_bulk_off.pass`, `results.half_gap_bulk_on.pass`).
- Produces: tracker ticks, the gate statement in the record, and git tag `rd-step-v1`.

- [ ] **Step 1: Check the gate from the saved results**

Run:

```bash
julia --project=research/baseline -e '
import JSON
r(p) = JSON.parsefile(joinpath("research/results/phase2", p, "summary.json"))
trace = r("p2_2_9_layer_trace")["pass"]
occ = r("p2_3_2_density_occupancy")["pass"]
ex = r("p2_3_exact_checks")["results"]
half = ex["half_gap_bulk_off"]["pass"] && ex["half_gap_bulk_on"]["pass"]
println("P2.2.9 trace: ", trace, "  P2.3.3 occupancy: ", occ, "  P2.3.4(a) h = w/2: ", half)
println(trace && occ && half ? "GATE PASS" : "GATE FAIL")'
```

Expected: `GATE PASS`. If it prints `GATE FAIL`, stop. Do not tick 2.3 or tag; report which item failed.

- [ ] **Step 2: Record the gate and tick the tracker**

Append to `research/phase2_record.md`:

```markdown
## Gate for section 2.4 (2026-10-03 or the run date)

P2.2.9 single-spin trace: PASS. P2.3.3 occupancy: PASS. P2.3.4(a) exact end point h = w/2: PASS. → Section 2.4 (other profiles) may start. All characterisation results above are at τ = 1e-2 ms and are provisional until the timestep study (P2.6).
```

In `research/project-progress-tracker.md`, set P2.3.1, P2.3.2, P2.3.3, P2.3.4, P2.3.5, P2.3.6, P2.3.7, P2.3.9, P2.3.11, P2.3.12 and P2.3.10 to ☑, each with its output folder. P2.3.8 stays ✖ (deferred). Increase the Phase 2 "Done" count by 11.

- [ ] **Step 3: Commit and tag**

```bash
git add research/project-progress-tracker.md research/phase2_record.md
git commit -m "Record step-profile verification and gate for section 2.4

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git tag -a rd-step-v1 -m "Step-profile R2(d) on walls verified (P2.2, P2.3); results at tau = 1e-2 ms provisional until P2.6"
git tag -l rd-step-v1
```

Expected: the last command prints `rd-step-v1`. Do not push the tag.
