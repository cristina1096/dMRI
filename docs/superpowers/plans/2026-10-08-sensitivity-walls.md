# Sensitivity to the Spatial Distribution of Near-Surface Relaxation (Walls) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Map when the free-decay signal between walls is sensitive to *how* near-surface relaxation is distributed, not only to how strong it is in total. The near-surface layer (step, linear, exponential) is compared with the conventional surface-relaxation model at the **same integrated relaxivity ρ**, while varying h, the compartment size w, D, ρ and the observation time. The signal difference is judged against the Monte Carlo SEM and against MRI noise at SNR 100 and 200.

**Architecture:**
- **Fixed ρ.** MCMR's conventional model removes a fraction θ·√τ of M⊥ at every wall collision. Gaussian steps hit a wall at a rate of n·√(D/(πτ)) per unit area. The surface relaxivity is therefore **ρ = θ·√(D/π)**, so the baseline is set with θ = ρ·√(π/D).
- **Checking that.** The baseline obeys the same mean-exponent identity as the layer, E[−ln M] = 2ρt/w. Task 2 confirms it at D = 1, 3 and 6 before any comparison is made.
- **Paired runs.** Neither the layer nor the conventional surface relaxation uses random numbers, so with the same seed both models follow identical trajectories. The per-seed difference ΔS(t) = S_layer − S_base is therefore measured with a very small paired SEM.
- **Sweeps.** One runner script does four one-at-a-time sweeps around a reference point (h, w, D, ρ). A summary script turns them into tables, figures and an observability statement. The timestep in each configuration is min(1e-2, 0.005/ΔR₂(0)), the adopted layer constraint, and both models use it.

**Tech Stack:** Julia 1.12, MCMRSimulator (this repo), `research/phase2/`, CairoMakie, `Test` stdlib.

**Spec:**
- The user's key investigation (2026-10-07): "Investigating when the simulated MRI signal becomes sensitive to the spatial distribution of near-surface relaxation, rather than only to its total strength … constant integrated near-surface relaxivity … varying h, compartment size, diffusion time, diffusion coefficient, and the shape and magnitude of the R₂(d) profile … relative to Monte Carlo uncertainty and representative MRI noise level."
- The user's decisions (2026-10-07):
  - Fixed ρ is the primary comparison; best fit is secondary.
  - Noise levels are SNR 100 and 200.
  - Profiles are step, linear and exponential, with step excluded for h > w/2.
  - Walls only, free decay.
- `research/project-progress-tracker.md` §2.5: P2.5.1–P2.5.5 and the scope decision.
- The user's framing: the baseline is a reference, not the truth; no fast-diffusion arguments.

## Global Constraints

- Work only on branch `cc/near-surface-r2`.
- **Reference point:** walls w = 2 µm, D = 3 µm²/ms, ρ = 0.01 µm/ms, h = 0.2 µm, layer on both faces, R₂_bulk = 0, no RF, readouts 0:10:200 ms. Seeds 1–10; spins per seed set by the pilot in Task 4 (default 50 000).
- **Baseline:** `Walls(repeats=w, surface_relaxation=θ)` with θ = ρ·√(π/D), and the same τ, seeds, spins and readouts as the layer run it is paired with.
- **Timestep:** τ = min(1e-2, 0.005/ΔR₂(0)), where ΔR₂(0) = (ρ/h)·g(0) for the profile. This is the adopted layer constraint (P2.6.5) with the 1e-2 ms reference cap.
- **Excluded:** the step profile with h > w/2 (decision 2026-10-07).
- **Observability:**
  - ΔS(t) is in units of S(0) = 1.
  - "MC-resolved" means |ΔS| ≥ 3 × paired SEM.
  - "Potentially observable at SNR X" means |ΔS| ≥ 2/X, i.e. 2σ of a single measurement with relative noise 1/X: 0.02 at SNR 100, 0.01 at SNR 200.
  - These are reported, not pass/fail. Only Task 2 has pass/fail.
- **Secondary (best fit):** the residual max |S_layer − e^(−R̂t)| after a mono-exponential fit. This is the best the conventional model could do, provided the baseline itself is mono-exponential within noise; that is checked per configuration through its curvature.
- **Long runs:** start detached (`nohup … & disown`) and wait with a background `until` loop, re-armed at the 2-hour limit.
- **Final review:** done by a fresh subagent reviewer (user's standing rule).
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `research/phase2/harness.jl` | Modify | `run_walls(…; D)`, `layer_walls(…; w)`, `rho_for(…; w)`, `profile_g0(shape, h; w)`, `theta_for`, `baseline_walls`, `rho_walls`, `tau_for` |
| `research/phase2/sensitivity.jl` | Create | `paired_difference`, `peak`, `monoexp_residual`, `centroid` |
| `research/phase2/test/test_harness.jl`, `test_sensitivity.jl`, `runtests.jl` | Modify / Create | research tests |
| `research/phase2/p2_5_1_theta_mapping.jl` | Create | mean-exponent identity for the baseline at D = 1, 3, 6 (and paired with the layer) |
| `research/phase2/p2_5_sensitivity.jl` | Create | one sweep (env `SWEEP` = h, w, D, rho), paired runs, per-configuration metrics |
| `research/phase2/p2_5_summary.jl` | Create | tables, figures, centroid replot, observability summary |
| `research/phase2_record.md`, `research/project-progress-tracker.md` | Modify | records, ticks |

## Review Focus

1. **Pairing really is paired.** The layer and baseline runs of a configuration must use the same seeds, spins, τ and readouts. Otherwise the paired SEM is meaningless. Test: Task 4's pilot checks that the per-seed positions are identical (a snapshot comparison).
2. **The θ mapping depends on D.** θ must be recomputed for each D. A fixed θ across the D sweep would compare different ρ. Test: Task 1 (`theta_for` at D = 1 and 6) and Task 2 (the identity at D = 1, 3, 6).
3. **Step exclusion.** No step configuration with h > w/2 may be run or reported. Test: Task 3 (the `configs` unit test).
4. **Timestep per configuration.** A strong layer (large ΔR₂(0)) must use the smaller τ, and its paired baseline must use the same τ. Test: Task 3 (`tau_for`), and the runner writes τ per configuration row.
5. **Exponential cutoff equals the w being swept.** In the w sweep the exponential's cutoff must follow w, which it does by default since cutoff = spacing. `rho_for`/`profile_g0` must use that w, not W_REF. Test: Task 1 (`profile_g0("exponential", h; w)`).

---

### Task 1: Harness generalisation and analysis helpers

**Files:**
- Modify: `research/phase2/harness.jl`, `research/phase2/test/test_harness.jl`, `research/phase2/test/runtests.jl`
- Create: `research/phase2/sensitivity.jl`, `research/phase2/test/test_sensitivity.jl`

**Interfaces:**
- Produces:
  - `run_walls(geometry; …, D=D_REF)`
  - `profile_g0(shape, h; w=W_REF)`
  - `rho_for(rate, h; shape="step", w=W_REF)`
  - `layer_walls(…; w=W_REF)`
  - `theta_for(ρ; D=D_REF)` = ρ·√(π/D)
  - `baseline_walls(; theta, w=W_REF)`
  - `rho_walls(; rho, h, shape, w=W_REF)` = `Walls(repeats=w, layer_rho=rho, layer_h=h, layer_shape=shape)`
  - `tau_for(rho, h, shape; w=W_REF)` = min(TAU_REF, DEFAULT_LAYER_SCALING/((ρ/h)·g(0)))
  - `paired_difference(Sa, Sb)` → (mean, sem) per readout
  - `peak(times, m)` → (index, t, value) of the largest |m|
  - `monoexp_residual(times, S)` → (R, residual)
  - `centroid(shape, h; cutoff)`

- [ ] **Step 1: Write the failing tests**

Append to `research/phase2/test/test_harness.jl`:

```julia
@testset "harness: sensitivity extensions" begin
    @test theta_for(0.01; D=3.0) ≈ 0.01 * sqrt(pi / 3)
    @test theta_for(0.01; D=1.0) ≈ 0.01 * sqrt(pi)
    @test profile_g0("exponential", 0.5; w=5.0) ≈ 1 / (1 - exp(-10))
    @test profile_g0("exponential", 0.5) ≈ 1 / (1 - exp(-4))                  # default w = W_REF
    @test rho_for(0.1, 0.5; shape="exponential", w=5.0) ≈ 0.05 * (1 - exp(-10))
    l = Simulation([]; geometry=layer_walls(h=0.5, rate=0.1, shape="exponential", w=5.0), verbose=false).geometry[1].layer[1]
    @test l.positive.cutoff == 5.0
    @test mr.Geometries.Internal.Layers.surface_rate(l.positive) ≈ 0.1
    b = Simulation([]; geometry=baseline_walls(theta=0.02, w=5.0), verbose=false)
    @test b.geometry[1].layer === nothing
    r = Simulation([]; geometry=rho_walls(rho=0.01, h=0.2, shape="linear", w=5.0), verbose=false).geometry[1].layer[1]
    @test (r.positive.rho, r.positive.h, r.positive.shape) == (0.01, 0.2, :linear)
    @test tau_for(0.01, 0.2, "step") == 1e-2                                    # ΔR2(0) = 0.05: constraint 0.1 > cap
    @test tau_for(0.1, 0.2, "linear") ≈ 0.005 / 1.0                             # ΔR2(0) = 2ρ/h = 1.0
    @test tau_for(0.1, 0.05, "step") ≈ 0.005 / 2.0
    # D keyword reaches the simulation
    @test run_walls(Walls(repeats=2.0); D=1.0, times=[0.0, 1.0], seeds=1:1, nspins=10).sim.diffusivity == 1.0
end
```

Create `research/phase2/test/test_sensitivity.jl`:

```julia
@testset "sensitivity tools" begin
    Sa = [1.0 0.9 0.8; 1.0 0.92 0.82]
    Sb = [1.0 0.95 0.9; 1.0 0.95 0.9]
    d = paired_difference(Sa, Sb)
    @test d.mean ≈ [0.0, -0.04, -0.09]
    @test d.sem ≈ [0.0, std([-0.05, -0.03]) / sqrt(2), std([-0.1, -0.08]) / sqrt(2)]
    p = peak([0.0, 10.0, 20.0], d.mean)
    @test (p.index, p.t, p.value) == (3, 20.0, d.mean[3])
    t = collect(0.0:10.0:100.0)
    f = monoexp_residual(t, exp.(-0.02 .* t))
    @test f.R ≈ 0.02 && f.residual < 1e-12
    g = monoexp_residual(t, 0.5 .* exp.(-0.01 .* t) .+ 0.5 .* exp.(-0.05 .* t))   # biexponential: residual > 0
    @test g.residual > 1e-3
    @test centroid("step", 0.2) ≈ 0.1
    @test centroid("linear", 0.3) ≈ 0.1
    @test centroid("exponential", 0.5; cutoff=2.0) ≈ 0.5 * (1 - 5 * exp(-4)) / (1 - exp(-4))
    @test centroid("exponential", 0.01; cutoff=2.0) ≈ 0.01 rtol=1e-10                # untruncated limit: h
end
```

In `research/phase2/test/runtests.jl`, after the `timestep.jl` include line, add:

```julia
isfile(joinpath(@__DIR__, "..", "sensitivity.jl")) && include(joinpath(@__DIR__, "..", "sensitivity.jl"))
```

- [ ] **Step 2: Run the research tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `UndefVarError: theta_for` and `UndefVarError: paired_difference`.

- [ ] **Step 3: Implement the harness changes**

In `research/phase2/harness.jl`:

- In `run_walls`, add the keyword `D=D_REF` and use `diffusivity=D` in the `Simulation` call.
- Replace `profile_g0`, `rho_for` and `layer_walls` with:

```julia
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
```

For the existing calls, `w=W_REF` reproduces the old geometry and ρ exactly.

- [ ] **Step 4: Create `research/phase2/sensitivity.jl`**

```julia
# Sensitivity-study helpers (P2.5): paired differences, peaks, mono-exponential residual, profile centroids.

"Per-readout mean and SEM over seeds of the paired difference Sa − Sb (rows = seeds, same seeds in both)."
function paired_difference(Sa::AbstractMatrix, Sb::AbstractMatrix)
    D = Sa .- Sb
    n = size(D, 1)
    return (mean=vec(mean(D, dims=1)), sem=n > 1 ? vec(std(D, dims=1)) ./ sqrt(n) : zeros(size(D, 2)), per_seed=D)
end

"Index, time and value of the largest |m|."
function peak(times, m)
    k = argmax(abs.(m))
    return (index=k, t=times[k], value=m[k])
end

"""
    monoexp_residual(times, S)

Best mono-exponential S ≈ e^(−R t) by least squares on ln S through the origin (t > 0), and the largest
|S − e^(−R t)| over the readouts. Stands for the best the conventional model can do when it is mono-exponential.
"""
function monoexp_residual(times, S)
    k = findall(t -> t > 0, times)
    R = -sum(times[k] .* log.(S[k])) / sum(times[k] .^ 2)
    return (R=R, residual=maximum(abs.(S .- exp.(-R .* times))))
end

"""
    centroid(shape, h; cutoff=Inf)

Mean distance from the wall of the relaxation, ∫ d·g(d/h) dd / ∫ g(d/h) dd: step h/2, linear h/3,
exponential truncated at `cutoff`: h·(1 − (c + 1)e^(−c))/(1 − e^(−c)), c = cutoff/h.
"""
function centroid(shape, h; cutoff=Inf)
    shape == "step" && return h / 2
    shape == "linear" && return h / 3
    shape == "exponential" || error("Unknown shape $shape")
    c = cutoff / h
    return isinf(c) ? h : h * (1 - (c + 1) * exp(-c)) / -expm1(-c)
end
```

- [ ] **Step 5: Run the research tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS. The earlier research tests are unchanged because `w=W_REF` and `D=D_REF` are the defaults.

- [ ] **Step 6: Commit**

```bash
git add research/phase2/harness.jl research/phase2/sensitivity.jl research/phase2/test/test_harness.jl research/phase2/test/test_sensitivity.jl research/phase2/test/runtests.jl
git commit -m "Harness: w and D parameters, theta mapping, baseline/rho walls; sensitivity helpers

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Confirm the ρ ↔ θ mapping (P2.5 prerequisite)

**Files:**
- Create: `research/phase2/p2_5_1_theta_mapping.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `theta_for`, `baseline_walls`, `rho_walls`, `run_walls(…; D, snapshots=true)` (Task 1).
- Produces: `results/phase2/p2_5_1_theta_mapping/summary.json` with `pass`.

- [ ] **Step 1: Create the script**

```julia
# P2.5 prerequisite — the conventional model's relaxivity is ρ = θ·√(D/π).
#
# MCMR multiplies M⊥ by e^(−θ√τ) at each wall collision; Gaussian steps hit a wall at n·√(D/(πτ)) per area and time,
# so E[−ln M_base(t)] = 2·θ·√(D/π)·t/w for a uniform start (the mean-exponent identity of the baseline). With
# θ = ρ·√(π/D) this equals the layer's 2ρt/w (P2.4.5). Checked at D = 1, 3, 6 and τ = 1e-2, 4e-2 (w = 2, ρ = 0.01,
# t = 50 ms, seeds 1–10). Pass: |z| ≤ 3 for the baseline at every (D, τ). Also reported: the paired difference of the
# mean exponents, baseline − exponential layer (h = 0.2) on identical trajectories.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_5_1_theta_mapping.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_5_1_theta_mapping")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const RHO, T = 0.01, 50.0
mean_expo(r) = [mean(-log(s.orientations[1].transverse) for s in snap.spins) for snap in r.snapshots]
positions(r) = [mr.position.(snap) for snap in r.snapshots]

rows = NamedTuple[]
for D in (1.0, 3.0, 6.0), τ in (1e-2, 4e-2)
    θ = theta_for(RHO; D=D)
    b = run_walls(baseline_walls(theta=θ); D=D, timestep=τ, times=[0.0, T], nspins=NSPINS, snapshots=true)
    l = run_walls(rho_walls(rho=RHO, h=0.2, shape="exponential"); D=D, timestep=τ, times=[0.0, T], nspins=NSPINS, snapshots=true)
    eb, el = mean_expo(b), mean_expo(l)
    expected = 2RHO * T / W_REF
    zb = (mean(eb) - expected) / (std(eb) / sqrt(length(eb)))
    dpair = eb .- el
    zp = mean(dpair) / (std(dpair) / sqrt(length(dpair)))
    same_paths = positions(b) == positions(l)
    push!(rows, (D=D, tau_ms=τ, theta=θ, expected=expected, base_mean=mean(eb), base_sem=std(eb) / sqrt(length(eb)), z_base=zb,
        layer_mean=mean(el), paired_diff=mean(dpair), paired_sem=std(dpair) / sqrt(length(dpair)), z_paired=zp,
        same_paths=same_paths, pass=abs(zb) <= 3))
    @printf("D = %.0f, τ = %.0e: θ = %.5f; baseline E[−ln M] = %.6f ± %.1e (2ρt/w = %.6f, z = %+.2f); paired base − layer = %+.1e ± %.0e (z = %+.1f); identical paths: %s → %s\n",
        D, τ, θ, mean(eb), std(eb) / sqrt(length(eb)), expected, zb, mean(dpair), std(dpair) / sqrt(length(dpair)), zp, same_paths, abs(zb) <= 3 ? "PASS" : "FAIL")
end
write_csv(joinpath(OUT, "mapping.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rho" => RHO, "T_ms" => T, "nspins" => NSPINS,
    "rows" => rows, "pass" => all(r.pass for r in rows), "paths_identical" => all(r.same_paths for r in rows)))
println("saved to $OUT")
```

- [ ] **Step 2: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_5_1_theta_mapping.jl`
Expected: 6 lines, all `→ PASS`, with `identical paths: true`. About 25 min.
- **If any line FAILs:** the mapping is off, for example because of multiple collisions per step at large τ. Record the measured ratio E[−ln M_base]/(2θt/w) as the empirical √(D/π) replacement, and use that in `theta_for` (a ruling). That keeps "same mean exponent" as the definition of fixed ρ.
- **If the paths are not identical:** the pairing assumption is wrong. Stop and investigate.

- [ ] **Step 3: Record and commit**

Add a summary row and a section "P2.5.0: conventional model at the same ρ — θ = ρ√(π/D)" to `research/phase2_record.md`. It should contain the derivation, the table (D, τ, θ, 2ρt/w, baseline mean ± SEM, z, paired difference ± SEM, identical paths, outcome), and one sentence stating that from here on, fixed ρ means θ = ρ√(π/D) for the baseline.

```bash
git add research/phase2/p2_5_1_theta_mapping.jl research/results/phase2/p2_5_1_theta_mapping research/phase2_record.md
git commit -m "Confirm the conventional model's relaxivity rho = theta*sqrt(D/pi)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: The sweep runner

**Files:**
- Create: `research/phase2/p2_5_sensitivity.jl`
- Modify: `research/phase2/sensitivity.jl`, `research/phase2/test/test_sensitivity.jl`

**Interfaces:**
- Consumes: Task 1 helpers.
- Produces:
  - `sweep_configs(sweep)::Vector{NamedTuple}`, with fields `(shape, h, w, D, rho)`, in `sensitivity.jl`
  - `results/phase2/p2_5_sensitivity_<sweep>/{rows.csv, curves.json, summary.json}`

- [ ] **Step 1: Write the failing test**

Append to `research/phase2/test/test_sensitivity.jl`:

```julia
@testset "sweep configurations" begin
    for s in ("h", "w", "D", "rho")
        c = sweep_configs(s)
        @test !isempty(c)
        @test all(x -> !(x.shape == "step" && x.h > x.w / 2), c)               # step excluded for h > w/2
        @test all(x -> x.shape in ("step", "linear", "exponential"), c)
    end
    h = sweep_configs("h")
    @test sort(unique([x.h for x in h])) == [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
    @test count(x -> x.shape == "step", h) == 5                                 # h = 2.0 excluded for step
    @test sort(unique([x.w for x in sweep_configs("w")])) == [1.0, 2.0, 5.0, 10.0]
    @test sort(unique([x.D for x in sweep_configs("D")])) == [1.0, 3.0, 6.0]
    @test sort(unique([x.rho for x in sweep_configs("rho")])) == [0.003, 0.01, 0.03, 0.1]
    @test_throws ErrorException sweep_configs("T")
end
```

- [ ] **Step 2: Run the research tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `UndefVarError: sweep_configs`.

- [ ] **Step 3: Add `sweep_configs` to `research/phase2/sensitivity.jl`**

```julia
"Reference point of the sensitivity study (Global Constraints of the plan)."
const SENS_REF = (h=0.2, w=2.0, D=3.0, rho=0.01)
const SENS_SHAPES = ("step", "linear", "exponential")

"""
    sweep_configs(sweep)

One-at-a-time configurations around SENS_REF for sweep "h", "w", "D" or "rho": every shape at every value,
except the step profile where h > w/2 (overlapping step layers relax the middle of the gap most; excluded).
"""
function sweep_configs(sweep::AbstractString)
    values = sweep == "h" ? [0.05, 0.1, 0.2, 0.5, 1.0, 2.0] :
             sweep == "w" ? [1.0, 2.0, 5.0, 10.0] :
             sweep == "D" ? [1.0, 3.0, 6.0] :
             sweep == "rho" ? [0.003, 0.01, 0.03, 0.1] : error("Unknown sweep $sweep; use h, w, D or rho.")
    key = Symbol(sweep)
    out = NamedTuple[]
    for v in values, shape in SENS_SHAPES
        c = merge(SENS_REF, NamedTuple{(key,)}((v,)))
        (shape == "step" && c.h > c.w / 2) && continue
        push!(out, (shape=shape, h=c.h, w=c.w, D=c.D, rho=c.rho))
    end
    return out
end
```

- [ ] **Step 4: Run the research tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS.

- [ ] **Step 5: Create the runner `research/phase2/p2_5_sensitivity.jl`**

```julia
# P2.5 — sensitivity to the spatial distribution of near-surface relaxation at fixed integrated relaxivity ρ.
#
# For every configuration of sweep SWEEP (h, w, D or rho; see sweep_configs): the layer (shape, h, ρ) and the
# conventional model with θ = ρ√(π/D) (P2.5.0), same walls, D, τ = tau_for(...), seeds 1–10, NSPINS spins, readouts
# 0:10:200 ms. Both models follow identical trajectories, so ΔS(t) = S_layer − S_base is a paired difference per seed.
# Per configuration: max |ΔS| and its time, paired SEM, |ΔS|/SEM, |ΔS| in units of 1/SNR for SNR 100 and 200,
# the best mono-exponential residual of the layer curve, and the baseline's curvature (mono-exponential check).
# Characterisation: nothing here is pass/fail.
#
# Run: SWEEP=h julia --project=research/baseline -t 8 research/phase2/p2_5_sensitivity.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "fit.jl"))
include(joinpath(@__DIR__, "sensitivity.jl"))

const SWEEP = get(ENV, "SWEEP", "h")
const NSPINS = parse(Int, get(ENV, "NSPINS", "50000"))
const TIMES = collect(0.0:10.0:200.0)
const OUT = phase2_outdir("p2_5_sensitivity_" * SWEEP)

base_cache = Dict{Tuple{Float64, Float64, Float64, Float64}, Any}()
function baseline_run(w, D, θ, τ)
    get!(base_cache, (w, D, θ, τ)) do
        run_walls(baseline_walls(theta=θ, w=w); D=D, timestep=τ, times=TIMES, nspins=NSPINS)
    end
end

rows = NamedTuple[]
curves = Dict{String, Any}[]
for c in sweep_configs(SWEEP)
    τ = tau_for(c.rho, c.h, c.shape; w=c.w)
    θ = theta_for(c.rho; D=c.D)
    lay = run_walls(rho_walls(rho=c.rho, h=c.h, shape=c.shape, w=c.w); D=c.D, timestep=τ, times=TIMES, nspins=NSPINS)
    base = baseline_run(c.w, c.D, θ, τ)
    d = paired_difference(lay.S, base.S)
    p = peak(TIMES, d.mean)
    fit = monoexp_residual(TIMES, lay.mean)
    cb = curvature(TIMES, base.mean)
    row = (sweep=SWEEP, shape=c.shape, h_um=c.h, w_um=c.w, D=c.D, rho=c.rho, theta=θ, tau_ms=τ,
        centroid_um=centroid(c.shape, c.h; cutoff=c.w),
        dS_peak=p.value, dS_peak_t_ms=p.t, dS_peak_sem=d.sem[p.index], dS_peak_z=p.value / max(d.sem[p.index], eps()),
        dS_in_snr100_units=abs(p.value) * 100, dS_in_snr200_units=abs(p.value) * 200,
        S_base_end=base.mean[end], S_layer_end=lay.mean[end],
        monoexp_R=fit.R, monoexp_residual=fit.residual, base_curvature=cb, runtime_s=lay.runtime)
    push!(rows, row)
    push!(curves, Dict("config" => Dict(pairs(c)), "tau_ms" => τ, "theta" => θ, "times" => TIMES,
        "S_layer" => lay.mean, "S_layer_sem" => lay.sem, "S_base" => base.mean, "S_base_sem" => base.sem,
        "dS" => d.mean, "dS_sem" => d.sem, "S_layer_per_seed" => [lay.S[i, :] for i in 1:size(lay.S, 1)]))
    @printf("%-11s h=%.3f w=%.1f D=%.0f ρ=%.3f τ=%.0e: max|ΔS| = %+.2e ± %.0e at %3.0f ms (%.0f SEM; ×SNR100 %.2f, ×SNR200 %.2f); mono-exp residual %.1e\n",
        c.shape, c.h, c.w, c.D, c.rho, τ, p.value, d.sem[p.index], p.t, abs(row.dS_peak_z), row.dS_in_snr100_units, row.dS_in_snr200_units, fit.residual)
    flush(stdout)
end
write_csv(joinpath(OUT, "rows.csv"), rows)
write_json(joinpath(OUT, "curves.json"), curves)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "sweep" => SWEEP, "nspins" => NSPINS,
    "times" => TIMES, "rows" => rows,
    "note" => "ΔS = S_layer − S_base on identical trajectories (same seeds); observable at SNR X if |ΔS| ≥ 2/X"))
println("saved to $OUT")
```

- [ ] **Step 6: Commit**

```bash
git add research/phase2/sensitivity.jl research/phase2/test/test_sensitivity.jl research/phase2/p2_5_sensitivity.jl
git commit -m "Add the P2.5 sensitivity sweep runner (paired layer vs conventional model at fixed rho)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Pilot and the four sweeps (P2.5.1, P2.5.2)

**Files:**
- Modify: `research/phase2_record.md` (pilot ruling only; the results are recorded in Task 5)

**Interfaces:**
- Consumes: the runner (Task 3).
- Produces: `results/phase2/p2_5_sensitivity_{h,w,D,rho}/`.

- [ ] **Step 1: Pilot**

Run: `SWEEP=rho NSPINS=50000 julia --project=research/baseline -t 8 -e 'include("research/phase2/harness.jl"); include("research/phase2/sensitivity.jl"); c = (shape="exponential", h=0.2, w=2.0, D=3.0, rho=0.01); τ = tau_for(c.rho, c.h, c.shape); a = run_walls(rho_walls(rho=c.rho, h=c.h, shape=c.shape); timestep=τ, times=0:10:200, nspins=50000, seeds=1:3, snapshots=true); b = run_walls(baseline_walls(theta=theta_for(c.rho)); timestep=τ, times=0:10:200, nspins=50000, seeds=1:3, snapshots=true); println("identical paths: ", [mr.position.(s) for s in a.snapshots] == [mr.position.(s) for s in b.snapshots]); d = paired_difference(a.S, b.S); p = peak(collect(0:10:200.0), d.mean); println("peak ΔS = ", p.value, " ± ", d.sem[p.index], "; runtime per run ", a.runtime, " s")'`
Expected:
- `identical paths: true`;
- a paired SEM at the peak of ≤ 1e-4, which is ≤ 1/50 of the SNR-200 threshold 0.01;
- the runtime per run, used to estimate the total.

If the SEM is larger than 1e-4, raise `NSPINS` for all sweeps and record a ruling. If the paths are not identical, stop: the pairing assumption fails.

- [ ] **Step 2: Run the four sweeps (detached chain)**

```bash
cat > /tmp/sens.sh <<'EOF'
#!/bin/zsh
set -e
cd /Users/ishimiyabisaki/Documents/mcmrsimulator.jl-main-cc-near-surface-r2
for s in h w D rho; do
  SWEEP=$s NSPINS=50000 julia --project=research/baseline -t 8 research/phase2/p2_5_sensitivity.jl > /tmp/sens_$s.log 2>&1; echo "$s done"
done
echo CHAIN DONE
EOF
chmod +x /tmp/sens.sh; nohup /tmp/sens.sh > /tmp/sens.out 2>&1 & disown
```

Wait for `CHAIN DONE` with a background `until` loop, re-armed at the 2-hour limit. Expected:
- one line per configuration: 17 (h), 12 (w), 9 (D) and 12 (rho);
- every sweep ends with `saved to`.

That is about 50 layer runs plus the cached baselines. Estimate the time from the pilot runtime; it is roughly 4–6 h.

- [ ] **Step 3: Commit the raw results**

```bash
git add research/results/phase2/p2_5_sensitivity_h research/results/phase2/p2_5_sensitivity_w research/results/phase2/p2_5_sensitivity_D research/results/phase2/p2_5_sensitivity_rho
git commit -m "Run the four P2.5 sensitivity sweeps (h, w, D, rho)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Summary, centroid replot and record (P2.5.2–P2.5.5)

**Files:**
- Create: `research/phase2/p2_5_summary.jl`
- Modify: `research/phase2_record.md`, `research/project-progress-tracker.md`

**Interfaces:**
- Consumes: `results/phase2/p2_5_sensitivity_{h,w,D,rho}/{rows.csv, curves.json}`.
- Produces: `results/phase2/p2_5_summary/{summary.json, observability.csv, p2_5_sweeps.png, p2_5_centroid.png, p2_5_curves.png}`.

- [ ] **Step 1: Create the summary script**

```julia
# P2.5 summary — observability of the profile (spatial distribution) at fixed ρ, from the four sweeps.
#
# (1) Per sweep: |ΔS| at its peak vs the swept parameter, per shape, with 3·SEM and the MRI thresholds 2/SNR
#     (SNR 100: 0.02, SNR 200: 0.01). (2) h sweep replotted against the profile centroid (P2.5.4). (3) ΔS(t) curves
#     for the h sweep. (4) observability.csv: one row per configuration.
#
# Run: julia --project=research/baseline research/phase2/p2_5_summary.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_5_summary")
const SWEEPS = ("h", "w", "D", "rho")
R = joinpath(REPO_ROOT, "research", "results", "phase2")
load_rows(s) = [Dict(zip(split(readlines(joinpath(R, "p2_5_sensitivity_" * s, "rows.csv"))[1], ","), split(l, ",")))
                for l in readlines(joinpath(R, "p2_5_sensitivity_" * s, "rows.csv"))[2:end]]
num(d, k) = parse(Float64, d[k])
classify(a, sem) = a >= 0.02 ? "SNR 100" : a >= 0.01 ? "SNR 200" : a >= 3sem ? "MC only" : "not resolved"

obs = NamedTuple[]
for s in SWEEPS, r in load_rows(s)
    a, se = abs(num(r, "dS_peak")), num(r, "dS_peak_sem")
    push!(obs, (sweep=s, shape=r["shape"], h_um=num(r, "h_um"), w_um=num(r, "w_um"), D=num(r, "D"), rho=num(r, "rho"),
        tau_ms=num(r, "tau_ms"), centroid_um=num(r, "centroid_um"), dS_peak=num(r, "dS_peak"), dS_peak_sem=se,
        dS_peak_t_ms=num(r, "dS_peak_t_ms"), level=classify(a, se), monoexp_residual=num(r, "monoexp_residual"),
        base_curvature=num(r, "base_curvature")))
    o = obs[end]
    @printf("%-3s %-11s h=%.3f w=%4.1f D=%.0f ρ=%.3f: ΔS = %+.2e ± %.0e at %3.0f ms → %s\n",
        s, o.shape, o.h_um, o.w_um, o.D, o.rho, o.dS_peak, se, o.dS_peak_t_ms, o.level)
end
write_csv(joinpath(OUT, "observability.csv"), obs)
counts = Dict(l => count(o -> o.level == l, obs) for l in ("SNR 100", "SNR 200", "MC only", "not resolved"))
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => obs, "counts" => counts,
    "thresholds" => Dict("SNR 100" => 0.02, "SNR 200" => 0.01, "MC" => "3·paired SEM")))
println("counts: ", counts)

using CairoMakie
param = Dict("h" => :h_um, "w" => :w_um, "D" => :D, "rho" => :rho)
fig = Figure(size=(1200, 330))
for (j, s) in enumerate(SWEEPS)
    ax = Axis(fig[1, j], xscale=log10, yscale=log10, xlabel=string(param[s]), ylabel="|ΔS| at peak", title="sweep $s")
    for shape in ("step", "linear", "exponential")
        sel = [o for o in obs if o.sweep == s && o.shape == shape]
        isempty(sel) && continue
        x = [getproperty(o, param[s]) for o in sel]
        scatterlines!(ax, x, max.(abs.([o.dS_peak for o in sel]), 1e-7), label=shape)
    end
    hlines!(ax, [0.02, 0.01], color=:gray, linestyle=:dash)
    j == 1 && axislegend(ax, position=:lt, labelsize=8)
end
save(joinpath(OUT, "p2_5_sweeps.png"), fig)

fig2 = Figure(size=(560, 380))
ax = Axis(fig2[1, 1], xscale=log10, yscale=log10, xlabel="centroid of the relaxation (µm)", ylabel="|ΔS| at peak", title="h sweep vs centroid (P2.5.4)")
for shape in ("step", "linear", "exponential")
    sel = [o for o in obs if o.sweep == "h" && o.shape == shape]
    scatterlines!(ax, [o.centroid_um for o in sel], max.(abs.([o.dS_peak for o in sel]), 1e-7), label=shape)
end
hlines!(ax, [0.02, 0.01], color=:gray, linestyle=:dash); axislegend(ax, position=:lt)
save(joinpath(OUT, "p2_5_centroid.png"), fig2)

curves = JSON.parsefile(joinpath(R, "p2_5_sensitivity_h", "curves.json"))
fig3 = Figure(size=(560, 380))
ax = Axis(fig3[1, 1], xlabel="t (ms)", ylabel="ΔS = S_layer − S_base", title="h sweep: ΔS(t) at fixed ρ = 0.01 µm/ms")
for c in curves
    cfg = c["config"]
    lines!(ax, Float64.(c["times"]), Float64.(c["dS"]), label="$(cfg["shape"]) h=$(cfg["h"])",
        linestyle=cfg["shape"] == "step" ? :dash : cfg["shape"] == "linear" ? :solid : :dot)
end
hlines!(ax, [-0.02, -0.01, 0.01, 0.02], color=:gray, linestyle=:dash)
Legend(fig3[1, 2], ax, labelsize=7, rowgap=0)
save(joinpath(OUT, "p2_5_curves.png"), fig3)
println("saved to $OUT")
```

- [ ] **Step 2: Run it**

Run: `julia --project=research/baseline research/phase2/p2_5_summary.jl`
Expected: 50 classification lines, the counts, and three figures.

- [ ] **Step 3: Record, tick the tracker, commit**

Add a section "P2.5: sensitivity to the spatial distribution at fixed ρ (walls, free decay)" to `research/phase2_record.md`. It should contain:
- the method paragraph (fixed ρ via θ = ρ√(π/D); paired trajectories; thresholds; τ rule; the step exclusion);
- one table per sweep from `observability.csv` (shape, swept value, τ, centroid, ΔS ± SEM, t of the peak, level, mono-exponential residual, baseline curvature);
- the three figures;
- **P2.5.5 summary paragraph (numbers only):**
  - in which configurations |ΔS| exceeds the SNR-100 and SNR-200 thresholds, and with which shape;
  - how |ΔS| scales with h, w, D and ρ;
  - whether plotting against the centroid brings the shapes together;
  - whether the layer curves stay mono-exponential (residual), i.e. whether a fitted conventional model would absorb the difference.

In `research/project-progress-tracker.md`, tick P2.5.1–P2.5.5 with their outputs. P2.5.1 uses the small-h rows of the h sweep (h ≤ 0.1, i.e. h/w ≤ 0.05). Increase the Phase 2 "Done" count by 5.

```bash
git add research/phase2/p2_5_summary.jl research/results/phase2/p2_5_summary research/phase2_record.md research/project-progress-tracker.md
git commit -m "P2.5: summarise sensitivity to the spatial distribution at fixed rho

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
