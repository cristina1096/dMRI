# Comparison with the conventional model inside cylinders (Phase 3.7a) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compare the near-surface layer (step, linear, exponential) with the conventional surface model at the same ρ, inside single cylinders, using the exact radial continuum. Monte Carlo is used only for spot checks.
- **Apparent T₂ against axon radius:** the true analogue of the paper's Fig. 9A, at the T₂-optimal echo time TE* = T₂.
- **Dependence on h (and centroid), D and ρ.**
- **A dimensionless map over ρR/D and h/R:** where the difference exceeds single-measurement noise at known ρ, and where it remains after fitting the conventional ρ freely.
- **The D-dependence signature.**

**Architecture:** Research code only.
- **Helpers:** tested in `research/phase3/analysis.jl` (apparent T₂, two-point T₂ noise, TE*, best-fit conventional, conventional ceiling).
- **Monte Carlo is primary** (user, 2026-11-08). One tested, paired MC runner computes every condition, with the comparison timestep rule. The exact radial continuum supports it: dense curves, a check of every MC point, TE*, and the conventional family for the best fit.
- **Runs:** detached overnight chains (≈ 25 h total).

**Tech Stack:** Julia 1.12, research environment `research/baseline` (MCMRSimulator, CairoMakie, SpecialFunctions), the radial continuum (`research/phase3/continuum_radial.jl`, validated for all profiles in P3.4.6/P3.4.6d).

**Spec:**
- `research/project-progress-tracker.md` 3.7: P3.7.1 (inside part), P3.7.3 (inside part), P3.7.4 (inside part).
- Phase 2 analogues in `research/phase2_record.md`:
  - P2.5.6/7, the Fig. 9/8A analogues;
  - the Review C1 follow-up map over ρw/D with best-fit conventional and the ceiling D·π²/w²;
  - P2.5.9, the D dependence of the apparent T₂.
- User decisions:
  - 3.7 is split into 3.7a (inside, this plan) and 3.7b (arrays and packing, annulus model plus MC);
  - all three profiles;
  - P3.6.2 option (b): comparison MC uses τ ≤ 0.135·c̄_min².

## Global Constraints

- **Branch.** Work only on branch `cc/near-surface-r2`. Never modify or commit on other branches. Do not push; the user pushes.
- **Commits.** Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Framing.**
  - The conventional model is a reference, not the truth.
  - No fast-diffusion or Brownstein–Tarr assumptions. The exact absorbing-wall eigenvalue D·j₀₁²/R² is used only as a bound, the conventional model's fastest possible decay.
- **Fixed ρ.**
  - ρ is per unit surface area (D3.1).
  - The conventional model uses the same ρ in the continuum (Robin) and θ = ρ√(π/D) in Monte Carlo (P3.1.6).
  - Bulk R₂ = 1/2000 ms⁻¹ wherever T₂ is reported (as in Phase 2 P2.5.6/7).
- **Physiological reference.** ρ_ref = 1/140 µm/ms, so that inside R = 1 the surface rate 2ρ/R = 1/70 ms⁻¹, which gives T₂ ≈ 70 ms (the paper's tissue T₂ scale). Radii R = 0.25, 0.5, 1, 2, 3 µm. D = 3 unless scanned.
- **Noise.**
  - Single-curve threshold |ΔS| ≥ 2/SNR (0.02 at SNR 100, 0.01 at SNR 200).
  - Apparent T₂ by two points: T₂_app = −TE/ln(S(TE)/S(0)), with σ_T₂/T₂ = √(1 + r²)/(SNR·r·|ln r|), r = S(TE)/S(0). At TE* = T₂_conv, 2σ = 0.029 at SNR 200 and 0.058 at SNR 100.
- **Comparison timestep (memory: comparison-timestep-rule).** Monte Carlo spot checks use τ = min(0.135·c̄_min², 0.005/ΔR₂(0)_max). The expected difference error is ≲ 10 %.
- **Commit scripts before running them** (provenance). Records: Design / Method / Result / Conclusion with formulas and figures. Signs and orderings are recorded as found.
- **Commands.**
  - Research tests: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
  - Scripts: `julia --project=research/baseline -t 8 research/phase3/<script>.jl`
- **Plan end.** A fresh subagent reviewer (most capable model) reviews the whole branch (user rule).

## Review Focus

1. **The conventional ceiling.** Inside a cylinder, a surface sink cannot decay faster than D·j₀₁²/R² (j₀₁ = 2.4048). The best-fit search must report when ρ_c runs to the edge or the ceiling, and must not quote a residual from an unconverged fit.
2. **TE* when T₂ is very short (large ρR/D).** The readout grid must resolve TE*, and S(TE*) must be computed directly, not interpolated from a coarse grid.
3. **Thin layers in the continuum.** h/R down to 0.01 must have enough rings. N = 2000 gives 20 rings at h/R = 0.01. Check grid convergence at the thinnest point of each map.
4. **Same centroid.** h is chosen per profile by its centroid (exponential support R, by bisection, as in P3.4.8). It must not be confused with the same h.
5. **Monte Carlo spot checks at strong relaxation.** Both timestep limits apply (rule and rate); readouts must start early enough to resolve the fast decay.

---

### Task 1: Comparison helpers with tests

**Files:**
- Modify: `research/phase3/analysis.jl`
- Modify: `research/phase3/test/test_analysis.jl`

**Interfaces — Produces:**
- `apparent_t2(S_TE, TE)` = −TE/ln(S_TE) (S normalised to S(0) = 1);
- `t2_noise_rel(r, snr)` = √(1 + r²)/(snr·r·|ln r|);
- `conventional_ceiling(R, D)` = D·(2.404825557695773/R)²;
- `te_star(Sfun; t_hi)`: the t with Sfun(t) = e^(−1), by bisection on [0, t_hi], where Sfun is monotone decreasing;
- `best_fit_conventional(S_layer, times, Sconv_of_rho; rho0, lo=exp(-4), hi=exp(8), ngrid=61, iters=60)` → `(rho_c, residual, at_edge::Bool)`. It minimises max_t |S_layer(t) − S_conv(t; ρ_c)| over ρ_c ∈ [lo·rho0, hi·rho0]: a log-grid scan, then golden-section search around the best grid point. `at_edge` is true when the best point is the first or last grid point.

- [ ] **Step 1: Write the failing tests** (append to `test_analysis.jl`)

```julia
@testset "comparison helpers (3.7a)" begin
    @test apparent_t2(exp(-1), 70.0) ≈ 70.0
    r = exp(-1)
    @test t2_noise_rel(r, 200) ≈ sqrt(1 + r^2) / (200 * r) rtol=1e-12          # |ln r| = 1
    @test 2 * t2_noise_rel(r, 200) ≈ 0.0289 atol=1e-4                            # the 2σ at TE* quoted in the plan
    @test conventional_ceiling(1.0, 3.0) ≈ 3 * 2.404825557695773^2
    @test te_star(t -> exp(-t / 70); t_hi=1000.0) ≈ 70.0 rtol=1e-9
    # best fit recovers a known rate for a mono-exponential family
    times = collect(0.0:5.0:100.0)
    Sl = exp.(-times ./ 40)
    fam(ρ) = exp.(-ρ .* times)
    (ρc, res, edge) = best_fit_conventional(Sl, times, fam; rho0=0.01)
    @test ρc ≈ 1 / 40 rtol=1e-6
    @test res < 1e-8
    @test !edge
    # a target faster than the family can reach → the fit runs to the edge and says so
    capped(ρ) = exp.(-min(ρ, 0.02) .* times)
    (_, res2, edge2) = best_fit_conventional(Sl, times, capped; rho0=0.01)
    @test edge2 || res2 > 1e-3
end
```

- [ ] **Step 2: Run to verify they fail** — `UndefVarError: apparent_t2`.

- [ ] **Step 3: Implement** (append to `analysis.jl`)

```julia
"Two-point apparent T₂ = −TE/ln(S(TE)/S(0)) for a signal normalised to S(0) = 1."
apparent_t2(S_TE, TE) = -TE / log(S_TE)

"Relative noise of the two-point T₂ at SNR `snr`: σ_T₂/T₂ = √(1 + r²)/(SNR·r·|ln r|), r = S(TE)/S(0)."
t2_noise_rel(r, snr) = sqrt(1 + r^2) / (snr * r * abs(log(r)))

"Fastest decay of any surface sink inside a cylinder (absorbing wall): D·j₀₁²/R², j₀₁ = 2.404825557695773."
conventional_ceiling(R, D) = D * (2.404825557695773 / R)^2

"The t in [0, t_hi] with Sfun(t) = e^(−1) for a monotonically decreasing Sfun (bisection, 200 steps)."
function te_star(Sfun; t_hi)
    lo, hi = 0.0, Float64(t_hi)
    Sfun(hi) > exp(-1) && error("te_star: S(t_hi) > e^(−1); increase t_hi")
    for _ in 1:200
        mid = (lo + hi) / 2
        Sfun(mid) > exp(-1) ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
    best_fit_conventional(S_layer, times, Sconv_of_rho; rho0, lo=exp(-4), hi=exp(8), ngrid=61, iters=60)

Minimax fit of the conventional model's ρ_c to a layer curve: ρ_c minimising max_t |S_layer(t) − S_conv(t; ρ_c)|.
A log-grid scan over [lo·rho0, hi·rho0], then golden-section search between the neighbours of the best grid point.
Returns (rho_c, residual, at_edge); at_edge = true when the best grid point is the first or last (the fit wants to
leave the searched range, e.g. beyond the conventional ceiling).
"""
function best_fit_conventional(S_layer, times, Sconv_of_rho; rho0, lo=exp(-4), hi=exp(8), ngrid=61, iters=60)
    obj(ρ) = maximum(abs.(S_layer .- Sconv_of_rho(ρ)))
    grid = rho0 .* exp.(range(log(lo), log(hi), length=ngrid))
    vals = obj.(grid)
    k = argmin(vals)
    at_edge = k == 1 || k == ngrid
    a, b = log(grid[max(k - 1, 1)]), log(grid[min(k + 1, ngrid)])
    g = (sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = obj(exp(c)), obj(exp(d))
    for _ in 1:iters
        if fc < fd
            b, d, fd = d, c, fc
            c = b - g * (b - a); fc = obj(exp(c))
        else
            a, c, fc = c, d, fd
            d = a + g * (b - a); fd = obj(exp(d))
        end
    end
    ρc = exp((a + b) / 2)
    return (ρc, obj(ρc), at_edge)
end
```

- [ ] **Step 4: Run to verify they pass.** **Step 5: Commit** ("Phase 3 comparison helpers: apparent T₂, two-point noise, TE*, conventional ceiling, minimax best-fit conventional").

---

## Revision (user, 2026-11-08): Monte Carlo is the primary method, the continuum is support

Tasks 2–5 are computed by **Monte Carlo**. The radial continuum supports them in three ways:
- **Dense backdrop curves** around the MC points;
- **An exact check** of every MC point: S and ΔS, |z| reported;
- **TE*:** the readout grid needs it before the MC run, so it is computed from the conventional continuum. The MC conventional T₂ is also reported.

The conventional **family** for the best-fit in Task 4 is the continuum Robin solution, which is exact: fitting by repeated MC runs at many θ would cost ≈ 60× more. The former Task 6 (MC spot checks) is superseded and removed.

**Budget rules.**
- **Timestep:** τ = min(0.135·c̄_min², 0.005/ΔR₂(0)_max) per condition (comparison timestep rule).
- **Centroids:** c̄ ≥ 0.1 µm wherever the question allows (τ ≤ 1.35e-3). c̄ = 0.05 only where it is the question.
- **Statistics:** seeds 1–10, 2·10⁴ spins in the box [−R, R]³. All models of a condition share τ, so ΔS and ΔT₂ are paired.
- **Readouts:** 0 to 1.5·TE*, in steps of TE*/10 (Tasks 2, 3, 5); to S_conv ≈ 0.01 (Task 4).
- **Execution:** long runs go in detached chains, and scripts are committed before they run.

### Task 2: The shared paired MC runner (with tests)

**Files:** Create `research/phase3/p3_7_mc.jl` (library: `run_condition`) and `research/phase3/p3_7_run.jl` (driver: `TASK` = fig9a | params | map | ddep → `research/results/phase3/p3_7_<task>/conditions.csv`, `curves.csv`, `summary.json`). Test: `research/phase3/test/test_p3_7_mc.jl`.

**Interfaces — Produces** `run_condition(; R, D, rho, R2_bulk, models, times, tau, nspins=20_000, seeds=1:10)`:
- `models` is a vector of (name, shape, h), with shape "conventional" for MCMR θ = ρ√(π/D);
- it returns a NamedTuple with, per model:
  - S(t): mean and SEM over seeds;
  - the paired ΔS(t) against "conventional": mean and SEM;
  - the continuum S(t) for the same model (step/linear/exponential layer, or Robin);
  - `identical_trajectories`.

- [ ] **Step 1: Failing test** (`test_p3_7_mc.jl`). A tiny condition: R = 1, D = 3, ρ = 0.02, conventional plus step h = 0.2, times [0, 10], τ = 1e-2, 500 spins, 2 seeds. Check:
  - the shapes of the outputs;
  - `identical_trajectories == true`;
  - S(0) = 1;
  - the continuum column equals `radial_signal` at the same times;
  - ΔS(0) = 0.
- [ ] **Step 2: Implement.** Build on `run_cylinders` (harness) for each model with `timestep = tau`. The pairing check uses the final positions; ΔS is computed per seed. Include `continuum_radial.jl` for the continuum columns.
- [ ] **Step 3: The driver** `p3_7_run.jl`. It holds the condition lists of Tasks 3–6 below and writes one row per (condition, model, t). Commit before any run.

### Task 3: Fig. 9A analogue — apparent T₂ against axon radius (P3.7.1)

**MC conditions** (`TASK=fig9a`): ρ = 1/140, D = 3, bulk 1/2000.
- R ∈ {0.5, 1, 2}, same centroid c̄ = 0.1 (step h = 0.2, linear h = 0.3, exponential by bisection with support R), plus the conventional model, τ = 1.35e-3;
- the same-h matching at R = 1 only: h = 0.2 for all, so c̄_min = 0.067 (linear) and τ = 6.0e-4.

**Quantities.** TE* = continuum conventional T₂ at each R. At TE*, interpolated between the readouts TE*·{0.9, 1.0, 1.1}, of which 1.0 is exact:
- T₂_app = −TE*/ln S(TE*) for every model;
- ΔT₂/T₂ (paired);
- 2σ two-point noise at SNR 100 and 200;
- the MC − continuum z.

**Support.** The continuum dense curve over R ∈ [0.25, 3] for all profiles in both matchings, overlaid on the figure.

**Cost** ≈ 8 h detached (R = 2 is the longest: TE* ≈ 129 ms, 1.4·10⁵ steps).

### Task 4: Parameter sweeps at R = 1 (P3.7.1)

**MC conditions** (`TASK=params`), one parameter at a time from the reference (ρ = 1/140, D = 3, c̄ = 0.1, the same-centroid matching, reused from Task 3 R = 1):
- **D** ∈ {0.5, 1, 2};
- **ρ** ∈ {1/280, 1/70, 1/35, 0.1} (T₂_conv ≈ 128, 34, 17, 5 ms);
- **c̄** ∈ {0.05, 0.25} (τ = 3.4e-4 and 8.4e-3).

The quantities are those of Task 3. **Support:** continuum dense curves over each sweep, plus ρ = 1/560 and c̄ = 0.01–0.025, which are continuum-only and labelled as such. **Cost** ≈ 12 h detached; the c̄ = 0.05 and ρ = 1/280 conditions dominate.

### Task 5: Dimensionless map over ρR/D and h/R (P3.7.3 inside)

**MC conditions** (`TASK=map`): R = 1, D = 3, no bulk; ρR/D ∈ {0.1, 1, 10} × h/R ∈ {0.1, 0.25, 0.5} × the three profiles, plus the conventional model per ρ. For each ρ, the τ is shared by all models (rule and rate limits). Readouts are 41 points to t_end = ln 100/λ₁,conv. At strong ρ the decay is fast, so runs are short.

**Quantities.**
- **Fixed ρ:** max_t |ΔS| from paired MC, with SEM.
- **Best fit:** `best_fit_conventional(S_layer_MC, times, ρc -> continuum Robin S(t; ρc))`. Report ρ_c/ρ, the residual (with the MC SEM of S_layer as its noise floor) and `at_edge`.
- **Ceiling flag:** whether the layer's lowest continuum rate exceeds D·j₀₁²/R².
- **Thresholds:** the first ρR/D where the fixed-ρ |ΔS| and the best-fit residual reach 0.01.

**Support.** The continuum dense map (ρR/D ∈ {0.003, …, 30} × h/R ∈ {0.01, …, 1}) as the heatmap behind the MC points. Grid convergence at h/R = 0.01, ρR/D = 1 (N = 2000 against 4000). **Cost** ≈ 3 h detached.

### Task 6: D-dependence signature (P3.7.3 inside)

**MC conditions** (`TASK=ddep`): R = 1, bulk 1/2000, ρ ∈ {0.1, 1}, D ∈ {0.5, 1, 2, 3}, conventional plus the three profiles at h/R ∈ {0.25, 0.5, 1.0} (linear ≤ 1; exponential support R). TE per ρ is the continuum conventional T₂ at D = 3.

**Quantities.** T₂_app(D) for every model and T₂(0.5)/T₂(3) − 1. The two-measurement band is √2 × the single 2σ.

**Support:** continuum values for every point and for ρ = 0.01 (continuum-only). **Cost** ≈ 2 h detached.

### Task 7: Record and tracker

- [ ] Write the record section "Comparison with the conventional model inside cylinders (P3.7a)" with Design / Method / Result / Conclusion and the figures:
  - the T₂(R) curve and ΔT₂/T₂ against R at the two matchings;
  - the parameter sweeps;
  - the map with thresholds, best-fit residuals and ceiling flags;
  - the D-dependence;
  - MC against the continuum for every point;
  - **a comparison with walls** (Phase 2 C1 map and P2.5.9): where the distribution becomes observable on cylinders compared with walls, in ρR/D against ρw/D.
- [ ] Every MC result is shown against the continuum (z of the difference).
- [ ] Tracker:
  - P3.7.1 ◐ (inside done; arrays in 3.7b);
  - P3.7.3 ◐ (inside done);
  - P3.7.4 ◐ (inside paragraph written).
- [ ] Research suite PASS; commit.

---

## After the last task

Final whole-branch review by a fresh subagent on the most capable model (user rule). Fix Critical and Important findings, and ledger the minors.
