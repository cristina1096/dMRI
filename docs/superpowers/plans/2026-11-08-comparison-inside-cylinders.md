# Comparison with the conventional model inside cylinders (Phase 3.7a) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compare the near-surface layer (step, linear, exponential) with the conventional surface model at the same ρ, inside single cylinders, using the exact radial continuum. Monte Carlo is used only for spot checks.
- **Apparent T₂ against axon radius:** the true analogue of the paper's Fig. 9A, at the T₂-optimal echo time TE* = T₂.
- **Dependence on h (and centroid), D and ρ.**
- **A dimensionless map over ρR/D and h/R:** where the difference exceeds single-measurement noise at known ρ, and where it remains after fitting the conventional ρ freely.
- **The D-dependence signature.**

**Architecture:** Research code only. Tested helpers in `research/phase3/analysis.jl` (apparent T₂, two-point T₂ noise, TE*, best-fit conventional, conventional ceiling). Five continuum scripts take seconds to minutes. One Monte Carlo spot-check script uses the comparison timestep rule.

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

### Task 2: Apparent T₂ against radius — the Fig. 9A analogue (P3.7.1)

**Files:** Create `research/phase3/p3_7_1_fig9a.jl` → `research/results/phase3/p3_7_1_fig9a/` (`t2_vs_R.csv`, `summary.json`, `p3_7_1_fig9a.png`)

**Method.** ρ = ρ_ref = 1/140 µm/ms (surface property fixed while R varies), D = 3, bulk R₂ = 1/2000. For R ∈ {0.25, 0.5, 0.75, 1, 1.5, 2, 3}:
- **conventional:** T₂_conv = TE* from the continuum (Robin, ρ), by `te_star` on `radial_signal(...)`;
- **every profile,** in two matchings:
  - same h: h = 0.1 µm, all R; linear needs h ≤ R, satisfied;
  - same centroid: c̄ = 0.05 µm (step h = 0.1, linear h = 0.15, exponential by bisection with support R);
- **at TE*:** S_layer(TE*) directly, T₂_app = apparent_t2(S, TE*), ΔT₂/T₂ = T₂_app/T₂_conv − 1;
- **noise:** 2σ at SNR 100 and 200 at TE*;
- **levels:** "observable at SNR X" if |ΔT₂/T₂| ≥ 2σ_X.
- **Also:** the T₂_conv(R) curve itself, the paper's Fig. 9A trend (T₂ ↑ with R through S/V = 2/R), and its ratio to 1/(R₂,bulk + 2ρ/R).

- [ ] **Step 1: Write the script**

```julia
# P3.7.1 (inside, Fig. 9A analogue): apparent T₂ against axon radius at fixed ρ, conventional vs step/linear/exponential,
# from the exact radial continuum. ρ = 1/140 µm/ms, D = 3, bulk R2 = 1/2000; TE* = T₂ of the conventional model at each R.
# Matchings: same h = 0.1 µm, and same centroid c̄ = 0.05 µm. Noise: two-point T₂ at TE*, 2σ at SNR 100 / 200.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_7_1_fig9a.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_7_1_fig9a")
const RHO, D, R2B = 1 / 140, 3.0, 1 / 2000
const RS = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
const SHAPES = ("step", "linear", "exponential")
cbar(shape, h, R) = shape == "step" ? h / 2 : shape == "linear" ? h / 3 : (c = R / h; h * (1 - (c + 1) * exp(-c)) / -expm1(-c))
function h_for_cbar(shape, cb, R)
    shape == "step" && return 2cb
    shape == "linear" && return 3cb
    lo, hi = 1e-6, 100R
    for _ in 1:200
        mid = sqrt(lo * hi); cbar("exponential", mid, R) < cb ? (lo = mid) : (hi = mid)
    end
    return sqrt(lo * hi)
end
S(R, kw...) = t -> radial_signal(; R=R, D=D, times=[t], R2_bulk=R2B, kw...)[1]
rows = NamedTuple[]
for R in RS
    Sc = S(R, :surface => RHO)
    T2c = te_star(Sc; t_hi=20 * 1 / (R2B + 2RHO / R))
    n200, n100 = 2t2_noise_rel(exp(-1), 200), 2t2_noise_rel(exp(-1), 100)
    for shape in SHAPES, (match, h) in (("same_h", 0.1), ("same_centroid", h_for_cbar(shape, 0.05, R)))
        (shape == "linear" && h > R) && continue
        Sl = radial_signal(R=R, D=D, times=[T2c], R2_bulk=R2B, layer=(shape, RHO, h))[1]
        T2l = apparent_t2(Sl, T2c)
        d = T2l / T2c - 1
        push!(rows, (R=R, shape=shape, match=match, h=h, centroid=cbar(shape, h, R), T2_conv=T2c, T2_layer=T2l, dT2_rel=d,
            twosigma_snr200=n200, twosigma_snr100=n100, observable_snr200=abs(d) >= n200, observable_snr100=abs(d) >= n100,
            T2_fast_diffusion_ref=1 / (R2B + 2RHO / R)))
        @printf("R = %.2f %-12s %-13s h = %.4f: T₂ conv %.2f, layer %.2f ms, ΔT₂/T₂ %+.2e (2σ@200 %.3f)\n", R, shape, match, h, T2c, T2l, d, n200)
    end
end
write_csv(joinpath(OUT, "t2_vs_R.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rho" => RHO, "D" => D, "R2_bulk" => R2B,
    "max_abs_dT2_rel" => maximum(abs(r.dT2_rel) for r in rows), "n_observable_snr200" => count(r -> r.observable_snr200, rows)))
fig = Figure(size=(1100, 420))
ax1 = Axis(fig[1, 1], xlabel="R (µm)", ylabel="T₂ (ms)", title="apparent T₂ vs radius, ρ = 1/140 µm/ms")
q = filter(r -> r.shape == "step" && r.match == "same_h", rows)
scatterlines!(ax1, [r.R for r in q], [r.T2_conv for r in q], color=:black, label="conventional")
for (sh, col) in zip(SHAPES, (:steelblue, :darkorange, :seagreen))
    q = filter(r -> r.shape == sh && r.match == "same_h", rows)
    scatterlines!(ax1, [r.R for r in q], [r.T2_layer for r in q], color=col, label="$sh (h = 0.1)")
end
axislegend(ax1, position=:lt, labelsize=9)
ax2 = Axis(fig[1, 2], xlabel="R (µm)", ylabel="ΔT₂/T₂ (layer − conventional)", title="difference vs 2σ two-point noise at SNR 200")
for (sh, col) in zip(SHAPES, (:steelblue, :darkorange, :seagreen)), (m, ls) in (("same_h", :solid), ("same_centroid", :dash))
    q = filter(r -> r.shape == sh && r.match == m, rows)
    lines!(ax2, [r.R for r in q], [r.dT2_rel for r in q], color=col, linestyle=ls, label="$sh, $m")
end
hlines!(ax2, [-2t2_noise_rel(exp(-1), 200), 2t2_noise_rel(exp(-1), 200)], color=:gray, linestyle=:dot)
axislegend(ax2, position=:lb, labelsize=8)
save(joinpath(OUT, "p3_7_1_fig9a.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 2: Commit the script, run it** (≈ 1–2 min), inspect the figure, commit the results.

---

### Task 3: Dependence on h (and centroid), D and ρ at R = 1 (P3.7.1)

**Files:** Create `research/phase3/p3_7_1_params.jl` → `research/results/phase3/p3_7_1_params/` (`sweeps.csv`, `summary.json`, `p3_7_1_params.png`)

**Method.** R = 1, bulk 1/2000, ΔT₂/T₂ at TE* = T₂_conv for each condition, with the 2σ thresholds as in Task 2. One parameter is varied at a time from the reference (ρ = 1/140, D = 3, h = 0.1 or c̄ = 0.05):
- h ∈ {0.02, 0.05, 0.1, 0.2, 0.5, 1.0} (linear h ≤ R) and centroid c̄ ∈ {0.01, 0.025, 0.05, 0.1, 0.25};
- D ∈ {0.5, 1, 2, 3};
- ρ ∈ {1/560, 1/280, 1/140, 1/70, 1/35, 0.1}, with T₂_conv from about 280 ms down to about 5 ms.

The script follows the Task 2 pattern: a loop over sweeps, writing one row per (sweep, value, shape, match) with the same columns plus `sweep` and `value`. The figure has four panels (h, centroid, D, ρ) of ΔT₂/T₂ with the SNR-200 band.

- [ ] **Step 1: Write the script** (same structure and helpers as Task 2; the sweep definitions above are its only inputs).
- [ ] **Step 2: Commit the script, run it, commit the results.**

---

### Task 4: The dimensionless map — fixed ρ and best-fit conventional (P3.7.3 inside)

**Files:** Create `research/phase3/p3_7_3_map.jl` → `research/results/phase3/p3_7_3_map/` (`map.csv`, `summary.json`, `p3_7_3_map.png`)

**Method.** This is the cylinder analogue of the Phase 2 Review C1 map.
- **Grid:** R = 1, D = 3, no bulk; ρR/D ∈ {0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30}; h/R ∈ {0.01, 0.05, 0.1, 0.25, 0.5} (plus 1.0 for step); all three profiles.
- **Readouts:** 41, from 0 until the conventional S ≈ 0.01 (t_end = ln 100/λ₁,conv, with λ₁,conv the lowest continuum rate of the Robin problem from `radial_rates`).
- **Fixed ρ:** max_t |S_layer − S_conv(ρ)|.
- **Best-fit conventional:** `best_fit_conventional` with `Sconv_of_rho = ρc -> radial_signal(R, D, times, surface = ρc)`, rho0 = ρ. Report ρ_c/ρ, the residual and `at_edge`.
- **Ceiling:** the conventional rate cannot exceed D·j₀₁²/R² = 17.35 ms⁻¹. Flag rows where the layer's lowest rate (`radial_rates(layer)[1]`) exceeds the ceiling.
- **Thresholds:** the first ρR/D where each quantity reaches 0.01 (SNR 200), per profile and h/R.
- **Grid convergence (Review Focus 3):** h/R = 0.01 at ρR/D = 1, N = 2000 against 4000; report the difference.

- [ ] **Step 1: Write the script.** Its loops are profile × h/R × ρR/D. Per point it calls one continuum solve for the layer, one for the conventional model at ρ, and about 61 + 60 solves for the best fit. In total about 9 × 16 × 125 ≈ 18 000 solves at about 0.05 s each, ≈ 15 min with 8 threads (`Threads.@threads` over the ρR/D values). The figure shows, per profile, heatmaps of max |ΔS| at fixed ρ and of the best-fit residual against (ρR/D, h/R), with the 0.01 contour.
- [ ] **Step 2: Commit the script, run it, commit the results.**

---

### Task 5: The D-dependence signature (P3.7.3 inside)

**Files:** Create `research/phase3/p3_7_3_ddep.jl` → `research/results/phase3/p3_7_3_ddep/`

**Method.** This is the analogue of P2.5.9 §4. R = 1, bulk 1/2000, ρ ∈ {0.01, 0.1, 1}, one TE per ρ (the conventional T₂ at D = 3), D ∈ {0.5, 1, 2, 3}. For the conventional model and each profile at h/R ∈ {0.05, 0.1, 0.25, 0.5, 1.0} (linear ≤ 1, exponential with support R): T₂_app(D) and T₂(0.5)/T₂(3) − 1. The two-measurement comparison band is √2 × the single 2σ. The figure plots the ratio against h/R per ρ, with the conventional value as a line.

- [ ] Write, commit, run, commit results (≈ 1 min).

---

### Task 6: Monte Carlo spot checks (P3.7.1 "continuum + MC")

**Files:** Create `research/phase3/p3_7_6_mc_spots.jl` → `research/results/phase3/p3_7_6_mc_spots/`

**Method.** Three points chosen from Tasks 2–4, paired against the conventional model in Monte Carlo (θ = ρ√(π/D)), with the timestep τ = min(0.135·c̄², 0.005/ΔR₂(0)) and readouts every TE*/10 up to 2·TE*:
1. **Fig. 9A point:** R = 0.5, exponential at same centroid c̄ = 0.05, ρ = 1/140, bulk 1/2000;
2. **Strong-relaxation point:** R = 1, ρR/D = 1, exponential h/R = 0.25, no bulk (where the best-fit residual of Task 4 is largest within the scanned non-edge rows; use the actual argmax from `map.csv`);
3. **D-dependence point:** R = 1, ρ = 0.1, D = 0.5, linear h/R = 0.5.

Criterion: S and ΔS against the continuum, max |z| ≤ 3 over readouts. Seeds 1–10, 2·10⁴ spins in the box. Run detached (point 2 has a large ΔR₂(0) and τ ≈ 1e-4, so it is the most expensive).

- [ ] Write, commit, run detached, commit results.

---

### Task 7: Record and tracker

- [ ] Write the record section "Comparison with the conventional model inside cylinders (P3.7a)" with Design / Method / Result / Conclusion and the figures:
  - the T₂(R) curve and ΔT₂/T₂ against R at the two matchings;
  - the parameter sweeps;
  - the map with thresholds, best-fit residuals and ceiling flags;
  - the D-dependence;
  - the MC spot checks;
  - **a comparison with walls** (Phase 2 C1 map and P2.5.9): where the distribution becomes observable on cylinders compared with walls, in ρR/D against ρw/D.
- [ ] Tracker:
  - P3.7.1 ◐ (inside done; arrays in 3.7b);
  - P3.7.3 ◐ (inside done);
  - P3.7.4 ◐ (inside paragraph written).
- [ ] Research suite PASS; commit.

---

## After the last task

Final whole-branch review by a fresh subagent on the most capable model (user rule). Fix Critical and Important findings, and ledger the minors.
