# Overlapping layers in arrays (Phase 3.5) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Characterise outer layers that overlap between neighbouring cylinders in arrays (P3.5.1–P3.5.4):
- which configurations overlap;
- whether the summed rate stays realistic (no interior maximum above the surface value);
- that the package sums multiple nearby surfaces exactly;
- how the signal changes with the overlap degree compared with the conventional model at the same ρ.

**Architecture:** Research code plus package tests.
- **Configurations and realism** (Tasks 1–2) are continuum and geometry only: a brute-force summed rate on a fine 2D grid of the unit cell.
- **Exact summation** (Task 3) is package tests against brute force.
- **The signal study** (Tasks 4–5) is paired Monte Carlo with the 3.7b runner (`run_condition_outside`, extended to rectangular arrays). No continuum exists for overlapping arrays; the annulus only sets readout windows and TE.
- **Two groups:**
  - A: Tasks 1–3, then review;
  - B: Tasks 4–6, then the final review.

**Tech Stack:** Julia 1.12; package `src/geometries/internal/layers.jl` (cylinder layers, image loop, `max_layer_rate`), package tests `test/test_layer_cylinders.jl`; research env `research/baseline`; `research/phase3/harness.jl` (`brute_rate_profile`), `p3_7b_mc.jl`, `continuum_annulus.jl`, `analysis.jl`.

**Spec:**
- `research/project-progress-tracker.md` §3.5 (P3.5.1–P3.5.4);
- Phase 2 analogue P2.5.9 (overlapping layers on walls: profile realism, comparison against overlap, D dependence) in `research/phase2_record.md`;
- 3.7b results for s = 2.5 (overlapping layers) in `research/phase3_record.md`.
- User decisions (2026-11-11):
  - vary h at fixed spacing;
  - ρR/D = 1 only for the signal study;
  - a small D-dependence check.

## Global Constraints

- **Branch.** Work only on `cc/near-surface-r2`. Never modify or commit other branches. Do not push.
- **Commits.** Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Framing.**
  - The conventional model is a reference, not the truth. No fast-diffusion or Brownstein–Tarr arguments.
  - Overlapping layers ADD, which is the package convention: the rate at a point is the sum over every cylinder within its support (P3.4 / D3.4).
- **Geometry.** R = 1 µm, outer layers only (impermeable cylinders; the inside layer is irrelevant outside), same ρ for every cylinder.
  - **Exponential support** is the default `layer_cutoff`: the smallest surface gap of each cylinder.
  - **Step and linear** support h must be ≤ that gap, or the package errors (`check_outer_reach`).
- **Configurations** (P3.5.1), classified per direction by gap g_x, g_y against the support L:
  - separated: L ≤ g/2 in both directions;
  - one-direction overlap;
  - two-direction overlap;
  - plus k-fold interstitial overlap where the four diagonal neighbours' supports meet (s/√2 − R < L in a square array).
- **Comparison timestep** (memory: comparison-timestep-rule): τ = min(0.135·c̄_min², 0.005/ΔR₂,max). Here ΔR₂,max is the **summed** peak rate over the unit cell (overlap included; 3.7 final-review lesson), taken from `max_layer_rate` of the package. c̄ < 0.05 µm is flagged as extrapolated.
- **Statistics and wording (3.7 lessons).**
  - Seeds 1000·k + (1:10) per condition (independent conditions); models within a condition are paired.
  - One threshold convention (`threshold_status`).
  - State what the data exclude, name pooling methods and correlated measurements, label post-hoc decisions.
- **Provenance.** Commit scripts before running them. Long runs go in detached `nohup` chains with retries and resume. Verify kills with `ps`.
- **Cost (user: avoid too-thin layers).** h ≥ 0.15 µm in MC. Report the smoke-projected runtime before each MC launch, and ask if it exceeds 1.5× the estimate.
- **Commands.**
  - Package tests: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders"])'` for the touched file, then the full `Pkg.test()` before the commit.
  - Research tests: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`.
- **Records.** Design / Method / Result / Conclusion with formulas and figures; every number checked against the committed CSVs.
- **Group ends.** A fresh subagent reviewer (most capable model) after each group.

## Review Focus

1. **k-fold interstitial overlap.** At wide spacing with long support (exponential at s = 4: support 2 µm, interstitial distance s/√2 − R = 1.83 µm) four layers meet. The brute force and the realism grid must include those points, and so must a package test.
2. **Rectangular arrays.** One-direction overlap (g_x ≠ g_y): the exponential support is the SMALLEST gap, in every direction. The brute force must use (sx, sy), not a single s.
3. **Summed-rate τ bound.** With overlap, the peak rate can exceed a single layer's ΔR₂(0) (step: ×2). τ must use the package's summed `max_layer_rate`, and the check is that ΔR₂,max·τ ≤ 0.005 for every condition.
4. **Realism criterion at the surface itself.** "No interior maximum above the surface value" compares max over the fluid region with the single-cylinder rate at d = 0. Points ON a surface where another layer also reaches (narrow gap) must be included.
5. **h = gap for step/linear.** Allowed within the package tolerance (support = gap). The brute force must treat d = L inclusively for step, consistently with the package.

---

## Group A: configurations, realism, exact summation

### Task 1: Configuration table and a rectangular brute force (P3.5.1)

**Files:**
- Modify: `research/phase3/harness.jl` (`brute_rate_profile` takes `sx, sy`, default `sx = sy = s`)
- Create: `research/phase3/p3_5_configs.jl`
- Test: `research/phase3/test/test_p3_5.jl`

**Interfaces — Produces:**
- `brute_rate_profile(x, y; R, s=Inf, sx=s, sy=s, inner, outer)`: summed rate over images within ±3 cells;
- `classify_overlap(sx, sy, R, L)`: returns (x_overlap::Bool, y_overlap::Bool, interstitial_fold::Int) with `x_overlap = L > (sx − 2R)/2` and `interstitial_fold = 4` if √((sx/2)² + (sy/2)²) − R < L, else 2 or 1 from the axes;
- `max_overlap_count(sx, sy, R, L; n=401)`: the maximum number of supports covering a point on a grid of the fluid region.

- [ ] **Step 1: Failing tests** (`test_p3_5.jl`):
  - for a square array (`sx = sy = s`), `brute_rate_profile` with `sx, sy` equals the old call with `s` at 100 random points;
  - in a rectangular array (2.5, 4.0) a point midway between x-neighbours gets two contributions, and one midway between y-neighbours (gap 2, L = 0.5) gets none;
  - `classify_overlap(2.5, 2.5, 1, 0.3) == (true, true, 2)`;
  - `classify_overlap(2.5, 4.0, 1, 0.3) == (true, false, 2)`;
  - `classify_overlap(4.0, 4.0, 1, 2.0)` gives interstitial fold 4;
  - `max_overlap_count(4.0, 4.0, 1, 2.0) == 4`.
- [ ] **Step 2:** Run the research tests. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** `p3_5_configs.jl` writes `results/phase3/p3_5_configs/configs.csv` for:
  - geometries: square 2.5 and 4.0, rectangular (2.5, 4.0);
  - supports: step/linear h ∈ {0.15, 0.25, 0.35, 0.5}, exponential with L = smallest gap;
  - columns: configuration (separated / one-direction / two-direction / interstitial k-fold), max overlap count, and the fraction of the fluid area covered by ≥ 2 supports.
- [ ] **Step 6:** Commit.

### Task 2: Realism of the summed rate in 2D (P3.5.2)

**Files:**
- Create: `research/phase3/p3_5_realism.jl` → `results/phase3/p3_5_realism/`

**Criterion:** realistic ⟺ max over the fluid region of the summed rate ≤ (1 + 1e-9) × the single-cylinder surface value ΔR₂(0) = (ρ/(hκ_out))·g(0), as on walls (P2.5.9).
- The grid is 801 × 801 over the unit cell plus exact evaluation at the gap midpoints, the interstitial point and the surface points facing a neighbour.
- Configurations: Task 1's table.
- **Expected and recorded as found:**
  - step: unrealistic wherever two supports overlap (×2), the counter-example;
  - linear: on walls the sum across a gap is exactly flat at h = gap. On cylinders the curvature (κ_out, 1/r spreading) changes this, so record the max/surface ratio;
  - exponential: record the ratio, including the 4-fold interstitial case at s = 4.
- **Figure:** maps of summed rate / surface value for the overlapping cases, with the maximum marked.

- [ ] **Step 1:** Commit the script, then run it (seconds to minutes).
- [ ] **Step 2:** Check the pass/fail table against Task 1's configurations, then commit the results.

### Task 3: Exact summation of multiple nearby surfaces in the package (P3.5.3)

**Files:**
- Modify: `test/test_layer_cylinders.jl` (package)

Existing coverage (keep): two and four overlapping non-repeating cylinders, a square array with images, two cylinders per cell in a rectangular array (step), and linear/exponential in a square array with overlap — all against brute force.

- [ ] **Step 1: New failing-first tests** (they fail only if the summation is wrong; if they pass at once, record that the behaviour already holds and why each test would catch a regression):
  - (a) **rectangular one-direction overlap** (repeats (2.5, 4.0)), linear h = 0.4 and exponential (default cutoff = 0.5): layer exponent over 200 random straight segments (length ≤ 0.5, both ends in the fluid) against the brute-force line integral (Gauss–Legendre 200 points per segment, profile sum over images), |Δ| ≤ 1e-10 relative to the exponent;
  - (b) **4-fold interstitial**: square 4.0, exponential cutoff 2.0 (default), 200 random segments through a 0.5 µm disc around the interstitial point, ≤ 1e-10;
  - (c) **h = gap**: step and linear h = 0.5 at square 2.5, points on the gap midline and on the facing surfaces: the rate equals the brute force exactly (Review Focus 5);
  - (d) **`max_layer_rate`** ≥ the brute-force grid maximum of the summed rate for (a), (b) and the step at h = gap (Review Focus 3).
- [ ] **Step 2:** Run the package tests for the file, then the full suite. Expected: PASS, the full suite count increased by the new tests.
- [ ] **Step 3:** Commit.
- [ ] **Group A review** (fresh subagent, most capable model), with Review Focus 1–5.

---

## Group B: signal against overlap degree, D check, record

### Task 4: Signal against the overlap degree at fixed ρ (P3.5.4)

**Files:**
- Modify: `research/phase3/p3_7b_mc.jl` (`run_condition_outside` takes `sx, sy`; default `sx = sy = s`), and its test
- Create: `research/phase3/p3_5_run.jl` (driver, crash-safe I/O from 3.7b), `research/phase3/p3_5_overlap_analysis.jl`

**MC conditions** (`TASK=overlap`):
- geometries: square s = 2.5 (two-direction) and rectangular (2.5, 4.0) (one-direction); R = 1, D = 3, no bulk, ρR/D = 1 (ρ = 3);
- models: conventional + step/linear/exponential at h ∈ {0.15, 0.25, 0.35, 0.5}, i.e. h/gap_x = 0.3 (separated), 0.5 (touching), 0.7, 1.0. The exponential support is fixed at the gap, so only its decay length changes;
- τ = min(0.135·c̄_min², 0.005/`max_layer_rate`) over the condition's models (the summed peak);
- 41 readouts to ln 100/λ₁, with λ₁ from the equal-area annulus (window only);
- seeds 1000·k + (1:10).
- **Estimated cost** ≈ 1.5 h (smoke-projected before launch).

**Quantities:**
- max_t |ΔS| (paired MC) against h/gap per profile and geometry;
- the same against the overlap fraction from Task 1;
- the comparison with the separated reference: the h/gap = 0.3 point and the 3.7b s = 4 map at ρR/D = 1;
- monotonicity in h, recorded as found (on walls overlap saturated the fixed-ρ difference, P2.5.9);
- conventional θ√τ; ΔR₂,max·τ for every condition.

- [ ] **Step 1: Failing test:** `run_condition_outside(sx=2.5, sy=4.0, …)` on a tiny condition: identical trajectories, S(0) = 1, rectangular geometry used (spins' periodic mask by (sx, sy)).
- [ ] **Step 2:** Implement and pass the test.
- [ ] **Step 3:** Write the driver and analysis. Smoke them (200 spins, 1 seed), report the projection, commit, run detached, analyse, commit the results.

### Task 5: D-dependence for overlapping layers (P3.5.4, small)

**MC conditions** (`TASK=ddep`):
- square s = 2.5, ρ = 1, bulk 1/2000, D ∈ {0.5, 3};
- conventional + three profiles at h ∈ {0.15, 0.35} (new overlap degrees);
- TE = the annulus conventional T₂ at D = 3, the same TE convention and geometry as 3.7b ddep. The 3.7b ddep points at h ∈ {0.25, 0.5} complete the sweep (independent seeds; stated).
- **Estimated cost** ≈ 0.5 h.

**Quantities:**
- q = T₂(0.5)/T₂(3) − 1 and Δq = q_layer − q_conv with paired σ, against the two-measurement band at SNR 200, as a function of h/gap;
- comparison with walls P2.5.9: there thick/overlapping layers were nearly D-independent.

- [ ] **Step 1:** Run in the same detached chain as Task 4, after it. Analyse, commit.
- [ ] **Group B review** (fresh subagent): Tasks 4–5.

### Task 6: Record and tracker

- [ ] Write the record section "Overlapping layers in arrays (P3.5)" with Design / Method / Result / Conclusion:
  - the configuration table;
  - the realism table and maps (step counter-example);
  - the exact-summation tests;
  - signal against overlap degree;
  - D-dependence;
  - the comparison with walls (P2.5.9);
  - links to the 3.7b overlap statements.
- [ ] Check every number against the CSVs.
- [ ] Tracker: P3.5.1–P3.5.4 → ☑ with one-line results; the Phase 3 counts row.
- [ ] Research suite and package suite PASS; commit.

## After the last task

Final review of the plan's range by a fresh subagent on the most capable model (user rule). Fix Critical and Important findings and ledger the minors. Report the rulings and deferred minors.
