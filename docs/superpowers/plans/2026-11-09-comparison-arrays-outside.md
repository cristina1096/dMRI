# Comparison with the conventional model in arrays and outside cylinders (Phase 3.7b) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compare the near-surface layer (step, linear, exponential) with the conventional surface model at the same ρ **outside** cylinders in square arrays, and in the **summed** intra + extra signal, completing P3.7.1–P3.7.4. The three questions are the same as inside (3.7a):
- **Apparent T₂ against packing density:** the analogue of the paper's Fig. 9B, at TE*.
- **A dimensionless outside map over ρR/D and h/R:** at two spacings, at fixed ρ and with the conventional ρ fitted freely.
- **The outside D-dependence signature.**

**Architecture:** Research code only.
- **Monte Carlo is primary** (user, 2026-11-08), with the paired runner design of 3.7a, extended to the outside compartment of square arrays.
- **Support model.** No exact continuum exists outside. The support is an **annulus continuum**: a Krogh-type cell, R ≤ r ≤ R_out with π R_out² = s², reflecting at R_out. Group A validates it against paired MC, and a stated rule decides where it may be used.
- **Summed signal.** Separate inside and outside runs are combined by volume fraction. The cylinders are impermeable, so the inside layer is irrelevant to outside spins: outside runs set only the outer layer, and inside runs only the inner layer (user, 2026-11-09).
- **Execution in four groups**, each ended by a fresh-subagent review (user, 2026-11-09):
  - A: Tasks 1–3, tools + annulus + validation;
  - B: Task 4, packing density, the long run, alone;
  - C: Tasks 5–6, outside map + D signature;
  - D: Task 7, record + final review.

**Tech Stack:** Julia 1.12, research environment `research/baseline` (MCMRSimulator, CairoMakie, SpecialFunctions). Existing research code: `research/phase3/harness.jl` (`run_cylinders`), `continuum_radial.jl`, `analysis.jl`, `p3_7_mc.jl` (`run_condition`, `read_done`, `write_merged`), `p3_7_shared.jl`.

**Spec:**
- `research/project-progress-tracker.md` 3.7:
  - P3.7.1, arrays part;
  - P3.7.2, packing density / volume fraction, intra and extra separately and summed;
  - P3.7.3 and P3.7.4, outside parts.
- `research/phase3_record.md`:
  - "Comparison with the conventional model inside cylinders (P3.7a)", whose methods, quantities and lessons are reused here;
  - P3.4.8 / P3.4.8e–h, the outside verification and centroid comparison in arrays;
  - P3.6.2, the comparison timestep rule.
- User decisions (2026-11-09):
  - four groups;
  - all three profiles outside;
  - the inside layer is off in outside runs;
  - the outside map and D signature are their own task at the reference spacing s = 2.5 plus one wider spacing;
  - the annulus model is the support.

## Global Constraints

- **Branch.** Work only on branch `cc/near-surface-r2`. Never modify or commit on other branches. Do not push; the user pushes.
- **Commits.** Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Framing.**
  - The conventional model is a reference, not the truth. No fast-diffusion or Brownstein–Tarr assumptions.
  - The annulus is an **approximation**, used only within the validated range (Task 3).
- **Fixed ρ.**
  - ρ is per unit surface area (D3.1). The outer layer is ΔR₂(d) = (ρ/h)·g(d/h)/κ_out with κ_out = 1 + c̄/R (package `kappa`).
  - The exponential's outside support is the cutoff, by default the smallest surface gap s − 2R (D3.5).
  - The conventional model uses θ = ρ√(π/D) in MC (P3.1.6) and Robin ρ in the annulus.
  - Bulk R₂ = 1/2000 ms⁻¹ wherever T₂ is reported.
- **Geometry.** R = 1 µm for every cylinder, square arrays with centre spacing s, gap s − 2R. Intra volume fraction f = πR²/s².
- **Reference values:** ρ_ref = 1/140 µm/ms, D = 3 µm²/ms, same-centroid matching c̄ = 0.1 µm (step h = 0.2, linear h = 0.3, exponential h by bisection with support = gap).
- **Noise and apparent T₂** as in 3.7a:
  - T₂_app = −TE/ln S(TE);
  - σ_T₂/T₂ = √(1 + r²)/(SNR·r·|ln r|), r = S(TE);
  - single-curve threshold 0.01 (2σ at SNR 200).
- **Comparison timestep (memory: comparison-timestep-rule).** τ = min(0.135·c̄_min², 0.005/ΔR₂(0)_max) over all layer models of a condition, at the D = 3 value for every D (ruling 3.7a). c̄_min < 0.05 µm is extrapolated and must be flagged.
- **Seeds (lesson, 3.7a Group B review I1).** Every condition gets its own seed block, seeds = 1000·k + (1:10) for condition index k, so different conditions are statistically independent. Within a condition all models share the block (paired). Exception: the single inside run of Task 4 is shared by every density's summed signal, and the record says so.
- **One threshold convention (lesson, 3.7a final review C1).** Every threshold is `threshold_status` (the log-log interpolated first crossing) on the quantity being compared, never "first grid point", including any number taken from earlier records.
- **A ceiling flag is not separation (lesson, 3.7a final review I1).** A fit flagged as beyond the conventional model's fastest decay reports no ρ_c and no residual. Whether such a curve is distinguishable is judged from the edge residual (an upper bound), never from the flag.
- **Statistics wording (lesson, I2).** Report what the data exclude ("no bias at the ±x % level; up to 2σ not excluded"), name the pooling method, and treat shared-seed conditions as correlated.
- **Commit scripts before running them** (provenance). Long runs go in detached `nohup` chains with done-markers and automatic resume (up to 3 attempts). Verify any killed process with `ps`; use `kill -9` for Julia smoke runs. Do not restart a waiter that hit the 2 h tool limit; check manually.
- **Cost (user, 2026-10-25: avoid too-thin layers).**
  - h ≥ 0.1 µm in MC.
  - At the start of every group, report the smoke-projected runtime.
  - If a group's projection exceeds 1.5× the estimate given here, stop and ask the user before launching.
- **Records.** Design / Method / Result / Conclusion with formulas and figures. Signs and orderings are recorded as found. Every number in the record is checked against the committed CSVs before commit.
- **Commands.**
  - Research tests: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
  - Scripts: `julia --project=research/baseline -t 8 research/phase3/<script>.jl`
- **Group ends.** A fresh subagent reviewer (most capable model) after each group (user rule). Fix Critical/Important findings with a failing test first where code is involved, ledger the minors.

## Review Focus

1. **Layer support vs the gap and the annulus width.**
   - The package rejects an outer support larger than the gap. The annulus is narrower than the gap (R_out − R = 0.41 µm at s = 2.5, where the gap is 0.5).
   - The annulus must truncate the layer at R_out and report the fraction of ∫ΔR₂ dA lost. It must never silently renormalise.
   - Test in Task 1.
2. **Exponential centroid outside.** h is found by bisection with support = gap. A target c̄ ≥ gap/2 is unreachable and must raise an error, not return a huge h (3.7a Group A minor M4).
   - Test in Task 2.
3. **Summed signal.**
   - Inside and outside runs must use identical readout times; this is asserted.
   - Their SEMs add in quadrature (independent runs). The paired ΔS is formed per compartment, then combined.
   - Test in Task 2.
4. **The annulus ceiling.** The conventional model's fastest decay outside is the Dirichlet-at-R, Neumann-at-R_out eigenvalue. It depends on s, not only R.
   - It must come from the exact Bessel root, with the large-ρ operator limit as a check.
   - Test in Task 1.
5. **Resume safety for the long runs (3.7a Group B minor M1).**
   - curves.csv is written before conditions.csv, each via a temporary file and rename.
   - Resume only labels present in both files; check the header before the first new condition.
   - Test in Task 2.

---

## Group A: tools, annulus model and its validation

### Task 1: Annulus continuum (support model outside)

**Files:**
- Create: `research/phase3/continuum_annulus.jl`
- Test: `research/phase3/test/test_continuum_annulus.jl`
- Modify: `research/phase3/test/runtests.jl` (include `continuum_annulus.jl`)

**Model.** For R ≤ r ≤ R_out:

  ∂m/∂t = D (1/r) ∂_r(r ∂_r m) − [R2_bulk + ΔR₂(r)] m,  uniform start m = 1,  ∂_r m = 0 at R_out;

- **layer:** ∂_r m = 0 at R, ΔR₂(r) = (ρ/h)·g((r − R)/h)/κ_out for 0 ≤ r − R ≤ L, where L = h (step, linear) or the cutoff (exponential), κ_out = 1 + c̄/R, and c̄ is the centroid with that support;
- **conventional:** D ∂_r m = ρ m at r = R (Robin, sink on the inner face);
- finite volumes on N rings of width (R_out − R)/N, exactly as `continuum_radial.jl`:
  - exact ring integrals of g;
  - face-centred Robin conductance 2πRρ/(1 + ρ dr/(2D));
  - a symmetric tridiagonal operator.
- **Truncation.** If R + L > R_out, the layer is cut at R_out. `annulus_layer_kept(shape, h, L, R, R_out)` = ∫_R^{min(R+L, R_out)} ΔR₂ r dr / ∫_R^{R+L} ΔR₂ r dr is reported, and the profile is never renormalised.

**Interfaces — Produces:**
- `krogh_radius(s) = s/√π`
- `outer_centroid(shape, h, L)`: step h/2, linear h/3, exponential h(1 − (c + 1)e^(−c))/(1 − e^(−c)) with c = L/h
- `outer_ring_integral(shape, h, R, L, lo, hi)`: ∫_lo^hi g((r − R)/h) r dr exactly (antiderivatives in r):
  - step r²/2;
  - linear 2[(1 + R/h) r²/2 − r³/(3h)];
  - exponential −h e^(−(r − R)/h)(r + h)/(1 − e^(−c));
  - clipped to [R, R + L].
- `annulus_layer_kept(shape, h, L, R, R_out)::Float64`
- `annulus_operator(; R, R_out, D, layer=nothing, surface=0.0, N=2000)`: returns (M, A), where `layer = (shape, ρ, h, L)`
- `annulus_signal(; R, R_out, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=2000)::Vector{Float64}`
- `annulus_rates(; kw...)`: ascending eigenvalues
- `annulus_signal_fun(; R, R_out, D, layer=nothing, surface=0.0, bulk=0.0, N=2000)`: a closure t → S(t) from one eigendecomposition
- `exact_rate_annulus_robin(R, R_out, D, ρ)`: lowest k² D of the Robin problem by bisection. With m = J₀(kr) + B Y₀(kr), B = −J₁(kR_out)/Y₁(kR_out), solve −D k[J₁(kR) + B Y₁(kR)] = ρ[J₀(kR) + B Y₀(kR)].
- `annulus_ceiling(R, R_out, D)`: the Dirichlet limit, lowest root of J₀(kR)Y₁(kR_out) − Y₀(kR)J₁(kR_out) = 0, as D k².

- [ ] **Step 1: Failing tests** (`test_continuum_annulus.jl`):
  - no sink, no bulk: `annulus_signal` ≡ 1 within 1e-12 at t ∈ {0, 10, 100};
  - `outer_ring_integral` over the full support equals h R κ_out with κ_out = 1 + `outer_centroid`/R within 1e-12, and matches a 10⁶-point midpoint quadrature within 1e-8, for step h = 0.2, linear h = 0.3, exponential (h = 0.2, L = 0.5);
  - **uniform limit:** step with h = R_out − R (s = 2.5) gives S(t) = exp(−ρ t/(hκ_out)) within 1e-9 at t = 10, 50 (ρ = 0.02);
  - **mean rate at t = 0:** Σ_k c_k λ_k (layer, no bulk) = ρ·2πR·kept/(π(R_out² − R²)) within 1e-6 relative. This holds both untruncated (step h = 0.2 at s = 2.5) and truncated (exponential h = 0.2, L = 0.5, at s = 2.5, kept < 1, which is asserted);
  - **Robin:** `annulus_rates(surface=ρ)[1]` against `exact_rate_annulus_robin` within 1e-6 relative, for ρ ∈ {0.007, 0.3, 30} at s ∈ {2.5, 4.0};
  - **ceiling:** `annulus_ceiling` against `annulus_rates(surface=1e9)[1]` within 1e-4 relative, and greater than the Robin rate at ρ = 30;
  - **grid:** step h = 0.1 at s = 2.5, ρ = 0.3, |S(N = 2000) − S(N = 4000)| ≤ 1e-6 on t ∈ 0:1:10.
- [ ] **Step 2:** Run the research tests. Expected: FAIL (`annulus_signal` not defined).
- [ ] **Step 3:** Implement `continuum_annulus.jl` with the header formulas, written independently of the package (like `continuum_radial.jl`). Use SpecialFunctions `besselj0/1`, `bessely0/1`. Bisection brackets: (0, k_D) for Robin, where k_D is the ceiling root; scan k from 1e-6 in steps of 1e-3·π/(R_out − R) for the first sign change for the ceiling.
- [ ] **Step 4:** Run the research tests. Expected: PASS (previous 204 + the new ones).
- [ ] **Step 5:** Commit `continuum_annulus.jl`, the test and `runtests.jl`.

### Task 2: Outside paired runner, summed-signal helpers, driver with safe resume

**Files:**
- Create: `research/phase3/p3_7b_mc.jl` (library) and `research/phase3/p3_7b_run.jl` (driver)
- Test: `research/phase3/test/test_p3_7b_mc.jl`
- Modify: `research/phase3/test/runtests.jl` (include `continuum_annulus.jl`, `p3_7b_mc.jl`)

**Interfaces — Consumes:** Task 1 (`krogh_radius`, `annulus_signal`, `outer_centroid`), `run_cylinders` (harness), `run_condition` (p3_7_mc.jl, for inside runs), `cbar`/`h_for_cbar` (p3_7_shared.jl, inside).

**Interfaces — Produces:**
- `h_for_cbar_out(shape, cb, gap)`:
  - step 2cb, linear 3cb;
  - exponential by bisection on `outer_centroid(…, L = gap)`;
  - errors if cb ≥ gap/2 (exponential unreachable) or if the support (2cb or 3cb) exceeds the gap.
- `run_condition_outside(; R=1.0, s, D, rho, R2_bulk, models, times, tau, nspins=20_000, seeds)`:
  - `models` holds (name, shape, h) with shape "conventional" (θ = ρ√(π/D)) or a profile on the **outer side only**: `Cylinders(radius=R, repeats=[s, s], layer_rho_outer=ρ, layer_h_outer=h, layer_shape=shape)`;
  - `region=:outside`, `sx = sy = s`, bounding box as in P3.4.8e;
  - returns the same NamedTuple fields as `run_condition`, with `S_ann` / `dS_ann` (annulus at R_out = `krogh_radius(s)`, exponential L = gap) in place of `S_cont` / `dS_cont`, plus `kept` per model.
- `sum_compartments(f, inside, outside, name)`:
  - asserts `inside.times == outside.times`;
  - returns S = f S_in + (1 − f) S_out;
  - SEM = √(f² SEM_in² + (1 − f)² SEM_out²);
  - ΔS = f ΔS_in + (1 − f) ΔS_out with the same quadrature for its SEM.
- `write_atomic(path, rows)`: `write_csv` to `path * ".tmp"`, then `mv(…; force=true)`.
- `read_done_both(cfile, kfile)`: labels present in both files. It errors on a header that differs from the expected one, and reports partial labels (in conditions but not curves) as not done.
- Driver `p3_7b_run.jl`:
  - `TASK` = validate | inside | density | map | ddep, writing `results/phase3/p3_7b_<task>/` (conditions.csv, curves.csv, summary.json);
  - writes curves before conditions, both via `write_atomic`;
  - seeds 1000·k + (1:10);
  - `NSPINS` and `RESUME` environment variables.

- [ ] **Step 1: Failing tests** (`test_p3_7b_mc.jl`):
  - `h_for_cbar_out("exponential", 0.1, 0.5)` reproduces c̄ = 0.1 within 1e-10; `h_for_cbar_out("exponential", 0.25, 0.5)` throws; `h_for_cbar_out("linear", 0.1, 0.25)` throws (support 0.3 > 0.25);
  - `run_condition_outside` on a tiny condition (s = 2.5, ρ = 0.02, conventional + step h = 0.2, times [0, 10], τ = 1e-2, 500 spins, seeds 1:2) gives `identical_trajectories == true`, S(0) = 1, ΔS(0) = 0, and an annulus column equal to `annulus_signal`;
  - `sum_compartments` combines the means and the quadrature SEM correctly on hand-made inputs, and throws on different times;
  - `write_atomic` leaves no `.tmp` file;
  - `read_done_both` treats a label present only in conditions.csv as not done, and throws on a header mismatch.
- [ ] **Step 2:** Run. Expected: FAIL (functions not defined).
- [ ] **Step 3:** Implement `p3_7b_mc.jl`.
- [ ] **Step 4:** Run. Expected: PASS.
- [ ] **Step 5: Driver and condition lists** (Tasks 3–6 below). Smoke with `NSPINS=200`, one seed, for every TASK:
  - record τ, readout windows, `identical_trajectories` and wall time per condition;
  - extrapolate to 2·10⁴ spins × 10 seeds;
  - ledger the projected cost per task;
  - delete the smoke outputs.
- [ ] **Step 6:** Commit library, driver, test and `runtests.jl` before any full run.

### Task 3: Annulus validation against paired Monte Carlo (decides how the annulus is used)

**MC conditions** (`TASK=validate`):
- s ∈ {2.5, 4.0} (gap 0.5 and 2.0; f = 0.503 and 0.196), D = 3, no bulk, ρ ∈ {0.02, 0.2};
- conventional + step / linear / exponential at the same centroid c̄ = 0.1 (step h = 0.2, linear h = 0.3, exponential h by `h_for_cbar_out` with L = gap);
- τ = 1.35e-3; readouts 0:5:50 for ρ = 0.02, and 41 readouts to ln 100/λ₁,ann for ρ = 0.2.
- Estimated cost ≈ 3–4 h.

**Quantities** (script `research/phase3/p3_7b_validate_analysis.jl` → `results/phase3/p3_7b_validate_analysis/`):
- per model and readout: e_S = (S_mc − S_ann)/S_ann with z = (S_mc − S_ann)/SEM_S;
- per layer model: max_t |ΔS| ratio ΔS_mc/ΔS_ann at the readout of max |ΔS_ann|, with σ = SEM_ΔS/|ΔS_ann|;
- the conventional rate: λ_mc from a log-linear fit over the readouts against λ₁,ann;
- `kept` per model.

**Decision rule** (stated before the run; ledger the outcome):
- **Full support.** The annulus is accepted as support where, at both spacings and both ρ, max |e_S| ≤ 1 % for all models and every ΔS ratio lies in [0.8, 1.2] within 2σ. It then provides the dense curves, TE*, the readout windows, the conventional fit family and the ceiling.
- **Partial support.** If only the conventional model passes, the annulus provides TE*, the readout windows and the ceiling (an approximate bound), but no layer dense curves and no best-fit family. The outside best fit is then reported only for MC layer curves against the annulus conventional family, with the conventional model's own MC − annulus difference quoted as the floor.
- **Failure.** If the conventional model fails too, the annulus is used only to choose TE and readout windows, and every outside result is MC-only.

- [ ] **Step 1:** Commit the analysis script, then run `TASK=validate` detached.
- [ ] **Step 2:** Run the analysis, write the decision into the ledger, and commit the results.
- [ ] **Group A review** (fresh subagent, most capable model): Tasks 1–3, with Review Focus 1, 2, 3, 4, 5. Fix and ledger.

---

## Group B: packing density (Fig. 9B analogue), alone

### Task 4: Apparent T₂ against packing density, intra, extra and summed (P3.7.1 arrays, P3.7.2)

**Conditions:**
- Spacings s ∈ {2.3, 2.5, 3.0, 4.0}: f = 0.594, 0.503, 0.349, 0.196; gap 0.3, 0.5, 1.0, 2.0.
- ρ = 1/140, D = 3, bulk 1/2000, same centroid c̄ = 0.1 on both sides:
  - inside: step 0.2, linear 0.3, exponential by `h_for_cbar` with support R;
  - outside: step 0.2, linear 0.3, exponential by `h_for_cbar_out` with L = gap.
- At s = 2.3 the linear support 0.3 equals the gap; this is allowed within the package tolerance and is checked in the smoke.
- **TE per density.** TE_sum*(s) is where the summed conventional signal f·S_in,cont(t) + (1 − f)·S_out,ann(t) = e⁻¹ (with bulk). Each density uses readouts 0:TE_sum*/10:1.5·TE_sum*.
- **Inside (`TASK=inside`):** one inside run at R = 1 (conventional + the three inside profiles, τ = 1.35e-3), with readouts = the union of the four density grids. It is shared by every density's summed signal.
- **Outside (`TASK=density`):** one outside run per s, conventional + three outer profiles, τ = 1.35e-3, readouts = that density's grid.
- **Estimated cost** ≈ 8–10 h: the inside run ≈ 2 h, outside 1.5–3 h per density. The smoke projection from Task 2 decides; above 15 h, ask the user.

**Quantities** (script `research/phase3/p3_7b_density_analysis.jl` → `results/phase3/p3_7b_density_analysis/`), per density, for intra, extra and summed, at TE_sum* (exact readout):
- T₂_app per model;
- the paired ΔT₂/T₂ per compartment and summed (σ via `sum_compartments`);
- 2σ two-point noise at SNR 100/200, with r = S_conv(TE_sum*) of the compartment in question;
- MC against support: intra z against `radial_signal`; extra z against the annulus only if Task 3 accepted it, otherwise not shown;
- the figure, Fig. 9B analogue:
  - T₂_app against f (conventional + profiles; intra, extra, summed);
  - |ΔT₂/T₂| on a log axis against the noise lines;
  - annulus/continuum dense curves over s ∈ [2.3, 6] if accepted.

- [ ] **Step 1:** Commit the analysis script. Run `TASK=inside` then `TASK=density` in one detached chain. Test the analysis on the first finished density and fix figure issues.
- [ ] **Step 2:** After the chain: full analysis, check the numbers, commit the results.
- [ ] **Group B review** (fresh subagent): Task 4, with Review Focus 3 and the statistics wording.

---

## Group C: outside map and D-dependence signature

### Task 5: Outside dimensionless map (P3.7.3 outside)

**MC conditions** (`TASK=map`):
- s ∈ {2.5, 4.0}; R = 1, D = 3, no bulk;
- ρR/D ∈ {0.1, 1, 10} × h/R ∈ {0.1, 0.25, 0.5} × three outer profiles, plus the conventional model, per (s, ρ).
- At s = 2.5, h = 0.5 equals the gap (allowed); the exponential support = gap.
- τ shared per condition (rule and rate limits). 41 readouts to t_end = ln 100/λ₁,conv, with λ₁ from the annulus.
- Estimated cost ≈ 8–10 h (6 conditions).

**Quantities** (script `research/phase3/p3_7b_map_analysis.jl`, structure as `p3_7_map_analysis.jl`):
- **fixed ρ:** max_t |ΔS| (paired MC, SEM at the arg-max readout);
- **best fit** (only under full support): `best_fit_conventional(S_layer_MC, times, ρc -> annulus conventional S)`. It reports ρ_c/ρ, the residual against the noise floor, and the flags `at_edge` and `ceiling` (layer λ₁,ann > `annulus_ceiling(R, R_out, D)`). Flagged fits get no ρ_c or residual; their edge residuals are reported separately as upper bounds (lesson I1);
- **thresholds** in ρR/D with `threshold_status`, for MC and, under full support, for the annulus dense map (ρR/D 0.003–30 × h/R 0.1–0.5 at each spacing);
- **checks:** conventional θ√τ and max|z_S| over readouts; annulus grid convergence at the thinnest map point (h/R = 0.1) at ρR/D = 1 and 30; the family check at the most demanding **valid** fit (lesson I4).

### Task 6: Outside D-dependence signature (P3.7.3 outside)

**MC conditions** (`TASK=ddep`):
- s = 2.5, bulk 1/2000, ρ ∈ {0.1, 1} × D ∈ {0.5, 1, 2, 3};
- conventional + three outer profiles at h/R ∈ {0.1, 0.25, 0.5};
- TE per ρ = the annulus conventional T₂ at D = 3 (exact readout 11 of 0:TE/10:1.5 TE).
- Estimated cost ≈ 2–3 h.

**Quantities** (script `research/phase3/p3_7b_ddep_analysis.jl`, structure as `p3_7_ddep_analysis.jl`):
- q = T₂(0.5)/T₂(3) − 1 per model;
- q_layer − q_conv against the two-measurement band 2·√(n(S₀.₅)² + n(S₃)²) at SNR 200;
- z_q = NaN when σ_q = 0;
- the band-edge cases are named explicitly (lesson I3), and the same-ρ assumption is stated.

- [ ] **Step 1:** Commit both analysis scripts. Run `TASK=map` then `TASK=ddep` in one detached chain. Test each analysis on its first condition.
- [ ] **Step 2:** After the chain: full analyses, check the numbers, commit the results.
- [ ] **Group C review** (fresh subagent): Tasks 5–6, with Review Focus 4 and the ceiling wording.

---

## Group D: record, tracker, final review

### Task 7: Record and tracker

- [ ] Write the record section "Comparison with the conventional model in arrays and outside cylinders (P3.7b)" in `research/phase3_record.md`, with Design / Method / Result / Conclusion and the figures:
  - the annulus model, its validation and the decision taken;
  - the Fig. 9B analogue (intra, extra, summed);
  - the outside map at two spacings, with thresholds, best fit and ceiling;
  - the outside D signature;
  - MC against support for every point;
  - **inside vs outside vs walls**, every threshold on one convention: `threshold_status`, re-applied to the 3.7a and Phase 2 maps by `p3_7_walls_comparison.jl`, extended to the outside map.
- [ ] Check every number against the committed CSVs before commit (Global Constraints).
- [ ] Tracker: P3.7.1–P3.7.4 → ☑ with one-line results (P3.7.2 newly), and the Phase 3 counts row.
- [ ] Research suite PASS; commit.

## After the last task

Final whole-branch review (3.7a + 3.7b) by a fresh subagent on the most capable model (user rule). Fix Critical and Important findings and ledger the minors. Report every ruling and deferred minor to the user.
