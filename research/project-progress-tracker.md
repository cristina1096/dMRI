# Project Progress Tracker — Distance-Dependent R₂(d) in MCMRSimulator

Working tracker for implementation. Tick a box when its pass criterion is met, not when the code merely runs.

**Status key:** ☐ not started · ◐ in progress · ☑ done · ✖ blocked (add a note)

**Calendar:** week 1 = 21–25 Sep 2026. Dates are targets; move them if the assessment calendar changes.

---

## Progress summary

| Phase | Weeks | Scope | Tasks | Done | Status |
|---|---|---|---|---|---|
| **Phase 1** | 1–2 (21 Sep – 2 Oct) | Verified baseline, unmodified code | 22 | 19 | ☑ exit criterion met (open: P1.1.3, P1.1.5, optional P1.3.5) |
| **Phase 2** | 3–11 (5 Oct – 4 Dec) | Implement and verify R₂(d) on wall geometry, no RF | 52 | 33 | ◐ |
| Phase 3 | 12–13 | Sensitivity studies (C1–C7), cylinder and sphere geometry | — | — | not yet planned |
| Phase 4 | 14–15 | Trade-off audit, write-up | — | — | not yet planned |

Update the "Done" column at the end of each week.

---

# PHASE 1 — Verified Baseline (Weeks 1–2)

**Goal.** Understand the existing code, reproduce published results with the **unmodified** simulator, and record every setting, so that any later change in signal can be attributed to the modification rather than to the setup.

**Exit criterion.** Rep1–Rep4 pass, the baseline record is complete, and the `baseline-v1` commit is tagged. Phase 2 does not start until this is met.

## 1.1 Environment

| ID | ☐ | Task | Pass criterion | Output |
|---|---|---|---|---|
| P1.1.1 | ☑ | Install Julia; install MCMRSimulator.jl **pinned to the paper's version (v1.0)** and MRIBuilder.jl | Package loads, bundled tutorial runs | Julia 1.12.7, MRIBuilder 0.4.1. **Decision 2026-10-02:** keep MCMRSimulator 1.1.0 (this repo, the code Phase 2 modifies) instead of v1.0.0; Rep1–Rep4 agree with the v1.0.0 notebooks |
| P1.1.2 | ☑ | Create git repository; commit `Project.toml` and `Manifest.toml` | Environment reproducible from a clean clone | branch `cc/near-surface-r2`; `research/baseline/Project.toml` + `Manifest.toml` committed (the root `Manifest.toml` is git-ignored by the project's rule) |
| P1.1.3 | ☐ | Set up Revise.jl and a persistent Julia session | Code edits reload without restart | — |
| P1.1.4 | ☑ | Clone the figure notebooks (`git.fmrib.ox.ac.uk/ndcn0236/mcmr_paper_figures`) | Notebooks open; note which are usable | used for Rep2 (`Figure_8_9/gradient_spin_echo.ipynb`), Rep3 (`Figure_S1/turtuosity.ipynb`), Rep4 (`Figure_S2/toy_model.ipynb`); notebook @ `c7b86f9`; two figure-axis issues recorded in `baseline_record.md` |
| P1.1.5 | ☐ | Apply for cluster access | Application submitted (lead time for Phase 3) | confirmation |

## 1.2 Code study

Trace one simulation end to end. Target output is a written map, not general familiarity.

| ID | ☐ | Task | What to locate | Output |
|---|---|---|---|---|
| P1.2.1 | ☑ | **Code structure** | Module layout; how geometry, sequence and simulation modules depend on each other; how MRIBuilder sequences are passed in | code map §1 — see `research/phase1_code_log.md` |
| P1.2.2 | ☑ | **Simulation workflow** | Main loop: setup → step → readout; where τ_max is computed from the constraint list (Supp. S1) | code map §2 |
| P1.2.3 | ☑ | **Geometry representation** | Storage of walls, cylinders, spheres, annuli, meshes; collision grid; repeating geometry; how surface vs volume parameters attach (Table 1) | code map §3 |
| P1.2.4 | ☑ | **Spin trajectories** | Gaussian step draw; collision detection; reflection and re-test; handling of perfectly permeable surfaces (§2.10) | code map §4 |
| P1.2.5 | ☑ | **Magnetisation update** | Bloch update between events; update-to-collision-time (§2.4); RF chunking (§2.5); **exact lines where `θ_relax` is applied** | code map §5 |
| P1.2.6 | ☑ | **Single-spin trace** — one cylinder, spin echo, ~100 spins; log one spin's position and magnetisation every step through at least one collision | Trace matches the reading of the code | trace log |
| P1.2.7 | ☑ | Mark insertion points for R₂(d): parameter definition, distance evaluation, relaxation application | Each insertion point named with file and function | code map §6 |

## 1.3 Reproduce baseline results

Use notebook values where available; values read off figures carry a few percent uncertainty.

| ID | ☐ | Result | Reference | Pass criterion | Output |
|---|---|---|---|---|---|
| P1.3.1 | ☑ | **Rep1** — surface relaxation, walls (w = 2 µm) and single cylinder: R₂ = θ_relax·√(D/π)·(S/V) | Exact Brownstein–Tarr rate; Eq. 17 recorded | Within 2 SEM of the exact rate (decided 2026-10-02; originally "within 2σ of Eq. 17") | fig + numbers — PASS vs exact: walls +1.7 SEM, cylinder −0.5 SEM. Eq. 17 is 0.08–0.11% above the exact rate (resolved: −5.9 / −10 SEM); see `baseline_record.md` |
| P1.3.2 | ☑ | **Rep2** — compartment T₂: ~80 ms intra-, ~60 ms extra-axonal, ~45 ms myelin | §3.3, Fig. 8 | Within 2σ, or 5% if read off figure | fig + numbers — within 5% of digitised notebook Fig. 8 at all TEs |
| P1.3.3 | ☑ | **Rep3** — timestep plateau, ADC vs τ, cylinder packings | Supp. Fig. S1 | Plateau onset near τ ≈ 0.01 ms | fig — onset 6e-3–2.5e-2 ms; notebook S1 x-axis mislabelled |
| P1.3.4 | ☑ | **Rep4** — surface-relaxation correction schemes, walls (w = 1 µm) | Supp. Fig. S2C | Logarithmic scheme flat to τ ≈ 10⁻³ ms; confirm θ_relax used | fig + θ_relax value — θ_relax = 1/T (θ = 10 for the τ≈1e-3 row) |
| P1.3.5 | ☐ | *(optional)* **Rep5** — diffusion vs Mitra / van Gelderen | Fig. 5 | Qualitative agreement | fig |

**If a result fails:** find the cause before continuing. It is usually a parameter mismatch, and the cause is worth recording.

## 1.4 Record

| ID | ☐ | Task | Pass criterion | Output |
|---|---|---|---|---|
| P1.4.1 | ☑ | Baseline record for every run (software versions + commit hashes, hardware, geometry with exact S/V, physics, sequence, τ_max and binding constraint, N_spins, seeds, signal, apparent T₂, σ, runtime, reference, pass/fail) | Every run has a complete entry | `baseline_record.md` — `research/baseline_record.md` |
| P1.4.2 | ☑ | Noise floor σ for each baseline configuration, ≥ 10 seeds | σ reported per configuration | table |
| P1.4.3 | ☑ | Save raw outputs with generating script and seed | Any result regenerable from the repo | `research/results/baseline/` + scripts in `research/baseline/`, committed |
| P1.4.4 | ☑ | Tag commit `baseline-v1` | Tag exists | git tag `baseline-v1` (2026-09-30). Rep1 cylinder criterion (Eq. 17 vs exact) left open: not used by Phase 2 |
| P1.4.5 | ☑ | Wall-geometry baselines for Phase 2 (B0–B5): test suite, null run, bulk R₂, density/occupancy, θ_relax reference curves, runtime vs τ | Each recorded with settings, σ and pass/fail or reference role | `baseline_record.md` B-series; B6 deferred |

## Phase 1 deliverables

- ☑ Code map with insertion points
- ☑ Rep1–Rep4 reproduced with pass/fail
- ☑ Baseline record, `baseline-v1` tag, pinned environment (MCMRSimulator 1.1.0 of this repo by decision, see P1.1.1)
- ☑ Noise floor σ per baseline configuration

---

# PHASE 2 — Implement and Verify R₂(d) on Wall Geometry (Weeks 3–11)

**Goal.** Implement the distance-dependent model and verify it on the simplest geometry available: planar walls, no RF pulses. Start with the step profile to verify each mechanical component separately, then add the remaining profiles, compare them under identical conditions, and establish timestep and integration-tolerance requirements.

**Why walls and no RF.** A planar wall makes every quantity exact: the distance is a single coordinate, the offset surface is another plane, the layer volume fraction is *f* = λ/*w* with no curvature, and along a straight step the distance changes **linearly in time**, so the relaxation integral has a closed form for every profile. Removing RF removes the second chunking driver, so any error is attributable to the R₂(d) code alone. Cylinder, sphere and finite-RF checks (V5 curved, V13) move to Phase 3.

**Exit criterion.** Section 2.3 passes; all four profiles pass Section 2.4; *c*_h is calibrated and ε_R₂ fixed (Section 2.6). Tag `rd-walls-v1`.

## Reference configuration

Use this for every Phase 2 test unless a task says otherwise, so results are comparable.

| Parameter | Value | Note |
|---|---|---|
| Geometry | Repeating planar walls, spacing *w* = 2 µm (`Walls(repeats=2)`) | Layer on **both faces** of every gap (symmetric setting `layer_rho`, `layer_h`), matching the baseline θ_relax reference (B4). For *h* > *w*/2 the two layers overlap and their contributions add (proposal Eq. 10). Side-specific settings are used only in P2.3.12 |
| D | 3 µm²·ms⁻¹ | |
| R₂^bulk | 0 | Isolates the layer term. Bulk is switched on only in the additivity test (P2.3.11) and in one V4 case |
| R₁ | off | |
| θ_relax, θ_perm, susceptibility, MT | 0 / off | Recorded explicitly |
| Sequence | No RF; spins start fully transverse; free decay | No refocusing needed since off-resonance is zero |
| Readout times | 0, 5, 10, …, 50 ms | Same as the B4 reference curves |
| Surface excess rate ΔR₂(0) | 0.1 ms⁻¹ | = (ρ/*h*)·*g*(0) (proposal Eq. 9). Must be ≥ 0.046 ms⁻¹ so every B4 θ is reachable (at *h* = *w*/2 the step profile gives *S*(50) = e^(−ΔR₂(0)·50)). Physiological values are too weak to test |
| Timestep | τ = 1e-2 ms for comparisons with B4; the converged τ from P2.6 once known | B4 reference is converged at 1e-2 ms; the default 0.04 ms biases the baseline by up to 0.055% |
| N_spins | 10⁵ | Adjust after P2.1.4 |
| Seeds | 1–10 | Same seeds as B1–B4, so runs are paired with the baseline |

**Exact expectations** (hold for any diffusion regime; compute before running):

| Case | Expected signal | Why it is exact |
|---|---|---|
| *h* = 0 (layer off) | identical to B1 (bulk off) or B2 (bulk on), bit for bit | no layer code runs |
| Step profile, *h* = *w*/2, both faces | *S*(*t*) = exp(−ΔR₂(0)·*t*) per spin, to roundoff | every point of the gap lies in exactly one layer, so every path sees ΔR₂(0) all the time |
| Any profile, R₂_bulk on | *S*_bulk+layer(*t*) = exp(−R₂_bulk·*t*)·*S*_layer(*t*), same seed, to roundoff | bulk is uniform and does not change trajectories (B2) |

Everything else is **measured**, not predicted, and compared with the B4 reference curves (baseline θ_relax model), which are a reference rather than the true answer.

## 2.1 Design and test harness (week 3)

**Implementation plan for P2.1.3–P2.1.5 and P2.2.9:** [docs/superpowers/plans/2026-10-03-phase2-harness-and-trace.md](../docs/superpowers/plans/2026-10-03-phase2-harness-and-trace.md).

| ID | ☐ | Task | Pass criterion | Output |
|---|---|---|---|---|
| P2.1.1 | ☑ | Define the parameter interface (proposal §3.2 III), **side-specific**: on `Walls`, `layer_rho` (ρ, µm/ms) and `layer_h` (*h*, µm) apply to both sides; `layer_rho_positive`, `layer_h_positive`, `layer_rho_negative`, `layer_h_negative` override one side (positive = side where the wall's local coordinate is larger, +x for `rotation=:x`). ΔR₂(0) = (ρ/*h*)·*g*(0) derived internally, never user-settable. Shape *g*: step only for now; a `layer_shape` field is added with the other profiles (P2.4) | Interface implemented and documented in field descriptions; both-sides, one-side and two-different-sides settings tested | plan Task 2 |
| P2.1.2 | ☑ | Construction-time validation: reject ρ < 0, *h* < 0, and *h* = 0 with ρ > 0; reject *h* larger than the wall spacing (the layer would pass through the neighbouring wall); overlapping layers (*w*/2 < *h* ≤ *w*) are allowed and summed (Eq. 10); warn that the layer is not applied while a spin is stuck (surface density > 0). *Not yet:* warning when ΔR₂(0) exceeds the rigid-lattice ceiling (needs a physical value) | Invalid inputs rejected with clear messages | plan Task 2 |
| P2.1.3 | ☑ | Test harness: runs a configuration over seeds, returns mean, σ, runtime, layer-visit count; computes ρ from a target ΔR₂(0) and *h* (fixed-ΔR₂(0) mode used by V4/V5); loads the B4 reference curves | Harness reproduces B1, B2 and one B4 curve exactly (same seeds) | `harness.jl` — `research/phase2/harness.jl`, tests in `research/phase2/test/` |
| P2.1.4 | ☑ | Noise floor of the model at ΔR₂(0) = 0.1 ms⁻¹, *h* = 0.1 and 0.5 µm (at *h* = 0 σ is exactly 0, B1) | σ of *S*(*t*) recorded per readout; N_spins adjusted so neighbouring *h* in the P2.3.4 sweep differ by > 10σ | σ table — `research/results/phase2/p2_1_4_noise_floor/` |
| P2.1.5 | ☑ | Expectations script: the exact cases above, plus loading and interpolating the B4 reference curves (θ ↔ *S*(*t*)) | Exact tables and reference interpolation available before any run | `expected.jl` — `research/phase2/expected.jl` |

## 2.2 Step profile — component implementation (weeks 3–5)

Implement and test each component separately, in this order. Each has its own unit test, so a failure points to one component.

**Implementation plan:** [docs/superpowers/plans/2026-10-01-near-surface-layer-walls.md](../docs/superpowers/plans/2026-10-01-near-surface-layer-walls.md) covers P2.1.1, P2.1.2 and P2.2.1–P2.2.8 step by step (code, tests, commands). Tick the rows here when the plan's task passes.

**Design note (step profile, no RF).** Along a straight free segment the distance to a wall changes linearly, so the entry and exit points of the layer and the time spent inside it have an exact closed form (`segment_overlap`). The existing `draw_step!` already splits each step at every reflection; within each straight piece no further splitting is needed. The layer factor exp(−∫ΔR₂ d*t*) multiplies the bulk factor from `relax!`, which is R₂_total = R₂_bulk + ΔR₂(*d*). Finite RF (where the order of sub-segments matters) is Phase 3.

| ID | ☐ | Component | What it does | Unit test | Pass criterion |
|---|---|---|---|---|---|
| P2.2.1 | ☑ | **Distance evaluation** | Signed distance from a point to each wall face carrying a layer, per side (*x* − *c* on the positive side, *c* − *x* on the negative side), over all repeated wall copies a segment can reach | Hand-computed segments near, between and far from walls; rotated and shifted walls; wall copies | Exact to floating-point precision (plan Task 3) |
| P2.2.2 | ☑ | **Crossing detection — entry** | Entry point of a straight segment into the layer 0 ≤ *d* ≤ *h* (`segment_overlap`) | Constructed segments with known crossings, including ones with **neither endpoint inside the layer**; points exactly at *d* = 0 and *d* = *h* | All crossings found, times exact (plan Task 1) |
| P2.2.3 | ☑ | **Crossing detection — exit** | Exit point, same function | As above | All crossings found (plan Task 1) |
| P2.2.4 | ☑ | **Crossing detection — with reflection** | Segment enters layer, reflects off wall, exits, within one step | Single spin moved 0.8 → −0.4 µm onto a wall at 0 with layer [0, 0.5]: ends at 0.4, time in layer 0.9 of 1.2 ms | *M*ₓᵧ = exp(−ΔR₂(0)·0.9 ms) to 10⁻¹² (plan Task 4) |
| P2.2.5 | ☑ | **Trajectory segmentation** | Existing `draw_step!` splits each step at every reflection; the layer is evaluated on each straight piece with the same start, end and duration as `relax!` | Durations of the pieces sum to τ; bulk factor bit-identical with and without layer (B2 comparison); bulk × layer = exact product per spin | Exact (plan Task 4) |
| P2.2.6 | ☑ | **Local relaxation evaluation** | ∫ΔR₂(*d*(*t*)) d*t* over each straight piece; for the step profile this is ΔR₂(0)·(time inside); overlapping layers add | Constructed segments vs hand calculation, one side, both sides, different sides, overlap (*h* > *w*/2) | Exact (plan Tasks 1, 3) |
| P2.2.7 | ☑ | **Transverse magnetisation update** | *M*ₓᵧ ← *M*ₓᵧ·exp(−∫ΔR₂ d*t*) after `relax!` for each piece; phase and *M*_z untouched; finite RF pulses rejected with an error | *h* = *w*/2, both sides: every spin decays at exactly ΔR₂(0); bulk on: per-spin product exp(−R₂_bulk·*t*)·*S*_layer | Exact to roundoff (plan Task 4) |
| P2.2.8 | ☑ | ***h* = 0 / ρ = 0 short-circuit** | No layer stored when every ρ is 0; `apply_layer!` returns before any work | Positions and magnetisation bit-identical to the unmodified code (bulk on); runtime per spin-step vs B1 (0.013 µs) | Bit-identical; < 10% slowdown (plan Tasks 4, 5); measured 0.0127 µs per spin-step with no layer (B1: 0.0130) |
| P2.2.9 | ☑ | End-to-end single-spin trace, as in P1.2.6, with the layer on | Logged decay matches manual integration of the logged path | Exact — `research/results/phase2/p2_2_9_layer_trace/` (max error 5e-8) |

## 2.3 Step profile — verification (weeks 5–6)

**Implementation plan for section 2.3 (P2.3.1–P2.3.12, tag `rd-step-v1`):** [docs/superpowers/plans/2026-10-03-phase2-step-verification.md](../docs/superpowers/plans/2026-10-03-phase2-step-verification.md). Depends on the harness plan above.

| ID | ☐ | Test | Method | Pass criterion | Output |
|---|---|---|---|---|---|
| P2.3.1 | ☑ | **V0b** baseline reproduction | *h* = 0 vs B1 (bulk off) and B2 (bulk on), same seeds | Bit-identical | table — `research/results/phase2/p2_3_exact_checks/` — bit-identical, bulk off and on |
| P2.3.2 | ☑ | **V9** equilibrium density | Relaxation off; histogram spin density vs *d* | Uniform within Poisson error (cf. RMS App. A for the failure mode) | histogram — `research/results/phase2/p2_3_2_density_occupancy/` — positions identical to no layer; χ²/dof 0.948 |
| P2.3.3 | ☑ | **Layer occupancy check** | Relaxation off; count spins inside the layer | Fraction within *h* of either face = 2*h*/*w* (*h* ≤ *w*/2) within Poisson error, as in B3. **Validates geometry independent of relaxation code** | table — `research/results/phase2/p2_3_2_density_occupancy/` — |z| ≤ 1.5 (gate PASS) |
| P2.3.4 | ☑ | **V5 attenuation vs *h*** (step profile, fixed surface excess rate ΔR₂(0), both faces) | (a) **Exact end point:** *h* = *w*/2, where every point of the gap lies in exactly one layer, so every spin decays at exactly ΔR₂(0) whatever its path. (b) **Characterisation:** sweep *h* = 0.05 – 1.0 µm; record *S*(*t*) vs *h* at all readouts, mean ± σ, seeds 1–10 | (a) *S*(*t*) = exp(−ΔR₂(0)·*t*) to the roundoff bound, per spin and ensemble. (b) Curve recorded with σ; **no functional form assumed** (not a linearity test) | table + fig — (a) `p2_3_exact_checks/` per spin ≤ 1.4e-13 (gate PASS); (b) `p2_3_4_h_sweep/` |
| P2.3.5 | ☑ | **Monotonicity** (proposal §3.2 VI) | Attenuation vs *h* from P2.3.4(b), fixed profile and fixed ΔR₂(0) | Attenuation strictly increasing with *h*: each step up in *h* increases 1 − *S*(50) by more than 2σ of the difference | fig — `p2_3_4_h_sweep/` — 13/13 steps |
| P2.3.6 | ☑ | **Decay shape** | ln *S*(*t*) over all readouts, for each *h* in P2.3.4(b) | Characterisation, no assumed form: record the curvature of ln *S*(*t*) vs *h* with σ, and compare with the curvature of the matched B4 curve | fig — `p2_3_4_h_sweep/` — curvature 0 within 1.3 SEM at every h |
| P2.3.7 | ☑ | **V11** intersection vs endpoint detection | Run both methods at τ = 1 µs | Endpoint method under-counts; bias quantified | table — `p2_3_7_endpoint/` — bias measured, sign changes with τ (−0.1 / −1.4 / +5.2 SEM at τ = 1e-3 / 1e-2 / 4e-2); "under-counts" not confirmed |
| P2.3.8 | ✖ | **V10** Route A cross-check — *deferred to future work (decision 2026-09-30)* | Same step profile built from perfectly permeable walls with per-compartment R₂ | Within 2σ of Route B | table |
| P2.3.9 | ☑ | **V4 fit *h*\*(θ) to the baseline** | For each θ in the B4 reference set (τ = 1e-2 ms curves): step profile, fixed ΔR₂(0) (≥ 0.046 ms⁻¹ so every θ is reachable, e.g. 0.1 ms⁻¹), same seeds 1–10 and 10⁵ spins as B4; find *h*\* minimising Σₜ[(*S*_model − *S*_θ)/σ]² over the readouts; bulk off, then one θ repeated with R₂_bulk = 1/80 | *h*\*(θ) found for every θ in the grid, with its uncertainty and the residual at every readout. **The baseline is a reference, not the truth:** residuals are reported as the difference between the two models, not as pass/fail | fig + table — `p2_3_9_fit_h/` — 12/12 θ fitted; whole curve matched within 4e-5; h*/θ ≈ 9.76 µm |
| P2.3.11 | ☑ | **Additivity** R₂_total = R₂_bulk + ΔR₂(*d*) | Same seeds, *h* = 0.2 µm, R₂_bulk = 0 and 1/80 ms⁻¹ | *S*_bulk+layer(*t*) = exp(−R₂_bulk·*t*)·*S*_layer(*t*) to the roundoff bound (as in B2/B4) | table — `p2_3_exact_checks/` — per spin 6.6e-14 |
| P2.3.12 | ☑ | **Side-specific layer** (characterisation) | (a) Mirror check: positive-only vs negative-only layer, same ρ and *h*. (b) One-sided and two-different-sides runs in the reference configuration, *S*(*t*) with σ | (a) Equal within noise (a plan Task 4 unit test checks this at small size). (b) Recorded; no baseline exists for asymmetric surfaces (MCMR's θ_relax is the same on both sides) | table + fig — `p2_3_12_sides/` — exact flipped-normal mirror test bit-identical (criterion changed from the statistical |z| ≤ 2, see `phase2_record.md`) |
| P2.3.10 | ☑ | **Tag** `rd-step-v1` | — | Tag exists | git tag — tag `rd-step-v1` (2026-10-03) |

**Gate:** do not start 2.4 until P2.3.3, P2.3.4(a) (exact *h* = *w*/2 end point) and P2.2.9 (single-spin trace) pass. Failure there means distance or segmentation is wrong, and every later profile would inherit it.

## 2.4 Additional profiles (weeks 6–8)

**Implementation plan for the linear profile (P2.4.1–P2.4.3, P2.4.7–P2.4.8, linear part):** [docs/superpowers/plans/2026-10-06-linear-profile.md](../docs/superpowers/plans/2026-10-06-linear-profile.md).

All profiles normalised so ∫₀^∞ *g*(*u*) d*u* = 1.

| Profile | *g*(*u*) | ∫*g* | Wall value *g*(0) | Note |
|---|---|---|---|---|
| Step | 1 for *u* < 1 | 1 | 1 | Done in 2.2 |
| Linear | 2(1 − *u*) for *u* < 1 | 1 | 2 | Special case *n* = 1 of polynomial |
| Polynomial, order *n* | (*n*+1)(1 − *u*)ⁿ for *u* < 1 | 1 | *n*+1 | Test *n* = 2, 3 |
| Exponential | *e*^(−*u*) | 1 | 1 | Infinite tail — see truncation note |
| *(optional)* Power law | (*m*−1)(1+*u*)^(−*m*) | 1 | *m*−1 | *m* = 4; needs *d*₀ |
| *(optional)* Half-Gaussian | √(2/π)·*e*^(−*u*²/2) | 1 | √(2/π) | |

**Truncation note.** With walls every *w*, the exponential tail is cut at *d* = *w*, losing a fraction *e*^(−*w*/λ) of ρ. Keep λ ≤ *w*/5 (< 1% loss) for exponential tests, or record the loss.

**Analytic segment integral.** On a planar wall, *d*(*t*) = *d*₀ + *v*ₙ*t* along a straight sub-segment, so ∫*g*(*d*(*t*)/λ) d*t* has a closed form for every profile above. Implement the closed form as the primary path.

| ID | ☐ | Task | Pass criterion | Output |
|---|---|---|---|---|
| P2.4.1 | ☐ | Shape-function module: *g*, its closed-form integral along a linear *d*(*t*), and *g*(0) for each profile | Unit tests vs symbolic or high-order numerical integration | `profiles.jl` |
| P2.4.2 | ☐ | Normalisation check: numerically integrate each *g* | ∫*g* = 1 to 10⁻⁸ | table |
| P2.4.3 | ☐ | Linear profile end-to-end | Single-spin trace exact (P2.2.9); P2.3.4(b)–P2.3.6 and P2.3.11 equivalents (the exact *h* = *w*/2 end point applies to the step profile only) | fig |
| P2.4.4 | ☐ | Polynomial *n* = 2, 3 end-to-end | As above | fig |
| P2.4.5 | ☐ | Exponential end-to-end, with truncation recorded | As above | fig |
| P2.4.6 | ☐ | Numerical-quadrature path (midpoint vs trapezoid, adaptive, tolerance ε_R₂) for profiles without a closed form | Agrees with closed form at small ε_R₂ | code |
| P2.4.7 | ☐ | **V4 per profile** — fit *h*\*(θ) to the B4 reference for each profile, fixed ΔR₂(0) | *h*\*(θ) with uncertainty and per-readout residuals for every profile; residuals reported, not pass/fail | fig + table |
| P2.4.8 | ☐ | Profile comparison at fixed ΔR₂(0): overlay the P2.3.4(b)-type curves *S*(*t*; *h*) for all profiles | Differences between profiles reported with σ; no expected equality assumed | fig |
| P2.4.9 | ☐ | *(optional)* Power law and half-Gaussian | As P2.4.3 | fig |
| P2.4.10 | ☐ | Tag `rd-profiles-v1` | Tag exists | git tag |

## 2.5 Shape comparison under consistent conditions (week 8)

Characterisation mode: **ρ fixed**, same reference configuration for every profile.

| ID | ☐ | Task | Pass criterion | Output |
|---|---|---|---|---|
| P2.5.1 | ☐ | Fixed ρ, λ/*w* ≪ 1 (e.g. 0.01–0.05), all profiles | Fixed ρ **and** fixed *h* (proposal §3.2 aim). Record the level of agreement between profiles with σ; no expected outcome assumed | table |
| P2.5.2 | ☐ | Fixed ρ, sweep λ/*w* from 0.01 to 0.5, all profiles | Curves of *S*(*t*) (and a descriptive rate) vs *h*/*w* per profile | fig |
| P2.5.3 | ☐ | Report centroid λ*ū* alongside λ for each profile (step 0.5λ, linear 0.33λ, poly *n* λ/(*n*+2), exponential 1.0λ) | Centroid column in every shape table | table |
| P2.5.4 | ☐ | Replot P2.5.2 against centroid instead of λ | Note whether shape differences shrink | fig |
| P2.5.5 | ☐ | Summary: where (if anywhere) profiles separate beyond 3σ on walls | Written paragraph with numbers | results note |

This is preliminary — the full sensitivity study in Phase 3 uses curved geometry and narrowing gaps, where differences are expected to be larger.

## 2.6 Timestep, tolerance and computational cost (weeks 9–11)

**Implementation plan (P2.6.1–P2.6.3, P2.6.5, P2.6.8, P2.6.10–P2.6.12; step profile, walls):** [docs/superpowers/plans/2026-10-06-timestep-study.md](../docs/superpowers/plans/2026-10-06-timestep-study.md). P2.6.4 follows the linear plan; P2.6.6, P2.6.7 and P2.6.9 have no numerical-integration path yet.

Converge one control at a time. Use the **same seed** across each sweep, then confirm with other seeds.

### Timestep

| ID | ☐ | Test | Method | Pass criterion | Output |
|---|---|---|---|---|---|
| P2.6.1 | ☑ | **V1** timestep plateau | Step profile, *h* = 0.1 µm; τ from 10⁻⁵ to 10⁻¹ ms, 3 per decade (B5: 10⁻⁵ ms costs ≈ 2 h per seed at 10⁵ spins × 50 ms, so use 10⁴ spins or a shorter sequence below 10⁻⁴ ms). If no plateau by 10⁻⁵ ms, record that and restrict *h* (see risks) | Plateau found; halving τ changes signal < 1σ; large-τ bias has the sign expected in proposal §3.2 V (too little attenuation) | toy scan + MCMR check: plateau to 0.1 ms at 1 % (τ_conv ≥ 0.1), 0.086 ms at 0.1 %; bias −1.1e-3 at 0.1 ms (less attenuation, as expected); MCMR = toy (\|z\| ≤ 1.3). `phase2_record.md` P2.6.1 |
| P2.6.2 | ☑ | **V2** *h*² scaling | Repeat at *h* = 0.05, 0.1, 0.2 µm; plot τ_conv vs *h*²/D (proposal Eq. 11) | Line through origin; slope = *c*_h | **FAIL**: τ_conv not ∝ *h*² (falls slightly with *h*). Decision point measured: ΔR₂(0)·τ_conv ≈ 0.05 (1 %), 0.003 (0.1 %); bias ∝ ΔR₂(0); ΔR₂(0) not ρ. `phase2_record.md` P2.6.2 |
| P2.6.3 | ☑ | **V3** 1/D scaling | λ = 0.1 µm, D = 1, 3, 6 | τ_conv ∝ 1/D | **FAIL**: τ_conv(0.1 %) = 0.043/0.086/0.089 ms for D = 1/3/6 (rises with D). `phase2_record.md` P2.6.3 |
| P2.6.4 | ☐ | Profile dependence of τ_conv | Repeat V1 for linear, polynomial, exponential | Record whether *c*_h depends on shape; if so, report per shape or use the most demanding | after the linear plan: rerun `p2_6_scan.jl` with a shape option; check ΔR₂(0)·τ_conv per shape |
| P2.6.5 | ☑ | Adopt constraint τ ≤ *c*_h·*h*²/D with one-decade margin; add to τ_max selection | Constraint active; report which term binds | **Done, changed form** (user decision 2026-10-06, since V2/V3 failed): τ ≤ c/ΔR₂(0)_max, option `layer`, default c = 0.005 (10× margin on 1 %); binds for ΔR₂(0) > 0.125 at w = 2 µm; verbose message names it. `phase2_record.md` P2.6.5 |

**Decision point:** if V2 does not give λ² scaling, test τ_conv against λ/ρ (the timescale 1/ΔR₂). Re-derive the constraint before adopting it.

### Integration tolerance

| ID | ☐ | Test | Method | Pass criterion | Output |
|---|---|---|---|---|---|
| P2.6.6 | ☐ | **V12** tolerance sweep | Quadrature path, ε_R₂ from 10⁻¹ to 10⁻⁶, each profile; compare to closed form | Signal error vs ε_R₂; fix ε_R₂ where error < 1σ | not applicable yet: every implemented profile has a closed form; revisit when one needs numerical integration |
| P2.6.7 | ☐ | Sub-division count vs ε_R₂ | Log sub-segments per step | Cost of tightening quantified | not applicable yet (see P2.6.6) |

### Computational cost

| ID | ☐ | Test | Method | Pass criterion | Output |
|---|---|---|---|---|---|
| P2.6.8 | ☑ | Runtime vs τ | Wall-clock per spin per ms of sequence, for the V1 sweep; ratio to the old code at the same τ (B5) | Cost curve alongside accuracy curve | overhead ×1.2–1.4 vs B5 at τ = 0.1–1e-4. `phase2_record.md` P2.6.8 |
| P2.6.9 | ☐ | Runtime: closed form vs quadrature | Same configuration, both paths | Speed-up quantified | not applicable yet (see P2.6.6) |
| P2.6.10 | ☑ | Runtime vs *h* at converged τ | *h* = 0.01–0.5 µm | Confirms feasible working range (proposal Table 1, risk "excessive timestep cost") | ≈ 0.03–0.05 h per full run for h = 0.01–0.5 (τ not h-limited); vs ΔR₂(0): 0.01 h (0.1) to 0.63 h (10). `phase2_record.md` |
| P2.6.11 | ☑ | Noise cost at fixed budget | Fewer spins required by smaller τ → σ increase (√ ratio) | σ penalty tabulated per *h* | σ penalty √(0.04/τ): 1.0 / 2.0 / 2.8 / 4.0 / 8.9 for ΔR₂(0) = 0.1 / 0.5 / 1 / 2 / 10 (independent of h). `phase2_record.md` |
| P2.6.12 | ☑ | Accuracy–cost summary | One figure: error vs runtime, points labelled by τ and ε_R₂ | Recommended settings stated | `results/phase2/p2_6_8_cost/cost.png`; recommended: default c = 0.005 (≈ 0.1 %), c = 0.05 for ≈ 1 %. Earlier τ = 1e-2 results: bias ≤ 6.4e-4, no reruns |
| P2.6.13 | ☐ | Tag `rd-walls-v1` | Tag exists | git tag |

## Phase 2 deliverables

- ☐ R₂(d) implementation on walls, all components unit-tested
- ☐ Step profile verified (V0b, V5 exact end point, V9, V11, occupancy) and characterised (V5 attenuation vs *h*, monotonicity, V4 *h*\*(θ) fit to the baseline). V10 deferred to future work
- ☐ Linear, polynomial (*n* = 2, 3) and exponential profiles implemented and verified; *h*\*(θ) fitted for every profile
- ☐ Preliminary shape comparison at fixed ρ on walls
- ☐ Calibrated *c*_h, chosen ε_R₂, accuracy–cost figure
- ☐ Tags `rd-step-v1`, `rd-profiles-v1`, `rd-walls-v1`

## Phase 2 risks

| Risk | Early warning | Response |
|---|---|---|
| Endpoint-based detection slips in | P2.2.2 test with neither endpoint inside fails | Use segment–plane intersection only |
| Density artefact near wall mimics relaxation | P2.3.2 non-uniform | Fix reflection handling before any relaxation test |
| Signal differences below noise | P2.1.4 shows 10σ not reachable | Raise ΔR₂(0) (verification mode) or N_spins |
| τ_conv too small to be practical | P2.6.1 no plateau by 10⁻⁵ ms at *h* = 0.1 µm (B5: ≈ 2 h per seed there) | Restrict *h* range; consider escape-from-a-layer approach (proposal Table 1) |
| Exponential tail truncation contaminates comparison | P2.4.5 loss > 1% | Keep λ ≤ *w*/5 or increase *w* |
| Shape-dependent *c*_h | P2.6.4 | Use the most demanding shape's constraint |

---

## Weekly log

| Week | Dates | Planned | Completed | Blockers / notes |
|---|---|---|---|---|
| 1 | 21–25 Sep | P1.1, P1.2.1–P1.2.3 | | |
| 2 | 28 Sep – 2 Oct | P1.2.4–P1.2.7, P1.3, P1.4 | P1.2.4–P1.2.7, Rep1–Rep4, P1.4.1–P1.4.4, wall baselines B0–B5 | Rep1 judged vs exact rate within 2 SEM (2026-10-02); V4/V5 and fast-diffusion-based Phase 2 tests rewritten as fitting/characterisation (2026-09-30) |
| 3 | 5–9 Oct | P2.1, P2.2.1 | | |
| 4 | 12–16 Oct | P2.2.2–P2.2.6 | | |
| 5 | 19–23 Oct | P2.2.7–P2.2.9, P2.3.1–P2.3.4 | | |
| 6 | 26–30 Oct | P2.3.5–P2.3.10, P2.4.1–P2.4.2 | | |
| 7 | 2–6 Nov | P2.4.3–P2.4.6 | | |
| 8 | 9–13 Nov | P2.4.7–P2.4.10, P2.5 | | |
| 9 | 16–20 Nov | P2.6.1–P2.6.3 | | |
| 10 | 23–27 Nov | P2.6.4–P2.6.7 | | |
| 11 | 30 Nov – 4 Dec | P2.6.8–P2.6.13 | | |
