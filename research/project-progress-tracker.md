# Project Progress Tracker — Distance-Dependent R₂(d) in MCMRSimulator

Working tracker for implementation. Tick a box when its pass criterion is met, not when the code merely runs.

**Status key:** ☐ not started · ◐ in progress · ☑ done · ✖ blocked (add a note)

**Calendar:** week 1 = 21–25 Sep 2026. Dates are targets; move them if the assessment calendar changes.

---

## Progress summary

| Phase | Weeks | Scope | Tasks | Done | Status |
|---|---|---|---|---|---|
| **Phase 1** | 1–2 (21 Sep – 2 Oct) | Verified baseline, unmodified code | 21 | 0 | ☐ |
| **Phase 2** | 3–11 (5 Oct – 4 Dec) | Implement and verify R₂(d) on wall geometry, no RF | 52 | 0 | ☐ |
| Phase 3 | 12–13 | Sensitivity studies (C1–C7), cylinder and sphere geometry | — | — | not yet planned |
| Phase 4 | 14–15 | Trade-off audit, write-up | — | — | not yet planned |

Update the "Done" column at the end of each week.

---

# PHASE 1 — Verified Baseline (Weeks 1–2)

**Goal.** Understand the existing code, reproduce published results with the **unmodified** simulator, and record every setting, so that any later change in signal can be attributed to the modification rather than to the setup.

**Exit criterion.** R1–R4 pass, the baseline record is complete, and the `baseline-v1` commit is tagged. Phase 2 does not start until this is met.

## 1.1 Environment

| ID | ☐ | Task | Pass criterion | Output |
|---|---|---|---|---|
| P1.1.1 | ☐ | Install Julia; install MCMRSimulator.jl **pinned to the paper's version (v1.0)** and MRIBuilder.jl | Package loads, bundled tutorial runs | — |
| P1.1.2 | ☐ | Create git repository; commit `Project.toml` and `Manifest.toml` | Environment reproducible from a clean clone | repo |
| P1.1.3 | ☐ | Set up Revise.jl and a persistent Julia session | Code edits reload without restart | — |
| P1.1.4 | ☐ | Clone the figure notebooks (`git.fmrib.ox.ac.uk/ndcn0236/mcmr_paper_figures`) | Notebooks open; note which are usable | notes |
| P1.1.5 | ☐ | Apply for cluster access | Application submitted (lead time for Phase 3) | confirmation |

## 1.2 Code study

Trace one simulation end to end. Target output is a written map, not general familiarity.

| ID | ☐ | Task | What to locate | Output |
|---|---|---|---|---|
| P1.2.1 | ☐ | **Code structure** | Module layout; how geometry, sequence and simulation modules depend on each other; how MRIBuilder sequences are passed in | code map §1 |
| P1.2.2 | ☐ | **Simulation workflow** | Main loop: setup → step → readout; where τ_max is computed from the constraint list (Supp. S1) | code map §2 |
| P1.2.3 | ☐ | **Geometry representation** | Storage of walls, cylinders, spheres, annuli, meshes; collision grid; repeating geometry; how surface vs volume parameters attach (Table 1) | code map §3 |
| P1.2.4 | ☐ | **Spin trajectories** | Gaussian step draw; collision detection; reflection and re-test; handling of perfectly permeable surfaces (§2.10) | code map §4 |
| P1.2.5 | ☐ | **Magnetisation update** | Bloch update between events; update-to-collision-time (§2.4); RF chunking (§2.5); **exact lines where `θ_relax` is applied** | code map §5 |
| P1.2.6 | ☐ | **Single-spin trace** — one cylinder, spin echo, ~100 spins; log one spin's position and magnetisation every step through at least one collision | Trace matches the reading of the code | trace log |
| P1.2.7 | ☐ | Mark insertion points for R₂(d): parameter definition, distance evaluation, relaxation application | Each insertion point named with file and function | code map §6 |

## 1.3 Reproduce baseline results

Use notebook values where available; values read off figures carry a few percent uncertainty.

| ID | ☐ | Result | Reference | Pass criterion | Output |
|---|---|---|---|---|---|
| P1.3.1 | ☐ | **R1** — surface relaxation, slab and single cylinder: R₂ = θ_relax·√(D/π)·(S/V) | Eq. 17, analytic | Within 2σ | fig + numbers |
| P1.3.2 | ☐ | **R2** — compartment T₂: ~80 ms intra-, ~60 ms extra-axonal, ~45 ms myelin | §3.3, Fig. 8 | Within 2σ, or 5% if read off figure | fig + numbers |
| P1.3.3 | ☐ | **R3** — timestep plateau, ADC vs τ, cylinder packings | Supp. Fig. S1 | Plateau onset near τ ≈ 0.01 ms | fig |
| P1.3.4 | ☐ | **R4** — surface-relaxation correction schemes, 1D slab | Supp. Fig. S2C | Logarithmic scheme flat to τ ≈ 10⁻³ ms; confirm θ_relax used | fig + θ_relax value |
| P1.3.5 | ☐ | *(optional)* **R5** — diffusion vs Mitra / van Gelderen | Fig. 5 | Qualitative agreement | fig |

**If a result fails:** find the cause before continuing. It is usually a parameter mismatch, and the cause is worth recording.

## 1.4 Record

| ID | ☐ | Task | Pass criterion | Output |
|---|---|---|---|---|
| P1.4.1 | ☐ | Baseline record for every run (software versions + commit hashes, hardware, geometry with exact S/V, physics, sequence, τ_max and binding constraint, N_spins, seeds, signal, apparent T₂, σ, runtime, reference, pass/fail) | Every run has a complete entry | `baseline_record.md` |
| P1.4.2 | ☐ | Noise floor σ for each baseline configuration, ≥ 10 seeds | σ reported per configuration | table |
| P1.4.3 | ☐ | Save raw outputs with generating script and seed | Any result regenerable from the repo | `results/baseline/` |
| P1.4.4 | ☐ | Tag commit `baseline-v1` | Tag exists | git tag |

## Phase 1 deliverables

- ☐ Code map with insertion points
- ☐ R1–R4 reproduced with pass/fail
- ☐ Baseline record, pinned environment, `baseline-v1` tag
- ☐ Noise floor σ per baseline configuration

---

# PHASE 2 — Implement and Verify R₂(d) on Wall Geometry (Weeks 3–11)

**Goal.** Implement the distance-dependent model and verify it on the simplest geometry available: planar walls, no RF pulses. Start with the step profile to verify each mechanical component separately, then add the remaining profiles, compare them under identical conditions, and establish timestep and integration-tolerance requirements.

**Why walls and no RF.** A planar wall makes every quantity exact: the distance is a single coordinate, the offset surface is another plane, the layer volume fraction is *f* = λ/*w* with no curvature, and along a straight step the distance changes **linearly in time**, so the relaxation integral has a closed form for every profile. Removing RF removes the second chunking driver, so any error is attributable to the R₂(d) code alone. Cylinder, sphere and finite-RF checks (V5 curved, V13) move to Phase 3.

**Exit criterion.** Section 2.3 passes; all four profiles pass Section 2.4; *c*_λ is calibrated and ε_R₂ fixed (Section 2.6). Tag `rd-walls-v1`.

## Reference configuration

Use this for every Phase 2 test unless a task says otherwise, so results are comparable.

| Parameter | Value | Note |
|---|---|---|
| Geometry | Repeating planar walls, spacing *w* = 2 µm | Layer on **one side** of each wall |
| D | 3 µm²·ms⁻¹ | |
| R₂^bulk | 0 | Isolates the layer term |
| R₁ | off | |
| θ_relax, θ_perm, susceptibility, MT | 0 / off | Recorded explicitly |
| Sequence | No RF; spins start fully transverse; free decay | No refocusing needed since off-resonance is zero |
| Readout times | 0, 10, 20, 30, 40, 50 ms | |
| Debug amplitude *A* | Chosen so attenuation at 50 ms is 0.2–0.8 | Physiological values are too weak to test |
| N_spins | 10⁵ | Adjust after P2.1.4 |
| Seeds | 10 | |

**Expected values, step profile, *A* = 0.5 ms⁻¹, *w* = 2 µm, t = 50 ms** — compute these before running:

| λ (µm) | *f* = λ/*w* | R₂ (ms⁻¹) | *M*ₓᵧ(50)/*M*ₓᵧ(0) |
|---|---|---|---|
| 0 | 0 | 0 | 1.000 |
| 0.1 | 0.05 | 0.025 | 0.287 |
| 0.2 | 0.10 | 0.050 | 0.082 |
| 0.4 | 0.20 | 0.100 | 0.007 |

(Reduce *A* or *t* if the λ = 0.4 point is below the noise floor.)

## 2.1 Design and test harness (week 3)

| ID | ☐ | Task | Pass criterion | Output |
|---|---|---|---|---|
| P2.1.1 | ☐ | Define the parameter interface: `rho`, `lambda`, `shape` as **surface** parameters; `A = rho/lambda` derived internally, never user-settable | Interface agreed and documented | design note |
| P2.1.2 | ☐ | Construction-time validation: reject λ < 0, ρ < 0; warn if *A*·*g*(0) exceeds the rigid-lattice ceiling; warn if λ > object size | Invalid inputs rejected with clear messages | code + tests |
| P2.1.3 | ☐ | Test harness: runs a configuration over seeds, returns mean, σ, runtime, layer-visit count; computes ρ from a target amplitude for verification mode | Harness reproduces a Phase 1 result | `harness.jl` |
| P2.1.4 | ☐ | Noise floor for the reference configuration at λ = 0 | σ recorded; N_spins adjusted so the λ = 0.1 → 0.2 difference exceeds 10σ | σ table |
| P2.1.5 | ☐ | Analytic expectations script: *f*(λ), R₂, attenuation for every test below | Tables generated before any run | `expected.jl` |

## 2.2 Step profile — component implementation (weeks 3–5)

Implement and test each component separately, in this order. Each has its own unit test, so a failure points to one component.

| ID | ☐ | Component | What it does | Unit test | Pass criterion |
|---|---|---|---|---|---|
| P2.2.1 | ☐ | **Distance evaluation** | Signed distance from a point to the nearest wall face carrying a layer | Random points in a cell vs hand-computed distance | Exact to floating-point precision |
| P2.2.2 | ☐ | **Crossing detection — entry** | Detect when a step segment crosses the offset plane *d* = λ inward | Constructed segments with known crossings, including ones with **neither endpoint inside the layer** | All crossings found, times exact |
| P2.2.3 | ☐ | **Crossing detection — exit** | Same, outward | As above | All crossings found |
| P2.2.4 | ☐ | **Crossing detection — with reflection** | Segment enters layer, reflects off wall, exits, within one step | Constructed cases | Correct sequence of sub-segments |
| P2.2.5 | ☐ | **Trajectory segmentation** | Split each step into sub-segments at every crossing and reflection, each labelled inside/outside | Sub-segment durations sum to τ; labels correct | Exact |
| P2.2.6 | ☐ | **Local relaxation evaluation** | ∫ΔR₂(*d*(*t*)) d*t* over each sub-segment; for the step profile this is *A*·(time inside) | Constructed segments vs hand calculation | Exact |
| P2.2.7 | ☐ | **Transverse magnetisation update** | *M*ₓᵧ ← *M*ₓᵧ·exp(−∫ΔR₂ d*t*), applied per sub-segment in time order; phase untouched | Single spin held inside the layer decays at exactly *A* | Exact |
| P2.2.8 | ☐ | **λ = 0 short-circuit** | No layer code executes when λ = 0 | Code-path check; runtime equal to baseline | No measurable slowdown |
| P2.2.9 | ☐ | End-to-end single-spin trace, as in P1.2.6, with the layer on | Logged decay matches manual integration of the logged path | Exact |

## 2.3 Step profile — verification (weeks 5–6)

| ID | ☐ | Test | Method | Pass criterion | Output |
|---|---|---|---|---|---|
| P2.3.1 | ☐ | **V0b** baseline reproduction | λ = 0 vs `baseline-v1` | Within 1σ | table |
| P2.3.2 | ☐ | **V9** equilibrium density | Relaxation off; histogram spin density vs *d* | Uniform within Poisson error (cf. RMS App. A for the failure mode) | histogram |
| P2.3.3 | ☐ | **Layer occupancy check** | Relaxation off; count spins inside the layer | Fraction = λ/*w* within Poisson error. **Validates geometry independent of relaxation code** | table |
| P2.3.4 | ☐ | **V5 linearity** (fixed amplitude, sweep λ = 0.05–0.4 µm) | R₂ vs λ | Straight line through origin; slope = *A*/*w* within 2σ; **no curvature** | fig |
| P2.3.5 | ☐ | **Monotonicity** | Attenuation vs λ | Strictly increasing | fig |
| P2.3.6 | ☐ | **Mono-exponential decay** | *M*ₓᵧ(t) at six readout times, thin layer | Log-linear in *t* (fast-exchange limit) | fig |
| P2.3.7 | ☐ | **V11** intersection vs endpoint detection | Run both methods at τ = 1 µs | Endpoint method under-counts; bias quantified | table |
| P2.3.8 | ☐ | **V10** Route A cross-check | Same step profile built from perfectly permeable walls with per-compartment R₂ | Within 2σ of Route B | table |
| P2.3.9 | ☐ | **V4** λ → 0 limit | λ decreasing, ρ fixed, vs baseline θ_relax result (ρ = θ_relax√(D/π)) | Converges within 2σ | fig |
| P2.3.10 | ☐ | **Tag** `rd-step-v1` | — | Tag exists | git tag |

**Gate:** do not start 2.4 until P2.3.3 and P2.3.4 pass. Failure there means distance or segmentation is wrong, and every later profile would inherit it.

## 2.4 Additional profiles (weeks 6–8)

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
| P2.4.3 | ☐ | Linear profile end-to-end | Passes P2.3.3–P2.3.6 equivalents | fig |
| P2.4.4 | ☐ | Polynomial *n* = 2, 3 end-to-end | As above | fig |
| P2.4.5 | ☐ | Exponential end-to-end, with truncation recorded | As above | fig |
| P2.4.6 | ☐ | Numerical-quadrature path (midpoint vs trapezoid, adaptive, tolerance ε_R₂) for profiles without a closed form | Agrees with closed form at small ε_R₂ | code |
| P2.4.7 | ☐ | **Slope-ratio test** — fixed wall value *A*, sweep λ; each profile's R₂ vs λ slope relative to step | Ratios = *k*ₕ = ∫*h*: step 1, linear 1/2, poly *n* 1/(*n*+1), exponential 1 — within 2σ | fig + table |
| P2.4.8 | ☐ | Step vs exponential cross-check (same *k*ₕ = 1) | Identical within 2σ in the thin-layer limit | table |
| P2.4.9 | ☐ | *(optional)* Power law and half-Gaussian | As P2.4.3 | fig |
| P2.4.10 | ☐ | Tag `rd-profiles-v1` | Tag exists | git tag |

## 2.5 Shape comparison under consistent conditions (week 8)

Characterisation mode: **ρ fixed**, same reference configuration for every profile.

| ID | ☐ | Task | Pass criterion | Output |
|---|---|---|---|---|
| P2.5.1 | ☐ | Fixed ρ, λ/*w* ≪ 1 (e.g. 0.01–0.05), all profiles | Record the level of agreement; expected to be within noise (only the zeroth moment survives) | table |
| P2.5.2 | ☐ | Fixed ρ, sweep λ/*w* from 0.01 to 0.5, all profiles | Curves of apparent R₂ vs λ/*w* per profile | fig |
| P2.5.3 | ☐ | Report centroid λ*ū* alongside λ for each profile (step 0.5λ, linear 0.33λ, poly *n* λ/(*n*+2), exponential 1.0λ) | Centroid column in every shape table | table |
| P2.5.4 | ☐ | Replot P2.5.2 against centroid instead of λ | Note whether shape differences shrink | fig |
| P2.5.5 | ☐ | Summary: where (if anywhere) profiles separate beyond 3σ on walls | Written paragraph with numbers | results note |

This is preliminary — the full sensitivity study in Phase 3 uses curved geometry and narrowing gaps, where differences are expected to be larger.

## 2.6 Timestep, tolerance and computational cost (weeks 9–11)

Converge one control at a time. Use the **same seed** across each sweep, then confirm with other seeds.

### Timestep

| ID | ☐ | Test | Method | Pass criterion | Output |
|---|---|---|---|---|---|
| P2.6.1 | ☐ | **V1** timestep plateau | Step profile, λ = 0.1 µm; τ from 10⁻⁷ to 10⁻¹ ms, 3 per decade | Plateau found; halving τ changes signal < 1σ; **large-τ bias has predicted sign (too little attenuation)** | fig |
| P2.6.2 | ☐ | **V2** λ² scaling | Repeat at λ = 0.05, 0.1, 0.2 µm; plot τ_conv vs λ²/D | Line through origin; slope = *c*_λ | fig + *c*_λ |
| P2.6.3 | ☐ | **V3** 1/D scaling | λ = 0.1 µm, D = 1, 3, 6 | τ_conv ∝ 1/D | fig |
| P2.6.4 | ☐ | Profile dependence of τ_conv | Repeat V1 for linear, polynomial, exponential | Record whether *c*_λ depends on shape; if so, report per shape or use the most demanding | table |
| P2.6.5 | ☐ | Adopt constraint τ ≤ *c*_λλ²/D with one-decade margin; add to τ_max selection | Constraint active; report which term binds | code |

**Decision point:** if V2 does not give λ² scaling, test τ_conv against λ/ρ (the timescale 1/ΔR₂). Re-derive the constraint before adopting it.

### Integration tolerance

| ID | ☐ | Test | Method | Pass criterion | Output |
|---|---|---|---|---|---|
| P2.6.6 | ☐ | **V12** tolerance sweep | Quadrature path, ε_R₂ from 10⁻¹ to 10⁻⁶, each profile; compare to closed form | Signal error vs ε_R₂; fix ε_R₂ where error < 1σ | fig + chosen ε_R₂ |
| P2.6.7 | ☐ | Sub-division count vs ε_R₂ | Log sub-segments per step | Cost of tightening quantified | table |

### Computational cost

| ID | ☐ | Test | Method | Pass criterion | Output |
|---|---|---|---|---|---|
| P2.6.8 | ☐ | Runtime vs τ | Wall-clock per spin per ms of sequence, for the V1 sweep | Cost curve alongside accuracy curve | fig |
| P2.6.9 | ☐ | Runtime: closed form vs quadrature | Same configuration, both paths | Speed-up quantified | table |
| P2.6.10 | ☐ | Runtime vs λ at converged τ | λ = 0.01–0.5 µm | Confirms feasible working range (proposal §5.4) | table |
| P2.6.11 | ☐ | Noise cost at fixed budget | Fewer spins required by smaller τ → σ increase (√ ratio) | σ penalty tabulated per λ | table |
| P2.6.12 | ☐ | Accuracy–cost summary | One figure: error vs runtime, points labelled by τ and ε_R₂ | Recommended settings stated | fig + note |
| P2.6.13 | ☐ | Tag `rd-walls-v1` | Tag exists | git tag |

## Phase 2 deliverables

- ☐ R₂(d) implementation on walls, all components unit-tested
- ☐ Step profile verified (V0b, V4, V5, V9, V10, V11, occupancy)
- ☐ Linear, polynomial (*n* = 2, 3) and exponential profiles implemented and verified; slope-ratio test passed
- ☐ Preliminary shape comparison at fixed ρ on walls
- ☐ Calibrated *c*_λ, chosen ε_R₂, accuracy–cost figure
- ☐ Tags `rd-step-v1`, `rd-profiles-v1`, `rd-walls-v1`

## Phase 2 risks

| Risk | Early warning | Response |
|---|---|---|
| Endpoint-based detection slips in | P2.2.2 test with neither endpoint inside fails | Use segment–plane intersection only |
| Density artefact near wall mimics relaxation | P2.3.2 non-uniform | Fix reflection handling before any relaxation test |
| Signal differences below noise | P2.1.4 shows 10σ not reachable | Raise *A* (verification mode) or N_spins |
| τ_conv too small to be practical | P2.6.1 plateau below 10⁻⁶ ms at λ = 0.1 µm | Restrict λ range; consider escape-from-a-layer approach |
| Exponential tail truncation contaminates comparison | P2.4.5 loss > 1% | Keep λ ≤ *w*/5 or increase *w* |
| Shape-dependent *c*_λ | P2.6.4 | Use the most demanding shape's constraint |

---

## Weekly log

| Week | Dates | Planned | Completed | Blockers / notes |
|---|---|---|---|---|
| 1 | 21–25 Sep | P1.1, P1.2.1–P1.2.3 | | |
| 2 | 28 Sep – 2 Oct | P1.2.4–P1.2.7, P1.3, P1.4 | | |
| 3 | 5–9 Oct | P2.1, P2.2.1 | | |
| 4 | 12–16 Oct | P2.2.2–P2.2.6 | | |
| 5 | 19–23 Oct | P2.2.7–P2.2.9, P2.3.1–P2.3.4 | | |
| 6 | 26–30 Oct | P2.3.5–P2.3.10, P2.4.1–P2.4.2 | | |
| 7 | 2–6 Nov | P2.4.3–P2.4.6 | | |
| 8 | 9–13 Nov | P2.4.7–P2.4.10, P2.5 | | |
| 9 | 16–20 Nov | P2.6.1–P2.6.3 | | |
| 10 | 23–27 Nov | P2.6.4–P2.6.7 | | |
| 11 | 30 Nov – 4 Dec | P2.6.8–P2.6.13 | | |
