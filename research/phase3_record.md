# Phase 3 · Record (R₂(d) on cylinders, no RF)

Settings, methods and results of Phase 3 (tracker: `research/project-progress-tracker.md`, Phase 3). The model is the user's own:

  R₂,total(x) = R₂,bulk + ΔR₂(d(x)),  ΔR₂(d) = (ρ/h)·g(d/h)/κ,

where d is the distance from a cylinder surface, measured into the compartment. The conventional surface-relaxation model of MCMR (per-collision loss θ√τ) is a reference, not the truth. No fast-diffusion assumption is used anywhere.

Each section follows the same pattern: **design** (what was chosen, the alternatives, why), **method** (the equations computed, the definition of every quantity, how uncertainty is estimated), **result**, and **conclusion** (numbers tied back to the formulas).

## Summary

| ID | Item | Criterion | Outcome |
|---|---|---|---|
| D3.1–D3.6 | Design decisions for cylinders | user decision | **Decided** (2026-10-24), see below |
| D3.1 check | Curvature normalisation κ = 1 ∓ c̄/R | quadrature = closed form | **PASS**: max deviation 1.6e-9 |
| P3.2.1–P3.2.4, P3.2.6 | Step layer on cylinders: exact shell fraction, images, overlap sum, fields, validation, timestep bound | package tests | **PASS**: 100 tests; full suite 3497/3497 |
| P3.1.1, P3.1.2 | C-B1 null, C-B2 bulk (single, array inside/outside) | exact / ≤ 1e-10 per spin, ≤ roundoff bound | **PASS** (3/3): per spin 5.4e-14; τ_default = 0.01 ms |
| P3.1.3 | C-B3 density and occupancy within h | \|z_χ\| ≤ 3; occupancy \|z\| ≤ 2 | χ² **PASS** (z +0.33 / +0.78); occupancy **FAIL** by criterion (outside h = 0.1/0.2: z +2.36/+2.89, nested shells); follow-up vs τ: no trend; identity P3.4.3 confirms uniform density to 1e-4 |
| P3.4.5 | Positions with the layer on = no layer | identical | **PASS** |
| P3.1.4 | C-B4 conventional θ reference curves (5 configurations × 3 θ) | τ_a vs τ_a/4 within 2 SEM; bulk factorisation | 14/15 (+ 1 confirmed at τ_a/16, z +0.83); bulk 5.9e-14 → **reference set** |
| P3.1.5 | C-B5 runtime | cost table | layer ×1.25–1.67 per spin-step (1 cylinder per cell); ×6.7 at 400 per cell (cost ∝ cylinders per cell) |
| P3.1.6 | Conventional mapping θ√(D/π)(S/V) on curved surfaces | \|z\| ≤ 3 | **PASS** (18/18, max \|z\| 2.17); ratio 0.99976–1.00013 for √(Dτ)/R ≤ 0.35 |
| P3.4.1 | Single-spin trace, square array, overlapping outer layers | every step ≤ 1e-6 | **PASS**: max 1.0e-7; 9808 reflected steps with a layer exponent, 2978 overlap steps |
| P3.2.5, P3.4.2, P3.4.4, P3.4.7 | V0b, uniform h = R, additivity, side-specific (full size) | exact / ≤ 1e-10 / leak ≤ 1e-9 | **PASS** (8/8) |
| P3.4.3 | Mean-exponent identity E[−ln M] = ρ(S/V)t, step (single, square, overlap, rectangular) | \|z\| ≤ 3 | **PASS** (12/12, \|z\| ≤ 1.25; ≤ 6e-5 relative) |
| P3.4.6 | Radial continuum benchmark (single cylinder, inside): solver vs exact Bessel eigenvalues; MC vs continuum at converged τ | solver tests pass; MC within 3 SEM | solver **PASS** (25 tests); MC 5/6 **PASS**; h = 0.05 paired ΔS **FAIL** (z 4.1 at τ = 6.25e-4, converging: Richardson −4.6e-5 vs −5.75e-5); C1 timestep bias on cylinders: thin-layer ΔS sign flips at τ = 1e-2 |
| P3.4.8 | Attenuation vs h at fixed ρ (step): inside continuum (R = 0.5, 1, 2), outside paired MC (s = 2.5, 3.0) | monotonicity, sign recorded; curvature of ln S | S decreases with h: inside 27/27 steps resolved; outside 10/13 steps beyond 3σ (all layer steps from h = 0.1 on), conventional → 0.05 not resolved (z −2.4, −1.8; residual τ bias); J ≈ f(h/R) inside (0.3 %); effects ≤ 0.33 % of S at ρ = 0.02 |
| P3.3.1–P3.3.3 | Linear and exponential on cylinders: quadrature between exact roots vs linear closed form, wall formulas, brute force; supports (D3.5), validation, simulation | package tests | **PASS**; step results bit-identical to before |
| P3.3.4, P3.6.3 | Quadrature tolerance and cost | error vs tolerance; runtime | linear vs closed form ≤ 9e-11; exponential 3.2e-11 at the default 1e-10 (kept); cost ×1.2 (linear), ×2.5 (exponential) vs step; no allocation |
| P3.4.1, P3.4.3 (linear, exponential) | Single-spin trace; mean-exponent identity | trace ≤ 1e-6; \|z\| ≤ 3 | **PASS**: trace 6.7e-12 (linear), 2.4e-8 (exponential); identity 12/12 each (\|z\| ≤ 1.27) |
| P3.4.6 (linear, exponential) | MC vs radial continuum inside (continuum now with exact ring averages for all profiles) | max \|z\| ≤ 3 at τ_b for S and paired ΔS | **PASS** (5/5; ΔS \|z\| ≤ 2.09); \|ΔS\| low by 4–38 % at τ = 1e-2 and still 2.5–12 % low (same sign, not resolved) at τ = 2.5e-3 (C1) |
| P3.4.8 (all profiles) | Attenuation vs h and same-centroid comparison: inside (continuum, R = 0.5, 1, 2), outside (paired MC, s = 3.0, thin end at τ = 6.25e-4) | monotonicity, sign recorded; centroid spread | S decreases with h: inside 84/84 steps (all profiles); outside every step resolved for linear and exponential (thin end at τ_c), step thin end unresolved; centroid collapse: 3–9 % inside, 1–22 % outside at the same c̄ (outside not pinned down: 3–6 % or 13–15 % at c̄ = 0.15–0.2 depending on τ) vs 34–66 % and 7–72 % at the same h |

---

## Design decisions (D3.1–D3.6)

Decided by the user on 2026-10-24, following the recommendations below. Notation:
- cylinder radius R; axis along z;
- r = distance of a point from the axis;
- inside, d = R − r ∈ [0, R]; outside, d = r − R ≥ 0;
- layer thickness h; profile g(u) with ∫g(u) du = 1 over its support;
- profile centroid c̄ = h∫u g(u) du (step h/2, linear h/3, exponential truncated at support c·h: h(1 − (c + 1)e^(−c))/(1 − e^(−c))).

### D3.1 — ρ is defined per unit surface area

**Problem.** On a wall the volume element at distance d is the same as at the surface, so the integral along the normal,

  ρ_n = ∫ ΔR₂(d) dd,

is also the relaxation per unit surface area. On a cylinder the shell at distance d has circumference 2π(R ∓ d), so per unit cylinder length

  ∫ ΔR₂ dA = 2πR ∫ ΔR₂(d) (1 ∓ d/R) dd,

with − inside and + outside. The two notions of "ρ" differ by the curvature factor.

**Options.**
- (a) ρ along the normal: ΔR₂(d) = (ρ/h) g(d/h), the same formula as on walls.
- (b) ρ per unit surface area: ρ = ∫ ΔR₂(d)(1 ∓ d/R) dd.

**Decision: (b).** Inserting ΔR₂ = (ρ/h) g(d/h)/κ gives the normalisation

  κ = ∫ g(u) (1 ∓ u·h/R) du = 1 ∓ c̄/R.

**Why.**
1. With (b), the volume mean of ΔR₂ over a compartment of volume V and surface A is exactly ρ·A/V. That is the same expression as for the conventional model with relaxivity ρ (rate ρ·S/V). So "same ρ" keeps meaning "same total relaxation per unit surface", which every comparison in Phase 3 relies on.
2. The mean-exponent identity used for verification keeps its simple form, E[−ln M] = ρ·(A/V)·t, because the uniform density is stationary.
3. On walls (R → ∞), κ → 1 and (b) reduces to the Phase 2 definition.

**Cost.**
- The surface rate becomes ΔR₂(0) = (ρ/h) g(0)/κ, which depends on curvature.
- Inside a cylinder κ < 1, so the profile is scaled up; outside, κ > 1 scales it down (figure, right panel).
- Results quoted per ρ are not directly comparable with option (a). When ΔR₂(0) is the quantity held fixed, both conventions are given.

**Check.** κ computed by midpoint quadrature (2·10⁵ points) agrees with 1 ∓ c̄/R to 1.6e-9 for step, linear and exponential, inside and outside, h/R = 0.05–1 (script [p3_0_design.jl](phase3/p3_0_design.jl), table `results/phase3/p3_0_design/kappa.csv`). Examples at h/R = 0.5: step inside 0.750 / outside 1.250; linear 0.833 / 1.167; exponential (inside support R, outside 2 µm) 0.657 / 1.463.

![κ vs h/R and the two definitions of ρ](results/phase3/p3_0_design/p3_0_kappa.png)

### D3.2 — Profiles: step first, then linear and exponential

**Decision.** The step profile is implemented and verified first. Linear and exponential follow once the integration path of D3.3 exists.

**Why.** For the step profile the exponent over a straight piece is ΔR₂(0)·dt·(fraction of the piece inside the shell). The fraction follows from exact line–circle intersections (D3.3). So the step profile on cylinders keeps an exact closed form and can be verified to roundoff, as on walls, before any numerical integration is introduced.

### D3.3 — Integration along a straight piece: quadrature for all profiles, linear closed form as its test

**Problem.** A straight piece from p₀ to p₁, with transverse displacement Δ⊥, has squared distance from the axis

  r(s)² = |p₀⊥ − c + s Δ⊥|² = a s² + b s + q,  s ∈ [0, 1],

where a = |Δ⊥|², b = 2(p₀⊥ − c)·Δ⊥ and q = |p₀⊥ − c|². So d(s) = R − √(as² + bs + q) (inside) is not linear in s. The exponent of the piece, of duration dt, is

  E = dt ∫₀¹ ΔR₂(d(s)) ds.

**Method.**
- The roots of r(s) = R − h and r(s) = R (or R and R + h outside) split [0, 1] into sub-intervals that lie entirely inside or entirely outside the shell. These are exact quadratic roots, with tangent and chord-only cases handled explicitly.
- **Step:** E = ΔR₂(0)·dt·(total length of the inside sub-intervals), exactly.
- **Linear and exponential:** an adaptive quadrature on each inside sub-interval, with tolerance ε_R₂ on E, fixed in P3.6.3.
- **Independent check:** the linear profile also has a closed form, because ∫√(as² + bs + q) ds is elementary (a logarithm/asinh form). The quadrature is tested against it (P3.3.2) and against the Phase 2 wall closed forms (P3.3.3).

**Why.** A single integration path for all profiles keeps the code small, and the closed forms make it testable. This activates the Phase 2 items deferred for walls: P2.4.6 (quadrature path), P2.6.6 (V12 tolerance), P2.6.7 and P2.6.9 (cost).

### D3.4 — Both sides, side-specific parameters

**Decision.** Each cylinder surface has an inner (inside) and an outer (outside) layer, each with its own (ρ, h, g). Both sides are symmetric by default. This is the walls' positive/negative pair transferred to cylinders.

**Why.** The intra- and extra-cylinder compartments can have different surface chemistry (the user's model, as in Phase 2). Setting one side to zero recovers a one-sided layer.

### D3.5 — Exponential support

**Decision.**
- **Inside:** support d ≤ R, so it reaches the axis, with the inside κ above.
- **Outside:** support d ≤ L, where L is the `layer_cutoff` and defaults to the smallest surface-to-surface gap between neighbouring cylinders, L = s_min − 2R, with s_min the smallest centre spacing. A single, non-repeating cylinder must set L explicitly.

**Why.** On walls the exponential is cut at the distance to the opposite surface (cutoff = wall spacing), so it reaches exactly the neighbouring surface. The cylinder analogue is the surface gap s − 2R.

**Clarification (2026-10-24).** The decision table first said "half the smallest centre spacing". That was a slip, corrected here to the gap, which matches the wall convention. Renormalising g over its support keeps ∫g = 1, so ρ stays the relaxation per unit area (D3.1).

### D3.6 — Arrays: square and rectangular

**Decision.** Square arrays (centre spacing s in x and y) and rectangular arrays (s_x ≠ s_y), as MCMR `Cylinders(repeats = [s_x, s_y])`. Random packings are an optional extension in 3.7.

**Why.** These give the three overlap classes of the proposal (Milestone 5) by construction. With the surface gaps G_x = s_x − 2R and G_y = s_y − 2R:
- **separated:** G_x > 2h and G_y > 2h;
- **one-direction overlap:** G_x ≤ 2h < G_y;
- **two-direction overlap:** G_x ≤ 2h and G_y ≤ 2h.

At overlap the contributions of all cylinders within reach are summed (Eq. 10 of the proposal), as on walls.

---

## Implementation of the step layer on cylinders (P3.2.1–P3.2.4, P3.2.6)

Plan [2026-10-24-cylinders-step.md](../docs/superpowers/plans/2026-10-24-cylinders-step.md). Code:
- [layers.jl](../src/geometries/internal/layers.jl): `CylinderLayer`, `kappa`, `cylinder_surface_rate`, `disk_length`, `shell_fraction`, `cylinder_exponent`, `image_centres`, `layer_exponent`, `max_layer_rate`;
- [fix.jl](../src/geometries/user/fix.jl): `cylinder_layers`, `check_outer_reach`, `obstruction_shifts`.

Tests: [test_layer_cylinders.jl](../test/test_layer_cylinders.jl), 100 tests. The full suite without plots passes, 3497/3497.

**Design.**
- **Exact shell time in closed form.** Along a straight piece p(s) = p₀ + sΔ, s ∈ [0, 1], the squared distance from the axis (in the transverse plane of the cylinder) is

  r(s)² = a s² + b s + q,  a = |Δ⊥|²,  b = 2(p₀⊥ − c)·Δ⊥,  q = |p₀⊥ − c|².

  It is convex in s, so {s : r(s)² ≤ X} is one interval of length L(X), bounded by the roots. The fraction of the piece in the shell r_lo ≤ r ≤ r_hi is

  f = L(r_hi²) − L(r_lo²),

  with no side bookkeeping: the inner shell is [R − h_in, R] and the outer shell is [R, R + h_out]. The roots use the cancellation-free form t = −(b + sign(b)√disc)/2, s = t/a and (q − X)/t. A piece with a = 0 (purely axial motion) is a stationary point.
- **Exponent of a piece.** For the step profile, E = ΔR₂(0)·dt·f exactly. With D3.1,

  ΔR₂(0) = (ρ/h)/κ,  κ = 1 − h/(2R) inside,  κ = 1 + h/(2R) outside.

- **Arrays.** Every image of every cylinder within reach R + h_out of the piece's bounding box is summed (proposal Eq. 10), so overlapping outer layers add. Rotated axes come from the group's rotation: motion along the axis does not change r.
- **Alternatives rejected.**
  - Sampling the piece: not exact, and slow.
  - Tracking which side a piece is on: unnecessary, because the two shells are disjoint radial ranges.
- **Validation.**
  - h_inner ≤ R.
  - h_outer ≤ the smallest surface gap |c_i − c_j| − R_i − R_j to any other cylinder or image (centres wrapped into the cell, images within ±1 cell). This is the wall rule h ≤ spacing, so an outer layer never relaxes spins inside another cylinder.
  - Only the step profile is accepted on cylinders so far (D3.2).
  - Fields: `layer_rho`, `layer_h`, `layer_shape`, each with `_inner` and `_outer` versions; JSON and CLI supported.
- **Timestep (P3.2.6).** τ ≤ c/ΔR₂,max with c = 0.005. Here ΔR₂,max is the larger of:
  - the inner surface rate;
  - the outer surface rate plus the outer rates of every cylinder or image whose outer shell intersects this one's (|c_i − c_j| < R_i + h_i + R_j + h_j).

  This equals the largest rate when shells do not overlap and is an upper bound otherwise. Example: a square array s = 2.5, R = 1, h = 0.3 has 5 overlapping outer rates at the 4 nearest neighbours.

**Verification (package tests).**
- Shell fraction against brute-force sampling (10⁵ points, 3000 pieces): max error < 3e-5, the sampling bound. Exact cases:
  - radial;
  - a chord;
  - tangent to the outer circle (fraction 0) and to the inner circle;
  - h = R;
  - stationary points.
- Review Focus cases:
  - nearly axial motion (|Δ⊥| = 1e-15);
  - a piece starting on r = R.
- Layer exponent against an independent brute-force sum of rates (2·10⁴ samples, 200 pieces, max error < 1e-4):
  - two and four overlapping cylinders;
  - a square array with images;
  - two cylinders per cell in a rectangular array;
  - a cylinder along x (r measured in the (y, z) plane).
- A hand-computed radial step with a reflection, inside and outside, to 1e-12.
- ρ = 0 bit-identical to no layer; h_inner = R uniform (every inside spin at e^(−0.1t)); bulk factorisation.
- Side-specific layers; a permeable cylinder; finite RF rejected.
- Mutation check: scaling the surface rate by 1.01 is caught by 22 tests.

**Roundoff at the surface (ruling).** A reflection point can lie a roundoff distance across r = R, inside the other side's shell (its boundary has measure zero). An outer-only layer then leaves inside spins at 1 − M⊥ ≈ 1e-15 per reflection rather than exactly 1, accumulating over 50 ms to at most 1.2e-11 at full size (10⁶ spins). The tests bound this leak at 1e-12 (package test, 2000 spins, 50 ms; max seen 1.3e-15) and 1e-9 (full size), eight orders below the layer's own effect (1 − S ≈ 0.87).

## C-B1, C-B2: null run and bulk R₂ on cylinders (P3.1.1, P3.1.2)

Script [p3_1_1_null_bulk.jl](phase3/p3_1_1_null_bulk.jl) → [results/phase3/p3_1_1_null_bulk/](results/phase3/p3_1_1_null_bulk/).

**Design.** These are the exact baselines of the old code on the geometries used in Phase 3:
- a single cylinder R = 1 µm (inside);
- a square array s = 2.5 µm, R = 1 (inside and outside; surface gap 0.5 µm).

Default timestep, seeds 1–10, 10⁵ spins per seed, readouts 0:5:50 ms. Spins for the single cylinder start uniformly in the box [−R, R]³; the region is set by the distance from the nearest axis at t = 0.

**Method.**
- **C-B1:** every spin keeps M⊥ = 1 exactly.
- **C-B2 (R₂ = 1/80 ms⁻¹, same seeds):**
  - positions identical to C-B1;
  - per spin, |M⊥/e^(−t/80) − 1| ≤ 1e-10;
  - ensemble within B = (t/τ + N)·ε = 2.3e-11.

**Result.**

| configuration | τ_default (ms) | spins in region (mean per seed) | C-B1 | C-B2 per spin | C-B2 ensemble | outcome |
|---|---|---|---|---|---|---|
| single, inside | 0.01 | 78 501 | exact | 5.4e-14 | 4.4e-14 | **PASS** |
| array, inside | 0.01 | 50 247 | exact | 5.4e-14 | 4.4e-14 | **PASS** |
| array, outside | 0.01 | 49 753 | exact | 5.4e-14 | 4.4e-14 | **PASS** |

**Conclusion.** τ_default = 0.03·R²/D = 0.01 ms, set by the tortuosity constraint with size scale R. The bulk factor is exact to roundoff on cylinders.

## C-B3: density and shell occupancy (P3.1.3), and with the layer on (P3.4.5)

Scripts [p3_1_3_density.jl](phase3/p3_1_3_density.jl) and [p3_1_3b_density_tau.jl](phase3/p3_1_3b_density_tau.jl) → [results/phase3/p3_1_3_density/](results/phase3/p3_1_3_density/), [results/phase3/p3_1_3b_density_tau/](results/phase3/p3_1_3b_density_tau/).

**Design.** By (M2) of Phase 2, the mean layer exponent is ΔR₂ times the time fraction in the shell. So the density near a curved surface must be uniform, including at finite τ, where finite steps are reflected on a curved surface. The layer must also not move spins.

**Method.** Positions at t = 50 ms, pooled over 10 seeds (785 012 spins inside, 497 526 outside).
- **Inside:** for a uniform density, r²/R² is uniform on [0, 1]. Test: 100-bin χ², z_χ = (χ² − dof)/√(2 dof).
- **Outside:** for R ≤ r ≤ s/2 whole annuli lie in the cell, so P(r_k ≤ r < r_k+1 | r ≤ s/2) = (r_k+1² − r_k²)/((s/2)² − R²). Test: 50 equal-area bins.
- **Occupancy within h:** binomial z = (n − Np)/√(Np(1 − p)), with
  - p = (R² − (R − h)²)/R² inside;
  - p = π((R + h)² − R²)/(s² − πR²) outside.
- **Criteria:** |z_χ| ≤ 3, |z| ≤ 2 per h.

**Result.**

| test | no layer | layer on (ρ = 0.02, h = 0.2) |
|---|---|---|
| inside χ²/dof (z) | 1.047 (+0.33) | 1.047 (+0.33) |
| outside χ²/dof (z) | 1.159 (+0.78) | 1.159 (+0.78) |
| occupancy z inside, h = 0.05 / 0.1 / 0.2 / 0.5 | −1.90 / −1.12 / −1.69 / −1.56 | same |
| occupancy z outside, h = 0.05 / 0.1 / 0.2 | +1.60 / **+2.36** / **+2.89** | same |
| positions identical to the no-layer run | — | **true** (P3.4.5 **PASS**) |

The criterion outcome is **FAIL** (outside h = 0.1 and 0.2), not loosened. The h = R row inside (p = 1) is an exact check: the fraction is 1.

**Follow-up: occupancy vs τ** (deviation δ = fraction/p − 1, σ = √((1 − p)/(Np)); t = 0 control, where the positions are uniform by construction):

| region, h | t = 0 (control) | τ = 4e-2 | τ = 1e-2 | τ = 2.5e-3 |
|---|---|---|---|---|
| inside, 0.05 | −0.02 % (z −0.05) | −0.32 % (−0.93) | −0.65 % (−1.90) | −0.24 % (−0.68) |
| inside, 0.2 | −0.02 % (−0.16) | +0.16 % (+1.07) | −0.26 % (−1.69) | +0.15 % (+0.98) |
| outside, 0.05 | −0.68 % (−1.64) | −0.26 % (−0.63) | +0.67 % (+1.60) | +0.18 % (+0.44) |
| outside, 0.2 | +0.26 % (+1.61) | −0.19 % (−1.19) | +0.46 % (+2.89) | +0.48 % (+3.05) |

![occupancy vs τ](results/phase3/p3_1_3b_density_tau/p3_1_3b_density_tau.png)

**Conclusion.**
- There is no monotone trend with τ, and the t = 0 control fluctuates by the same amount.
- The nested shells (h = 0.05 ⊂ 0.1 ⊂ 0.2) reuse the same spins, so one fluctuation appears in several z values. The χ² over disjoint bins passes in both regions.
- The decisive check is the mean-exponent identity (P3.4.3 below). A real +0.46 % excess within 0.2 µm of the outer surface would shift E[−ln M⊥] by ≈ +0.5 %, about 20σ. The measured value agrees to +5e-6 relative.
- The C-B3 deviations are therefore fluctuations: the density is uniform near curved surfaces at this precision. The layer moves no spin.

## C-B4: conventional θ reference curves on cylinders (P3.1.4)

Scripts [p3_1_4_theta_reference.jl](phase3/p3_1_4_theta_reference.jl) and [p3_1_4b_tau16.jl](phase3/p3_1_4b_tau16.jl) → [results/phase3/p3_1_4_theta_reference/](results/phase3/p3_1_4_theta_reference/) (`reference_curves.csv` is the reference set).

**Design.** These curves are the conventional model's reference on cylinders (a reference, not the truth). Configurations:
- single cylinder inside, R = 0.5, 1, 2;
- square array s = 2.5, R = 1, inside and outside.

θ = 0.003, 0.01, 0.03; D = 3; seeds 1–10, 2·10⁴ spins per seed.

**Method.**
- **Convergence:** τ_a = min(1e-2, 0.03R²/D) against τ_a/4, z = (S_a(50) − S_b(50))/√(SEM_a² + SEM_b²), one test per curve; the 11 readouts are correlated. Pass if |z| ≤ 2.
- **Bulk factorisation:** same seeds with R₂ = 1/80; the per-spin ratio must be e^(−t/80).

**Result.** 14 of 15 curves have |z| ≤ 2 (range −0.98 … +1.19). The exception is array outside θ = 0.03, z = +2.06. Rerun at τ_a/16 = 6.25e-4 ms, it agrees with τ_a/4: 0.05194 against 0.05192 ± 2.2e-5, z = +0.83. τ_a/4 is the reference for that curve (ruling). Bulk factorisation: max per-spin deviation 5.9e-14.

| configuration | S(50), θ = 0.003 / 0.01 / 0.03 (τ_a/4) |
|---|---|
| single R = 0.5 | 0.55643 / 0.14178 / 0.00286 |
| single R = 1 | 0.74596 / 0.37666 / 0.05368 |
| single R = 2 | 0.86372 / 0.61398 / 0.23254 |
| array inside | 0.74596 / 0.37666 / 0.05369 |
| array outside | 0.74359 / 0.37265 / 0.05194 |

**Conclusion.** The conventional reference is τ-converged at these settings. A single cylinder and the array give the same inside curve (the compartments are identical), which is a consistency check.

## Conventional mapping on curved surfaces (P3.1.6)

Script [p3_1_6_mapping.jl](phase3/p3_1_6_mapping.jl) → [results/phase3/p3_1_6_mapping/](results/phase3/p3_1_6_mapping/).

**Design.** MCMR's conventional model multiplies M⊥ by e^(−θ√τ) at each collision, so for every path −ln M⊥ = θ√τ·(number of collisions). On walls the collision rate per unit area is n√(D/(πτ)), which gives the relaxivity ρ = θ√(D/π) (Phase 2, P2.5.0). On a curved surface, finite steps √(2Dτ) could change the collision rate by O(√(Dτ)/R). The fixed-ρ comparisons of 3.7 need the mapping to hold, or a measured correction.

**Method.**

  E[−ln M⊥(t)] = θ√(D/π)·(S/V)·t,  S/V = 2/R inside,  2πR/(s² − πR²) outside.

- Measured at t = 50 ms for 5 configurations, D = 1, 3, 6 (D scan on single R = 1 and array outside), τ = 1e-2 and 2.5e-3, θ = 0.01; seeds 1–10, 2·10⁴ spins.
- z = (mean − expected)/SEM. Pass if |z| ≤ 3; otherwise the measured ratio becomes the mapping.

**Result.** 18 of 18 have |z| ≤ 3 (max |z| = 2.17). The ratio measured/expected lies in 0.99976–1.00013 over √(Dτ)/R = 0.04–0.35.

![mapping vs √(Dτ)/R](results/phase3/p3_1_6_mapping/p3_1_6_mapping.png)

**Conclusion.** The conventional relaxivity ρ = θ√(D/π) holds on these curved surfaces to ≈ 1e-4, with no resolved trend in √(Dτ)/R up to 0.35. Phase 3 comparisons at fixed ρ use θ = ρ√(π/D) unchanged.

## C-B5: runtime on cylinders (P3.1.5)

Script [p3_1_5_runtime.jl](phase3/p3_1_5_runtime.jl) → [results/phase3/p3_1_5_runtime/](results/phase3/p3_1_5_runtime/). Single runs of 10⁴ spins × 1000 steps, 8 threads; ratios carry a load uncertainty of 10–20 %.

| τ (ms) | single, no layer (µs/spin-step) | single, layer | ratio | array, no layer | array, layer | ratio |
|---|---|---|---|---|---|---|
| 0.1 | 0.0738 | 0.0924 | 1.25 | 0.0824 | 0.1312 | 1.59 |
| 1e-2 | 0.0442 | 0.0666 | 1.51 | 0.0495 | 0.0806 | 1.63 |
| 1e-3 | 0.0370 | 0.0557 | 1.50 | 0.0393 | 0.0657 | 1.67 |
| 1e-4 | 0.0350 | 0.0538 | 1.54 | 0.0380 | 0.0603 | 1.59 |

**Many cylinders per cell** (final-review finding; script [p3_1_5b_runtime_many.jl](phase3/p3_1_5b_runtime_many.jl), `runtime_many.csv`).
- Set-up: n × n cylinders (R = 0.7, spacing 2) per cell of side 2n, the same packing density for every n; layer ρ = 0.02, h = 0.1; 2000 spins × 1000 steps.
- The layer exponent loops over every cylinder of the group on every straight piece, at ≈ 1 ns per cylinder even when none of its images is within reach. So the cost grows with the total number of cylinders per cell.

| cylinders per cell | no layer (µs/spin-step) | layer | ratio |
|---|---|---|---|
| 1 | 0.119 | 0.123 | 1.0 |
| 9 | 0.100 | 0.136 | 1.4 |
| 100 | 0.100 | 0.271 | 2.7 |
| 400 | 0.107 | 0.721 | 6.7 |

The final reviewer measured ×17 at 400 cylinders per cell in another configuration (R = 0.7, h = 0.1, a 40 × 40 cell).

**Conclusion.**
- With one cylinder per cell, the step layer costs ×1.25–1.67 per spin-step on cylinders, against ×1.2–1.4 on walls. A reference run (10⁵ spins × 50 ms × 10 seeds at τ = 0.01) takes ≈ 5–7 min.
- The cost grows with the total number of cylinders per cell (×6.7 at 400 here). Packings of hundreds of axons are therefore correct but slow.
- A cell list of layer indices (shell bounding boxes R + h_out, like MCMR's hit grid) would restore a cost proportional to the cylinders within reach. It is deferred to the plan that introduces packings (ruling; tracker P3.6.4).

## Verification of the step layer on cylinders (P3.4.1, P3.4.3)

### P3.4.1: single-spin trace

Script [p3_4_1_trace.jl](phase3/p3_4_1_trace.jl) → [results/phase3/p3_4_1_trace/](results/phase3/p3_4_1_trace/).

**Design.** The configuration is the hardest one for the bookkeeping: a square array s = 2.5, R = 1 with h = 0.3 on both sides, so the outer shells overlap (gap 0.5 < 2h). It also has reflections inside shells and image cylinders.

**Method.** Each step's prediction uses the simulator's reflected path and an independent brute-force summed rate (5 × 5 images, written from D3.1 in the harness):

  M⊥,pred = M⊥,prev · exp(−R₂,bulk·τ − Σ_j dt_j · (1/K)Σ_k ΔR₂(p_j(s_k))),  dt_j = τ·ℓ_j/Σℓ,  K = 10⁴.

Parameters: ρ = 0.03, so ΔR₂ = 0.1176 (inside) and 0.0870 (outside); R₂,bulk = 1/80; τ = 0.01; 100 spins × 500 steps. Midpoint sampling error per boundary crossing ≤ ΔR₂·τ/(2K) ≈ 6e-8.

**Result.** Max |M⊥,sim/M⊥,pred − 1| = 1.0e-7 (tolerance 1e-6) → **PASS**. 9808 steps had a reflection and a non-zero layer exponent; this counts reflected steps that touched a shell, not reflections located inside a shell. 2978 steps passed through an overlap of two outer shells, i.e. had a sample with summed rate 2·ΔR₂,out.

**Conclusion.** The error is at the sampling resolution (≈ 2 boundary crossings × 6e-8). The exact shell fractions, κ, the image loop and the overlap sum agree with brute force on real reflected paths.

### P3.2.5, P3.4.2, P3.4.4, P3.4.7: exact cases at full size

Script [p3_4_exact.jl](phase3/p3_4_exact.jl) → [results/phase3/p3_4_exact/](results/phase3/p3_4_exact/). Seeds 1–10, 10⁵ spins per seed, readouts 0:5:50, τ = 1e-2 (V0b: default τ).

**Design.** These cases are exact in any diffusion regime, so a failure means the code is wrong, not that there is noise.
- **V0b:** ρ = 0 must leave the run bit-identical to no layer.
- **Uniform inner layer:** step, h = R, so κ = 1 − 1/2 = 1/2. Then ΔR₂ = (ρ/R)/κ = 2ρ/R everywhere inside: with ρ = 0.05 and R = 1, every inside spin decays as M⊥ = e^(−(0.1 + R₂,bulk)t) whatever its path. This is the cylinder analogue of the wall cases h = w/2 (step) and h = w (linear), and it checks κ directly.
- **Additivity:** the bulk factor multiplies every step, so the per-spin ratio is e^(−t/80).
- **Side-specific:** an inner-only layer must not touch outside spins, and an outer-only layer must not touch inside spins.

**Method.**
- Per spin: ≤ 1e-10.
- Ensemble: ≤ B = (t/τ + N)·ε = 2.3e-11.
- Side leak: max(1 − M⊥) ≤ 1e-9 (roundoff ruling above).

**Result.**

| check | configuration | result | outcome |
|---|---|---|---|
| V0b inside / outside | square array, ρ = 0, h = 0.3 vs no layer, R₂ = 1/80 | positions and every M⊥ identical | **PASS** / **PASS** |
| uniform, bulk off | single R = 1, inner h = 1, ρ = 0.05 | per spin 9.7e-14, ensemble 6.6e-14 | **PASS** |
| uniform, bulk on (1/80) | same | per spin 1.3e-13, ensemble 1.1e-13 | **PASS** |
| additivity inside / outside | square array ρ = 0.02, h = 0.2 | positions identical; per spin 6.6e-14 / 6.7e-14 | **PASS** / **PASS** |
| inner-only layer | square array | outside spins: max leak 2.4e-12; inside S(50) = 0.1356 | **PASS** |
| outer-only layer | square array | inside spins: max leak 1.2e-11; outside S(50) = 0.1327 | **PASS** |

**Conclusion.** The κ normalisation gives the exact uniform rate 2ρ/R, so ρ is per unit surface area, as D3.1 requires. ρ = 0 runs the unmodified code. Bulk and layer factorise exactly. Each side acts only in its own compartment, up to roundoff at reflection points.

### P3.4.3: mean-exponent identity on cylinders

Script [p3_4_3_identity.jl](phase3/p3_4_3_identity.jl) → [results/phase3/p3_4_3_identity/](results/phase3/p3_4_3_identity/).

**Design.** With ρ defined per unit surface area (D3.1) and a uniform stationary density, every compartment satisfies

  E[−ln M⊥(t)] = t·⟨ΔR₂⟩_compartment = ρ·(S/V)·t

for any h and τ. Each shell integrates to exactly ρ per unit area of its surface:

  ∫ shell ΔR₂ dA = (ρ/(hκ))·π((R ± h)² − R²) = (ρ/(hκ))·2πRh(1 ± h/2R) = 2πRρ.

Overlapping outer layers add their integrals, and validation keeps them out of other cylinders. The identity therefore tests κ, the shell fractions, the image loop, the overlap sum and the density end to end. It does not constrain S = ⟨M⊥⟩, which is where profiles and the conventional model differ.

**Method.** ρ = 0.02, step, t = 50 ms, τ = 1e-2 and 4e-2, seeds 1–10, 10⁵ spins per seed. z = (mean − ρ(S/V)t)/SEM; pass if |z| ≤ 3.

**Result.**

| configuration | τ (ms) | ρ(S/V)t | E[−ln M⊥] ± SEM | z |
|---|---|---|---|---|
| single R = 1, h = 0.2 | 1e-2 | 2.000000 | 2.000019 ± 9.2e-5 | +0.21 |
| | 4e-2 | | 1.999988 ± 8.1e-5 | −0.15 |
| single R = 0.5, h = 0.2 | 1e-2 | 4.000000 | 3.999918 ± 8.3e-5 | −1.00 |
| | 4e-2 | | 4.000063 ± 1.1e-4 | +0.55 |
| square inside, h = 0.2 | 1e-2 | 2.000000 | 2.000105 ± 1.5e-4 | +0.68 |
| | 4e-2 | | 2.000003 ± 1.2e-4 | +0.03 |
| square outside, separated (h = 0.2) | 1e-2 | 2.021352 | 2.021362 ± 6.5e-5 | +0.16 |
| | 4e-2 | | 2.021373 ± 9.4e-5 | +0.23 |
| square outside, overlap (h_out = 0.4 > gap/2) | 1e-2 | 2.021352 | 2.021363 ± 4.8e-5 | +0.23 |
| | 4e-2 | | 2.021311 ± 5.4e-5 | −0.76 |
| rectangular 2.4 × 3.0 outside, overlap in x only (h_out = 0.3) | 1e-2 | 1.548190 | 1.548257 ± 9.6e-5 | +0.70 |
| | 4e-2 | | 1.548305 ± 9.2e-5 | +1.25 |

**Conclusion.**
- 12 of 12 pass with |z| ≤ 1.25. The runs share seeds, so the z values are correlated.
- The agreement is ≤ 6e-5 relative, inside and outside, with and without overlap, for both τ. It holds even at R = 0.5, where κ = 0.8 differs most from 1.
- This confirms the D3.1 normalisation and the overlap sum. It also confirms that the density near the surfaces is uniform at the 1e-4 level (see C-B3).

## P3.4.6: radial continuum benchmark (single cylinder, inside)

Plan [2026-10-25-radial-continuum-h-sweep.md](../docs/superpowers/plans/2026-10-25-radial-continuum-h-sweep.md).
- Solver: [continuum_radial.jl](phase3/continuum_radial.jl); tests: [test_continuum_radial.jl](phase3/test/test_continuum_radial.jl).
- Monte Carlo check: [p3_4_6_continuum_check.jl](phase3/p3_4_6_continuum_check.jl), [p3_4_6b_tau16.jl](phase3/p3_4_6b_tau16.jl), [p3_4_6c_verdict.jl](phase3/p3_4_6c_verdict.jl) → [results/phase3/p3_4_6_continuum_check/](results/phase3/p3_4_6_continuum_check/).

**Design.**
- **Why a continuum.** At fixed ρ every model has the same mean exponent, so the differences between models are of order 1e-4 in S. That is the size of the Monte Carlo timestep error (Review C1 on walls). An exact solution without Monte Carlo and without a timestep is therefore needed as the reference for differences.
- **Why inside only.** Inside a single cylinder the problem is radially symmetric and one-dimensional. Outside, in an array, it is two-dimensional (2D solver optional, not built).

**Method.**

  ∂m/∂t = D (1/r) ∂_r(r ∂_r m) − [R₂,bulk + ΔR₂(r)] m,  0 ≤ r ≤ R,  m(r, 0) = 1,

- **Layer:** ∂_r m(R) = 0, with ΔR₂(r) = (ρ/h)·g((R − r)/h)/κ and κ = ∫g(u)(1 − uh/R)du (D3.1). κ is computed by quadrature, independently of the package's 1 − h/2R.
- **Conventional model:** D ∂_r m(R) = −ρ m(R) (Robin), with ρ = θ√(D/π) (valid on cylinders, P3.1.6).
- **Discretisation:** finite volumes, N = 2000 rings of width dr.
  - ring areas A_i = π(r_{i+1/2}² − r_{i−1/2}²);
  - conductances 2πD r_{i+1/2}/dr;
  - the Robin sink 2πRρ on the outer ring;
  - the step profile averaged exactly over each ring (area overlap).
- **Solution:** with y = √A·m the system is ẏ = −My, where M = A^(−1/2)(K + AV)A^(−1/2) is symmetric tridiagonal with eigenpairs (λ_k, u_k). Then

  S(t) = Σ_k (u_kᵀ√A)² e^(−λ_k t) / Σ_i A_i.

**Exact tests (25 assertions, all pass).**
- Free and bulk decay.
- h = R: the uniform rate 2ρ/R.
- κ for step and linear profiles.
- The initial slope −ρ(S/V) = −2ρ/R for both models.
- **Lowest rate, conventional model:** λ₁ = Dk², with D k J₁(kR) = ρ J₀(kR). Agreement rtol 1e-3, the first-order error of the Robin cell.
- **Lowest rate, step layer:** core m = J₀(kr), λ = Dk². Shell m = I₀(βr) + (I₁(βR)/K₁(βR))K₀(βr), with β = √((q − λ)/D) and q = ρ/(h(1 − h/2R)). The two are matched by continuity of m′/m at r = R − h. Agreement rtol 1e-5, including a strong case ρ = 0.5 (q ≈ 5 ms⁻¹).
- The Jensen bound S ≥ e^(−ρ(S/V)t), including h = 1e-4, below the grid spacing.
- Grid convergence: N = 1000 against 4000 within 1e-6 (layer) and 1e-4 (conventional).

**Ruling.** The exact-case tolerance is 2e-6, not the planned 1e-12. The eigenvalue roundoff is ≈ ε‖M‖ ≈ ε·4DN²/R² ≈ 1e-8 ms⁻¹, so S is off by ≈ t·1e-8. The measured error is 4e-8 at t = 50.

**Monte Carlo against the continuum.**
- **Set-up:** R = 1, D = 3, ρ = 0.02, the conventional model and the step layer h = 0.05, 0.1, 0.2, 0.5, 1.0; seeds 1–10, 2·10⁴ spins in the box (π/4 inside).
- **Pairing:** all models share τ, so their trajectories are identical (checked from the final positions and recorded: true in every pass), and ΔS = S_layer − S_conv is paired per seed.
- **Timesteps:** τ_a = 1e-2, τ_b = 2.5e-3, and τ_c = 6.25e-4 for the conventional model and h = 0.05, 0.1.
- **Criterion:** max over readouts |z| ≤ 3 for S and for ΔS, at the smallest τ run for each model.

| model | ΔS(50) MC, τ = 1e-2 | τ = 2.5e-3 | τ = 6.25e-4 | continuum | max \|z\| (S / ΔS) at smallest τ | outcome |
|---|---|---|---|---|---|---|
| conventional | — | — | — | — | 0.84 / — | **PASS** |
| step h = 0.05 | +1.95e-4 ± 1.8e-5 | +2.54e-5 ± 1.3e-5 | −2.81e-5 ± 7.1e-6 | −5.75e-5 | 2.49 / **4.14** | **FAIL** |
| step h = 0.1 | +3.78e-5 ± 1.6e-5 | −6.66e-5 ± 9.8e-6 | −9.77e-5 ± 8.6e-6 | −1.109e-4 | 1.95 / 1.80 | **PASS** |
| step h = 0.2 | −1.33e-4 ± 1.6e-5 | −1.88e-4 ± 1.6e-5 | — | −2.046e-4 | 1.17 / 1.02 | **PASS** |
| step h = 0.5 | −3.87e-4 ± 1.8e-5 | −3.81e-4 ± 2.7e-5 | — | −3.870e-4 | 0.90 / 0.92 | **PASS** |
| step h = 1.0 (uniform) | −4.54e-4 ± 1.9e-5 | −4.42e-4 ± 3.3e-5 | — | −4.508e-4 | S exact to 4.9e-8 / 0.72 | **PASS** |

**Notes on this table.**
- **Conventional continuum.** The values were computed with the first version of the conventional continuum (Robin sink on the outer ring's centre). The corrected face-centred scheme (final review I3) shifts S_conv(50) by +4.6e-7, which moves the ΔS z values by ≤ 0.07; not rerun.
- **Which summary to read.** The authoritative verdict is `summary_final.json` (p3_4_6c_verdict.jl). `summary.json` of the first pass says FAIL because it predates the τ_c pass and the exact h = 1 comparison.
- **Script provenance.** `p3_4_6b_tau16.jl` was trimmed after it ran: its original verdict section crashed after writing `comparison_tau16.csv`, and the verdict moved to p3_4_6c.
- **Uniform case.** The h = 1.0 S comparison uses t = 50 only. The run is deterministic, and t = 50 is where the roundoff is largest.

The h = 1.0 S comparison is exact rather than a z (ruling): every spin decays at exactly 0.1 ms⁻¹, so its SEM is ≈ 1e-17 and a z would only measure the continuum's roundoff.

![MC vs continuum](results/phase3/p3_4_6_continuum_check/p3_4_6_continuum_check.png)

**Conclusion.**
- **The continuum is validated.** It reproduces the Monte Carlo S(t) of every model within 3 SEM at a converged τ. For h ≥ 0.1 it also reproduces the small paired difference ΔS.
- **Review C1 recurs on cylinders.** At fixed ρ, the default τ = 1e-2 biases the thin-layer difference so much that its sign flips (h = 0.05: +1.95e-4 against −5.75e-5). For h ≤ 0.1 it still matters at τ = 2.5e-3.
- **h = 0.05 still fails at τ = 6.25e-4** (z = 4.1), but it converges towards the continuum: the error from the continuum shrinks from 2.5e-4 to 8.3e-5 to 2.9e-5 for each fourfold reduction in τ (∝ τ^0.8). A linear-in-τ Richardson estimate from the last two passes gives −4.6e-5 (h = 0.05) and −1.081e-4 (h = 0.1, continuum −1.109e-4).
- **Consequence.** Inside single cylinders, fixed-ρ differences use the continuum. Monte Carlo differences for h ≲ 0.1 need τ ≲ 1e-3 or an extrapolation in τ. This feeds P3.6.2: for comparison studies, set τ from the error in the difference, not in R.

## P3.4.8: attenuation vs h at fixed ρ (step)

**Design.**
- **Why fixed ρ.** With ρ fixed, the mean exponent E[−ln M⊥] = ρ(S/V)t is the same for every h (P3.4.3). The h dependence of S is then entirely the Jensen term of Phase 2 (M3):

  S(t) = e^(−ρ(S/V)t)·e^(J(t)),  J = ln E[e^(−X)] + E[X] ≥ 0,  J ≈ Var[X]/2.

- **Limits.** The two ends of the sweep are known: h → 0 tends to the conventional surface model, and h = R inside is uniform (J = 0).
- **Sign.** The direction of the change with h is recorded, not assumed.
- **Range.** h ≥ 0.05 µm for Monte Carlo (user request 2026-10-25: thin layers force τ ∝ h). The continuum covers h/R = 0.005–1 at no extra cost.

### Inside single cylinders (continuum)

Script [p3_4_8_h_sweep_inside.jl](phase3/p3_4_8_h_sweep_inside.jl) → [results/phase3/p3_4_8_h_sweep_inside/](results/phase3/p3_4_8_h_sweep_inside/) (`sweep.csv`, `steps.csv`).
- R = 0.5, 1, 2; D = 3; ρ = 0.02; t = 50 ms. E[X] = 0.02·(2/R)·50 = 4, 2, 1 for R = 0.5, 1, 2.
- Curvature of ln S from a quadratic fit over 0:5:50.
- Grid error: |S(N = 2000) − S(N = 4000)| ≤ 5.2e-08 for the layers. For the conventional model, |S(N = 1000) − S(N = 2000)| ≤ 7e-09, with the Robin condition on the face.
- Correction (final review I3): the first version put the Robin sink on the outer ring's centre, a first-order error of ≈ 5e-7 in S_conv(50). That changed ΔS vs conventional by 10–20 % at h/R ≤ 0.01; the numbers below are corrected.

| h/R | S(50), R = 1 | J(50), R = 0.5 / 1 / 2 | ΔS(50) vs conventional, R = 0.5 / 1 / 2 |
|---|---|---|---|
| 0.005 | 0.135781 | 3.287e-03 / 3.285e-03 / 3.278e-03 | -8.1e-7 / -6.0e-6 / -1.6e-5 |
| 0.02 | 0.135763 | 3.156e-03 / 3.154e-03 / 3.148e-03 | -3.2e-6 / -2.4e-5 / -6.4e-5 |
| 0.05 | 0.135729 | 2.904e-03 / 2.902e-03 / 2.897e-03 | -7.9e-6 / -5.8e-5 / -1.6e-4 |
| 0.1 | 0.135675 | 2.510e-03 / 2.508e-03 / 2.504e-03 | -1.5e-5 / -1.1e-4 / -3.0e-4 |
| 0.2 | 0.135581 | 1.818e-03 / 1.817e-03 / 1.815e-03 | -2.8e-5 / -2.1e-4 / -5.6e-4 |
| 0.3 | 0.135505 | 1.251e-03 / 1.251e-03 / 1.250e-03 | -3.8e-5 / -2.8e-4 / -7.6e-4 |
| 0.5 | 0.135399 | 4.712e-04 / 4.714e-04 / 4.710e-04 | -5.2e-5 / -3.9e-4 / -1.1e-3 |
| 0.7 | 0.135349 | 9.790e-05 / 9.775e-05 / 9.761e-05 | -5.9e-5 / -4.4e-4 / -1.2e-3 |
| 1.0 (uniform) | 0.135335 | ≈ 0 (−2e-7 … −7e-9, roundoff) | -6.1e-5 / -4.5e-4 / -1.2e-3 |
| conventional | 0.135787 | 3.331e-03 / 3.329e-03 / 3.322e-03 | 0 |

![h sweep inside](results/phase3/p3_4_8_h_sweep_inside/p3_4_8_h_sweep_inside.png)

**Result.**
- **Monotone.** All 27 neighbouring steps are resolved (more than 10× the grid error), and in every one S(50) falls as h grows. Every layer has S below the conventional model.
- **Single-exponential.** The curvature of ln S is ≤ 2.6e-9 ms⁻² everywhere.

**Conclusion.**
- **Ordering.** At fixed ρ, spreading the relaxation deeper into the cylinder lowers the spin-to-spin variance of the exponent and so lowers S. In order: the conventional surface sink has the largest J and the highest S; the uniform layer (h = R) has J = 0 and the lowest S, e^(−2ρt/R). This is the same ordering as on walls.
- **J depends on h/R almost only.** At the same ρ, t and D (only R was varied), J is the same function of h/R for R = 0.5, 1 and 2 to within 0.3 %, although E[X] differs by a factor of 4. A small systematic trend remains: J falls slightly as R grows at every h/R, e.g. 3.287 / 3.285 / 3.278e-3 at h/R = 0.005. It is far above the grid error, and is presumably a finite-t or higher-cumulant term. This is consistent with Var[X] ∝ (ρ²t/D)·f(h/R): the rate fluctuation of a spin scales as ρ/R and its correlation time as R²/D, so R cancels. Consequences:
  - ΔS ≈ S·(J − J_conv) grows with R here only because S = e^(−2ρt/R) is larger for larger R;
  - the relative difference ΔS/S, ≈ −3.3e-3 between the uniform layer and the conventional model, is the same for every R at the dimensionless ρ²t/D = 6.7e-3.
- **Observability.** The largest difference is |ΔS(50)| = 1.2e-3 (R = 2, h = R), against 0.01 at SNR 200. Inside single cylinders at this ρ the h dependence is not observable from a single curve. Where it might become observable is P3.7.

### Outside, square arrays (paired Monte Carlo)

Scripts [p3_4_8_h_sweep_outside.jl](phase3/p3_4_8_h_sweep_outside.jl) and [p3_4_8b_outside_tau16.jl](phase3/p3_4_8b_outside_tau16.jl) → [results/phase3/p3_4_8_h_sweep_outside/](results/phase3/p3_4_8_h_sweep_outside/) (`sweep.csv`, `steps.csv`, `steps_tau16.csv`).

**Method.**
- **Geometry:** outer layer only, R = 1, ρ = 0.02, D = 3, in two arrays:
  - s = 2.5 (gap 0.5): h = 0.05–0.5; layers overlap for h > 0.25;
  - s = 3.0 (gap 1.0): h = 0.05–1.0; layers overlap for h > 0.5.
- **Why Monte Carlo.** There is no continuum outside.
- **Pairing.** All h and the conventional model of one array share τ, so their trajectories are identical (checked from the final positions and recorded: true). Neighbouring h are compared by paired per-seed differences d = S(50; h_{k+1}) − S(50; h_k).
- **Timesteps:**
  - τ_a = 1e-2 and τ_b = 2.5e-3: seeds 1–10, 2·10⁴ spins per seed;
  - τ_c = 6.25e-4 for the conventional model and h ≤ 0.2: 10⁴ spins per seed.

  The thin end was rerun at τ_c because P3.4.6 showed a Review C1 bias for h ≤ 0.1.
- **Resolved step:** |z| > 3. The sign is recorded, not assumed.

**Result.**
- E[−ln M⊥] = ρ(S/V)t within noise for every h and τ: 2.0213 at s = 2.5 and 1.0725 at s = 3.0, a re-check of P3.4.3.
- Each neighbour step below is reported at the smallest τ run for it: τ_c for the steps involving h ≤ 0.2, τ_b otherwise.
  - Thin end: the shift from τ_b to τ_c is ≤ 2.1σ. That test has little power, because τ_c used half the spins and its shift z treats the two runs as independent although they share initial positions.
  - Inside, at the same τ_c and h = 0.05, about half of the true difference was still timestep bias (P3.4.6). The thin-end steps outside are therefore probably biased towards zero, by up to tens of percent.
  - Thicker steps: the τ_a → τ_b shift is ≤ 1.5σ, except s = 2.5, 0.4 → 0.5 (2.6σ).

| array | step | paired ΔS(50) ± SEM | z | τ |
|---|---|---|---|---|
| s = 2.5 | conventional → h 0.05 | −2.34e-5 ± 9.8e-6 | −2.4 | τ_c |
| | 0.05 → 0.1 | −4.77e-5 ± 1.2e-5 | −3.9 | τ_c |
| | 0.1 → 0.2 | −9.32e-5 ± 1.7e-5 | −5.4 | τ_c |
| | 0.2 → 0.3 (overlap begins) | −4.69e-5 ± 8.9e-6 | −5.3 | τ_b |
| | 0.3 → 0.4 | −2.75e-5 ± 4.5e-6 | −6.1 | τ_b |
| | 0.4 → 0.5 (= gap) | −2.29e-5 ± 4.4e-6 | −5.2 | τ_b |
| s = 3.0 | conventional → h 0.05 | −4.36e-5 ± 2.4e-5 | −1.8 | τ_c |
| | 0.05 → 0.1 | −5.49e-5 ± 2.4e-5 | −2.3 | τ_c |
| | 0.1 → 0.2 | −1.19e-4 ± 1.4e-5 | −8.5 | τ_c |
| | 0.2 → 0.3 | −1.37e-4 ± 1.0e-5 | −13.3 | τ_b |
| | 0.3 → 0.5 | −1.67e-4 ± 2.2e-5 | −7.8 | τ_b |
| | 0.5 → 0.7 (overlap) | −8.24e-5 ± 1.4e-5 | −6.0 | τ_b |
| | 0.7 → 1.0 (= gap) | −4.27e-5 ± 6.8e-6 | −6.3 | τ_b |

Total difference from the conventional model at τ_b: S(50) − S_conv(50) = −2.57e-4 at s = 2.5, h = 0.5 (−0.19 % of S) and −6.52e-4 at s = 3.0, h = 1.0 (−0.19 %). These τ_b totals contain the thin-end step at τ_b. With the τ_c thin-end steps they are ≈ −0.20 to −0.21 %, and the residual bias at the thin end would make them slightly larger again.

**Curvature of ln S outside** ([p3_4_8c_outside_curvature.jl](phase3/p3_4_8c_outside_curvature.jl), `curvature.csv`; added after the final review, I2).
- Method: per seed, a fit ln S = a + bt + ct² over 0:5:50 at τ = 1e-2, 2·10⁴ spins.
- Result: c is zero within noise for the conventional model and every h in both arrays: |c| ≤ 1.9e-7 ms⁻², |z| ≤ 1.91.
- Conclusion: outside, as inside, the curves are single-exponential at this precision.

![h sweep outside](results/phase3/p3_4_8_h_sweep_outside/p3_4_8_h_sweep_outside.png)

**Conclusion.**
- **Same ordering as inside, resolved from h = 0.05 on.** Outside, S falls with h through every layer step. 10 of the 13 steps are resolved beyond 3σ, and every layer step from 0.1 → 0.2 on is resolved.
- **The outside position of the conventional model is consistent but not resolved.** conventional → 0.05 has z = −2.4 (s = 2.5) and −1.8 (s = 3.0), and at τ_b the s = 3.0 step was even positive (+2.6e-5). The step 0.05 → 0.1 at s = 3.0 has z = −2.3. That the conventional model has the highest S outside rests on the inside result (continuum) and on these unresolved steps, whose residual timestep bias works towards zero.
- **Overlap** changes nothing qualitatively: the decrease simply continues more slowly once the layers overlap.
- **The relative effect is smaller than inside.** −0.19 % of S at h = gap, against −0.33 % for the uniform layer inside. Outside, even h = gap leaves the cell corners far from every surface, so the rate never becomes uniform.
- **Timestep.** At the default τ = 1e-2 the thin-end steps are biased, reversing the sign of conventional → 0.05: +1.03e-4 and +1.75e-4. At τ_c = 6.25e-4 they agree with τ_b within 2.1σ, but by the inside evidence a residual bias of the order of the step itself may remain.

**Overall conclusion (P3.4.8).**
- At fixed ρ on cylinders, S(t) decreases monotonically with h: resolved at every step inside (continuum), and from h = 0.05 on outside (paired MC). It lies between two limits:
  - the conventional surface model, the h → 0 limit, has the largest J and the highest S;
  - the most spread-out layer has the lowest S.
- All the differences are ≤ 0.33 % of S at ρ = 0.02 (ρ²t/D = 6.7e-3), far below single-measurement noise.
- Inside, J depends on h/R almost alone (to 0.3 %, at fixed ρ, t and D). The comparison over R, D and ρ, and the apparent T₂, belong to P3.7.

## Linear and exponential profiles on cylinders (P3.3.1–P3.3.4, D3.3, D3.5; P3.6.3)

Plan [2026-11-06-cylinders-linear-exponential.md](../docs/superpowers/plans/2026-11-06-cylinders-linear-exponential.md), Group 1 (package core). Code: [layers.jl](../src/geometries/internal/layers.jl) (`quad_adaptive`, `disk_interval`, `shell_profile_integral`, `cylinder_side_exponent`) and [fix.jl](../src/geometries/user/fix.jl) (`cylinder_layers`, `surface_gaps`). Tests: [test_layer_cylinders.jl](../test/test_layer_cylinders.jl).

**Design.**
- **No closed form along a piece.** For the step profile, the exponent of a piece follows exactly from the time spent in the shell. For non-step profiles, d(s) = R ∓ √(as² + bs + q) is not linear in s, so g(d(s)/h) has no closed form for the exponential profile. The linear profile does have one, but writing a closed form per profile does not scale.
- **D3.3 (decided):** one numerical path for all non-step profiles, with the linear closed form as an independent test. The step profile keeps its exact closed form.
- **Why Gauss–Legendre.** The support of a piece is cut at the exact shell roots (from the step code's quadratic) and at the closest approach to the axis, s* = −b/2a. On each part r(s) is smooth and monotone. Gauss–Legendre therefore converges very fast and needs no special handling of the kink at r = 0, which only the inner exponential reaches.
- **Exponential support (D3.5).**
  - Inside, the exponential reaches the axis (support R, c = R/h).
  - Outside, its support is the new field `layer_cutoff`. The default is the smallest surface gap to any neighbouring cylinder or image, the wall analogue (cutoff = spacing).
  - A cylinder with no neighbours and an outer exponential must set `layer_cutoff`; validation errors otherwise.
- **Support, not h.** Image reach, the outer-reach validation and the timestep bound all use the support L: h for step and linear, the cutoff for the exponential. An exponential with a short decay length but a long support is therefore summed wherever it reaches.

**Method.** For one side of a cylinder, over a straight piece of duration dt,

  E = (ρ/h)/κ · dt · ∫₀¹ g(d(s)/h) ds,  d(s) = R − r(s) inside,  r(s) − R outside,  r(s)² = as² + bs + q,

with:
- κ = 1 ∓ c̄/R, where c̄ = h/2 (step), h/3 (linear), h(1 − (c + 1)e^(−c))/(1 − e^(−c)) (exponential, c = L/h);
- the integral over the support taken as the sum of its parts {s : r_lo ≤ r(s) ≤ r_hi}, bounded by the roots of r(s)² = r_lo² and r(s)² = r_hi², each further split at s*;
- each part integrated by 8-point Gauss–Legendre, the nodes computed with Golub–Welsch (no new dependency), halving until the whole-interval and two-halves estimates agree within LAYER_QUAD_TOL·(length). The default LAYER_QUAD_TOL is 1e-10.

**Verification (package tests, all pass; full suite without plots 3553/3553).**
- **Linear against its closed form:** ∫g ds = ∫(2 − 2(R − r)/h) ds, with ∫√(as² + bs + q) ds = (2as + b)√Q/(4a) + (4aq − b²)/(8a^(3/2))·ln|2√a√Q + 2as + b|. Agreement within 1e-10 on chords, pieces through the core, and pieces inside the shell.
- **Wall regression:** along a radial piece, d(s) is exactly linear, so the cylinder integral must equal the wall closed form (`side_exponent`) for linear and exponential. Agreement within rtol 1e-9.
- **Brute force:** 300 random pieces × 4 cases (linear and exponential, inside and outside), 2·10⁵ samples each; max error < 2e-4, the sampling bound.
- **Through the axis:** the inner exponential reaches r = 0, and a radial piece through it gives exactly h, since ∫₀¹ g ds = h∫g du.
- **Tangent and grazing pieces:** tangent to r = R, to r = R − h, and to the exponential's support edge. All finite and correct.
- **Simulation:** a hand-computed radial step with a reflection, inside, for linear (κ = 1 − h/3R) and exponential (∫g dt from e^(−d/h) on each straight piece), to rtol 1e-9 and 1e-8.
- **Layer exponent against an independent brute-force summed rate:** square array with overlapping outer layers, both profiles.
- **Support-based reach:** an outer exponential with h = 0.05 and support 0.5 contributes at 0.4 from the surface. The timestep bound counts the 4 neighbours whose supports intersect.
- **Side-specific and bulk:** both hold for the new profiles.
- **Validation:** shapes, default and explicit cutoff, cutoff larger than the gap, an isolated cylinder without a cutoff, h_inner > R for linear. JSON and CLI round trips.
- **Step results bit-identical:** in the reference run (square array, 1000 spins, 20 ms) the summed M⊥ is 447.8892123934866 both at the base commit (3bcd2c6) and after this plan.

**Tolerance and cost (P3.3.4, P3.6.3).** Script [p3_3_quadrature.jl](phase3/p3_3_quadrature.jl) → [results/phase3/p3_3_quadrature/](results/phase3/p3_3_quadrature/).
- Accuracy: 2·10⁴ random pieces, compared with a tolerance-1e-14 reference.
- Runtime: square array s = 2.5, ρ = 0.02, h = 0.2 on both sides, τ = 1e-2, 10⁴ spins × 1000 steps.

Accuracy is measured against the same algorithm at tol 1e-14 (self-convergence). The linear profile is also checked against its independent closed form ([p3_3b_linear_closed_form.jl](phase3/p3_3b_linear_closed_form.jl), 2·10⁴ pieces): max error 8.2e-14 for h = 0.3 and 9.3e-11 for h = R, where pieces pass through the axis and the shell reaches r = 0.

The values below were measured after the allocation fix (Group 1 review); the first version allocated 200–400 B per piece and was ≈ 1.2–2× slower. Times are single-shot wall times on a shared machine (±20 %).

| case | max error vs tol-1e-14 reference at tol = 1e-4 / 1e-6 / 1e-8 / 1e-10 / 1e-12 | ns per piece at 1e-6 / 1e-10 |
|---|---|---|
| linear, inside (h = 0.3) | ≤ 6.1e-16 at every tolerance | 243 / 240 |
| linear, outside | ≤ 6.7e-16 at every tolerance | 244 / 244 |
| exponential, inside (h = 0.2, support R) | 2.3e-7 / 2.3e-7 / 2.1e-8 / 3.2e-11 / 8.5e-13 | 613 / 1206 |
| exponential, outside (h = 0.15, support 0.5) | ≤ 2.2e-14 / … / 3.9e-16 | 580 / 622 |

| profile | µs per spin-step (tol 1e-10) | at tol 1e-6 | relative to step |
|---|---|---|---|
| step | 0.089 | — | 1 |
| linear | 0.106 | 0.107 | ×1.2 |
| exponential | 0.221 | 0.222 | ×2.5 |

**Conclusion.**
- **Linear is converged at the first level.** Its integrand, a linear function of √(quadratic), is so smooth on each part that 8-point Gauss–Legendre needs no bisection. The result agrees with the independent closed form to 8e-14, or 9e-11 for pieces through the axis when h = R.
- **The exponential's error is controlled by the tolerance.** At the default 1e-10 the worst error, 3.2e-11, occurs for inner pieces through the steep part near the surface. Pieces outside are exact to roundoff.
- **What 1e-10 means for the signal.** The error of a piece's exponent is ≤ (ρ/h)/κ·dt·3.2e-11, i.e. ≲ 1e-10 of the per-step exponent ΔR₂(0)·τ. That is far below the timestep error of R (≈ 1e-3) and the Monte Carlo noise.
- **No reason to relax it.** In the simulation a looser tolerance saves no measurable runtime (0.221 against 0.222 µs per spin-step). Diffusion-length pieces are short and converge at the first level, so the cost is dominated by the fixed 3 × 8 evaluations per part. Only long random pieces through the steep inner exponential get cheaper at 1e-6 (1206 → 613 ns), and such pieces are rare in the simulation. So 1e-10 stays (P3.6.3).
- **Cost.** The non-step profiles cost ×1.2 (linear) and ×2.5 (exponential) per spin-step relative to step in this array, on top of the step layer's own ×1.6 over no layer (C-B5). After the fix the quadrature path allocates nothing (tested).

**Rulings (Group 1).**
- Two old assertions from the step-only stage ("linear must be rejected", "centroid(:linear) throws") were replaced by the new behaviour.
- **Unrelated test failure.** One full test-suite run had a single failure in `test_known_sequences.jl`, "PGSE in realistic scanner with free diffusion": 37677.5 against 36787.9 ± 2 %.
  - That test is not seeded and uses no geometry, so the layer code is never entered.
  - It passed 3/3 times when run on its own, and in every earlier full run on this branch.
  - A second full run passed 3553/3553, including that test.

## Verification of linear and exponential on cylinders (Group 2: P3.4.1, P3.4.3, P3.4.6 per profile)

### Radial continuum for all profiles (exact ring averages)

Solver [continuum_radial.jl](phase3/continuum_radial.jl); tests [test_continuum_radial.jl](phase3/test/test_continuum_radial.jl).

**Design.** The first version averaged non-step profiles over each ring by sub-sampling. Layers thinner than a ring then lost part of their integral, ∫ΔR₂ dA = 0.996·2πRρ at h = 1e-4 (a deferred minor of the P3.4.6 review). Every profile has an elementary antiderivative of g((R − r)/h)·r, so the ring averages can be exact.

**Method.** The ring integral ∫_lo^hi g((R − r)/h)·r dr is clipped to the support [R − L, R]:
- step: r²/2;
- linear: 2[(1 − R/h)r²/2 + r³/(3h)];
- exponential: h·e^((r − R)/h)(r − h)/(1 − e^(−c)), with c = L/h and support L = R inside (D3.5).

The ring rate is V_i = (ρ/(hκ))·2π·∫/A_i, with κ = ∫_{R−L}^R g·r dr/(hR) from the same exact integral.

**Ruling.** κ was first computed by a midpoint rule over u ∈ [0, c]. That failed for thin exponential layers, where c = R/h = 10⁴ and a step of 0.1 in u gives an error of ≈ 1e-4. κ now comes from the exact ring integral. The κ tests against the closed forms (step 1 − h/2R, linear 1 − h/3R, exponential 1 − c̄/R) remain the independent check.

**Tests (156/156 in the research suite after the Group 2 review):**
- **ring_integral against an independent fine midpoint rule** of g((R − r)/h)·r on 90 random sub-intervals inside the support (linear h = 0.3; exponential h = 0.2 and 1.5), rtol 1e-9. This test was added after the Group 2 review: since κ comes from the same ring integral, ∫ΔR₂ dA = 2πRρ and the initial slope −2ρ/R hold by construction for any antiderivative. A deliberately wrong exponential antiderivative passed them and failed this test (60 failures).
- ∫ΔR₂ dA = 2πRρ for linear and exponential, including h = 1e-4 (below the grid spacing). This is by construction now, so it only guards the ring summation. The tolerance is 1e-8 because the test extracts V as a difference of diagonals ≈ D/dr².
- κ of the exponential with support R.
- The initial slope −2ρ/R and the Jensen bound.
- Grid convergence: N = 1000 against 2000 within 2e-8.
- The step results are unchanged, including the exact Bessel eigenvalue tests.

### Harness for any profile

[harness.jl](phase3/harness.jl): `g_profile`, `side_prefactor` (ρ/h)/κ with κ by quadrature, and `brute_rate_profile` (summed rate over 5 × 5 images). All are written independently of the package.

### P3.4.1 and P3.4.3 for linear and exponential

Scripts [p3_4_1_trace.jl](phase3/p3_4_1_trace.jl) and [p3_4_3_identity.jl](phase3/p3_4_3_identity.jl) with `SHAPE` → `results/phase3/p3_4_1_trace_<shape>/`, `results/phase3/p3_4_3_identity_<shape>/`. Configurations are as for the step profile:
- **trace:** square array s = 2.5, R = 1, h = 0.3, ρ = 0.03. For the exponential, h is the decay length, the support is R inside and the default cutoff 0.5 (the gap) outside;
- **identity:** ρ = 0.02, t = 50 ms, τ = 1e-2 and 4e-2. For the exponential the outer cutoff is the gap: 0.5 in the square array, 0.4 in the rectangular one. Single cylinders set `layer_cutoff = 1.0` for their irrelevant outer side (ruling).

**Method.**
- **Trace:** each step's M⊥ is predicted from the simulator's reflected path. The exponent of each piece is τ·ℓ_j/Σℓ times the mean of the brute-force summed rate over 10⁴ points.
- **Identity:** E[−ln M⊥] = ρ(S/V)t, which is exact for any profile because ρ is per unit surface area (D3.1) and every shell integrates to 2πRρ.

| profile | trace: max \|M⊥,sim/M⊥,pred − 1\| | reflected steps with a layer exponent / steps through an overlap (gates, see note) | identity: rows passing, max \|z\| |
|---|---|---|---|
| step (P3.4.1/3, for reference) | 1.0e-7 | 9808 / 2978 | 12/12, 1.25 |
| linear | 6.7e-12 | 9808 / 2946 | 12/12, 1.17 |
| exponential | 2.4e-8 | 9808 / 12599 | 12/12, 1.27 |

**Note on the count columns.** Both counts are gates (> 0), not comparable statistics.
- "Reflected steps with a layer exponent" counts every reflecting step with a non-zero exponent: 9808 for every shape. For the exponential (support R) this is every reflection.
- The overlap count for linear and exponential samples every 100th point. Linear and step share paths and supports, so its 2946 against step's 2978 is a detector difference.

**Conclusion.**
- **Trace.** Both profiles agree with brute force on real reflected paths; the tolerance is 1e-6.
  - **Linear, 7e-12:** far below the step's 1e-7, because the linear profile has no jump, so midpoint sampling is nearly exact.
  - **Exponential, 2.4e-8:** comes from the jump at the outer truncation r = R + 0.5 (≈ 7e-6 relative per crossing at 10⁴ samples). That edge coincides with a neighbour's surface only along the axes; elsewhere spins cross it. Smooth near-surface pieces sample to ≈ 3e-10.
- **Overlap.** The exponential passes through overlapping outer layers in 4× more steps than the others, because its default support reaches the neighbouring surfaces.
- **Identity.** It holds for both profiles in every configuration, inside and outside, at both τ.
  - For linear, the outside rows include separated (h = 0.2 < gap/2) and overlapping layers.
  - For the exponential every outside row overlaps, since its default support is the full gap. The row named "square_sep_out" is therefore not separated for the exponential.
  - This confirms the κ normalisation, the supports and the quadrature path end to end.

### P3.4.6 for linear and exponential: Monte Carlo against the continuum (inside)

Script [p3_4_6d_profiles.jl](phase3/p3_4_6d_profiles.jl) → [results/phase3/p3_4_6d_profiles/](results/phase3/p3_4_6d_profiles/).

**Design.** The comparison is the same as for the step profile (P3.4.6): single cylinder R = 1, D = 3, ρ = 0.02, with the conventional model at the same ρ.
- **Layers:** linear h = 0.3 and 0.6 (centroids 0.1 and 0.2); exponential h = 0.1 and 0.2 (support R; centroids ≈ 0.1 and 0.2).
- **No thin layers** (user, 2026-10-25): the largest surface rate, 0.23 ms⁻¹, allows τ_a = 1e-2.
- **Pairing:** all models share τ (trajectories identical, checked), and ΔS = S_layer − S_conv is paired per seed.
- **Timesteps:** τ_a = 1e-2 and τ_b = 2.5e-3; seeds 1–10, 2·10⁴ spins in the box.
- **Criterion:** at τ_b, max over readouts |z| ≤ 3 for S and for ΔS.

**Result.**

| model | S(50) MC at τ_b ± SEM | continuum | max \|z\| S | ΔS(50) MC τ_a | ΔS(50) MC τ_b | ΔS continuum | max \|z\| ΔS at τ_b |
|---|---|---|---|---|---|---|---|
| conventional | 0.135777 ± 3.3e-5 | 0.135787 | 0.73 | — | — | — | — |
| linear h = 0.3 | 0.135585 ± 2.2e-5 | 0.135576 | 1.02 | −1.43e-4 ± 1.5e-5 | −1.92e-4 ± 1.5e-5 | −2.11e-4 | 1.27 |
| linear h = 0.6 | 0.135443 ± 1.3e-5 | 0.135444 | 0.83 | −3.28e-4 ± 1.6e-5 | −3.34e-4 ± 2.3e-5 | −3.43e-4 | 0.81 |
| exponential h = 0.1 | 0.135593 ± 2.5e-5 | 0.135578 | 0.90 | −1.28e-4 ± 1.5e-5 | −1.84e-4 ± 1.2e-5 | −2.08e-4 | 2.09 |
| exponential h = 0.2 | 0.135473 ± 1.7e-5 | 0.135470 | 0.93 | −2.85e-4 ± 1.5e-5 | −3.04e-4 ± 1.9e-5 | −3.16e-4 | 0.79 |

→ **PASS** at τ_b for every model.

**Conclusion.**
- **The continuum is validated for linear and exponential** inside a cylinder, for S and for the small difference from the conventional model. It stands in for Monte Carlo inside single cylinders for every profile.
- **The Review C1 timestep effect recurs.**
  - **At the default τ = 1e-2:** |ΔS| is underestimated by 32 % (linear h = 0.3), 4 % (linear h = 0.6, not significant, z 1.9), 38 % (exponential h = 0.1) and 10 % (exponential h = 0.2). S itself also fails there for every model, including the conventional one (max |z| 2.05–5.37).
  - **At τ = 2.5e-3:** the criterion passes, but every layer's |ΔS| is still below the continuum, by 9 %, 2.5 %, 12 % and 4 % (z 1.27, 0.37, 2.09, 0.65; all the same sign). The shrinkage from τ_a is roughly proportional to τ.
  - **What this means.** The residual bias is not resolved at τ = 2.5e-3, consistent with the step profile at h = 0.05 still failing at 6.25e-4 (P3.4.6). Monte Carlo differences of surface-concentrated profiles at τ = 2.5e-3 should be taken as ≈ 10 % low. The continuum is the reference inside.
- **Centroid comparison.** The two profiles have similar ΔS at similar centroids, but this is only an indication:
  - c̄ ≈ 0.1: −2.11e-4 (linear) and −2.08e-4 (exponential);
  - c̄ ≈ 0.2: −3.43e-4 and −3.16e-4, 8 % apart. The exponential h = 0.2 with support R has c̄ = 0.193, not 0.2.

  It rests on two continuum points per profile; the systematic comparison is P3.4.8 (Group 3).

## Attenuation vs h and the centroid comparison for all profiles (Group 3: P3.4.8 per profile)

### Inside single cylinders (continuum)

Script [p3_4_8d_profiles_inside.jl](phase3/p3_4_8d_profiles_inside.jl) → [results/phase3/p3_4_8d_profiles_inside/](results/phase3/p3_4_8d_profiles_inside/) (`sweep.csv`, `steps.csv`, `centroid.csv`).

**Design.**
- **Question.** On walls (Phase 2, P2.5.6/7 and Review C1), the three profiles differed by up to 69 % of the effect at the same h, but agreed within 11 % at the same centroid c̄ = h∫u g du. Does that hold on cylinders, where κ and the shell geometry depend on curvature?
- **Why the continuum.** Inside single cylinders it gives exact S(t) for every profile (validated in P3.4.6), so the comparison has neither noise nor timestep bias.

**Method.**
- Fixed ρ = 0.02, D = 3, t = 50 ms, R = 0.5, 1 and 2.
- h sweep: step and linear h/R = 0.005–1; exponential h/R = 0.005–5 (support R).
- ΔS(50) = S − S_conv. Neighbouring h steps count as resolved when they exceed 10× the grid error.
- **Same centroid:** for c̄/R = 0.02, 0.05, 0.1 and 0.2, h is chosen per profile so that it has that centroid: step 2c̄, linear 3c̄, exponential by bisection with support R.
- **Spread:** max over profile pairs of |ΔS_a − ΔS_b|/max(|ΔS_a|, |ΔS_b|). For comparison the same is computed at the same h (h = 2c̄).

**Result.**
- **Monotone:** all 84 neighbouring h steps (3 profiles × 3 radii) are resolved, and S(50) falls with h in every one.
- **Thick limits (R = 1):**
  - step at h = R is uniform: ΔS = −4.51e-4;
  - linear at h = R is not uniform (unlike step at h = R): ΔS = −4.21e-4;
  - exponential approaches the uniform limit as h grows: −4.40e-4 at h = R, −4.48e-4 at 2R, −4.51e-4 at 5R.

| R | c̄/R | ΔS(50) step / linear / exponential | spread, same centroid | spread, same h |
|---|---|---|---|---|
| 1 | 0.02 | −4.67e-5 / −4.88e-5 / −5.15e-5 | 0.093 | 0.663 |
| 1 | 0.05 | −1.11e-4 / −1.16e-4 / −1.19e-4 | 0.065 | 0.618 |
| 1 | 0.10 | −2.05e-4 / −2.11e-4 / −2.08e-4 | 0.026 | 0.527 |
| 1 | 0.20 | −3.42e-4 / −3.43e-4 / −3.22e-4 | 0.060 | 0.339 |
| 0.5 | 0.02–0.20 | (scaled by S, same pattern) | 0.094 / 0.065 / 0.026 / 0.059 | 0.663 / 0.618 / 0.527 / 0.339 |
| 2 | 0.02–0.20 | (scaled by S, same pattern) | 0.093 / 0.065 / 0.026 / 0.060 | 0.663 / 0.618 / 0.528 / 0.339 |

![profiles inside](results/phase3/p3_4_8d_profiles_inside/p3_4_8d_profiles_inside.png)

**Conclusion.**
- **The centroid collapse of Phase 2 holds on cylinders.** At fixed ρ, the three profiles' differences from the conventional model agree within 2.6–9.4 % at the same centroid, against 34–66 % at the same h.
  - To first order, the centroid is the quantity that S responds to: how far, on average, the relaxation sits from the surface.
  - The remaining spread is largest for thin layers (9 % at c̄/R = 0.02, where the exponential is the most negative) and at c̄/R = 0.2 (6 %, where the exponential is the least negative). It reflects the profiles' different spread around the centroid.
- **The spreads are the same for R = 0.5, 1 and 2,** consistent with the h/R scaling of J found for the step profile (P3.4.8). The relative comparison depends on c̄/R, not on R itself.
- **S falls monotonically with h at fixed ρ for every profile,** as for the step profile. All lie between the conventional model (h → 0) and the uniform limit.

### Outside, square array s = 3.0 (paired Monte Carlo)

Scripts:
- [p3_4_8e_profiles_outside.jl](phase3/p3_4_8e_profiles_outside.jl) → `results/phase3/p3_4_8e_profiles_outside_<shape>/`;
- thin-end pass [p3_4_8g_outside_tau16.jl](phase3/p3_4_8g_outside_tau16.jl) (`steps_tau16.csv`, `vs_conv_tau16.csv`);
- centroid analysis [p3_4_8h_outside_centroid.jl](phase3/p3_4_8h_outside_centroid.jl) → [results/phase3/p3_4_8h_outside_centroid/](results/phase3/p3_4_8h_outside_centroid/).

**Design.**
- **Geometry and h.** Outer layer only, R = 1, gap 1.0, ρ = 0.02, the same seeds, spins and τ as the step run at s = 3.0 (P3.4.8), so all three profiles can be compared:
  - linear h = 0.1–1.0, starting at 0.1 because its surface rate 2ρ/(hκ) at 0.05 would force τ < 1e-2;
  - exponential h = 0.05–1.0, with support equal to the gap, so the outer layers always overlap.
- **Timesteps.** Each profile is paired with its own conventional run at τ = 1e-2 and 2.5e-3. There is no continuum outside, and inside, |ΔS| of surface-concentrated profiles was still ≈ 10 % low at τ = 2.5e-3 (Group 2). So a third timestep, τ = 6.25e-4, was run at the thin end: the conventional model plus the 3 thinnest h per profile, 10⁴ spins per seed (option B, user decision 2026-11-07).

**Result: monotonicity.**
- **The identity holds.** E[−ln M⊥] = ρ(S/V)t = 1.0725 within noise for every h at τ = 1e-2 and 2.5e-3. The τ = 6.25e-4 pass did not record it.
- **Linear.** All layer steps from 0.1 → 0.2 on are resolved and negative at τ = 2.5e-3 (|z| 7.7–15.3). Conventional → 0.1 is not resolved at τ = 2.5e-3 (z −1.4; at τ = 1e-2 it even has the wrong sign, +8.9e-5), but is resolved at τ = 6.25e-4: −7.0e-5 ± 2.2e-5 (z −3.2).
- **Exponential.** All 6 steps are resolved and negative at τ = 2.5e-3, including conventional → 0.05 (z −5.5). The thin-end steps are confirmed at τ = 6.25e-4 (z −5.3, −19.1, −9.7).
- **Timestep shifts at the thin end.**
  - **Against the conventional model:** the shifts of ΔS from τ = 2.5e-3 to 6.25e-4 are within 1.6σ, and the SEM at τ = 6.25e-4 is 10–31 % of the value.
  - **Between neighbouring h (paired):** the exponential's thin-end steps shift by about 3σ: 0.05 → 0.1 from −1.64e-4 to −1.22e-4 (z −2.99), and 0.1 → 0.2 from −1.97e-4 to −1.32e-4 (z −3.14). That is 26–33 % of each step, so a τ dependence of the exponential thin end is resolved at ≈ 3σ.
  - **Direction:** |ΔS| gets smaller as τ decreases, the opposite of the inside C1 bias. Whether τ = 6.25e-4 is converged is not established.

**Result: the centroid comparison outside.** Centroids: step h/2, linear h/3, exponential h(1 − (c + 1)e^(−c))/(1 − e^(−c)) with c = 1/h. ΔS(50) against the conventional model uses τ = 6.25e-4 at the thin end and τ = 2.5e-3 otherwise, interpolated in c̄.

| c̄ (µm) | ΔS step | ΔS linear | ΔS exponential | spread (same c̄) |
|---|---|---|---|---|
| 0.05 | −9.85e-5 | −1.12e-4 | −1.26e-4 | 0.22 |
| 0.10 | −2.18e-4 | −2.37e-4 | −2.48e-4 | 0.12 |
| 0.15 | −3.59e-4 | −3.68e-4 | −3.19e-4 | 0.13 |
| 0.20 | −4.43e-4 | −4.65e-4 | −3.96e-4 | 0.15 |
| 0.30 | −5.68e-4 | −5.75e-4 | −5.74e-4 | 0.01 |

**Sensitivity to the timestep.** With the exponential's thin-end values taken at τ = 2.5e-3 instead of 6.25e-4, the spread at c̄ = 0.05 / 0.1 / 0.15 / 0.2 becomes 0.12 / 0.18 / 0.03 / 0.06, against 0.22 / 0.12 / 0.13 / 0.15 above.

**Same h.** At the same h the spread is much larger: 0.72, 0.60, 0.57, 0.32 and 0.07 at h = 0.1, 0.2, 0.3, 0.5 and 1.0 (`same_h_outside.csv`).

![outside centroid](results/phase3/p3_4_8h_outside_centroid/p3_4_8h_outside_centroid.png)

**Conclusion (outside).**
- **Ordering.** For linear and exponential, S decreases with h from the conventional model on, and every step is resolved (the thin end at τ = 6.25e-4). For the step profile outside, the thin end stays unresolved even at τ = 6.25e-4 (conventional → 0.05: z −1.8; 0.05 → 0.1: z −2.3; see P3.4.8).
- **The centroid collapse is approximate outside, and its size is not pinned down.**
  - At the same h the profiles differ by 7–72 %; at the same centroid by 1–22 %.
  - At c̄ = 0.15–0.2 the same-centroid spread is 13–15 % with the exponential's τ = 6.25e-4 thin-end values, but 3–6 % with its τ = 2.5e-3 values. The exponential thin end shows a ≈ 3σ τ dependence (above), and no SEM of the spread was computed, so which value is right is undetermined.
  - So the collapse outside is consistent with the inside result (3–9 %) but not established at that precision.
- **Possible real reasons for a looser collapse.**
  - Outside, the exponential's support always overlaps its neighbours', so its profile is not a single-surface profile.
  - The planar centroid h∫u g du ignores the r/R area weighting of a curved shell.
- **What would resolve it.** Per-seed values and a paired cross-profile SEM (the profiles share trajectories at each τ), more spins at τ = 6.25e-4 or smaller, or an outside 2D continuum.
- **Timestep.** Outside fixed-ρ differences at the default τ = 1e-2 are biased enough to flip a sign (linear thin end). For comparison studies, τ ≤ 2.5e-3 is needed, and τ ≤ 6.25e-4 for thin layers (P3.6.2).
