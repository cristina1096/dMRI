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
